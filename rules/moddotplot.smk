"""
Adopted from liftover_acro_parm.smk created by Jiadong Lin
"""

import os
import pandas as pd
import numpy as np
import re
import glob


configfile: 'config/config_asm_qc.yaml'

MANIFEST = config.get('MANIFEST', 'config/manifest.tab')

raw_manifest_df = pd.read_csv(MANIFEST, sep='\t')

## Universial conversion of manifest df

add_haps = {"H2":"hap2", "UNASSIGNED":"unassigned"}
df_transform = list()
for idx, row in raw_manifest_df.iterrows():
    df_transform.append({"SAMPLE": f"%s_hap1"%row["SAMPLE"], "ASM":row["H1"]}) # required

    for add_hap in add_haps:
        if (add_hap in raw_manifest_df.columns) and (not pd.isna(row[add_hap])):
            df_transform.append({"SAMPLE": f"%s_%s"%(row["SAMPLE"], add_haps[add_hap]), "ASM": row[add_hap]})
        
manifest_df = pd.DataFrame(df_transform)
manifest_df.set_index("SAMPLE", inplace=True)

BED = config.get('BED', 'config/moddot_regions.bed')
bed_df = pd.read_csv(
    BED, sep="\t", header=None, names=["chr", "start", "end"], dtype=str
)

bed_df["NAME"] = bed_df["chr"] + "_" + (bed_df["start"].astype(int) + 1).astype(str) + "_" + bed_df["end"]
bed_df.set_index("NAME", drop=True, inplace=True)

acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']

try:
    REF = config["REF"]["CHM13"]["PATH"]
except KeyError:
    REF = "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta"

# rustybam/0.1.16

def get_region_name(wildcards):
    return bed_df[bed_df["chr"] == wildcards.region].index[0]

def cigar_tuple(cigar):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    opVals = re.findall(r'(\d+)([\w=])', cigar)
    lengths = [int(opVals[i][0]) for i in range(0, len(opVals))]
    ops = [opVals[i][1] for i in range(0, len(opVals))]
    tuples = []
    for i in range(len(lengths)):
        tuples.append([int(lengths[i]), ops[i]])
    return tuples

def find_tigs(wildcards):
    tig_dict = {}
    for sample in manifest_df.index:
        for region in bed_df.index:
            tig_dict[f'{sample}_{region}_pq'] = rules.tag_contigs.output.pq.format(sample=sample,region=region)
            tig_dict[f'{sample}_{region}_p'] = rules.tag_contigs.output.p.format(sample=sample,region=region)
    return tig_dict

def get_all_flag(wildcards):
    return [ f"moddotplot/results/{wildcards.sample}/.{region}.done" for region in acros]

rule summarize_moddot_results:
    input:
        get_all_flag
    output:
        tsv = "moddotplot/results/{sample}.generated_acros.tsv"
    resources:
        mem=4,
        hrs=1,
    threads: 1
    run:
        header = ["Haplotype"]+acros
        called = [[wildcards.sample]]
        plot_dir=f"moddotplot/results/{wildcards.sample}"
        for chrom in acros:
            region_name = bed_df[bed_df["chr"] == chrom].index[0]
            try:
                pdf = glob.glob(f"{plot_dir}/{region_name}-*_FULL.pdf")[0]
                called[0].append(os.path.basename(pdf).replace(f"{region_name}-","").replace("_FULL.pdf",""))
            except IndexError:
                called[0].append("NA")
        result_df = pd.DataFrame(called, columns = header)
        result_df.to_csv(output.tsv, sep="\t", index=False)
                

checkpoint subset_target_region:
    input:
        bed=BED,
    output:
        bed="moddotplot/target_beds/{region}.bed",
    resources:
        mem=4,
        hrs=2,
    threads: 1
    run:
        out_df = bed_df[bed_df["chr"] == wildcards.region]
        out_df[["chr", "start", "end"]].to_csv(output.bed, sep="\t", header=False, index=False)

rule liftover:
    input:
        bed=rules.subset_target_region.output.bed,
        paf='saffire/CHM13/results/{sample}/alignments/{sample}.minimap2.paf',
    output:
        paf= temp("moddotplot/tmp/liftover/{sample}.{region}.paf"),
    resources:
        mem=24,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        rustybam liftover --bed {input.bed} {input.paf} > {output.paf}
        """


rule trim_paf:
    input:
        paf="moddotplot/tmp/liftover/{sample}.{region}.paf",
    output:
        paf="moddotplot/liftover/paf/{sample}/{region}.trimmed.paf",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        rustybam trim-paf {input.paf} > {output.paf}
        """

rule paf_stats:
    input:
        paf="moddotplot/liftover/paf/{sample}/{region}.trimmed.paf",
    output:
        stats="moddotplot/liftover/stats/{sample}/{region}.trimmed.stats",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell:
        """
        rustybam stats --paf --qbed {input.paf} > {output.stats}
        """

rule group_tigs:
    input:
        tigs = unpack(find_tigs),
        paf = "saffire/CHM13/results/{sample}/alignments/{sample}.minimap2.paf",
    output:
        tigs = "moddotplot/liftover/stats/{sample}/lifted_contigs.tsv",
    threads: 1
    resources:
        mem=10,
        hrs=24,
        disk_free=1

    run:
        acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']

        fout = open(output.tigs, 'w')
        print('Sample\tArm\tTig\tTig_start\tTig_end\tSize\tChrom', file=fout)
        for sample in manifest_df.index:
            qry_chrom_tracker = {}
            qrylen_dict = {}
            with open(input.paf) as f:
                for line in f:
                    query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]
                    qrylen_dict[query] = query_len
                    if int(alignment_len) < 200000:
                        continue
                    if query not in qry_chrom_tracker:
                        qry_chrom_tracker[query] = []
                    qry_chrom_tracker[query].append(target)

            for region in bed_df.index:
                for line in open(input.tigs[f'{sample}_{region}_pq']):
                    entries = line.strip().split('\t')
                    qlen = qrylen_dict[entries[0]]
                    print(f"{sample}\tpq\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)

                for line in open(input.tigs[f'{sample}_{region}_p']):
                    entries = line.strip().split('\t')
                    qlen = qrylen_dict[entries[0]]
                    not_acro = False
                    for ele in qry_chrom_tracker[entries[0]]:
                        if ele not in acros:
                            not_acro = True
                            break
                    if not not_acro:
                        print(f"{sample}\tp\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)


rule tag_contigs:
    input:
        paf = "moddotplot/liftover/paf/{sample}/{region}.trimmed.paf"
    output:
        pq = "moddotplot/contigs/{sample}/{region}_pq_contig.bed",
        p = "moddotplot/contigs/{sample}/{region}_p_contig.bed",
        aln = "moddotplot/contigs/{sample}/{region}_contig_aln.paf"
    threads: 1,
    resources:
        mem=10,
        hrs=24,
        disk_free=1
    run:
        centro_dict = {
                'chr13': {
                    'p' : (15547593, 16522942),
                    'q' : (16522942, 17498291),
                },
                'chr14': {
                    'p' : (10092112, 11400261),
                    'q' : (11400261, 12708411),
                },
                'chr15': {
                    'p' : (16678794, 17186630),
                    'q' : (17186630, 17694466),
                },
                'chr21': {
                    'p' : (10962853, 11134529),
                    'q' : (11134529, 11306205),
                },
                'chr22': {
                    'p' : (12788180, 14249622),
                    'q' : (14249622, 15711065),
                }
            }

        seq2len_dict = {}
        query2target2info_dict = {}
        ref_tag = wildcards.region.split('_')[0]
        sample_name = wildcards.sample.split('_')[0]
        qry_chrom_tracker = {}
        with open(input.paf) as f:
            for line in f:
                query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]

                cg = [i.split(":")[-1] for i in line.strip().split("\t")[12:] if i[:2] == 'cg'][0]
                cg_tuple = cigar_tuple(cg)
                iden = round((sum([int(i[0]) for i in cg_tuple if i[1] == '=' or i[1] == 'M']) / sum(
                    [int(i[0]) for i in cg_tuple if i[1] in {'=', 'M', 'X', 'D', 'I'}])) * 100,2)

                query_len = int(query_len)
                query_start = int(query_start)
                query_end = int(query_end)
                target_len = int(target_len)
                target_start = int(target_start)
                target_end = int(target_end)

                seq2len_dict[query] = query_len
                seq2len_dict[target] = target_len

                if query not in qry_chrom_tracker:
                    qry_chrom_tracker[query] = []
                qry_chrom_tracker[query].append(target)


                if int(alignment_len) < 200000:
                    continue
                if query not in query2target2info_dict:
                    query2target2info_dict[query] = {'p': [], 'q': []}

                if target_end - 1000000 >= centro_dict[target]['q'][1]:
                    arm = 'q'
                    query2target2info_dict[query][arm].append((
                    (query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, iden,
                    line.strip()))
                elif target_end <= centro_dict[target]['p'][0]:
                    arm = 'p'
                    query2target2info_dict[query][arm].append((
                    (query_start, query_end), strand, (target_start, target_end), num_matches, alignment_len, iden,
                    line.strip()))


        query2arm2target2info_dict = {}

        fout_pq = open(output.pq,'w')
        fout_p = open(output.p, 'w')
        faln = open(output.aln, 'w')
        pq_contig_tracker = []
        for query, query2arm_dict in query2target2info_dict.items():
            query_len = seq2len_dict[query]
            ## contigs aligned to both p and q arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) > 0:
                distal_aln = sorted(query2arm_dict['q'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{query}',file=fout_pq)
                pq_contig_tracker.append(query)

            ## contig only aligned to p arm
            if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) == 0:
                distal_aln = sorted(query2arm_dict['p'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                qry_seq_start = 1
                qry_seq_end = distal_aln[0][1]
                if distal_aln[1] == '-':
                    qry_seq_start = distal_aln[0][0]
                    qry_seq_end = query_len
                print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{query}',file=fout_p)

            for arm, qry2arm_list in query2arm_dict.items():
                for aln_info in qry2arm_list:
                    print(aln_info[-1], file=faln)


rule get_pq_fa:
    input:
        bed="moddotplot/contigs/{sample}/{region}_pq_contig.bed",
        hap = "fcs_cleaned_fasta/{sample}/{sample}.fasta"
    output:
        fa="moddotplot/fasta/{sample}/{region}_pq_contig.fa"
    params:
        region_name = get_region_name        
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        if [ -s {input.bed} ]; then
            echo ">{params.region_name}-$(awk '{{print $1}}' {input.bed})" > {output.fa}
            bedtools getfasta -fi {input.hap} -bed {input.bed} | tail -n +2 >> {output.fa}
            samtools faidx {output.fa}
        else
            touch {output.fa}
        fi
        """

rule pq_selfplot:
    input:
        bed = "moddotplot/contigs/{sample}/{region}_pq_contig.bed",
        fa = "moddotplot/fasta/{sample}/{region}_pq_contig.fa"
    output:
        flag = "moddotplot/results/{sample}/.{region}.done"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1,
    singularity:
        "docker://eichlerlab/moddotplot:0.9.0"
    shell:
        """
        if [ -s {input.bed} ]; then
            moddotplot static -f {input.fa} --no-hist --no-bed -o $( dirname {output.flag} )
            touch {output.flag}
        else
            touch {output.flag}
        fi
        """