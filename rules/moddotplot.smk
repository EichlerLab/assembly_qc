"""
Adopted from acro_parm.smk created by Jiadong Lin
"""

import os
import pandas as pd
import numpy as np
import re
import glob
import pysam


# configfile: 'config/config_asm_qc.yaml'

# MANIFEST = config.get('MANIFEST', 'config/manifest.tab')

# raw_manifest_df = pd.read_csv(MANIFEST, sep='\t', comment='#', na_values=["","NA","na","N/A"])

# ## Universial conversion of manifest df

# add_haps = {"H2":"hap2", "UNASSIGNED":"unassigned"}
# df_transform = list()
# for idx, row in raw_manifest_df.iterrows():
#     df_transform.append({"SAMPLE": f"%s_hap1"%row["SAMPLE"], "ASM":row["H1"]}) # required

#     for add_hap in add_haps:
#         if (add_hap in raw_manifest_df.columns) and (not pd.isna(row[add_hap])):
#             df_transform.append({"SAMPLE": f"%s_%s"%(row["SAMPLE"], add_haps[add_hap]), "ASM": row[add_hap]})
        
# manifest_df = pd.DataFrame(df_transform)
# manifest_df.set_index("SAMPLE", inplace=True)

# BED = config.get('BED', 'config/moddot_regions.bed')
# bed_df = pd.read_csv(
#     BED, sep="\t", header=None, names=["chr", "start", "end"], dtype=str
# )

human_acro_bed = [
    ["chr13","0","22508596"],
    ["chr14","0","17708411"],
    ["chr15","0","22694466"],
    ["chr21","0","16306378"],
    ["chr22","0","20711065"]
]

bed_df = pd.DataFrame(human_acro_bed, columns = ["chr", "start", "end"], dtype=str)
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

def find_tigs_flag(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/work/find_tigs/flags/{wildcards.hap}.{region}_done" for region in bed_df.index ]

def find_pq_tigs(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/work/find_tigs/beds/{wildcards.hap}/{region}_pq_contig.bed" for region in bed_df.index ]

def find_p_tigs(wildcards): 
    return [ f"results/{wildcards.sample}/moddotplot/work/find_tigs/beds/{wildcards.hap}/{region}_p_contig.bed" for region in bed_df.index ]

def get_all_flag(wildcards):
    return [ f"results/{wildcards.sample}/moddotplot/work/selfplot/flags/{wildcards.hap}.{region}.done" for region in acros]

rule summarize_moddot_results:
    input:
        get_all_flag
    output:
        tsv = "results/{sample}/moddotplot/outputs/summary/{hap}.generated_acros.tsv"
    resources:
        mem=4,
        hrs=1,
    threads: 1
    run:
        header = ["Haplotype"]+acros
        called = [[wildcards.sample]]
        plot_dir=f"results/{wildcards.sample}/moddotplot/outputs/plots/{wildcards.hap}"
        for chrom in acros:
            region_name = bed_df[bed_df["chr"] == chrom].index[0]
            pdf_files = glob.glob(f"{plot_dir}/{region_name}_*_FULL.pdf")
            if len(pdf_files) > 0:
                contig_names = [os.path.basename(pdf).replace(f"{region_name}_","").replace("_FULL.pdf","") for pdf in pdf_files]
                called[0].append(",".join(contig_names))
            else:
                called[0].append("NA")
        result_df = pd.DataFrame(called, columns = header)
        result_df.to_csv(output.tsv, sep="\t", index=False)
                

checkpoint subset_target_region:
    output:
        bed="resources/acro_target_beds/{region}.bed",
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
        paf="results/{sample}/saffire/work/alignments/CHM13/pafs/{hap}.minimap2.paf"
    output:
        paf= "results/{sample}/moddotplot/work/liftover/CHM13/pafs/{hap}/{region}.paf",
    resources:
        mem=24,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell: """
        rustybam liftover --bed {input.bed} {input.paf} > {output.paf}
        """

rule trim_paf_moddot:
    input:
        paf= rules.liftover.output.paf
    output:
        paf = "results/{sample}/moddotplot/work/liftover/CHM13/trimmed_pafs/{hap}/{region}.paf"
    threads: 8
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        rustybam trim-paf {input.paf} > {output.paf}
        """

rule paf_stats:
    input:
        paf=rules.trim_paf_moddot.output.paf
    output:
        stats="results/{sample}/moddotplot/work/liftover/CHM13/paf_stats/{hap}/{region}.stats",
    resources:
        mem=12,
        hrs=2,
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    shell: """
        rustybam stats --paf --qbed {input.paf} > {output.stats}
        """

rule tag_contigs:
    input:
        paf = rules.trim_paf_moddot.output.paf,
        paf_stats = rules.paf_stats.output.stats
    output:
        flag = "results/{sample}/moddotplot/work/find_tigs/flags/{hap}.{region}.pq_contig.done"
    params:
        pq = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed",
        p = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_p_contig.bed",
        aln = "results/{sample}/moddotplot/work/find_tigs/pafs/{hap}/{region}_contig_aln.paf",
    threads: 1,
    resources:
        mem=10,
        hrs=24,
        disk_free=1
    run:
        os.makedirs(f"results/{wildcards.sample}/moddotplot/work/find_tigs/pafs/{wildcards.hap}", exist_ok=True)
        os.makedirs(f"results/{wildcards.sample}/moddotplot/work/find_tigs/beds/{wildcards.hap}", exist_ok=True)
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
        qry_chrom_tracker = {}
        with open(input.paf) as f:
            for line in f:
                try:
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

                    # if int(alignment_len) < 200000:
                    if int(alignment_len) <= 100000 or iden < 90:
                        continue
                    if query not in query2target2info_dict:
                        query2target2info_dict[query] = {'p': [], 'q': []}

                    if target_end - 1000000 >= centro_dict[target]['q'][1] and int(alignment_len) > 1000000:
                        arm = 'q'
                        query2target2info_dict[query][arm].append((
                        (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                        line.strip()))

                    if target_start <= centro_dict[target]['p'][1]:
                        arm = 'p'
                        query2target2info_dict[query][arm].append((
                        (query_start, query_end), strand, (target_start, target_end), target, num_matches, alignment_len, iden,
                        line.strip()))
                except:
                    continue


        query2arm2target2info_dict = {}

        fout_pq = open(params.pq,'w')
        fout_p = open(params.p, 'w')
        faln = open(params.aln, 'w')
        pq_contig_tracker = []
        for query, query2arm_dict in query2target2info_dict.items():
            try:
                query_len = seq2len_dict[query]
                ## contigs aligned to both p and q arm
                if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) > 0:
                    distal_aln = sorted(query2arm_dict['q'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                    qry_seq_start = 1
                    qry_seq_end = distal_aln[0][1]
                    target_chrom = distal_aln[3]
                    if distal_aln[1] == '-':
                        qry_seq_start = distal_aln[0][0]
                        qry_seq_end = query_len
                    # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{query}',file=fout_pq)
                    ## used for Arang's annotated assembly
                    new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                    print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{new_query}',file=fout_pq)
                    # pq_contig_tracker.append(query)

                ## contig only aligned to p arm
                if len(query2arm_dict['p']) > 0 and len(query2arm_dict['q']) == 0:
                    distal_aln = sorted(query2arm_dict['p'],key=lambda x: (x[2][0], x[2][1]), reverse=True)[0]
                    qry_seq_start = 1
                    qry_seq_end = distal_aln[0][1]
                    target_chrom = distal_aln[3]
                    if distal_aln[1] == '-':
                        qry_seq_start = distal_aln[0][0]
                        qry_seq_end = query_len
                    # print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{query}_{target_chrom}',file=fout_p)
                    ## used for Arang's annotated assembly
                    new_query = query.split('_')[1] + '_' + query.split('_')[0] if len(query.split('_')) > 1 else f'{query}_{target_chrom}'
                    print(f'{query}\t{qry_seq_start}\t{qry_seq_end}\t{wildcards.sample}_{new_query}',file=fout_p)

                for arm, qry2arm_list in query2arm_dict.items():
                    for aln_info in qry2arm_list:
                        print(aln_info[-1], file=faln)
            except:
                continue

        fout_pq.close()
        fout_p.close()
        faln.close()
        
        with open(str(output.flag), "w") as f_flag:
            print(f"{wildcards.sample}:{wildcards.region} Done.", file=f_flag)

rule get_pq_fa:
    input:
        flag = rules.tag_contigs.output.flag,
        hap = rules.rename_fasta.output.final_fasta
    output: 
        fa = "results/{sample}/moddotplot/work/get_pq_tigs/fasta/{hap}/{region}.pq_contig.fa",
    params:
        bed="results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed",
        region_name = get_region_name        
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    run:
        input_bed = params.bed
        input_hap = pysam.FastaFile(input.hap)
        output_fa = output.fa
        region_name = params.region_name
        output_fa_data = []
        if (os.path.isfile(input_bed)) and (os.path.getsize(input_bed) > 0):
            with open(input_bed) as finp_bed:
                finp_bed_lines = finp_bed.read().strip().split("\n")
            for bed_line in finp_bed_lines:
                token = bed_line.split("\t")
                contig_name, start, end = token[:3]
                start, end = int(start)+1, int(end)
                output_fa_seq_name = f">{region_name}_{contig_name}"
                subseq = input_hap.fetch(contig_name, start, end)
                output_fa_data += [output_fa_seq_name, subseq]
        with open(output_fa,"w") as fout:
            fout.write("\n".join(output_fa_data))
        if len(output_fa_data) > 0:
            shell(f"samtools faidx {output_fa}")

rule pq_selfplot:
    input:
        flag = rules.tag_contigs.output.flag,
        fa = rules.get_pq_fa.output.fa,
    output:
        flag = "results/{sample}/moddotplot/work/selfplot/flags/{hap}.{region}.done"
    params:
        bed = "results/{sample}/moddotplot/work/find_tigs/beds/{hap}/{region}_pq_contig.bed"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1,
    singularity:
        "docker://eichlerlab/moddotplot:0.9.0"
    shell: """
        outdir=$(dirname {output.flag} | sed "s/work\/selfplot\/flags/outputs\/plots\/{wildcards.hap}/g")
        if [[ ! -d $outdir ]];then
            mkdir -p $outdir
        fi
        if [[ -s {params.bed} ]]; then
            moddotplot static -f {input.fa} --no-hist --no-bed -o $outdir
            echo "{wildcards.sample} {wildcards.region}" > {output.flag}
        else
            echo "{wildcards.sample} {wildcards.region}" > {output.flag}
        fi
        """



# rule group_tigs:
#     input:
#         flag = find_tigs_flag,
#         paf = "results/{sample}/saffire/outputs/trimmed_pafs/CHM13/{hap}.minimap2.trimmed.paf",
#     output:
#         tigs = "results/{sample}/moddotplot/work/liftover/CHM13/contig_stats/{hap}.lifted_contigs.tsv",
#         tab = "results/{sample}/moddotplot/work/liftover/CHM13/contig_stats/{hap}.lifted_contigs.tab"
#     params:
#         p_tigs = find_p_tigs,
#         pq_tigs = find_pq_tigs,
#     threads: 1
#     resources:
#         mem=10,
#         hrs=24,
#         disk_free=1

#     run:
#         sample = wildcards.sample
#         input_pq_tigs = params.pq_tigs
#         input_p_tigs = params.p_tigs
#         input_paf = input.paf
        
#         acros = ['chr13', 'chr14', 'chr15', 'chr21', 'chr22']
#         fout = open(output.tigs, 'w')
#         print('Sample\tArm\tTig\tTig_start\tTig_end\tSize\tChrom',file=fout)
#         haps = open(output.tab, 'w')
#         print('SAMPLE\tHAP', file=haps)
        
#         qry_chrom_tracker = {}
#         qrylen_dict = {}
#         with open(input_paf) as f:
#             lines = f.read().strip().split("\n")
#         for line in lines:
#             query, query_len, query_start, query_end, strand, target, target_len, target_start, target_end, num_matches, alignment_len = line.strip().split("\t")[:11]
#             qrylen_dict[query] = query_len
#             if int(alignment_len) < 200000:
#                 continue
#             if query not in qry_chrom_tracker:
#                 qry_chrom_tracker[query] = []
#             if not target in qry_chrom_tracker[query]:
#                 qry_chrom_tracker[query].append(target)

#         for pq_tig in input_pq_tigs:
#             if not os.path.getsize(pq_tig)>0:
#                 continue
#             print (pq_tig)
#             region = os.path.basename(pq_tig).split("_")[0]
#             print (region)
#             with open(pq_tig) as finp_pq_tig:
#                 print (pq_tig)
#                 pq_tig_lines = finp_pq_tig.read().strip().split("\n")
#             for pq_tig_line in pq_tig_lines:
#                 print (pq_tig_line)
#                 entries = pq_tig_line.strip().split('\t')
#                 print (entries)
#                 print (["entry",entries[0]])
#                 print ("qrylen_dict",qrylen_dict)
#                 qlen = qrylen_dict[entries[0]]
#                 print (qlen)
#                 print(f"{sample}\tpq\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}")
#                 print (f'moddotplot/fasta/{sample}/{region}_pq_contig.fa')
#                 print(f"{sample}\tpq\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)
#                 print(pq_tig, file=haps)

#         for p_tig in input_p_tigs:
#             if not os.path.getsize(p_tig)>0:
#                 continue
#             with open(p_tig) as finp_p_tig:
#                 p_tig_lines = finp_p_tig.read().strip().split("\n")
#             for p_tig_line in p_tig_lines:
#                 entries = p_tig_line.strip().split('\t')
#                 qlen = qrylen_dict[entries[0]]
#                 not_acro = False
#                 for ele in qry_chrom_tracker[entries[0]]:
#                     if ele not in acros:
#                         not_acro = True
#                         break
#                 if not not_acro:
#                     print(f"{sample}\tp\t{entries[0]}\t{entries[1]}\t{entries[2]}\t{qlen}\t{region}", file=fout)

#         fout.close()
#         haps.close()
