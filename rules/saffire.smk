import pandas as pd
from pybedtools import BedTool
import numpy as np
import networkx as nx
import more_itertools as mit


configfile: 'config/config_asm_qc.yaml'

PARTS=config.get('PARTS', 30)
# MINIMAP_PARAMS = config.get('MINIMAP_PARAMS', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
MINIMAP_PARAMS = config.get('MINIMAP_PARAMS', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
MANIFEST = config.get('MANIFEST', 'config/manifest.tab')
SV_SIZE = config.get('SV_SIZE', '30000')

REF_DICT = config["REF"]
ALIGNER = config.get('ALIGNER',"minimap2")


raw_manifest_df = pd.read_csv(MANIFEST, sep='\t', comment='#', na_values=["","NA","na","N/A"])

## Universial conversion of manifest df

add_haps = {"H2":"hap2", "UNASSIGNED":"unassigned"}
df_transform = list()
for idx, row in raw_manifest_df.iterrows():
    df_transform.append({"SAMPLE": f"%s_hap1"%row["SAMPLE"], "ASM":row["H1"]}) # required

    for add_hap in add_haps:
        if (add_hap in raw_manifest_df.columns) and (not pd.isna(row[add_hap])):
            df_transform.append({"SAMPLE": f"%s_%s"%(row["SAMPLE"], add_haps[add_hap]), "ASM": row[add_hap]})
        
manifest_df = pd.DataFrame(df_transform)
manifest_df.set_index("SAMPLE",inplace=True)

#-----------------------------------------


scattergather:
    split=PARTS

def find_ref_fa(wildcards):
    return REF_DICT[wildcards.ref]["PATH"]

def find_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

def find_fasta_index(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta.fai"

def find_contigs(wildcards):
    return gather.split('saffire/{ref}/tmp/{sample}.{aligner}.{scatteritem}.broken.paf', ref=wildcards.ref, sample=wildcards.sample, aligner=wildcards.aligner)

def find_all_paf(wildcards):
    avail_pafs = []
    haps = ["hap1","hap2","unassigned"]
    values = raw_manifest_df.loc[raw_manifest_df["SAMPLE"] == wildcards.asm][["H1","H2","UNASSIGNED"]].iloc[0].tolist()
    for idx, hap in enumerate(haps):
        if not str(values[idx]) == "nan":
            avail_pafs.append (f'saffire/{wildcards.ref}/results/{wildcards.asm}_{hap}/alignments/{wildcards.asm}_{hap}.{wildcards.aligner}.paf')
    return avail_pafs


checkpoint copy_ref_fasta:
    input:
        raw_ref=lambda wildcards: config["REF"][wildcards.ref]["PATH"],
    output:
        ref = "saffire/{ref}/reference/{ref}.fa",
        fai = "saffire/{ref}/reference/{ref}.fa.fai",
        genome_index = "saffire/{ref}/reference/{ref}.genome.txt"
    threads: 1,
    resources:
        mem = 4,
        hrs = 1,
    shell:
        """
        cp {input.raw_ref} {output.ref}
        samtools faidx {output.ref}
        cut -f1 {output.fai}> {output.genome_index}
        """

rule get_meryl_db:
    input:
        ref="saffire/{ref}/reference/{ref}.fa",
    output:
        meryl_db_name = directory("saffire/{ref}/merylDB/{ref}"),
    threads: 1
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    resources:
        mem = 16,
        hrs = 72
    shell:
        '''
        meryl count k=19 output {output.meryl_db_name} {input.ref}
        '''

rule get_kmer_count:
    input:
        meryl_db_name = rules.get_meryl_db.output.meryl_db_name
    output:
        kmer_count = "saffire/{ref}/merylDB/repetitive_k19_{ref}.txt"
    threads: 1
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    resources:
        mem = 8,
        hrs = 72
    shell:
        '''
        meryl print greater-than distinct=0.9998 {input.meryl_db_name} > {output.kmer_count}
        '''

rule get_batch_ids:
    input:
        fai=find_fasta_index,
    output:
        batches=temp(
            scatter.split(
                "saffire/tmp/{{sample}}.{scatteritem}.batch",
            )
        ),
    resources:
        mem=4,
        hrs=5,
        disk_free=5,
    threads: 1
    run:
        NIDS = len(output.batches)

        batch_dict = {}

        for i in range(NIDS):
            batch_dict[i] = []

        with open(input.fai, "r") as infile:
            fai_list = [line.split("\t")[0] for line in infile]


        for j in range(len(fai_list)):
            batch_dict[j % NIDS].append(fai_list[j])

        outs = [open(f, "w+") for f in output.batches]

        for i in range(NIDS):
            outs[i].write("\n".join(batch_dict[i]) + "\n")
            outs[i].close()

rule split_fasta:
    input:
        fasta=find_fasta,
        batch_file="saffire/tmp/{sample}.{scatteritem}.batch",
    output:
        scatter_fasta=temp(
            "saffire/tmp/{sample}.{scatteritem}.fasta"
        ),
    resources:
        mem=4,
        hrs=24,
        disk_free=1,
    threads: 4
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    shell:
        """
        samtools faidx {input.fasta} -r {input.batch_file} > {output.scatter_fasta}
        """

rule make_minimap_paf:
    input:
        fa = find_fasta,
        ref = find_ref_fa
    output:
        paf = "saffire/{ref}/results/{sample}/alignments/{sample}.minimap2.paf",
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',
    benchmark: "saffire/{ref}/benchmarks/{sample}.minimap2_paf.benchmark.txt",
    threads: 12
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 120
    shell:
        '''
        minimap2 -c -t {threads} -K {resources.mem}G --cs {params.map_opts} --secondary=no --eqx -Y {input.ref} {input.fa} > {output.paf}
        '''

rule make_minimap_bam:
    input:
        fa = find_fasta,
        ref = find_ref_fa
    output:
        bam = "saffire/{ref}/results/{sample}/alignments/{sample}.minimap2.bam",
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',        
    benchmark: "saffire/{ref}/benchmarks/{sample}.minimap2_bam.benchmark.txt",
    threads: 12
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 120
    shell:
        '''
        minimap2 -c -t {threads} -K {resources.mem}G --cs -a {params.map_opts} --MD --secondary=no --eqx -Y {input.ref} {input.fa} | samtools view -S -b /dev/stdin | samtools sort -@ {threads} -o {output.bam}
        '''

rule make_bed:
    input:
        paf = "saffire/{ref}/results/{sample}/alignments/{sample}.{aligner}.paf",
        genome_index = "saffire/{ref}/reference/{ref}.genome.txt"
    output:
        bed = "saffire/{ref}/results/{sample}/beds/{sample}.{aligner}.bed",
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',              
    threads: 1
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 1
    shell:
        '''
        awk -vOFS="\t" '{{print $6,$8,$9,$1,$12}}' {input.paf} | bedtools sort -g {input.genome_index} > {output.bed}
        '''

rule split_paf:
    input:
        paf = "saffire/{ref}/results/{sample}/alignments/{sample}.{aligner}.paf",
    output:
        flag = temp(scatter.split("saffire/{{ref}}/tmp/{{sample}}.{{aligner}}.{scatteritem}.paf")),
        temp_paf = temp("saffire/{ref}/tmp/{sample}_uniform.{aligner}.paf")
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    run:
        with open(input.paf, 'r') as infile:
            for i, line in enumerate(infile):
                if i == 0:
                    all_tags = set([x.split(':')[0] for x in line.split('\t')[12:]])
                else:
                    all_tags = all_tags.intersection(set([x.split(':')[0] for x in line.split('\t')[12:]]))

        out_list = []

        with open(input.paf, 'r') as infile:
            for line in infile:
                out_list.append(line.split('\t')[0:12]+[x for x in line.split('\t')[12:] if x.split(':')[0] in all_tags])

        with open(output.temp_paf, 'w') as outfile:
            for item in out_list:
                outfile.write('\t'.join(item))

        df = pd.read_csv(output.temp_paf, sep='\t', low_memory=False, header=None)
        col_out = df.columns.values
        for i, contig in enumerate(df[0].unique()):
            out_num = (i % PARTS) + 1
            df.loc[df[0] == contig][col_out].to_csv(f'saffire/{wildcards.ref}/tmp/{wildcards.sample}.{ALIGNER}.{out_num}-of-{PARTS}.paf', sep='\t', index=False, header=False, mode='a+')

rule trim_break_orient_paf:
    input:
        paf = 'saffire/{ref}/tmp/{sample}.{aligner}.{scatteritem}.paf'
    output:
        contig = temp('saffire/{ref}/tmp/{sample}.{aligner}.{scatteritem}.orient.paf'),
        broken = temp('saffire/{ref}/tmp/{sample}.{aligner}.{scatteritem}.broken.paf')
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    params:
        sv_size = SV_SIZE,
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 24,
        hrs = 24
    shell:
        '''
        rustybam orient {input.paf} | rustybam trim-paf | rb filter --paired-len 1000000 > {output.contig}
        rustybam break-paf --max-size {params.sv_size} {output.contig} > {output.broken}
        '''
rule concat_pafs:
    input:
        paf = find_all_paf
    output:
        paf = 'saffire/{ref}/results/merged_paf/{aligner}/{asm}.concat.paf'
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads:1 
    resources:
        mem=4,
        hrs=1
    shell:
        '''
        cat {input.paf} > {output.paf}
        '''

rule combine_paf:
    input:
        paf = find_contigs,
        flag = rules.split_paf.output.flag
    output:
        paf = 'saffire/{ref}/tmp/{sample}.{aligner}.broken.paf'
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    shell:
        '''
        cat {input.paf} > {output.paf}
        '''

rule saff_out:
    input:
        paf = rules.combine_paf.output.paf
    output:
        saf = 'saffire/{ref}/results/{sample}/saff/{sample}.{aligner}.saf'
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 8,
        hrs = 24
    shell:
        '''
        rb stats --paf {input.paf} > {output.saf}
        '''

rule check_covered_chromosomes:
    input:
        paf = "saffire/{ref}/results/{sample}/alignments/{sample}.{aligner}.paf"
    output:
        tsv = "saffire/{ref}/results/{sample}/beds/{sample}.{aligner}.chrom_cov.tsv"
    wildcard_constraints:
        sample='|'.join(manifest_df.index),
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    run:
        from collections import defaultdict

        chrom_length_dict = dict()
        chrom_intervals = defaultdict(list)
        with open(input.paf) as finp:
            read_lines = finp.read().strip().split("\n")
        for line in read_lines:
            token = line.split("\t")[:9]
            chrom, start, end, contigname, chrom_length = token[5], int(token[7]), int(token[8]), token[0], int(token[6])
            if not chrom in chrom_length_dict:
                chrom_length_dict[chrom] = chrom_length
            chrom_intervals[chrom].append((start,end))
        
        merged_dict = dict()
        for chrom in sorted(chrom_intervals):
            intervals = sorted(chrom_intervals[chrom])
            merged = []
            current_start, current_end = intervals[0]
            for start, end in intervals[1:]:
                if start <= current_end:
                    current_end = max(current_end, end)
                else:
                    merged.append((current_start, current_end))
                    current_start, current_end = start, end
            merged.append((current_start, current_end))
            merged_dict[chrom] = merged
        chroms = [f"chr{i}" for i in range(1,23)]+["chrX","chrY"]
        with open(output.tsv, "w") as fout:
            print ("chrom\tcovered_pct", file=fout)
            for chrom in chroms:
                try:
                    covered_length = 0
                    chrom_length = chrom_length_dict[chrom]
                    for intervals in merged_dict[chrom]:
                        covered_length += (intervals[1]-intervals[0])
                    covered_rate = covered_length/chrom_length*100

                except:
                    covered_rate = 0.0
                print (f"{chrom}\t{covered_rate:.2f}", file=fout)
                
            

        