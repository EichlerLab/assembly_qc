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
GENOME_INDEX = "/net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/chrom_index.txt"
REF = config.get('REF')
REF_NAME = REF.split('/')[-1].split('.fa')[0]
ALIGNER = config.get('ALIGNER',"minimap2")

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
manifest_df.set_index("SAMPLE",inplace=True)

#-----------------------------------------


scattergather:
    split=PARTS

wildcard_constraints:
    sample='|'.join(manifest_df.index),


def find_fasta(wildcards):
    return f"QC_results/fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

def find_fasta_index(wildcards):
    return f"QC_results/fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta.fai"

def find_contigs(wildcards):
    return gather.split('QC_results/saffire/tmp/{sample}.{aligner}.{scatteritem}.broken.paf', sample=wildcards.sample, aligner=ALIGNER)

def find_scatter_sams(wildcards):
    return gather.split("QC_results/saffire/tmp/{sample}.winnowmap.{scatteritem}.sam", sample=wildcards.sample)



localrules: 
    all,

rule all:
    input:
        expand("QC_results/saffire/results/{sample}/{sample}.{aligner}.saf", sample=manifest_df.index, aligner=ALIGNER)

rule gather_bams:
    input:
        expand("QC_results/saffire/results/{sample}/alignments/{sample}.{aligner}.bam", sample=manifest_df.index, aligner=ALIGNER)

rule gather_beds:
    input:
        expand("QC_results/saffire/results/{sample}/beds/{sample}.{aligner}.bed", sample=manifest_df.index, aligner=ALIGNER)

rule gather_pafs:
    input:
        expand("QC_results/saffire/results/{sample}/alignments/{sample}.{aligner}.paf", sample=manifest_df.index, aligner=ALIGNER)

rule get_meryl_db:
    input:
        ref=REF
    output:
        meryl_db_name = directory(f"QC_results/saffire/merylDB/{REF_NAME}"),
    threads: 1
    singularity:
        "docker://eichlerlab/merqury/1.3.1"
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
        kmer_count = f"QC_results/saffire/merylDB/repetitive_k19_{REF_NAME}.txt"
    threads: 1
    singularity:
        "docker://eichlerlab/merqury/1.3.1"
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
                "QC_results/saffire/tmp/{{sample}}.{scatteritem}.batch",
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
        batch_file="QC_results/saffire/tmp/{sample}.{scatteritem}.batch",
    output:
        scatter_fasta=temp(
            "QC_results/saffire/tmp/{sample}.{scatteritem}.fasta"
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


rule make_winnowmap_scatter_paf:
    input:
        ref = REF,
        kmer_count = rules.get_kmer_count.output.kmer_count,
        fa = "QC_results/saffire/tmp/{sample}.{scatteritem}.fasta"
    output:
        scatter_paf=temp(
            "QC_results/saffire/tmp/{sample}.winnowmap.{scatteritem}.aln.paf"
        ),
    benchmark: "QC_results/saffire/benchmarks/{sample}.winnowmap_paf.{scatteritem}.benchmark.txt",
    resources:
        mem=lambda wildcards, attempt: max(attempt * 2, 16),
        hrs=168,
        disk_free=5,
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    threads: 12,
    shell:
        """
        winnowmap -W {input.kmer_count} -c -t {threads} -K {resources.mem}G --cs {params.map_opts} --secondary=no --eqx -Y {input.ref} {input.fa} > {output.scatter_paf}
        """


rule make_winnowmap_scatter_sam:
    input:
        ref = REF,
        kmer_count = rules.get_kmer_count.output.kmer_count,
        fa = "QC_results/saffire/tmp/{sample}.{scatteritem}.fasta"
    output:
        scatter_sam=temp(
            "QC_results/saffire/tmp/{sample}.winnowmap.{scatteritem}.sam"
        ),
    benchmark: "QC_results/saffire/benchmarks/{sample}.winnowmap_sam.{scatteritem}.benchmark.txt",        
    resources:
        mem=lambda wildcards, attempt: max(attempt * 2, 16),
        hrs=168,
        disk_free=5,
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    threads: 8
    shell:
        """
        winnowmap -W {input.kmer_count} -c -t {threads} -K {resources.mem}G --cs -a {params.map_opts} --MD --secondary=no --eqx -Y {input.ref} {input.fa} > {output.scatter_sam}
        """

rule combine_winnowmap_scatter_sam:
    input:
        scatter_sams=find_scatter_sams,
    output:
        bam="QC_results/saffire/results/{sample}/alignments/{sample}.winnowmap.bam",
    resources:
        mem=12,
        hrs=48,
        disk_free=5,
    singularity:
        "docker://eichlerlab/binf-basics:0.1"
    threads: 8
    shell:
        """
        awk '!seen[$0]++' {input.scatter_sams} | samtools view -bS - | samtools sort -@ {threads} - > {output.bam}
        """

rule combine_winnowmap_paf:
    input:
        scatter_pafs=gather.split(
            "QC_results/saffire/tmp/{sample}.winnowmap.{scatteritem}.aln.paf",
            sample=manifest_df.index
        ),
    output:
        paf="QC_results/saffire/results/{sample}/alignments/{sample}.winnowmap.paf"
    resources:
        mem=8,
        hrs=12,
        disk_free=1,
    threads: 1
    shell:
        """
        awk '!seen[$0]++' {input.scatter_pafs} > {output.paf}
        """	

rule make_minimap_paf:
    input:
        fa = find_fasta,
        ref = REF
    output:
        paf = "QC_results/saffire/results/{sample}/alignments/{sample}.minimap2.paf",
    benchmark: "QC_results/saffire/benchmarks/{sample}.minimap2_paf.benchmark.txt",
    threads: 8
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

rule make_minimap_reverse_paf:
    input:
        ref = find_fasta,
        fa = REF
    output:
        paf = "QC_results/saffire/results/{sample}/alignments/{sample}.minimap2.rev.paf",
    benchmark: "QC_results/saffire/benchmarks/{sample}.minimap2_rev_paf.benchmark.txt",
    threads: 8
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
        ref = REF
    output:
        bam = "QC_results/saffire/results/{sample}/alignments/{sample}.minimap2.bam"
    benchmark: "QC_results/saffire/benchmarks/{sample}.minimap2_bam.benchmark.txt",
    threads: 8
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources: 
        mem = 60,
        hrs = 120
    shell:
        '''
        minimap2 -c -t {threads} -K {resources.mem}G --cs -a {params.map_opts} --MD --secondary=no --eqx -Y {input.ref} {input.fa} | samtools view -S -b /dev/stdin | samtools sort -@ {threads} /dev/stdin > {output.bam}
        '''

rule make_bed:
    input:
        paf = "QC_results/saffire/results/{sample}/alignments/{sample}.{aligner}.paf",
    output:
        bed = "QC_results/saffire/results/{sample}/beds/{sample}.{aligner}.bed",
    params:
        genome_index = GENOME_INDEX
    threads: 1
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources: 
        mem = 12,
        hrs = 1
    shell:
        '''
        awk -vOFS="\t" '{{print $6,$8,$9,$1,$12}}' {input.paf} | bedtools sort -g {params.genome_index} -i - > {output.bed}
        '''


rule split_paf:
    input:
        paf = "QC_results/saffire/results/{sample}/alignments/{sample}.{aligner}.paf",
    output:
        flag = temp(scatter.split("QC_results/saffire/tmp/{{sample}}.{{aligner}}.{scatteritem}.paf")),
        temp_paf = temp('QC_results/saffire/tmp/{sample}_uniform.{aligner}.paf')
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
            df.loc[df[0] == contig][col_out].to_csv(f'QC_results/saffire/tmp/{wildcards.sample}.{ALIGNER}.{out_num}-of-{PARTS}.paf', sep='\t', index=False, header=False, mode='a+')

rule trim_break_orient_paf:
    input:
        paf = 'QC_results/saffire/tmp/{sample}.{aligner}.{scatteritem}.paf'
    output:
        contig = temp('QC_results/saffire/tmp/{sample}.{aligner}.{scatteritem}.orient.paf'),
        broken = temp('QC_results/saffire/tmp/{sample}.{aligner}.{scatteritem}.broken.paf')
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources: 
        mem = 24,
        hrs = 24
    shell:
        '''
        rustybam orient {input.paf} | rustybam trim-paf | rb filter --paired-len 1000000 > {output.contig}
        rustybam break-paf --max-size {SV_SIZE} {output.contig} > {output.broken}
        '''

rule combine_paf:
    input:
        paf = find_contigs,
        flag = rules.split_paf.output.flag
    output:
        paf = 'QC_results/saffire/tmp/{sample}.{aligner}.broken.paf'
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
        saf = 'QC_results/saffire/results/{sample}/{sample}.{aligner}.saf'
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

