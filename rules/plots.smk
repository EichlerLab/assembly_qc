import pandas as pd
import os
import sys


configfile: 'config/config_asm_qc.yaml'
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
IDEO_PLOT_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/ideo_plot.py"
N50_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/n50"
PLOIDY_PLOT_SCRIPT = f"{SNAKEMAKE_DIR}/../scripts/ploidy.R"
ALIGNER = config.get('ALIGNER',"minimap2")


def find_fasta(wildcards):
    return f"QC_results/fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

def find_contig_fasta(wildcards):
    return f"QC_results/fcs_cleaned_fasta/{wildcards.sample}/contig_fasta/{wildcards.sample}.fasta"

def get_asm_manifest_df(manifest_df):

    ## Universial conversion of manifest df

    add_haps = {"H2":"hap2", "UNASSIGNED":"unassigned"}
    df_transform = list()
    for idx, row in manifest_df.iterrows():
        df_transform.append({"SAMPLE": f"%s_hap1"%row["SAMPLE"], "ASM":row["H1"]}) # required

        for add_hap in add_haps:
            if (add_hap in manifest_df.columns) and (not pd.isna(row[add_hap])):
                df_transform.append({"SAMPLE": f"%s_%s"%(row["SAMPLE"], add_haps[add_hap]), "ASM": row[add_hap]})

    return pd.DataFrame(df_transform)

full_manifest_df = pd.read_csv(MANIFEST, header=0, sep='\t', comment='#').fillna('NA')
conv_manifest_df = get_asm_manifest_df(full_manifest_df)

full_manifest_df.set_index("SAMPLE",inplace=True) ## manifest df for merqury
conv_manifest_df.set_index("SAMPLE",inplace=True) ## shortened df using only SAMPLE / ASM fields


rule gather_plots:
    input:
        expand("QC_results/contig_stats/n50/{sample}.n50.stats", sample = conv_manifest_df.index),
        expand("QC_results/contig_stats/plots/{sample}.dist.log.png", sample = conv_manifest_df.index),
        expand("QC_results/contig_stats/plots/{sample}.scatter.log.png", sample = conv_manifest_df.index),
        expand("QC_results/ideo_plots/{aligner}/{asm}.ideoplot.pdf", asm = full_manifest_df.index, aligner = ALIGNER),
        expand("QC_results/ploidy_plots/{aligner}/{asm}.summary.txt", asm = full_manifest_df.index, aligner = ALIGNER),


rule make_contig_stats:
    input:
        asm_fasta = find_fasta,
    output:
        n50_stats = "QC_results/contig_stats/n50/{sample}.n50.stats",
        dist_plot = "QC_results/contig_stats/plots/{sample}.dist.log.png",
        scatter_plot = "QC_results/contig_stats/plots/{sample}.scatter.log.png",
        flag = "QC_results/contig_stats/{sample}.contig_stats.done"
    threads:
        1
    resources:
        mem = 8,
        hrs = 12,
    params:
        script_path = N50_SCRIPT
    shell:
        '''
        {params.script_path} {input.asm_fasta} -p {output.scatter_plot} -d {output.dist_plot} -t {wildcards.sample} -l > {output.n50_stats}
        touch {output.flag}
        '''

rule make_ideo_plot:
    input:
        hap_one_bed = "QC_results/saffire/results/{asm}_hap1/beds/{asm}_hap1.{aligner}.bed",
        hap_two_bed = "QC_results/saffire/results/{asm}_hap2/beds/{asm}_hap2.{aligner}.bed",
    output:
        pdf = "QC_results/ideo_plots/{aligner}/{asm}.ideoplot.pdf",
    threads:
        1
    resources:
        mem = 8,
        hrs = 1,
    params:
        script_path = IDEO_PLOT_SCRIPT,
    shell:
        '''
        {params.script_path} --asm1 {input.hap_one_bed} --asm2 {input.hap_two_bed} -r chm13 -s {wildcards.asm} -o {output.pdf}
        '''

rule make_ploidy_plot:
    input:
        hap_one_paf = "QC_results/saffire/results/{asm}_hap1/alignments/{asm}_hap1.{aligner}.paf",
        hap_two_paf = "QC_results/saffire/results/{asm}_hap2/alignments/{asm}_hap2.{aligner}.paf",

    output:
        pdf = "QC_results/ploidy_plots/{aligner}/{asm}.ploidy.pdf",
        summary = "QC_results/ploidy_plots/{aligner}/{asm}.summary.txt",
    threads:
        1
    resources:
        mem = 8,
        hrs = 4,
    params:
        script_path = PLOIDY_PLOT_SCRIPT,
        snakemake_dir = f"{SNAKEMAKE_DIR}/../scripts",
    singularity:
        "docker://eichlerlab/binf-rplot:0.1"
    shell:
        '''
        Rscript {params.script_path} {params.snakemake_dir} {wildcards.asm} {input.hap_one_paf} {input.hap_two_paf} {output.pdf} {output.summary}
        '''
