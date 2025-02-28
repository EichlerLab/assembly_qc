import pandas as pd
import os
import sys


configfile: 'config/config_asm_qc.yaml'
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

N50_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/n50"
ALIGNER = config.get('ALIGNER',"minimap2")


def find_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

def find_contig_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/contig_fasta/{wildcards.sample}.fasta"


def find_all_bed_set(wildcards):
    haplotype_set = [f"{wildcards.asm}_hap1"]
    if not pd.isna(full_manifest_df.at[wildcards.asm, "H2"]):
        haplotype_set.append(f"{wildcards.asm}_hap2")        
    all_bed_set = [f"saffire/{wildcards.ref}/results/{sample}/beds/{sample}.{wildcards.aligner}.bed" for sample in haplotype_set]
    return all_bed_set


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


rule get_contig_length_plot:
    input:
        asm_fasta = find_fasta,
    output:
        dist_plot = "plots/contigs/{sample}.dist.log.png",
        scatter_plot = "plots/contigs/{sample}.scatter.log.png",
    threads:
        1
    resources:
        mem = 8,
        hrs = 12,
    params:
        script_path = N50_SCRIPT
    shell:
        '''
        {params.script_path} {input.asm_fasta} -p {output.scatter_plot} -d {output.dist_plot} -t {wildcards.sample} -l -n
        '''

rule get_ideo_plot:
    input:
        hap_beds = find_all_bed_set
    output:
        pdf = "plots/ideo/{aligner}/{asm}_to_{ref}.ideoplot.pdf",
    threads:
        1
    resources:
        mem = 8,
        hrs = 1,
    params:
        ref = "{ref}",
        path = lambda wildcards: config["REF"][wildcards.ref]["PATH"],
        cyto = lambda wildcards: config["REF"][wildcards.ref].get("CYTO","Null"),
        chroms = lambda wildcards: config["REF"][wildcards.ref]["CHROMS"],
    singularity:
        "docker://eichlerlab/ideoplot:0.1",
    script:
        f"{SNAKEMAKE_DIR}/../scripts/ideo_plot.py"

rule get_ploidy_plot:
    input:
        hap_one_paf = "saffire/{ref}/results/{asm}_hap1/alignments/{asm}_hap1.{aligner}.paf",
        hap_two_paf = "saffire/{ref}/results/{asm}_hap2/alignments/{asm}_hap2.{aligner}.paf",
    output:
        pdf = "plots/ploidy/{ref}/{aligner}/{asm}.ploidy.pdf",
        summary = "plots/ploidy/{ref}/{aligner}/{asm}.summary.txt",
    threads:
        1
    resources:
        mem = 8,
        hrs = 4,
    params:
        script_dir = f"{SNAKEMAKE_DIR}/../scripts",
    singularity:
        "docker://eichlerlab/binf-rplot:0.1"
    shell:
        '''
        Rscript {params.script_dir}/ploidy.R {params.script_dir} {wildcards.asm} {input.hap_one_paf} {input.hap_two_paf} {output.pdf} {output.summary}
        '''
