import pandas as pd
import os
import sys

N50_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/n50"


def find_all_bed_set(wildcards):
    try:
        h2_val = full_manifest_df.at[wildcards.sample, "H2"]
    except KeyError:
        h2_val = None 

    has_h2 = isinstance(h2_val, str) and h2_val.strip() != "" and os.path.isfile(h2_val)
    haps = ["hap1", "hap2"] if has_h2 else ["hap1"]
    outputs = [
        f"results/{wildcards.sample}/saffire/work/alignments/{wildcards.ref}/beds/{hap}.minimap2.bed"
        for hap in haps
    ]
    return outputs[0] if len(outputs) == 1 else outputs


def find_all_saf_paf(wildcards):
    avail_pafs = []
    haps = ["hap1","hap2","un"]
    values = full_manifest_df.loc[wildcards.sample][["H1","H2","UNASSIGNED"]].tolist()
    for idx, hap in enumerate(haps):
        if not str(values[idx]) == "nan":
            avail_pafs.append (f"results/{wildcards.sample}/saffire/outputs/trimmed_pafs/{wildcards.ref}/{hap}.minimap2.trimmed.paf")
    return avail_pafs[0] if len(avail_pafs) == 1 else avail_pafs

rule get_scaffold_length_plot:
    input:
        hap_fasta = rules.rename_fasta.output.final_fasta,
    output:
        scatter_plot = "results/{sample}/plots/outputs/contig_length/{hap}.scaffold.scatter_logged.png",
    threads:
        1
    resources:
        mem = 8,
        hrs = 12,
    params:
        script_path = N50_SCRIPT
    shell: """
        {params.script_path} {input.hap_fasta} -p {output.scatter_plot} -t "{wildcards.sample}-{wildcards.hap}" -l -n
        """

rule get_contig_length_plot:
    input:
        hap_fasta = rules.split_scaffolds.output.contig_fasta,
    output:
        scatter_plot = "results/{sample}/plots/outputs/contig_length/{hap}.contig.scatter_logged.png",
    threads:
        1
    resources:
        mem = 8,
        hrs = 12,
    params:
        script_path = N50_SCRIPT
    shell: """
        {params.script_path} {input.hap_fasta} -p {output.scatter_plot} -t "{wildcards.sample}-{wildcards.hap}" -l -n
        """

rule get_ideo_plot:
    input:
        hap_beds = find_all_bed_set
    output:
        pdf = "results/{sample}/plots/outputs/ideo/{ref}/pdf/{sample}.minimap2.ideoplot.pdf",
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
        wide = False
    singularity:
        "docker://eichlerlab/ideoplot:0.1",
    script:
        f"{SNAKEMAKE_DIR}/scripts/ideo_wide_plot.py"


rule get_wide_ideo_plot:
    input:
        hap_beds = find_all_bed_set
    output:
        pdf = "results/{sample}/plots/outputs/ideo/{ref}/pdf/{sample}.minimap2.ideoplot_wide.pdf",
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
        wide = True
    singularity:
        "docker://eichlerlab/ideoplot:0.1",
    script:
        f"{SNAKEMAKE_DIR}/scripts/ideo_wide_plot.py"

rule get_ploidy_plot:
    input:
        paf_set = find_all_saf_paf
    output:
        pdf = "results/{sample}/plots/outputs/ploidy/{ref}/pdf/{sample}.minimap2.ploidy.pdf",
        summary = "results/{sample}/plots/outputs/ploidy/{ref}/summary/{sample}.minimap2.ploidy_summary.txt",
    threads:
        1
    resources:
        mem = 8,
        hrs = 4,
    params:
        script_dir = f"{SNAKEMAKE_DIR}/scripts",
    singularity:
        "docker://eichlerlab/binf-rplot:0.1"
    shell:
        '''
        Rscript {params.script_dir}/ploidy.R {params.script_dir} {wildcards.sample} {input.paf_set} {output.pdf} {output.summary}
        '''
