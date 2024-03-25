import pandas as pd
import os
from pathlib import Path

configfile: "config/config.yaml"
MANIFEST = config.get("MANIFEST", "config/manifest.tab")


MODE = config.get("MODE", "lite")
LINEAGE = config.get("LINEAGE", "primates")
MB_DOWNLOADS = config.get("MB_DOWNLOADS", "/net/eichler/vol28/software/modules-sw/compleasm/0.2.2/Linux/CentOS7/x86_64/compleasm_kit/mb_download/")

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

## Universial conversion of manifest df
manifest_df = pd.read_csv(MANIFEST, sep='\t')

if ("H1" in manifest_df.columns) and ("H2" in manifest_df):
    df_transform = list()
    for idx, row in manifest_df.iterrows():
        df_transform.append({"SAMPLE": "%s_hap1"%row["SAMPLE"], "ASM":row["H1"]})
        df_transform.append({"SAMPLE": "%s_hap2"%row["SAMPLE"], "ASM":row["H2"]})

    manifest_df = pd.DataFrame(df_transform)

manifest_df.set_index("SAMPLE",inplace=True)
##-------------------------------------


def get_fasta(wildcards):
    return f"QC_results/contamination_screening/results/{wildcards.sample}/fasta/{wildcards.sample}.fasta"

wildcard_constraints:
    sample="|".join(manifest_df.index),

localrules: compleasm_run, all,

rule all:
    input:
        expand("QC_results/compleasm/results/{sample}/summary.txt", sample=manifest_df.index.values),

rule compleasm_run:
    input:
        asm_fasta=get_fasta,
    output:
        summary = "QC_results/compleasm/results/{sample}/summary.txt",
    threads: 16
    resources:
        mem=16,
        hrs=2,
    params:
        mb_downloads=MB_DOWNLOADS,
        mode=MODE,
        lineage=LINEAGE,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.9.2",
        f"compleasm/0.2.2"
    shell:
        """
        compleasm.py run -a {input.asm_fasta} -o QC_results/compleasm/results/{wildcards.sample} -t {threads} -l {params.lineage} -m {params.mode} -L {params.mb_downloads}
        """

