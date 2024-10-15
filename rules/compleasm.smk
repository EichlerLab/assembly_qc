import pandas as pd
import os
from pathlib import Path

configfile: "config/config.yaml"
MANIFEST = config.get("MANIFEST", "config/manifest.tab")


MODE = config.get("MODE", "lite")
LINEAGE = config.get("LINEAGE", "primates")
MB_DOWNLOADS = config.get("MB_DOWNLOADS", "/net/eichler/vol28/software/modules-sw/compleasm/0.2.6/Linux/Ubuntu22.04/x86_64/mb_downloads/")

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

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



def get_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

wildcard_constraints:
    sample="|".join(manifest_df.index),

localrules: compleasm_run, all,

rule all:
    input:
        expand("compleasm/results/{sample}/summary.txt", sample=manifest_df.index.values),

rule compleasm_run:
    input:
        asm_fasta=get_fasta,
    output:
        summary = "compleasm/results/{sample}/summary.txt",
    threads: 16
    resources:
        mem=16,
        hrs=2,
    params:
        mb_downloads=MB_DOWNLOADS,
        mode=MODE,
        lineage=LINEAGE,
    singularity:
        "docker://eichlerlab/compleasm:0.2.6"
    shell:
        """
        compleasm.py run -a {input.asm_fasta} -o compleasm/results/{wildcards.sample} -t {threads} -l {params.lineage} -m {params.mode} -L {params.mb_downloads}
        """

