import pandas as pd
import os
import sys

# define your configuration via python
# or define a yaml e.g.

configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
IDEO_PLOT_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/ideo_plot.py"
N50_SCRIPT = "/net/eichler/vol28/7200/software/pipelines/compteam_tools/n50"
PLOIDY_PLOT_SCRIPT = f"{SNAKEMAKE_DIR}/scripts/ploidy.R"
ALIGNER = config.get('ALIGNER',"minimap2")


def get_asm_manifest_df(manifest_df):
    df_transform = list()
    for idx, row in manifest_df.iterrows():
        df_transform.append({"SAMPLE": "%s_hap1"%row["SAMPLE"], "ASM":row["H1"]})
        df_transform.append({"SAMPLE": "%s_hap2"%row["SAMPLE"], "ASM":row["H2"]})

    return pd.DataFrame(df_transform)


full_manifest_df = pd.read_csv(MANIFEST, header=0, sep='\t', comment='#').fillna('NA')
conv_manifest_df = get_asm_manifest_df(full_manifest_df)

full_manifest_df.set_index("SAMPLE",inplace=True) ## manifest df for merqury
conv_manifest_df.set_index("SAMPLE",inplace=True) ## shortened df using only SAMPLE / ASM fields

module merqury:
    snakefile:
        "rules/merqury.smk"
    config:
        config
use rule * from merqury as mer_*

module saffire:
    snakefile:
        "rules/saffire.smk"
    config:
        config
use rule * from saffire as saf_*

module fcs:
    snakefile:
        "rules/fcs_gx.smk"
    config:
        config
use rule * from fcs as fcs_*

module compleasm:
    snakefile:
        "rules/compleasm.smk"
    config:
        config
use rule * from compleasm as cpl_*

module plots:
    snakefile:
        "rules/plots.smk"
    config:
        config
use rule * from plots as plt_*

localrules:
    mer_agg_hapmers,
    fcs_all_gx,
    fcs_run_fcs,
    fcs_all_adaptor,
    fcs_clean_fasta,
    cpl_compleasm_run,


rule all:
    # fcs as default
    input:
        # Saffire
        expand(
            "QC_results/saffire/results/{sample}/alignments/{sample}.{aligner}.bam",
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),
        expand(
            rules.saf_make_bed.output.bed,
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),
        # Merqury
        expand(            
            "QC_results/merqury/results/{asm}/{asm}.qv",
            asm=full_manifest_df.index.values,
        ),
        # Compleasm
        expand(
            rules.cpl_compleasm_run.output.summary,
            sample=conv_manifest_df.index.values,
        ),
        # Contig stats
        expand(
            rules.plt_make_contig_stats.output.flag,
            sample=conv_manifest_df.index,
        ),
        # ideo plot for assembly
        expand(
            rules.plt_make_ideo_plot.output.pdf,
            asm=full_manifest_df.index,
            aligner = ALIGNER,
        ),
        # ploidy plot for assembly
        expand(
            rules.plt_make_ploidy_plot.output.summary,
            asm=full_manifest_df.index,
            aligner = ALIGNER,
        ),

    default_target: True

