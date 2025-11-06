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
TAXID = config.get("TAXID", "9606")
REF_DICT = config["REF"]
ALIGNER = config.get('ALIGNER',"minimap2")

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

full_manifest_df = pd.read_csv(MANIFEST, header=0, sep='\t', comment='#', na_values=["","NA","na","N/A"])
conv_manifest_df = get_asm_manifest_df(full_manifest_df)

full_manifest_df.set_index("SAMPLE",inplace=True) ## manifest df for merqury
conv_manifest_df.set_index("SAMPLE",inplace=True) ## shortened df using only SAMPLE / ASM fields

def get_all_inputs():
    inputs = [
        expand("stats/seq_stats/{sample}.scaffold.stats",
            sample=conv_manifest_df.index.values,
        ),
        expand("stats/seq_stats/{sample}.contig.stats",
            sample=conv_manifest_df.index.values,
        ),
        expand("merqury/results/{asm}/{asm}.qv",
            asm=full_manifest_df.index.values,
        ),
        expand("saffire/{ref}/results/{sample}/alignments/{sample}.{aligner}.paf",
            ref=REF_DICT,
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),
    	expand("saffire/{ref}/results/{sample}/alignments/{sample}.{aligner}.bam",
            ref=REF_DICT,
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),  
        expand("saffire/{ref}/results/{sample}/beds/{sample}.{aligner}.bed",
            ref=REF_DICT,
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),
        expand("compleasm/results/{sample}/summary.txt",
            sample=conv_manifest_df.index.values,
        ),
        expand("plots/contigs/{sample}.contig.scatter.log.png",
            sample=conv_manifest_df.index.values,
        ),
        expand("plots/ideo/{aligner}/{asm}_to_{ref}.ideoplot.pdf",
            ref=REF_DICT,
            asm=full_manifest_df.index.values,
            aligner=ALIGNER,            
        ),
        expand("saffire/{ref}/results/merged_paf/{aligner}/{asm}.concat.paf",
            ref=REF_DICT,
            asm=full_manifest_df.index.values,
            aligner=ALIGNER,
        ),
        expand("stats/acro_stats/{sample}.generated_acros.tsv",
            sample=conv_manifest_df.index.values,
        ),
        # expand("moddotplot/liftover/stats/{sample}/lifted_contigs.tsv",
        #     sample=conv_manifest_df.index.values,
        # )
    ]

    ploidy_inputs = expand("plots/ploidy/CHM13/{aligner}/{asm}.ploidy.pdf",
        asm=[
            asm for asm in full_manifest_df.index.values
            if not pd.isna(full_manifest_df.at[asm, "H2"])
        ],
        aligner=ALIGNER,
    )
    moddot_inputs = expand("stats/acro_stats/{sample}.generated_acros.tsv",
        sample=conv_manifest_df.index.values,
    )

    if int(TAXID) == 9606: # human
        inputs.extend(ploidy_inputs)
        inputs.extend(moddot_inputs)

    return inputs

def get_plot_inputs():
    inputs = [

        expand("plots/ideo/{aligner}/{asm}_to_{ref}.ideoplot.pdf",
            ref=REF_DICT,
            asm=full_manifest_df.index.values,
            aligner=ALIGNER,
        ),
        expand("plots/contigs/{sample}.contig.scatter.log.png",
            sample=conv_manifest_df.index.values,
        ),
    ]

    ploidy_inputs = expand("plots/ploidy/CHM13/{aligner}/{asm}.ploidy.pdf",
        asm=[
            asm for asm in full_manifest_df.index.values
            if not pd.isna(full_manifest_df.at[asm, "H2"])
        ],
        aligner=ALIGNER,
    )

    moddot_inputs = [
        expand("stats/acro_stats/{sample}.generated_acros.tsv",
            sample=conv_manifest_df.index.values,
        ),
        # expand("moddotplot/liftover/stats/{sample}/lifted_contigs.tsv",
        #     sample=conv_manifest_df.index.values,
        # )
    ]

    if int(TAXID) == 9606: # human
        inputs.extend(ploidy_inputs)    
        inputs.extend(moddot_inputs)

    return inputs

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

module stats:
    snakefile:
        "rules/contig_stats.smk"
    config:
        config
use rule * from stats as cst_*

if int(TAXID) == 9606:  # only for human
    module moddot:
        snakefile:
            "rules/moddotplot.smk"
        config:
            config
    use rule * from moddot as mdp_*

localrules:
    mer_agg_hapmers,
    fcs_all_gx,
    fcs_run_fcs,
    fcs_all_adaptor,
    fcs_clean_fasta,
    cpl_compleasm_run,

rule all:
    input:
        get_all_inputs()

    default_target: True

rule get_saf:
    input:
        expand("saffire/{ref}/results/{sample}/saff/{sample}.{aligner}.saf",
            ref=REF_DICT,
            sample=conv_manifest_df.index.values,
            aligner=ALIGNER,
        ),

rule get_compleasm:
    input:
        expand("compleasm/results/{sample}/summary.txt",
            sample=conv_manifest_df.index.values,
        ),

rule get_cleaned_fasta:
    input:
        expand("fcs_cleaned_fasta/{sample}/{sample}.fasta",
            sample=conv_manifest_df.index.values,
        ),
        expand("contamination_screening/results/{sample}/fasta/{sample}-mito.fasta",
            sample=conv_manifest_df.index.values,
        )

rule get_qv:
    input:
        expand("merqury/results/{asm}/{asm}.qv",
            asm=full_manifest_df.index.values,
        ),

rule get_plots:
    input:
        get_plot_inputs()

rule get_stats:
    input:
        expand("stats/seq_stats/{sample}.scaffold.stats",
            sample=conv_manifest_df.index.values,
        ),
        expand("stats/seq_stats/{sample}.contig.stats",
            sample=conv_manifest_df.index.values,
        ),