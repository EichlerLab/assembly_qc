import pandas as pd
import os
from pathlib import Path


configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

TAXID = config.get("TAXID", "9606")
SOURCE_DIR = "/net/eichler/vol28/7200/software/pipelines/foreign_contamination_screen"
MITO_DB = config.get("mito", f"{SOURCE_DIR}/db/mito_ebv.fa")
RDNA_DB = config.get("rdna", f"{SOURCE_DIR}/db/rdna.fa")

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
    return manifest_df.at[wildcards.sample, "ASM"]


def get_fai(wildcards):
    return manifest_df.at[wildcards.sample, "ASM"] + ".fai"


def find_gx_report(wildcards):
    (IDS,) = glob_wildcards(
        "QC_results/contamination_screening/results/{sample}/fcs_gx/{gx_name}.fcs_gx_report.txt".format(
            sample=wildcards.sample, gx_name="{gx_name}"
        )
    )

    return expand(
        "QC_results/contamination_screening/results/{sample}/fcs_gx/{gx_name}.fcs_gx_report.txt",
        sample=wildcards.sample,
        gx_name=IDS,
    )


def report_files(wildcards):
    asmfile = manifest_df.at[wildcards.sample, "ASM"]
    reportfile = Path(asmfile).stem + f".{TAXID}.fcs_gx_report.txt"
    return expand(rules.run_fcs.output, sample=manifest_df.index.values, TAXID=TAXID)


wildcard_constraints:
    sample="|".join(manifest_df.index),
    sub="|".join(["gx", "adaptor"]),


localrules:
    run_fcs,
    all,
    all_adaptor,
    all_gx,
    clean_fasta,


rule all:
    input:
        expand(
            "QC_results/contamination_screening/results/{sample}/trim.bed",
            sample=manifest_df.index,
        ),


rule all_adaptor:
    input:
        expand(
            "QC_results/contamination_screening/results/{sample}/fcs_adaptor/fcs_adaptor_report.txt",
            sample=manifest_df.index,
        ),


rule all_gx:
    input:
        expand(
            "QC_results/contamination_screening/results//{sample}/fcs_{sub}/.{sample}.done",
            sample=manifest_df.index,
            sub=["gx"],
        ),


rule clean_fasta:
    input:
        expand("QC_results/contamination_screening/results/{sample}/fasta/{sample}.fasta", sample=manifest_df.index),


checkpoint run_fcs:
    input:
        asm_fasta=get_fasta,
    output:
        flag=touch("QC_results/contamination_screening/results/{sample}/fcs_gx/.{sample}.done"),
    threads: 1
    resources:
        mem=16,
        hrs=12,
    benchmark: "QC_results/contamination_screening/benchmarks/fcs_gx/{sample}.tsv"
    params:
        taxid=TAXID,
        GXDB_LOC="/tmp/GXDB/gxdb/",
        fcs_img=f"{SOURCE_DIR}/images/fcs-gx.sif",
        fcs_script=f"{SOURCE_DIR}/fcs.py",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.9.2",
    shell:
        """
        python3 {params.fcs_script} --image {params.fcs_img} screen genome --fasta {input.asm_fasta} --out-dir $( dirname {output.flag} ) --gx-db {params.GXDB_LOC}  --tax-id {params.taxid}
        """


rule run_fcs_adapter:
    input:
        asm_fasta=get_fasta,
    output:
        report_txt="QC_results/contamination_screening/results/{sample}/fcs_adaptor/fcs_adaptor_report.txt",
    threads: 1
    log:
        "log/run_fcs_adaptor_{sample}.log",
    singularity:
        f"{SOURCE_DIR}/images/fcs-adaptor.sif"
    resources:
        mem=16,
        hrs=12,
    params:
        taxid=TAXID,
        fcs_adaptor_img=f"{SOURCE_DIR}/fcsadaptor/fcs-adaptor.sif",
        fcs_adaptor=f"{SOURCE_DIR}/fcsadaptor/run_fcsadaptor.sh",
    shell:
        """
        /app/fcs/bin/av_screen_x -o $( dirname {output.report_txt} ) --euk {input.asm_fasta}
        """


rule blast_mito:
    input:
        asm_fasta=get_fasta,
        mito_db=MITO_DB,
    output:
        blast_out="QC_results/contamination_screening/results/{sample}/fcs_mito/fcs_mito.txt",
    threads: 1
    log:
        "log/blast_mito_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "ncbi-blast/2.15.0",
    shell:
        """
        blastn -query {input.asm_fasta} -db {input.mito_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        """


rule filter_mito:
    input:
        blast_out=rules.blast_mito.output.blast_out,
        fai=get_fai,
    output:
        mito_bed="QC_results/contamination_screening/results/{sample}/fcs_mito/mito.bed",
    threads: 1
    log:
        "log/filter_mito_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "bedtools/2.29.0",
    shell:
        """
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.5) print $1,$2,$3}}' > {output.mito_bed}
        """


rule trim_bed:
    input:
        gx_report=find_gx_report,
        flag=rules.run_fcs.output.flag,
        adapt_report=rules.run_fcs_adapter.output.report_txt,
        mito_bed=rules.filter_mito.output.mito_bed,
    output:
        trim_file="QC_results/contamination_screening/results/{sample}/trim.bed",
    threads: 1
    log:
        "log/run_fcs_trim_{sample}.log",
    resources:
        mem=16,
        hrs=12,
    run:
        out_df = pd.DataFrame()
        mito_df = pd.read_csv(
            input.mito_bed,
            sep="\t",
            header=None,
            names=["#seq_id", "start_pos", "end_pos"],
        )
        mito_df["reason"] = "mito_ebv_rdna"
        out_df = pd.concat([out_df, mito_df])
        df_gx = pd.concat(
            [
                pd.read_csv(
                    x,
                    sep="\t",
                    header=None,
                    comment="#",
                    names=[
                        "#seq_id",
                        "start_pos",
                        "end_pos",
                        "seq_len",
                        "action",
                        "div",
                        "agg_cont_cov",
                        "top_tax_name",
                    ],
                )
                for x in input.gx_report
            ]
        )
        df_gx["start_pos"] = df_gx["start_pos"] - 1
        df_gx["reason"] = "foreign_contam"
        out_df = pd.concat(
            [out_df, df_gx[["#seq_id", "start_pos", "end_pos", "reason"]]]
        )
        df_adapt = pd.read_csv(input.adapt_report, sep="\t")
        df_adapt_trim = df_adapt.loc[df_adapt["action"] == "ACTION_TRIM"].copy()
        df_adapt_exc = df_adapt.loc[df_adapt["action"] == "ACTION_EXCLUDE"].copy()
        if len(df_adapt_trim) > 0:
            df_adapt_trim["#seq_id"] = df_adapt_trim["#accession"]
            df_adapt_trim["range"] = df_adapt_trim["range"].str.split(",")
            df_adapt_trim = df_adapt_trim.explode("range")
            df_adapt_trim["start_pos"] = (
                df_adapt_trim["range"]
                .str.split("..", expand=True, regex=False)[0]
                .astype(int)
                - 1
            )
            df_adapt_trim["end_pos"] = (
                df_adapt_trim["range"]
                .str.split("..", expand=True, regex=False)[1]
                .astype(int)
            )
            df_adapt_trim["reason"] = "adapter"
            out_df = pd.concat(
                [out_df, df_adapt_trim[["#seq_id", "start_pos", "end_pos", "reason"]]]
            )
        if len(df_adapt_exc) > 0:
            df_adapt_exc["#seq_id"] = df_adapt_exc["#accession"]
            df_adapt_exc["start_pos"] = 0
            df_adapt_exc["end_pos"] = df_adapt_exc["length"]
            df_adapt_exc["reason"] = "adapter"
            out_df = pd.concat(
                [out_df, df_adapt_exc[["#seq_id", "start_pos", "end_pos", "reason"]]]
            )
        out_df[["#seq_id", "start_pos", "end_pos", "reason"]].to_csv(
            output.trim_file, sep="\t", header=False, index=False
        )


rule coerce_bed:
    input:
        trim_file=rules.trim_bed.output.trim_file,
        fai=get_fai,
    output:
        regions_file="QC_results/contamination_screening/temp/{sample}/regions.out",
    threads: 1
    log:
        "log/coerce_bed_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "bedtools/2.29.0",
    shell:
        """
        bedtools complement -i <( sort -k1,1 -k2,2n {input.trim_file} ) -g <( sort -k1 {input.fai} ) | awk '{{if ( $3 - $2 > 1000 ) print $1":"$2+1"-"$3}}' | sort > {output.regions_file}
        """


rule trim_sequence:
    input:
        regions_file=rules.coerce_bed.output.regions_file,
        asm=get_fasta,
    output:
        cleaned_fasta="QC_results/contamination_screening/temp/{sample}/fasta/{sample}.fasta",
        cleaned_index="QC_results/contamination_screening/temp/{sample}/fasta/{sample}.fasta.fai",
    threads: 1
    log:
        "log/trim_sequence_{sample}.log",
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.12",
    resources:
        mem=4,
        hrs=12,
    shell:
        """
        if [[ $( wc -l {input.regions_file} | awk '{{print $1}}' ) == 0 ]]; then
            rsync -av {input.asm} {output.cleaned_fasta}
        else
            samtools faidx -r {input.regions_file} {input.asm} | sed 's/:/#/g' > {output.cleaned_fasta}
        fi 
        samtools faidx {output.cleaned_fasta} 
        """


rule blast_rdna:
    input:
        asm_fasta=rules.trim_sequence.output.cleaned_fasta,
        rdna_db=RDNA_DB,
    output:
        blast_out="QC_results/contamination_screening/results/{sample}/rdna/rdna.txt",
    threads: 1
    log:
        "log/blast_rdna_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "ncbi-blast/2.15.0",
    shell:
        """
        blastn -query {input.asm_fasta} -db {input.rdna_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        """


rule filter_rdna:
    input:
        blast_out=rules.blast_rdna.output.blast_out,
        fai=rules.trim_sequence.output.cleaned_index,
    output:
        rdna_ctg="QC_results/contamination_screening/results/{sample}/rdna/rdna.bed",
        other_ctg="QC_results/contamination_screening/results/{sample}/clean/clean.bed",
    threads: 1
    log:
        "log/filter_mito_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "bedtools/2.29.0",
    shell:
        """
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.95) print $1}}' > {output.rdna_ctg}
        if [[ $( cat {output.rdna_ctg} | wc -l ) == 0 ]]; then
            cp {input.fai} {output.other_ctg}
        else
            cat {output.rdna_ctg} | tr "\\n" "|" | sed 's/|$//' | xargs -i grep -vE {{}} {input.fai} |  cut -f 1 > {output.other_ctg}
        fi 
        """


rule split_rdna:
    input:
        fasta=rules.trim_sequence.output.cleaned_fasta,
        fai=rules.trim_sequence.output.cleaned_index,
        rdna=rules.filter_rdna.output.rdna_ctg,
        others=rules.filter_rdna.output.other_ctg,
    output:
        cleaned_fasta="QC_results/contamination_screening/results/{sample}/fasta/{sample}.fasta",
        cleaned_fai="QC_results/contamination_screening/results/{sample}/fasta/{sample}.fasta.fai",
        rdna_fasta="QC_results/contamination_screening/results/{sample}/fasta/{sample}-rdna.fasta",
        rdna_fai="QC_results/contamination_screening/results/{sample}/fasta/{sample}-rdna.fasta.fai",
    threads: 1
    log:
        "log/filter_rdna_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "samtools/1.19",
    shell:
        """
        if [[ $( wc -l {input.rdna} | awk '{{print $1}}' ) == 0 ]]; then
            cp -l {input.fasta} {output.cleaned_fasta}
            touch {output.rdna_fasta}
        else
            samtools faidx -r {input.rdna} {input.fasta} > {output.rdna_fasta}
            samtools faidx -r {input.others} {input.fasta} > {output.cleaned_fasta}
        fi 
        """
