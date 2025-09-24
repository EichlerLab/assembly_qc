import pandas as pd
import os
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
import pysam
import numpy as np
import re

configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')

TAXID = config.get("TAXID", "9606")
SOURCE_DIR = "/net/eichler/vol28/7200/software/pipelines/foreign_contamination_screen"
MITO_DB = config.get("mito", f"{SOURCE_DIR}/db/mito_ebv.fa")
RDNA_DB = config.get("rdna", f"{SOURCE_DIR}/db/rdna.fa")

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



def get_fasta(wildcards):
    return manifest_df.at[wildcards.sample, "ASM"]


def find_gx_report(wildcards):
    (IDS,) = glob_wildcards("contamination_screening/results/{sample}/fcs_gx/{gx_name}.fcs_gx_report.txt".format(
            sample=wildcards.sample, gx_name="{gx_name}"
        )
    )

    return expand("contamination_screening/results/{sample}/fcs_gx/{gx_name}.fcs_gx_report.txt",
        sample=wildcards.sample,
        gx_name=IDS,
    )


def report_files(wildcards):
    asmfile = manifest_df.at[wildcards.sample, "ASM"]
    reportfile = Path(asmfile).stem + f".{TAXID}.fcs_gx_report.txt"
    return expand(rules.run_fcs.output, sample=manifest_df.index.values, TAXID=TAXID)


wildcard_constraints:
    sample = "|".join(manifest_df.index),
    sub = "|".join(["gx", "adaptor"]),


localrules:
    run_fcs,
    all,
    all_adaptor,
    all_gx,
    clean_fasta,


rule all:
    input:
        expand("contamination_screening/results/{sample}/trim.bed",
            sample=manifest_df.index,
        ),


rule all_adaptor:
    input:
        expand("contamination_screening/results/{sample}/fcs_adaptor/fcs_adaptor_report.txt",
            sample=manifest_df.index,
        ),


rule all_gx:
    input:
        expand("contamination_screening/results/{sample}/fcs_{sub}/.{sample}.done",
            sample=manifest_df.index,
            sub=["gx"],
        ),


rule clean_fasta:
    input:
        expand("cleaned_fasta/{sample}/scaffold/{sample}.fasta", 
            sample=manifest_df.index
        ),


checkpoint run_fcs:
    input:
        asm_fasta = get_fasta,
    output:
        flag = touch("contamination_screening/results/{sample}/fcs_gx/.{sample}.done"),
    threads: 1
    resources:
        mem = 16,
        hrs = 12,
    benchmark: "contamination_screening/benchmarks/fcs_gx/{sample}.tsv"
    params:
        taxid=TAXID,
        GXDB_LOC="/data/scratch/GXDB/gxdb/",
        fcs_img=f"{SOURCE_DIR}/images/fcs-gx.sif",
        fcs_script=f"{SOURCE_DIR}/fcs.py",
    shell: """
        python3 {params.fcs_script} --image {params.fcs_img} screen genome --fasta {input.asm_fasta} --out-dir $( dirname {output.flag} ) --gx-db {params.GXDB_LOC}  --tax-id {params.taxid}
        """

rule index_asm_fasta:
    input:
        asm_fasta = get_fasta
    output:
        link_fasta = "contamination_screening/raw_fasta/{sample}.fasta",
        fai = "contamination_screening/raw_fasta/{sample}.fasta.fai"
    threads: 1
    resources:
        mem = 8,
        hrs = 4,
    shell: """
        ln -s {input.asm_fasta} {output.link_fasta}
        samtools faidx {output.link_fasta}
        """


rule remove_short_contigs:
    input:
        fasta = "contamination_screening/raw_fasta/{sample}.fasta"
    output:
        filtered_fasta = "contamination_screening/filtered_fasta/{sample}.filtered.fasta",
        filtered_fai = "contamination_screening/filtered_fasta/{sample}.filtered.fasta.fai"
    log:
        "log/contamination_screening/remove_short_contigs_{sample}.log",
    threads: 1
    resources:
        mem = 8,
        hrs = 2,
    run:
        filtered_fasta = output.filtered_fasta
        below_ten_fasta = filtered_fasta.replace(".filtered.fasta",".below_10bp.fasta")

        raw_fasta_records = list(SeqIO.parse(input.fasta, "fasta"))
        filtered_fasta_records = list()
        below_ten_fasta_records = list()

        for raw_fasta_record in raw_fasta_records:
            contig_name = str(raw_fasta_record.id)
            contig_seq = str(raw_fasta_record.seq)
            if len(contig_seq) < 10:
                below_ten_fasta_records.append(SeqRecord(Seq(contig_seq), id=contig_name, description=""))
            else:
                filtered_fasta_records.append(SeqRecord(Seq(contig_seq), id=contig_name, description=""))

        with open(filtered_fasta, "w") as fout_f:
            fasta_writer = FastaWriter(fout_f, wrap=None)
            fasta_writer.write_file(filtered_fasta_records)
        pysam.faidx(filtered_fasta) # indexing
        if len(below_ten_fasta_records) > 0:
            with open(below_ten_fasta, "w") as fout_b:
                fasta_writer = FastaWriter(fout_b, wrap=None)
                fasta_writer.write_file(below_ten_fasta_records)
            pysam.faidx(below_ten_fasta) # indexing


rule run_fcs_adapter:
    input:
        asm_fasta = "contamination_screening/filtered_fasta/{sample}.filtered.fasta",
    output:
        report_txt = "contamination_screening/results/{sample}/fcs_adaptor/fcs_adaptor_report.txt",
    threads: 1
    log:
        "log/contamination_screening/run_fcs_adaptor_{sample}.log",
    singularity:
        f"{SOURCE_DIR}/images/fcs-adaptor.sif"
    resources:
        mem = 16,
        hrs = 12,
    params:
        taxid=TAXID,
        fcs_adaptor_img = f"{SOURCE_DIR}/fcsadaptor/fcs-adaptor.sif",
        fcs_adaptor = f"{SOURCE_DIR}/fcsadaptor/run_fcsadaptor.sh",
    shell: """
        /app/fcs/bin/av_screen_x -o $( dirname {output.report_txt} ) --euk {input.asm_fasta}
        """

rule blast_mito:
    input:
        asm_fasta = "contamination_screening/filtered_fasta/{sample}.filtered.fasta",
        mito_db = MITO_DB,
    output:
        blast_out = "contamination_screening/results/{sample}/fcs_mito/fcs_mito.txt",
    threads: 1
    log:
        "log/contamination_screening/blast_mito_{sample}.log",
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/ncbi-tk:0.1"
    shell: """
        blastn -query {input.asm_fasta} -db {input.mito_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        """


rule filter_mito:
    input:
        blast_out = rules.blast_mito.output.blast_out,
        fai = "contamination_screening/filtered_fasta/{sample}.filtered.fasta.fai",
    output:
        mito_bed = "contamination_screening/results/{sample}/fcs_mito/mito.bed",
    threads: 1
    log:
        "log/contamination_screening/filter_mito_{sample}.log",
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.5) print $1,$2,$3}}' > {output.mito_bed}
        """

rule extract_mito:
    input:
        blast_out = rules.blast_mito.output.blast_out,
        fasta = "contamination_screening/raw_fasta/{sample}.fasta",
    output:
        mito_fasta = "contamination_screening/results/{sample}/fasta/{sample}-mito.fasta",
    params:
        mt_id = "NC_012920.1",
        mt_len = 16569,
        min_pid = 99.0,
        min_len = 16000
    threads: 1
    log:
        "log/contamination_screening/extract_mito_{sample}.log",
    resources:
        mem = 4,
        hrs = 4,
    run:
        blast_header = ["qseqid","sseqid","pid","length","mismatch","gapopen","qstart_blast","qend_blast","sstart","send","evalue","bitscore"]
        blast_df = pd.read_csv(input.blast_out, sep="\t", header=None, names=blast_header)
        full_mt = (
            (blast_df["sseqid"] == params.mt_id) & 
            (
                ((blast_df["sstart"]==1) & (blast_df["send"]==params.mt_len))|
                ((blast_df["sstart"]==params.mt_len) & (blast_df["send"]==1))
            )
        )

        mito_df = blast_df[full_mt & (blast_df["length"] >= params.min_len) & (blast_df["pid"] >= params.min_pid)].copy()

        if mito_df.empty:
            open(output.mito_fasta, "w").close()
            return

        mito_df.loc[:, "strand"] = mito_df.apply(lambda row: "+" if row["qstart_blast"] < row["qend_blast"] else "-", axis=1)
        mito_df["qstart"] = mito_df[["qstart_blast","qend_blast"]].min(axis=1)
        mito_df["qend"] = mito_df[["qstart_blast","qend_blast"]].max(axis=1)
        mito_df["qstart_bed"] = mito_df["qstart"] - 1
        mito_df = mito_df.sort_values(
            by=["qseqid","bitscore","pid","length","mismatch","gapopen"],
            ascending=[True, False, False, False, True, True]
        )
        best = mito_df.drop_duplicates(subset=["qseqid"], keep="first").copy()

        fasta_dict = SeqIO.to_dict(SeqIO.parse(input.fasta, "fasta"))

        out_records = []
        check_dup_seq = dict()
        num_id = 0
        num_rep_id = 0
        for idx, row in best.iterrows():
            contig = row["qseqid"]

            start = int(row["qstart_bed"])
            end = int(row["qend"])
            subrec = fasta_dict[contig][start:end]
            subrec_seq = str(subrec.seq).upper()
            if row["strand"] == "-":
                subrec = subrec.reverse_complement(id=True, name=True, description=True)

            new_id = f"{contig}:{int(row['qstart'])}-{int(row['qend'])}"
            num_id += 1
            subrec.id = new_id
            subrec.name = new_id
            subrec.description = ""

            if not subrec_seq in check_dup_seq:
                num_rep_id += 1
                check_dup_seq[subrec_seq] = []
                check_dup_seq[subrec_seq].append(new_id)
                out_records.append(subrec)
            else:
                check_dup_seq[subrec_seq].append(new_id)

        with open(output.mito_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(out_records)

rule trim_bed:
    input:
        gx_report = find_gx_report,
        flag = rules.run_fcs.output.flag,
        adapt_report = rules.run_fcs_adapter.output.report_txt,
        mito_bed = rules.filter_mito.output.mito_bed,
    output:
        trim_file = "contamination_screening/results/{sample}/trim.bed",
    threads: 1
    log:
        "log/contamination_screening/run_fcs_trim_{sample}.log",
    resources:
        mem = 16,
        hrs = 12,
    run:
        out_df = pd.DataFrame()
        mito_df = pd.read_csv(
            input.mito_bed,
            sep="\t",
            header=None,
            names=["#seq_id", "start_pos", "end_pos"],
        )
        mito_df["reason"] = "mito_ebv_rdna"
        print (mito_df)
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
        trim_file = rules.trim_bed.output.trim_file,
        fai = "contamination_screening/filtered_fasta/{sample}.filtered.fasta.fai",
    output:
        regions_file = "contamination_screening/temp/{sample}/regions.out",
    threads: 1
    log:
        "log/contamination_screening/coerce_bed_{sample}.log",
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        bedtools complement -i <( sort -k1,1 -k2,2n {input.trim_file} ) -g <( sort -k1 {input.fai} ) | awk '{{if ( $3 - $2 > 1000 ) print $1":"$2+1"-"$3}}' | sort > {output.regions_file}
        """


rule trim_sequence:
    input:
        regions_file=rules.coerce_bed.output.regions_file,
        asm="contamination_screening/filtered_fasta/{sample}.filtered.fasta",
    output:
        cleaned_fasta="contamination_screening/temp/{sample}/fasta/{sample}.gx_adapt_cleaned.fasta", ## contamination & adaptor cleaned
        cleaned_index="contamination_screening/temp/{sample}/fasta/{sample}.gx_adapt_cleaned.fasta.fai",
    threads: 1
    log:
        "log/contamination_screening/trim_sequence_{sample}.log",
    singularity: "docker://eichlerlab/binf-basics:0.1"
    resources:
        mem=4,
        hrs=12,
    ### removed sed 's/:/#/g'
    shell: """
        samtools faidx -r {input.regions_file} {input.asm} > {output.cleaned_fasta}
        samtools faidx {output.cleaned_fasta} 
        """


rule blast_rdna:
    input:
        asm_fasta=rules.trim_sequence.output.cleaned_fasta,
        rdna_db=RDNA_DB,
    output:
        blast_out="contamination_screening/results/{sample}/rdna/rdna.txt",
    threads: 1
    log:
        "log/contamination_screening/blast_rdna_{sample}.log",
    resources:
        mem=4,
        hrs=12,
    singularity: "docker://eichlerlab/ncbi-tk:0.1"
    shell: """
        blastn -query {input.asm_fasta} -db {input.rdna_db} -outfmt 6 -evalue 1e-30 > {output.blast_out}
        """


rule filter_rdna:
    input:
        blast_out = rules.blast_rdna.output.blast_out,
        fai = rules.trim_sequence.output.cleaned_index,
    output:
        rdna_ctg = "contamination_screening/results/{sample}/rdna/rdna.bed",
        other_ctg = "contamination_screening/results/{sample}/clean/clean.bed",
    threads: 1
    log:
        "log/contamination_screening/filter_mito_{sample}.log",
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        bedtools coverage -a <( awk -vOFS="\\t" '{{print $1,"0",$2}}' {input.fai} ) -b <( cut -f 1,7,8 {input.blast_out} ) | awk -vOFS="\\t" '{{ if ($NF >= 0.95) print $1}}' > {output.rdna_ctg}
        if [[ $( cat {output.rdna_ctg} | wc -l ) == 0 ]]; then
            cp {input.fai} {output.other_ctg}
        else
            cat {output.rdna_ctg} | tr "\\n" "|" | sed 's/|$//' | xargs -i grep -vE {{}} {input.fai} |  cut -f 1 > {output.other_ctg}
        fi 
        """


rule split_rdna:
    input:
        fasta = rules.trim_sequence.output.cleaned_fasta,
        fai = rules.trim_sequence.output.cleaned_index,
        rdna = rules.filter_rdna.output.rdna_ctg,
        others = rules.filter_rdna.output.other_ctg,
    output:
        cleaned_fasta = "contamination_screening/temp/{sample}/fasta/{sample}.gx_adapt_rdna_cleaned.fasta", # r-DNA-cleaned
        rdna_fasta = "contamination_screening/results/{sample}/fasta/{sample}-rdna.fasta",
        cleaned_fai = "contamination_screening/temp/{sample}/fasta/{sample}.gx_adapt_rdna_cleaned.fasta.fai",
    threads: 1
    log:
        "log/contamination_screening/filter_rdna_{sample}.log",
    resources:
        mem = 4,
        hrs = 12,
    singularity: "docker://eichlerlab/binf-basics:0.1"
    shell: """
        if [[ $( wc -l {input.rdna} | awk '{{print $1}}' ) == 0 ]]; then
            cp {input.fasta} {output.cleaned_fasta}
            touch {output.rdna_fasta}
        else
            samtools faidx -r {input.rdna} {input.fasta} > {output.rdna_fasta}
            samtools faidx -r {input.others} {input.fasta} > {output.cleaned_fasta}
        fi 
        samtools faidx {output.cleaned_fasta}
        """


rule rename_fasta:
    input:
        asm_fasta = "contamination_screening/filtered_fasta/{sample}.filtered.fasta",
        cleaned_fasta = rules.split_rdna.output.cleaned_fasta,
        mito_fasta = rules.extract_mito.output.mito_fasta,
    output:
        renamed_final_fasta = "fcs_cleaned_fasta/{sample}/{sample}.fasta",
        renamed_final_fai = "fcs_cleaned_fasta/{sample}/{sample}.fasta.fai"
    threads: 1
    log:
        "log/contamination_screening/rename_fasta_{sample}.log",
    resources:
        mem = 8,
        hrs = 12,
    run:
        renamed_final_fasta = output.renamed_final_fasta

        original_fasta = pysam.FastaFile(input.asm_fasta)
        cleaned_fasta_records = list(SeqIO.parse(input.cleaned_fasta, "fasta"))
        final_fasta_records = list()

        for record in cleaned_fasta_records:
            seq_name = str(record.id)
            #original_seq_name = seq_name.split(":")[0]
            original_seq_name = re.sub(r":[^:]+?$", "", seq_name)
            cleaned_sequence = str(record.seq)
            raw_sequence = original_fasta.fetch(original_seq_name)
            if cleaned_sequence == raw_sequence:
                cleaned_seq_name = original_seq_name
            else:
                cleaned_seq_name = seq_name.replace(":","_trim_")
            final_fasta_records.append(SeqRecord(Seq(cleaned_sequence), id=cleaned_seq_name, description=""))
        with open(renamed_final_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(final_fasta_records)

        pysam.faidx(renamed_final_fasta) # indexing
