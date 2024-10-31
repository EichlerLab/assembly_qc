import pandas as pd
import os
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import pysam

configfile: "config/config_asm_qc.yaml"
MANIFEST = config.get('MANIFEST', 'config/manifest_asm_qc.tab')
SNAKEMAKE_ROOT_DIR = os.path.dirname(workflow.snakefile).replace("/rules","")


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

raw_manifest_df.set_index("SAMPLE", inplace=True)

#-----------------------------------------

def find_scaftig_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/{wildcards.sample}.fasta"

def find_contig_fasta(wildcards):
    return f"fcs_cleaned_fasta/{wildcards.sample}/contig_fasta/{wildcards.sample}.fasta"

def find_all_fasta_set(wildcards):
    fasta_set = [f"{wildcards.asm}_hap1"]
    if not pd.isna(raw_manifest_df.at[wildcards.asm, "H2"]):
        fasta_set.append(f"{wildcards.asm}_hap2")
    if not pd.isna(raw_manifest_df.at[wildcards.asm, "UNASSIGNED"]):
        fasta_set.append(f"{wildcards.asm}_unassigned")
    fasta_set = [f"fcs_cleaned_fasta/{sample}/{sample}.fasta" for sample in fasta_set] + [f"fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta" for sample in fasta_set]
    return fasta_set

def find_all_telo_tsvs(wildcards):
    telo_haps = [f"{wildcards.asm}_hap1"]
    if not pd.isna(raw_manifest_df.at[wildcards.asm, "H2"]):
        telo_haps.append(f"{wildcards.asm}_hap2")
    telo_tsvs = [f"stats/telo/{sample}.scaftig.telo.tsv" for sample in telo_haps] + [f"stats/telo/{sample}.contig.telo.tsv" for sample in telo_haps]
    return telo_tsvs

def find_all_qv(wildcards):
    all_qv = []
    return all_qv

rule split_scaftigs:
    input:
        scaftig_fasta=find_scaftig_fasta,
    output:
        contig_fasta = "fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta"
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 8,
        hrs=12,
    run:

        contig_fasta = output.contig_fasta
        contig_fasta_records = list()
        gap_pattern = re.compile(r"N+")
        fasta_records = list(SeqIO.parse(input.scaftig_fasta, "fasta"))

        for record in fasta_records:

            seq_name = record.id
            raw_sequence = str(record.seq)
            gaps = list(gap_pattern.finditer(raw_sequence))

            non_gap_coords = []
            prev_gap_end = 0

            if len(gaps) == 0:
                contig_fasta_records.append(SeqRecord(Seq(raw_sequence), id=seq_name, description=""))
                continue

            for gap in gaps:
                g_start, g_end = gap.start(), gap.end() # g_start; 0-based
                g_start += 1 # converted 1-based
                if g_start == 1:
                    prev_gap_end = g_end
                    continue
                non_gap_coords.append((prev_gap_end+1, g_start-1)) # 1-based coords for non-gap segment
                prev_gap_end = g_end    
            if prev_gap_end < len(raw_sequence):
                non_gap_coords.append((prev_gap_end+1, len(raw_sequence)))

            for non_gap_coord in non_gap_coords:
                non_gap_start, non_gap_end = non_gap_coord
                contig_name = f"{seq_name}_sub_{non_gap_start}-{non_gap_end}"
                contig_sequence = raw_sequence[non_gap_start-1:non_gap_end]
                contig_fasta_records.append(SeqRecord(Seq(contig_sequence), id=contig_name, description=""))

        with open(contig_fasta, "w") as fout:
            fasta_writer = FastaWriter(fout, wrap=None)
            fasta_writer.write_file(contig_fasta_records)
        pysam.faidx(contig_fasta)

rule get_telo_stats:
    input:
        scaftig_fasta = "fcs_cleaned_fasta/{sample}/{sample}.fasta",
        contig_fasta = "fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta",
    output:
        scaftig_telo_tbl = "stats/telo/{sample}.scaftig.telo.tbl",
        contig_telo_tbl = "stats/telo/{sample}.contig.telo.tbl",
    singularity:
        "docker://eichlerlab/binf-basics:0.1",
    threads: 1,
    resources:
        hrs=4,
        mem=6,
    shell:
        """
        echo -e "seq_name\tstart\tend\tseq_length" > {output.scaftig_telo_tsv} && seqtk telo {input.scaftig_fasta} >> {output.scaftig_telo_tbl}
        echo -e "seq_name\tstart\tend\tseq_length" > {output.contig_telo_tsv} && seqtk telo {input.contig_fasta} >> {output.contig_telo_tbl}
        """

rule get_contig_stats:
    input:
        fasta = "fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta"
    output:
        stats = "stats/seq_stats/{sample}.contig.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_scaftig_stats:
    input:
        fasta = "fcs_cleaned_fasta/{sample}/{sample}.fasta"
    output:
        stats = "stats/seq_stats/{sample}.scaftig.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

# rule summarize_stats:
#     input:
#         qv_files = find_all_qv
#     output:
#         summary = "stats/summary/{asm}.summary.tsv"
#     threads: 1,
#     resources:
#         mem=8,
#         hrs=4,
#     run:

