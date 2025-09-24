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

raw_manifest_df.set_index("SAMPLE", inplace=True)

#-----------------------------------------

def find_scaffold_fasta(wildcards):
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
    telo_tsvs = [f"stats/telo/{sample}.scaffold.telo.tsv" for sample in telo_haps] + [f"stats/telo/{sample}.contig.telo.tsv" for sample in telo_haps]
    return telo_tsvs

def find_all_qv(wildcards):
    all_qv = []
    return all_qv

rule split_scaffolds:
    input:
        scaffold_fasta=find_scaffold_fasta,
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
        fasta_records = list(SeqIO.parse(input.scaffold_fasta, "fasta"))

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
        scaffold_fasta = "fcs_cleaned_fasta/{sample}/{sample}.fasta",
        contig_fasta = "fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta",
    output:
        scaffold_telo_tbl = "stats/telo/{sample}.scaffold.telo.tbl",
        contig_telo_tbl = "stats/telo/{sample}.contig.telo.tbl",
    singularity:
        "docker://eichlerlab/binf-basics:0.1",
    threads: 1,
    resources:
        hrs=4,
        mem=6,
    shell:
        """
        echo -e "seq_name\tstart\tend\tseq_length" > {output.scaffold_telo_tbl} && seqtk telo {input.scaffold_fasta} >> {output.scaffold_telo_tbl}
        echo -e "seq_name\tstart\tend\tseq_length" > {output.contig_telo_tbl} && seqtk telo {input.contig_fasta} >> {output.contig_telo_tbl}
        """

rule get_contig_stats:
    input:
        fasta = "fcs_cleaned_fasta/{sample}/contig_fasta/{sample}.fasta",
        telo_tbl = "stats/telo/{sample}.contig.telo.tbl",
    output:
        stats = "stats/seq_stats/{sample}.contig.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_scaffold_stats:
    input:
        fasta = "fcs_cleaned_fasta/{sample}/{sample}.fasta",
        chrom_cov = "saffire/CHM13/results/{sample}/beds/{sample}.minimap2.chrom_cov.tsv",
        telo_tbl = "stats/telo/{sample}.scaffold.telo.tbl",
    output:
        stats = "stats/seq_stats/{sample}.scaffold.stats",
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

#number_of_contigs	number_of_scaffolds	number_of_near_T2T_contigs	number_of_near_T2T_scaffolds	scaffold_l50	scaffold_n50	contig_l50	contig_n50	total_ungapped_length	largest_contig_size	gc_content  percent_single_copy	percent_multi_copy	percent_fragmented	percent_missing