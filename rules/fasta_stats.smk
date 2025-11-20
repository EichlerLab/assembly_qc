import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import SeqIO
import pysam

SNAKEMAKE_ROOT_DIR = os.path.dirname(workflow.snakefile).replace("/rules","")

#-----------------------------------------

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
        scaffold_fasta=rules.rename_fasta.output.final_fasta,
    output:
        contig_fasta = "results/{sample}/contamination_screening/outputs/contig_fasta/{sample}_{hap}.fasta"
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
        scaffold_fasta = rules.rename_fasta.output.final_fasta,
        contig_fasta = rules.split_scaffolds.output.contig_fasta
    output:
        scaffold_telo_tbl = "results/{sample}/stats/work/telo/{hap}.scaffold.telo.tbl",
        scaffold_flag = "results/{sample}/stats/work/telo/flags/scaffold_{hap}.done",
        contig_telo_tbl = "results/{sample}/stats/work/telo/{hap}.contig.telo.tbl",
        config_flag = "results/{sample}/stats/work/telo/flags/contig_{hap}.done",
    singularity:
        "docker://eichlerlab/binf-basics:0.1",
    threads: 1,
    resources:
        hrs=1,
        mem=6,
    shell:
        """
        echo -e "seq_name\tstart\tend\tseq_length" > {output.scaffold_telo_tbl} && seqtk telo {input.scaffold_fasta} >> {output.scaffold_telo_tbl}
        touch {output.scaffold_flag}
        echo -e "seq_name\tstart\tend\tseq_length" > {output.contig_telo_tbl} && seqtk telo {input.contig_fasta} >> {output.contig_telo_tbl}
        touch {output.config_flag}
        """

rule get_contig_stats:
    input:
        fasta = rules.split_scaffolds.output.contig_fasta,
        telo_tbl = rules.get_telo_stats.output.contig_telo_tbl
    output:
        stats = "results/{sample}/stats/work/fasta_stats/{hap}.contig.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule get_scaffold_stats:
    input:
        fasta = rules.rename_fasta.output.final_fasta,
        telo_tbl = rules.get_telo_stats.output.scaffold_telo_tbl,
    output:
        stats = "results/{sample}/stats/work/fasta_stats/{hap}.scaffold.stats",
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: attempt * 16,
        hrs=4,
    script:
        f"{SNAKEMAKE_ROOT_DIR}/scripts/fasta_stats.py"

rule summary_hap_stats:
    input:
        scaffold_stats = rules.get_scaffold_stats.output.stats,
        contig_stats = rules.get_scaffold_stats.output.stats,
        covered_chrom_tsv = "results/{sample}/saffire/outputs/chrom_cov/CHM13/{hap}.minimap2.chrom_cov.tsv"
    output:
        hap_summary = "results/{sample}/stats/outputs/summary_by_hap/{hap}.summary.stats"
    threads: 1,
    resources:
        hrs=1,
        mem=4,
    run:
        scaffold_stats = input.scaffold_stats
        contig_stats = input.contig_stats
        covered_chrom_tsv = input.covered_chrom_tsv
        df_scaffold = pd.read_csv(scaffold_stats, sep="\t").add_suffix("_scaffold")
        df_scaffold = df_scaffold.rename(columns = {"sample_scaffold":"sample", "haplotype_scaffold":"haplotype"})
        
        df_contig = pd.read_csv(contig_stats, sep="\t").add_suffix("_contig")
        df_contig = df_contig.rename(columns = {"sample_contig":"sample", "haplotype_contig":"haplotype"})

        df = pd.merge(df_scaffold, df_contig, on=["sample","haplotype"], how="outer")
        print (df.info())
        



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