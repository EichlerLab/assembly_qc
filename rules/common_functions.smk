
def find_cleaned_hap_scaffold_fasta(wildcards):
    return f"results/{wildcards.sample}/contamination_screening/outputs/final_fasta/{wildcards.sample}_{wildcards.hap}.fasta"

def find_cleaned_hap_contig_fasta(wildcards):
    return f"results/{wildcards.sample}/contamination_screening/outputs/contig_fasta/{wildcards.sample}_{wildcards.hap}.fasta"