import pandas as pd
from pybedtools import BedTool
import numpy as np
import networkx as nx
import more_itertools as mit

PARTS=config.get('PARTS', 15)
# MINIMAP_PARAMS = config.get('MINIMAP_PARAMS', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
MINIMAP_PARAMS = config.get('MINIMAP_PARAMS', '-x asm20 -m 10000 -z 10000,50 -r 50000 --end-bonus=100 -O 5,56 -E 4,1 -B 5')
MANIFEST = config.get('MANIFEST', 'config/manifest.tab')
SV_SIZE = config.get('SV_SIZE', '30000')

#-----------------------------------------


scattergather:
    split=PARTS


def find_ref_fa(wildcards):
    return REF_DICT[wildcards.ref]["PATH"]

def find_contigs(wildcards):
    return gather.split("results/{sample}/saffire/work/split_paf/{ref}/tmp/{hap}.minimap2.{scatteritem}.broken.paf", sample=wildcards.sample, ref=wildcards.ref, hap=wildcards.hap)

def find_all_trimmed_paf(wildcards):
    avail_pafs = []
    haps = ["hap1","hap2","un"]
    values = full_manifest_df.loc[wildcards.sample][["H1","H2","UNASSIGNED"]].tolist()
    for idx, hap in enumerate(haps):
        if not str(values[idx]) == "nan":
            avail_pafs.append (f"results/{wildcards.sample}/saffire/outputs/trimmed_pafs/{wildcards.ref}/{hap}.minimap2.trimmed.paf")
    return avail_pafs


rule copy_ref_fasta:
    input:
        raw_ref=find_ref_fa,
    output:
        ref = "resources/reference/{ref}/genome.fa",
        fai = "resources/reference/{ref}/genome.fa.fai",
        genome_index = "resources/reference/{ref}/genome_index.txt"
    threads: 1,
    resources:
        mem = 4,
        hrs = 1,
    shell: """
        cp {input.raw_ref} {output.ref}
        samtools faidx {output.ref}
        cut -f1 {output.fai}> {output.genome_index}
        """

rule make_minimap_paf:
    input:
        fa = rules.rename_fasta.output.final_fasta,
        ref = "resources/reference/{ref}/genome.fa"
    output:
        paf = "results/{sample}/saffire/work/alignments/{ref}/pafs/{hap}.minimap2.paf",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',
    threads: 12
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        minimap2 -c -t {threads} -K {resources.mem}G --cs {params.map_opts} --secondary=no --eqx -Y {input.ref} {input.fa} > {output.paf}
        """

rule trim_paf:
    input:
        paf = "results/{sample}/saffire/work/alignments/{ref}/pafs/{hap}.minimap2.paf",
    output:
        paf = "results/{sample}/saffire/outputs/trimmed_pafs/{ref}/{hap}.minimap2.trimmed.paf",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',
    threads: 12
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        rustybam trim-paf --remove-contained {input.paf} > {output.paf}
        """

# added -L param in case the CIGAR is > 65535 since older tools are unable to convert alignments with >65535 CIGAR ops to BAM. (https://lh3.github.io/minimap2/minimap2.html)
rule filter_paf:
    input:
        paf=rules.trim_paf.output.paf,
    output:
        filt_paf="results/{sample}/chain_files/work/filter_paf/{ref}/{hap}.filt.paf",
        filt_invert_paf="results/{sample}/chain_files/work/filter_paf/{ref}/{hap}.filt.invert.paf",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',
    threads: 12    
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 12,
        hrs = 48
    params:
        min_aln_len=100_000,
    shell: """
        rb filter -a {params.min_aln_len} {input.paf} > {output.filt_paf}
        rb invert {output.filt_paf} > {output.filt_invert_paf}
        """

rule make_chain_file:
    input:
        filt_paf = rules.filter_paf.output.filt_paf,
        filt_invert_paf = rules.filter_paf.output.filt_invert_paf
    output:
        chain="results/{sample}/chain_files/outputs/{sample}_{hap}_To_{ref}.chain",
        chain_invert="results/{sample}/chain_files/outputs/{sample}_{hap}_To_{ref}.invert.chain",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',
    threads: 12    
    singularity:
        "docker://eichlerlab/paf2chain:0.1.1"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        paf2chain --input {input.filt_paf} > {output.chain}
        paf2chain --input {input.filt_invert_paf} > {output.chain_invert}
        """

rule make_minimap_bam:
    input:
        fa = "results/{sample}/contamination_screening/outputs/final_fasta/{sample}_{hap}.fasta",
        ref = "resources/reference/{ref}/genome.fa"
    output:
        bam = "results/{sample}/saffire/work/alignments/{ref}/bams/{hap}.minimap2.bam",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',        
    threads: 12
    params:
        map_opts = MINIMAP_PARAMS,
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 48
    shell: """
        minimap2 -c -t {threads} -K {resources.mem}G --cs -a {params.map_opts} --MD --secondary=no --eqx -Y -L {input.ref} {input.fa} | samtools view -S -b /dev/stdin | samtools sort -@ {threads} -o {output.bam}
        """

rule make_bed:
    input:
        paf = rules.trim_paf.output.paf,
        genome_index = "resources/reference/{ref}/genome_index.txt"
    output:
        bed = "results/{sample}/saffire/work/alignments/{ref}/beds/{hap}.minimap2.bed",
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',              
    threads: 1
    singularity:
        "docker://eichlerlab/align-basics:0.1"
    resources:
        mem = 12,
        hrs = 1
    shell: """
        awk -vOFS="\t" '{{print $6,$8,$9,$1,$12}}' {input.paf} | bedtools sort -g {input.genome_index} > {output.bed}
        """

rule split_paf:
    input:
        paf = "results/{sample}/saffire/work/alignments/{ref}/pafs/{hap}.minimap2.paf",
    output:
        flag = temp(scatter.split("results/{{sample}}/saffire/work/split_paf/{{ref}}/tmp/{{hap}}.minimap2.{scatteritem}.paf")),
        temp_paf = temp("results/{sample}/saffire/work/split_paf/{ref}/tmp/{hap}_uniform.minimap2.paf")
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    run:
        with open(input.paf, 'r') as infile:
            for i, line in enumerate(infile):
                if i == 0:
                    all_tags = set([x.split(':')[0] for x in line.split('\t')[12:]])
                else:
                    all_tags = all_tags.intersection(set([x.split(':')[0] for x in line.split('\t')[12:]]))

        out_list = []

        with open(input.paf, 'r') as infile:
            for line in infile:
                out_list.append(line.split('\t')[0:12]+[x for x in line.split('\t')[12:] if x.split(':')[0] in all_tags])
        output_temp_paf = output.temp_paf
        with open(output_temp_paf, 'w') as outfile:
            for item in out_list:
                outfile.write('\t'.join(item))

        df = pd.read_csv(output_temp_paf, sep='\t', low_memory=False, header=None)
        col_out = df.columns.values
        for i, contig in enumerate(df[0].unique()):
            out_num = (i % PARTS) + 1
            out_paf = f'results/{wildcards.sample}/saffire/work/split_paf/{wildcards.ref}/tmp/{wildcards.hap}.minimap2.{out_num}-of-{PARTS}.paf'
            df.loc[df[0] == contig][col_out].to_csv(out_paf, sep='\t', index=False, header=False, mode='a+')
        for check_num in range(1,PARTS+1):
            check_out = f'results/{wildcards.sample}/saffire/work/split_paf/{wildcards.ref}/tmp/{wildcards.hap}.minimap2.{check_num}-of-{PARTS}.paf'
            if (not os.path.exists(check_out)):
                open(check_out,"w").close()

rule trim_break_orient_paf:
    input:
        paf = "results/{sample}/saffire/work/split_paf/{ref}/tmp/{hap}.minimap2.{scatteritem}.paf"
    output:
        contig = temp("results/{sample}/saffire/work/split_paf/{ref}/tmp/{hap}.minimap2.{scatteritem}.orient.paf"),
        broken = temp("results/{sample}/saffire/work/split_paf/{ref}/tmp/{hap}.minimap2.{scatteritem}.broken.paf")
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    params:
        sv_size = SV_SIZE,
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 24,
        hrs = 24
    shell: """
        if [[ ! -s {input}.paf ]];then
            : > {output.contig}
            : > {output.broken}
        else
            rustybam orient {input.paf} | rustybam trim-paf | rb filter --paired-len 1000000 > {output.contig}
            rustybam break-paf --max-size {params.sv_size} {output.contig} > {output.broken}
        fi
        """

rule concat_pafs:
    input:
        paf = find_all_trimmed_paf
    output:
        paf = "results/{sample}/saffire/outputs/merged_paf/{ref}/minimap2.concat.paf",
        flag = "results/{sample}/saffire/work/merged_paf/flags/{ref}.minimap2.done"
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',             
    threads:1 
    resources:
        mem=4,
        hrs=1
    shell: """
        cat {input.paf} > {output.paf}
        touch {output.flag}
        """

rule combine_paf:
    input:
        paf = find_contigs,
        flag = rules.split_paf.output.flag
    output:
        paf = "results/{sample}/saffire/work/combine_paf/{ref}/tmp/{hap}.minimap2.broken.paf",
        flag = "results/{sample}/saffire/work/combine_paf/flags/{ref}.{hap}.minimap2.done"
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    shell: """
        cat {input.paf} > {output.paf}
        touch {output.flag}
        """

rule saff_out:
    input:
        paf = rules.combine_paf.output.paf
    output:
        saf = "results/{sample}/saffire/outputs/safs/{ref}/{hap}.minimap2.saf",
        flag = "results/{sample}/saffire/work/make_saf/flags/{ref}.{hap}.minimap2.done"
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',             
    threads: 1
    singularity:
        "docker://eichlerlab/rustybam:0.1.33"
    resources:
        mem = 8,
        hrs = 24
    shell: """
        rb stats --paf {input.paf} > {output.saf}
        touch {output.flag}
        """

rule check_covered_chromosomes:
    input:
        paf = rules.trim_paf.output.paf
    output:
        tsv = "results/{sample}/saffire/outputs/chrom_cov/{ref}/{hap}.minimap2.chrom_cov.tsv",
        flag = "results/{sample}/saffire/work/chrom_cov/flags/{ref}.{hap}.minimap2.done"
    wildcard_constraints:
        ref='[A-Za-z0-9_-]+',
    threads: 1
    resources:
        mem = 8,
        hrs = 24
    run:
        from collections import defaultdict

        chrom_length_dict = dict()
        chrom_intervals = defaultdict(list)
        with open(input.paf) as finp:
            read_lines = finp.read().strip().split("\n")
        for line in read_lines:
            token = line.split("\t")[:9]
            chrom, start, end, contigname, chrom_length = token[5], int(token[7]), int(token[8]), token[0], int(token[6])
            if not chrom in chrom_length_dict:
                chrom_length_dict[chrom] = chrom_length
            chrom_intervals[chrom].append((start,end))
        
        merged_dict = dict()
        for chrom in sorted(chrom_intervals):
            intervals = sorted(chrom_intervals[chrom])
            merged = []
            current_start, current_end = intervals[0]
            for start, end in intervals[1:]:
                if start <= current_end:
                    current_end = max(current_end, end)
                else:
                    merged.append((current_start, current_end))
                    current_start, current_end = start, end
            merged.append((current_start, current_end))
            merged_dict[chrom] = merged
        chroms = [f"chr{i}" for i in range(1,23)]+["chrX","chrY"]
        with open(output.tsv, "w") as fout:
            print ("chrom\tcovered_pct", file=fout)
            for chrom in chroms:
                try:
                    covered_length = 0
                    chrom_length = chrom_length_dict[chrom]
                    for intervals in merged_dict[chrom]:
                        covered_length += (intervals[1]-intervals[0])
                    covered_rate = covered_length/chrom_length*100
                except:
                    covered_rate = 0.0
                print (f"{chrom}\t{covered_rate:.2f}", file=fout)
        
        open(output.flag,"w").close()