#!/bin/env python

import sys
import os
import pandas as pd
import numpy as np
import logging
import pysam
import re

#fasta = input.fasta
#fai = f"{fasta}.fai"
#tsv = output.tsv

fasta_set_files = snakemake.input.fasta_set
telo_tsvs = snakemake.input.telo_tsvs
output_files = snakemake.output

nt_df_header = ["sample","asm_type","haplotype","num_of_contigs","num_of_scaftigs","bases","ungapped_bases","cnt_A","cnt_T","cnt_G","cnt_C","cnt_N","fai"]

def get_telo_info(telo_tsv):
    if re.search(".scaftig.", telo_tsv):
        asm_type = "scaftig"
    else:
        asm_type = "contig"
    telo_df = pd.read_csv(telo_tsv, sep="\t", header=0)
    telo_df["width"] = telo_df["end"] - telo_df["start"]
    telo_df = telo_df[telo_df["width"]>=100]
    telo_df[["forward","backward"]] = 0
    telo_df.loc[(telo_df["start"] == 0), "forward"] = 1
    telo_df.loc[(telo_df["end"] == telo_df["backward"])] = 1
    print (telo_df)
    

def get_nt_length(fasta):
    
    if re.search("/contig_fasta/", fasta):
        asm_type = "contig"
    else:
        asm_type = "scaftig"
    fasta_name = os.path.basename(fasta)
    for hap_tag in ["_hap1.fasta","_hap2.fasta","_unassigned.fasta"]:
        if re.search(hap_tag, fasta_name):
            sample = fasta_name.replace(hap_tag,"")
            haplotype = hap_tag.split("_")[1].split(".")[0]
            break
    
    fai = f"{fasta}.fai"

    cnt_A = 0
    cnt_T = 0
    cnt_G = 0
    cnt_C = 0
    cnt_N = 0
    bases = 0
    ungapped_bases = 0
    num_of_contigs = 0
    num_of_scaftigs = 0

    with pysam.FastaFile(fasta) as asm_fasta:
        for ref in asm_fasta.references:
            ref_bases = len(asm_fasta.fetch(ref))
            ref_cnt_N = asm_fasta.fetch(ref).upper().count("N")
            cnt_A += asm_fasta.fetch(ref).upper().count("A")
            cnt_T += asm_fasta.fetch(ref).upper().count("T")
            cnt_G += asm_fasta.fetch(ref).upper().count("G")
            cnt_C += asm_fasta.fetch(ref).upper().count("C")
            cnt_N += ref_cnt_N
            bases += ref_bases
            ungapped_bases += (ref_bases-ref_cnt_N)
            if ref_cnt_N > 0:
                num_of_scaftigs += 1
            else:
                num_of_contigs += 1

    data_set = [sample, asm_type, haplotype, num_of_contigs, num_of_scaftigs, bases, ungapped_bases, cnt_A, cnt_T, cnt_G, cnt_C, cnt_N, fai]
    
    return pd.DataFrame([data_set],columns = nt_df_header)


def get_contig_stats(fai_files):
    len_list = pd.Series(dtype=int)
    for fai in fai_files.split():
        fai_df = pd.read_csv(fai, sep="\t", header=None, usecols=[0,1])
        len_list = pd.concat([len_list , pd.Series(fai_df[1].copy())])
    vals = len_list.sort_values(ascending=False)
    vals_csum = np.cumsum(vals)
    num_contigs = len(len_list)

    over_100k_bases = np.sum(len_list.loc[len_list >= 100000])
    n50 = vals.iloc[np.sum(vals_csum <= (vals_csum.iloc[-1] // 2))]
    l50 = np.sum(vals > n50)
    largest_size = max(len_list)
    aun = np.sum(len_list*len_list)/np.sum(len_list)

    return pd.Series([int(over_100k_bases), int(l50), int(n50), int(largest_size), int(aun)])

def main():
    output = output_files[0]
    nt_df = pd.DataFrame(columns = nt_df_header)
    for telo_tsv in telo_tsvs:
        get_telo_info(telo_tsv)

    for fasta in fasta_set_files:
        print (fasta)
        continue
        df = get_nt_length(fasta)
        nt_df = pd.concat([nt_df, df], ignore_index=True)

    sys.exit()
    nt_df = nt_df[~((nt_df["asm_type"] == "scaftig") & (nt_df["cnt_N"] == 0))].reset_index(drop=True)
    samples = nt_df["sample"].unique()
    asm_types = nt_df["asm_type"].unique()
    numeric_cols = ["num_of_contigs","num_of_scaftigs","bases","ungapped_bases","cnt_A","cnt_T","cnt_G","cnt_C","cnt_N"]
    for sample in samples:
        for asm_type in asm_types:
            print (sample, asm_type)
            all_assigned_df = nt_df[(nt_df["sample"] == sample) & (nt_df["asm_type"] == asm_type) & (nt_df["haplotype"] != "unassigned")]
            summed_nt_values = all_assigned_df[numeric_cols].sum()
            new_row_df = pd.DataFrame({
                "sample": [sample],
                "asm_type": [asm_type],
                "haplotype": ["all_assigned"],
                "num_of_contigs": [summed_nt_values["num_of_contigs"]],
                "num_of_scaftigs": [summed_nt_values["num_of_scaftigs"]],
                "bases": [summed_nt_values["bases"]],
                "ungapped_bases": [summed_nt_values["ungapped_bases"]],
                "cnt_A": [summed_nt_values["cnt_A"]],
                "cnt_T": [summed_nt_values["cnt_T"]],
                "cnt_G": [summed_nt_values["cnt_G"]],
                "cnt_C": [summed_nt_values["cnt_C"]],
                "cnt_N": [summed_nt_values["cnt_N"]],
                "fai": [" ".join(all_assigned_df["fai"].unique())],
            })
            print (new_row_df)
            print (nt_df)
            nt_df = pd.concat([nt_df, new_row_df])
            print ("MERGED")
            print (nt_df)
    nt_df["gc_content"] = ((nt_df["cnt_G"]+nt_df["cnt_C"])/nt_df["ungapped_bases"]) * 100
    contig_stats_header = ["over_100k_bases", "l50", "n50", "largest_size", "aun"]
    nt_df[contig_stats_header] = nt_df["fai"].apply(get_contig_stats)
    nt_df = nt_df[["sample","asm_type","haplotype","num_of_contigs","num_of_scaftigs","bases","gc_content","l50", "n50", "largest_size", "aun"]]
#    nt_df.to_csv(output, sep="\t", index=None, header=True)
    

if __name__=="__main__":
    main()