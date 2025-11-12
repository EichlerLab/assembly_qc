import glob


wildcard_constraints:
    read_sample="|".join(full_manifest_df["ILLUMINA"].unique()),
    sample="|".join(full_manifest_df.index),
    read="\d+",


def find_cleaned_hap_fasta(which_one):
    def inner(wildcards):
        if which_one == "all":
            if not os.path.isfile(full_manifest_df.at[wildcards.sample, "H2"]):
                haps = ["hap1"]
            else:
                haps = ["hap1","hap2"]
        elif which_one in ["hap1","hap2","un"]:
            haps = [which_one]
        else:
            print (f"input type {which_one} is not available.")
            sys.exit(1)
        outputs = [ f"results/{wildcards.sample}/contamination_screening/outputs/final_fasta/{wildcards.sample}_{hap}.fasta" for hap in haps ]
        if len(outputs) == 1:
            return outputs[0]
        else:
            return outputs
    return inner


def find_read_fastq(wildcards):
    fofn_df = pd.read_csv(
        list(full_manifest_df[full_manifest_df["ILLUMINA"] == wildcards.read_sample]["FOFN"].unique())[
            0
        ],
        sep="\t",
        header=None,
        names=["fastq"],
    )
    return fofn_df.at[int(wildcards.read), "fastq"]


def agg_reads(wildcards):
    fofn_df = pd.read_csv(
        list(full_manifest_df[full_manifest_df["ILLUMINA"] == wildcards.read_sample]["FOFN"].unique())[
            0
        ],
        sep="\t",
        header=None,
        names=["fastq"],
    )
    return directory(
        expand(
            rules.run_meryl.output.meryl, read_sample=wildcards.read_sample, read=fofn_df.index
        )
    )


def find_meryl(wildcards):
    return directory(
        rules.meryl_combine.output.meryl.format(
            read_sample=full_manifest_df.at[wildcards.sample, "ILLUMINA"],
        )
    )


def find_trios(wildcards):
    return expand("results/{sample}/merqury/outputs/trio/{sample}_trio.spectra-sample.st.png",
        sample=full_manifest_df[full_manifest_df["TRIO"] == "YES"].index,
    )


def find_mat_meryl(wildcards):
    mother = full_manifest_df.loc[
        (full_manifest_df["ILLUMINA"] == wildcards.read_sample) & (full_manifest_df["MO_ID"] != "NA")
    ].iloc[0]["MO_ID"]
    return directory(expand(rules.meryl_combine.output.meryl, read_sample=mother))


def find_pat_meryl(wildcards):
    father = full_manifest_df.loc[
        (full_manifest_df["ILLUMINA"] == wildcards.read_sample) & (full_manifest_df["FA_ID"] != "NA")
    ].iloc[0]["FA_ID"]
    return directory(expand(rules.meryl_combine.output.meryl, read_sample=father))


def find_hapmer_hist(wildcards):
    return expand(
        rules.hapmers.output.inherited_hist,
        read_sample=full_manifest_df.at[wildcards.sample, "ILLUMINA"],
    )


def find_hapmers(wildcards):
    read_sample = full_manifest_df.at[wildcards.sample, "ILLUMINA"]
    mother = full_manifest_df.at[wildcards.sample, "MO_ID"]
    father = full_manifest_df.at[wildcards.sample, "FA_ID"]
    return " ".join(
        [
            os.path.abspath(f"resources/meryl/{read_sample}/{read_sample}_all.meryl"),
            os.path.abspath(
                f"resources/meryl/{read_sample}/hapmers/{read_sample}_all_and_{mother}_all.only.meryl"
            ),
            os.path.abspath(
                f"resources/meryl/{read_sample}/hapmers/{read_sample}_all_and_{father}_all.only.meryl"
            ),
        ]
    )


localrules:
    agg_hapmers,


rule non_trio_only:
    input:
        expand("results/{sample}/merqury/outputs/{sample}.qv",
            sample=full_manifest_df.loc[full_manifest_df["H1"] != "NA"].index,
        ),


rule trio_only:
    input:
        find_trios,


rule run_meryl:
    input:
        fastq=find_read_fastq,
    output:
        meryl=temp(directory("resources/meryl/{read_sample}/{read}.meryl")),
    resources:
        mem=lambda wildcards, attempt: attempt * 120,
        hrs=48,
    threads: 1
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    shell:
        """
        meryl k=21 count memory={resources.mem} {input.fastq} output {output.meryl}
        """


rule meryl_combine:
    input:
        meryl=agg_reads,
    output:
        meryl=directory("resources/meryl/{read_sample}/{read_sample}_all.meryl"),
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 128,
        hrs=48,
    threads: 1
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    shell:
        """
        if [[ $( echo {input.meryl} | wc -w ) == 1 ]]; then
            cp -rl {input.meryl} {output.meryl}
        else
            meryl union-sum {input.meryl} output {output.meryl}
        fi
        """


rule merqury_script:
    input:
        meryl=find_meryl,
        cleaned_haps = find_cleaned_hap_fasta(which_one="all")
    output:
        run_script="results/{sample}/merqury/work/{sample}_run.sh",
    params:
        unassigned_contig=find_cleaned_hap_fasta(which_one="un"),
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 4,
        hrs=2,
    threads: 1
    run:
        meryl_abs = os.path.abspath(input.meryl)
        haps_abs = [ os.path.abspath(cleaned_hap) for cleaned_hap in input.cleaned_haps ]
        sample_all = " ".join(haps_abs)
        with open(output.run_script, "w") as outfile:
            outfile.write("#!/usr/bin/env bash \n")
            outfile.write(f"merqury.sh {meryl_abs} {sample_all} {wildcards.sample}\n")
            if os.path.isfile(params.unassigned_contig):
                unassigned_abs = os.path.abspath(params.unassigned_contig)
                outfile.write(f"merqury.sh {meryl_abs} {unassigned_abs} {wildcards.sample}_unassigned\n")
        os.chmod(output.run_script, 0o755)


rule merqury_run:
    input:
        run_script=rules.merqury_script.output.run_script,
    output:
        png="results/{sample}/merqury/outputs/{sample}.qv",
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 4,
        hrs=96,
    threads: 16
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    shell:
        """
        pushd merqury/results/{wildcards.sample}/; ./$( basename {input.run_script} ); popd
        """


rule hapmers:
    input:
        mat_meryl=find_mat_meryl,
        pat_meryl=find_pat_meryl,
        sample_meryl=rules.meryl_combine.output.meryl,
    output:
        inherited_hist="resources/meryl/{read_sample}/hapmers/outputs/inherited_hapmers.hist",
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 4,
        hrs=96,
    threads: 8
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    shell:
        """
        pushd $( dirname {output.inherited_hist} )
        $( dirname $( which merqury.sh ) )/trio/hapmers.sh ../../../{input.mat_meryl} ../../../{input.pat_meryl} ../../../{input.sample_meryl}
        popd
        """


rule agg_hapmers:
    input:
        find_hapmer_hist,
    output:
        hflag=touch(temp("results/{sample}/merqury/work/agg_hapmers/trio/flags/hapmers_done")),


rule merqury_trio_script:
    input:
        hapmer_flag=rules.agg_hapmers.output.hflag,
        sample_meryl=find_meryl,
        hap_one=find_cleaned_hap_fasta(which_one="hap1"),
        hap_two=find_cleaned_hap_fasta(which_one="hap2"),
    output:
        run_script="results/{sample}/merqury/work/trio/{sample}_run.sh",
    params:
        meryl_hapmers=find_hapmers,
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 4,
        hrs=2,
    threads: 1
    run:
        meryl_all = params.meryl_hapmers
        sample_all = " ".join(
            [os.path.abspath(input.hap_one), os.path.abspath(input.hap_two)]
        )
        out_sample_name = "_".join([wildcards.sample, "trio"])
        with open(output.run_script, "w") as outfile:
            outfile.write("#!/bin/bash \n")
            outfile.write(f"merqury.sh {meryl_all} {sample_all} {out_sample_name}\n")
        os.chmod(output.run_script, 0o755)


rule merqury_trio_run:
    input:
        run_script=rules.merqury_trio_script.output.run_script,
    output:
        png="results/{sample}/merqury/outputs/trio/{sample}_trio.spectra-sample.st.png",
    resources:
        mem=lambda wildcards, attempt: (2 ** (attempt-1)) * 4,
        hrs=96,
    threads: 16
    singularity:
        "docker://eichlerlab/merqury:1.3.1"
    shell:
        """
        pushd merqury/results/{wildcards.sample}/trio; ./$( basename {input.run_script} ); popd
        """
