from pathlib import Path

MODE = config.get("MODE", "lite")
LINEAGE = config.get("LINEAGE", "primates")
BUSCO_DB_PATH = config.get("DB_DIR", "/net/eichler/vol28/eee_shared/buscodb/")


def find_all_compleasm_results(wildcards):
    avail_results = []
    haps = ["hap1","hap2","un"]
    values = full_manifest_df.loc[wildcards.sample][["H1","H2","UNASSIGNED"]].tolist()
    for idx, hap in enumerate(haps):
        if not str(values[idx]) == "nan":
            avail_results.append (f"results/{wildcards.sample}/compleasm/work/summary/{hap}/summary.txt")
    return avail_results


localrules: compleasm_run

rule compleasm_run:
    input:
        asm_fasta=rules.rename_fasta.output.final_fasta,
    output:
        summary = "results/{sample}/compleasm/work/summary/{hap}/summary.txt",
    threads: 16
    resources:
        mem=16,
        hrs=2,
    params:
        db_path=BUSCO_DB_PATH,
        mode=MODE,
        lineage=LINEAGE,
    singularity:
        "docker://eichlerlab/compleasm:0.2.6"
    shell:
        """
        compleasm.py run -a {input.asm_fasta} -o results/{wildcards.sample}/compleasm/work/summary/{wildcards.hap} -t {threads} -l {params.lineage} -m {params.mode} -L {params.db_path}
        """


rule summarize_compleasm_results:
    input:
        all_results = find_all_compleasm_results
    output:
        summary = "results/{sample}/compleasm/outputs/summary/{sample}.summary.tsv"
    threads: 1,
    resources:
        mem=4,
        hrs=1,
    run:
        summary_data = []
        compleasm_cols = ["ASSEMBLY", "HAPLOTYPE", "S", "D", "F", "I", "M", "N"]
        for hap_result in input.all_results:
            path_token = hap_result.split("/")
            sample = path_token[1]
            hap = path_token[5]

            hap_record = {"ASSEMBLY":sample, "HAPLOTYPE":hap}
            with open(hap_result) as finp:
                for line in finp:
                    if not line.strip() or line.startswith("#"):
                        continue
                    key, remains = line.strip().split(":")
                    value = int(remains.split(",")[-1])
                    hap_record[key] = value
            summary_data.append(hap_record)
        df = pd.DataFrame(summary_data)
        df.reindex(columns=compleasm_cols)
        df.to_csv(output.summary, sep="\t", index=False)
        

        
            
                    

    