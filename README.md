# assembly_qc
The assembly QC pipeline is an integrated pipeline that begins with screening for foreign contamination, incorporating steps from previously conducted assembly QC pipelines. The following results can be obtained:

 - Basic Statistics: Includes N50 for both scaftigs and contigs
 - Telomere Region Annotation: BED files generated using seqtk telo
 - Saffire result includes assembly to reference PAF
 - Gene Completeness Assessment: Compleasm (miniBUSCO)
 - Assembly Quality Score (QV): Merqury
 - Projection Plot
 - ModDotPlot for acrocentric regions: Available for human assemblies only


## Getting Started
Example Config: config/config_asm_qc.yaml
```commandline
### general
MANIFEST: config/manifest_asm_qc.tab
REF:
  CHM13:
    PATH: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta
    CHROMS: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/genome.txt
    CYTO: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/anno/cyto.bed ## Optional
  GRCh38:
    PATH: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/hg38.no_alt.fa
    CHROMS: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/genome.txt
    CYTO: /net/eichler/vol28/eee_shared/assemblies/hg38/no_alt/anno/cyto.bed ## Optional

### compleasm
DB_PATH: /net/eichler/vol28/eee_shared/buscodb/
LINEAGE: primates
MODE: busco

### foreign contamination screening
TAXID: 9606
INCLUDE_MITO: true
```

Example Manifest ( tab-delemetered )
```commandline
SAMPLE    H1    H2    UNASSIGNED      ILLUMINA      FOFN    TRIO    MO_ID   FA_ID
SAMPLE  test_data/asm/verkko.hap1.fasta test_data/asm/verkko.hap2.fasta test_data/asm/verkko.unassigned.fasta SAMPLE  test_data/fofn/reads.fofn NO  NA  NA
```

## Anlaysis
Begin with a dry-run
```commandline
./runcluster 30 -np
```

If dry-run looks good, proceed with:
```commandline
./runcluster 30 
```

## Major Analysis Output


```commandline
.
├── resources
│   ├── acro_target_beds
│   │   └── *.bed
│   ├── meryl
│   │   └── SAMPLE
│   │       └── SAMPLE_all.meryl
│   └── reference
│       └── {REF}
│           ├── genome.fa
│           ├── genome.fa.fai
│           └── genome_index.txt
└── results
    └── SAMPLE
        ├── assembly_eval_config
        │   └── output
        │       └── config_file
        │           └── SAMPLE.config.yaml
        ├── chain_files
        │   ├── outputs
        │   │   └── {REF}
        │   │       ├── {hap}.query_to_target.chain
        │   │       └── {hap}.query_to_target.invert.chain
        │   └── work
        │       └── filter_paf
        ├── compleasm
        │   ├── outputs
        │   │   └── summary
        │   │       └── SAMPLE.summary.tsv
        │   └── work
        │       └── summary
        │           └── {hap}
        │               ├── primates_odb10
        │               └── summary.txt
        ├── complete_flag
        ├── contamination_screening
        │   ├── outputs
        │   │   ├── contig_fasta
        │   │   │   └── SAMPLE_{hap}.fasta
        │   │   ├── final_fasta
        │   │   │   └── SAMPLE_{hap}.fasta
        │   │   ├── mito_fasta
        │   │   │   └── {hap}-mt.fasta
        │   │   └── rdna_fasta
        │   │       └── {hap}-rdna.fasta
        │   └── work
        │       ├── blast
        │       │   ├── beds
        │       │   └── outputs
        │       ├── extract_mito
        │       │   └── flags
        │       ├── fastq_cleaning
        │       │   ├── gx_adapt_cleaning
        │       │   │   ├── beds
        │       │   │   └── cleaned_fasta
        │       │   ├── rdna_cleaning
        │       │   │   └── beds
        │       │   └── short_contig_flitering
        │       │       └── cleaned_fasta
        │       ├── fcs_adaptor
        │       │   └── outputs
        │       │       └── {hap}
        │       │           ├── cleaned_sequences
        │       │           │   └── {hap}.filtered.fasta
        │       │           ├── *.log
        │       │           ├── *.jsonl
        │       │           ├── fcs_adaptor_report.txt
        │       │           ├── pipeline_args.yaml
        │       │           └── validate_fasta.txt
        │       ├── fcs_gx
        │       │   ├── flags
        │       │   └── outputs
        │       │       └── {hap}
        │       │           ├── {hap}.filtered.9606.fcs_gx_report.txt
        │       │           └── {hap}.filtered.9606.taxonomy.rpt
        │       ├── hap2_{hap}_merge_for_nucfreq
        │       │   └── flags
        │       ├── rdna_cleaning
        │       │   └── cleaned_fasta
        │       └── temp
        ├── merqury
        │   ├── outputs
        │   │   ├── logs
        │   │   ├── SAMPLE.qv
        │   │   ├── SAMPLE_{hap}.completeness.stats
        │   │   ├── SAMPLE_{hap}.qv
        │   └── work
        │       └── SAMPLE_run{hap}.sh
        ├── moddotplot
        │   ├── outputs
        │   │   ├── contig_stats
        │   │   │   ├── {hap}.CHM13_lifted_contigs.tab
        │   │   │   └── {hap}.CHM13_lifted_contigs.tsv
        │   │   ├── plots
        │   │   │   └── {hap}
        │   │   │       ├── {acro_chr}_{position}_{contig_name}_FULL.pdf
        │   │   │       └── {acro_chr}_{position}_{contig_name}_FULL.png
        │   │   └── summary
        │   │       └── {hap}.generated_acros.tsv
        │   └── work
        │       ├── find_tigs
        │       │   ├── beds
        │       │   ├── flags
        │       │   └── pafs
        │       ├── get_pq_tigs
        │       │   ├── fasta
        │       │   └── flags
        │       ├── liftover
        │       │   └── CHM13 *fixed
        │       │       ├── pafs
        │       │       ├── paf_stats
        │       │       └── trimmed_pafs
        │       └── selfplot
        │           └── flags
        ├── plots
        │   └── outputs
        │       ├── contig_length
        │       │   └── {hap}.scaffold.scatter_logged.png
        │       ├── ideo
        │       │   └── {REF}
        │       │       └── pdf
        │       │           ├── SAMPLE.minimap2.ideoplot.pdf
        │       │           └── SAMPLE.minimap2.ideoplot_wide.pdf
        │       └── ploidy
        │           └── CHM13 *fixed
        │               ├── pdf
        │               │   └── SAMPLE.minimap2.ploidy.pdf
        │               └── summary
        │                   └── SAMPLE.minimap2.ploidy_summary.txt
        ├── saffire
        │   ├── outputs
        │   │   ├── chrom_cov
        │   │   │   └── {REF}
        │   │   │       └── {hap}.minimap2.chrom_cov.tsv
        │   │   ├── safs
        │   │   │   └── {REF}
        │   │   │       └── {hap}.minimap2.saf
        │   │   └── trimmed_pafs
        │   │       └── {REF}
        │   │           └── {hap}.minimap2.trimmed.paf
        │   └── work
        │       ├── alignments
        │       │   └── {REF}
        │       │       ├── beds
        │       │       └── pafs
        │       ├── chrom_cov
        │       │   └── flags
        │       ├── combine_paf
        │       │   ├── flags
        │       │   └── {REF}
        │       └── make_saf
        │           └── flags
        └── stats
            ├── outputs
            │   ├── summary
            │   │   └── SAMPLE.summary.stats
            │   └── summary_by_hap
            │       └── {hap}.summary.stats
            └── work
                ├── fasta_stats
                │   ├── full_genome.contig.stats
                │   ├── full_genome.stats
                │   ├── {hap}.contig.stats
                │   └── {hap}.scaffold.stats
                ├── full_genome
                │   ├── contig_fasta
                │   │   └── SAMPLE.fasta
                │   └── SAMPLE.fasta
                └── telo
                    ├── flags
                    │   └── *.done
                    ├── full_genome_contigs.telo.tbl
                    ├── full_genome.telo.tbl
                    ├── {hap}.contig.telo.tbl
                    └── {hap}.scaffold.telo.tbl
```
## Partial Reulsts Options
 - get_saf_only : saffire results including asm to ref paf
 - get_busco_only : BUSCO results (compleasm)
 - get_cleaned_fasta_only : Cleaned fasta
 - get_stats_only : basic statistics such as N50
 - get_qv_only: QV score
 - get_plots_only : ploidy, projection, 
 - get_moddotplots_only: ModDotPlot

```commandline
./runcluster 30 get_qv get_compleasm
```


## Overview
![pipeline vector](https://github.com/youngjun0827/assembly_qc/blob/main/rules/assembly_qc.workflow.png)
