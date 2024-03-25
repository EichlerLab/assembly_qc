# assembly_qc
The assembly QC pipeline begins by screening for foreign contamination.</br>
Once contaminants have been removed from the assembly in FASTA format, the pipeline generates QC results utilizing tools like Merqury and Compleam (a.k.a miniBUSCO). 
It also provides basic statistics, including the N50 metric, and generates plots that illustrate the characteristics of the assembly.


## Getting Started
Example Config: config/config_asm_qc.yaml
```commandline
### general
MANIFEST: config/manifest_asm_qc.tab
REF: /net/eichler/vol28/eee_shared/assemblies/CHM13/T2T/v2.0/T2T-CHM13v2.fasta

### compleasm
MB_DOWNLOADS: /net/eichler/vol28/software/modules-sw/compleasm/0.2.2/Linux/CentOS7/x86_64/compleasm_kit/mb_download/
LINEAGE: primates
MODE: busco

### foreign contamination screening
TAXID: 9606

### saffire; should be either or both "minimap2" and "winnowmap"
ALIGNER:
  - minimap2
```

Example Manifest ( tab-delemetered )
```commandline
SAMPLE    H1    H2    ILLUMINA    FOFN    TRIO    MO_ID   FA_ID
unique_sample_name    unique_sample_name/assembly.haplotype1.fasta   unique_sample_name/assembly.haplotype2.fasta    sample_name    sample_name_illumina.fofn   NO      NA      NA
```

## Anlaysis
Begin with a dry-run
```commandline
./runlocal 30 -np
```

If dry-run looks good, proceed with:
```commandline
./runlocal 30 
```

## Major Analysis Output


```commandline
QC_results
├── compleasm
│   └── results
│       ├── unique_sample_name_hap1
│       │   ├── primates_odb10/
│       │   └── summary.txt
│       └── unique_sample_name_hap2
│           ├── primates_odb10/
│           └── summary.txt
├── contamination_screening
│   ├── results
│   │   ├── unique_sample_name_hap1
│   │   │   ├── clean/
│   │   │   ├── fasta
│   │   │   │   ├── unique_sample_name_hap1.fasta
│   │   │   │   └── unique_sample_name_hap1-rdna.fasta
│   │   │   ├── fcs_adaptor/
│   │   │   │   ├── cleaned_sequences/
│   │   │   │   └── fcs_adaptor_report.txt
│   │   │   ├── fcs_gx
│   │   │   │   └── assembly.haplotype1.9606.fcs_gx_report.txt
│   │   │   ├── fcs_mito
│   │   │   │   └── fcs_mito.txt
│   │   │   ├── rdna
│   │   │   │   └── rdna.txt
│   │   │   └── trim.bed
│   │   └── unique_sample_name_hap2
│   │       ├── clean/
│   │       ├── fasta
│   │       │   ├── unique_sample_name_hap2.fasta
│   │       │   └── unique_sample_name_hap2-rdna.fasta
│   │       ├── fcs_adaptor/
│   │       │   ├── cleaned_sequences/
│   │       │   └── fcs_adaptor_report.txt
│   │       ├── fcs_gx
│   │       │   └── assembly.haplotype2.9606.fcs_gx_report.txt
│   │       ├── fcs_mito
│   │       │   └── fcs_mito.txt
│   │       ├── rdna
│   │       │   └── rdna.txt
│   │       └── trim.bed
│   └── temp
│       ├── unique_sample_name_hap1
│       │   ├── fasta
│       │   │   └── unique_sample_name_hap1.fasta
│       │   └── regions.out
│       └── unique_sample_name_hap2
│           ├── fasta
│           │   └── unique_sample_name_hap2.fasta
│           └── regions.out
├── contig_stats
│   ├── n50
│   │   ├── unique_sample_name_hap1.n50.stats
│   │   └── unique_sample_name_hap2.n50.stats
│   └── plots
│       ├── unique_sample_name_hap1.dist.log.png
│       ├── unique_sample_name_hap1.scatter.log.png
│       ├── unique_sample_name_hap2.dist.log.png
│       └── unique_sample_name_hap2.scatter.log.png
├── ideo_plots
│   └── {aligner}
│       └── unique_sample_name.ideoplot.pdf
├── merqury
│   ├── meryl/
│   └── results
│       └── unique_sample_name
│           ├── unique_sample_name.unique_sample_name_hap1.qv
│           ├── unique_sample_name.unique_sample_name_hap2.qv
│           ├── unique_sample_name_hap1.meryl/
│           ├── unique_sample_name_hap2.meryl/
│           ├── unique_sample_name.qv
│           └── completeness.stats
├── ploidy_plots
│   └── minimap2
│       ├── unique_sample_name.ploidy.pdf
│       └── unique_sample_name.summary.txt
└── saffire
    ├── benchmarks/
    └── results
        ├── unique_sample_name_hap1
        │   ├── alignments
        │   │   ├── unique_sample_name_hap1.minimap2.bam
        │   │   └── unique_sample_name_hap1.minimap2.paf
        │   └── beds
        │       └── unique_sample_name_hap1.minimap2.bed
        └── unique_sample_name_hap2
            ├── alignments
            │   ├── unique_sample_name_hap2.minimap2.bam
            │   └── unique_sample_name_hap2.minimap2.paf
            └── beds
                └── unique_sample_name_hap2.minimap2.bed
```
## Overview
![pipeline vector](https://github.com/youngjun0827/assembly_qc/blob/main/rules/assembly_qc.workflow.png)
