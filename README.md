# assembly_qc
assembly QC pipeline starts with foreign contamination screening


## Major Analysis Output

QC_results/
├── compleasm
│   └── results
│       ├── unique_sample_name_hap1
│       │   ├── primates_odb10/
│       │   └── **summary.txt**
│       └── unique_sample_name_hap2
│           ├── primates_odb10/
│           └── **summary.txt**
├── contamination_screening
│   ├── results
│   │   ├── unique_sample_name_hap1
│   │   │   ├── clean/
│   │   │   ├── fasta
│   │   │   │   ├── **unique_sample_name_hap1.fasta**
│   │   │   │   └── unique_sample_name_hap1-rdna.fasta
│   │   │   ├── fcs_adaptor/
│   │   │   │   ├── cleaned_sequences/
│   │   │   │   └── **fcs_adaptor_report.txt**
│   │   │   ├── fcs_gx
│   │   │   │   └── **assembly.haplotype1.9606.fcs_gx_report.txt**
│   │   │   ├── fcs_mito
│   │   │   │   └── **fcs_mito.txt**
│   │   │   ├── rdna
│   │   │   │   └── **rdna.txt**
│   │   │   └── trim.bed
│   │   └── unique_sample_name_hap2
│   │       ├── clean/
│   │       ├── fasta
│   │       │   ├── **unique_sample_name_hap2.fasta**
│   │       │   └── unique_sample_name_hap2-rdna.fasta
│   │       ├── fcs_adaptor/
│   │       │   ├── cleaned_sequences/
│   │       │   └── **fcs_adaptor_report.txt**
│   │       ├── fcs_gx
│   │       │   └── **assembly.haplotype2.9606.fcs_gx_report.txt**
│   │       ├── fcs_mito
│   │       │   └── **fcs_mito.txt**
│   │       ├── rdna
│   │       │   └── **rdna.txt**
│   │       └── trim.bed
│   └── temp
│       ├── unique_sample_name_hap1
│       │   ├── fasta
│       │   │   └── **unique_sample_name_hap1.fasta**
│       │   └── regions.out
│       └── unique_sample_name_hap2
│           ├── fasta
│           │   └── **unique_sample_name_hap2.fasta**
│           └── regions.out
├── contig_stats
│   ├── n50
│   │   ├── **unique_sample_name_hap1.n50.stats**
│   │   └── **unique_sample_name_hap2.n50.stats**
│   └── plots
│       ├── **unique_sample_name_hap1.dist.log.png**
│       ├── **unique_sample_name_hap1.scatter.log.png**
│       ├── **unique_sample_name_hap2.dist.log.png**
│       └── **unique_sample_name_hap2.scatter.log.png**
├── ideo_plots
│   └── {aligner}
│       └── **unique_sample_name.ideoplot.pdf**
├── merqury
│   ├── meryl/
│   └── results
│       └── unique_sample_name
│           ├── unique_sample_name.unique_sample_name_hap1.qv
│           ├── unique_sample_name.unique_sample_name_hap2.qv
│           ├── unique_sample_name_hap1.meryl/
│           ├── unique_sample_name_hap2.meryl/
│           ├── **unique_sample_name.qv**
│           └── completeness.stats
├── ploidy_plots
│   └── minimap2
│       ├── **unique_sample_name.ploidy.pdf**
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

## Overview
![pipeline vector](https://github.com/youngjun0827/assembly_qc/blob/main/rules/assembly_qc.workflow.png)
