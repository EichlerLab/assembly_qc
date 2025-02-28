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

### ModDotPlot
BED: config/moddot_regions.bed

### compleasm
DB_PATH: /net/eichler/vol28/eee_shared/buscodb/
LINEAGE: primates
MODE: busco

### foreign contamination screening
TAXID: 9606
```

Example Manifest ( tab-delemetered )
```commandline
SAMPLE    H1    H2    UNASSIGNED      ILLUMINA      FOFN    TRIO    MO_ID   FA_ID
unique_sample_name    unique_sample_name/assembly.haplotype1.fasta   unique_sample_name/assembly.haplotype2.fasta    NA    sample_name    sample_name_illumina.fofn   NO      NA      NA
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
.
├── assembly_stats
│   └── unique_sample_name
│       └── unique_sample_name.tsv
├── compleasm
│   └── results
│       ├── unique_sample_name_hap1
│       │   ├── {database}
│       │   └── summary.txt
│       └── unique_sample_name_hap2
│           ├── {database}
│           └── summary.txt
├── contamination_screening
│   ├── filtered_fasta
│   │   ├── unique_sample_name_hap1.filtered.fasta
│   │   ├── unique_sample_name_hap2.filtered.fasta
│   ├── raw_fasta
│   │   ├── unique_sample_name_hap1.fasta (link to the input hap1 fasta)
│   │   └── unique_sample_name_hap2.fasta (link to the input hap2 fasta)
│   └── results
│       ├── unique_sample_name_hap1
│       └── unique_sample_name_hap2
├── fcs_cleaned_fasta
│   ├── unique_sample_name_hap1
│   │   ├── unique_sample_name_hap1.fasta
│   │   └── contig_fasta
│   │       └── unique_sample_name_hap1.fasta
│   └── unique_sample_name_hap2
│       ├── unique_sample_name_hap2.fasta
│       └── contig_fasta
│           └── unique_sample_name_hap2.fasta
├── merqury
│   ├── meryl
│   │   ├── sample_name
│   │   │   └── sample_name_all.meryl
│   └── results
│       └── unique_sample_name
│           ├── unique_sample_name.unique_sample_name_hap1.qv
│           ├── unique_sample_name.unique_sample_name_hap2.qv
│           ├── unique_sample_name.completeness.stats
│           ├── unique_sample_name.qv
│           ├── unique_sample_name_hap1.meryl
│           └── unique_sample_name_hap2.meryl
├── moddotplot (only for human ; TAXID 9606)
│   ├── contigs
│   │   ├── unique_sample_name_hap1
│   │   └── unique_sample_name_hap2
│   ├── fasta
│   │   ├── unique_sample_name_hap1
│   │   └── unique_sample_name_hap2
│   ├── liftover
│   │   └── paf
│   │       ├── unique_sample_name_hap1
│   │       └── unique_sample_name_hap2
│   ├── results
│   │   ├── unique_sample_name_hap1.generated_acros.tsv
│   │   ├── unique_sample_name_hap2.generated_acros.tsv
│   │   ├── unique_sample_name_hap1 (some can be missing)
│   │   │   ├── CHROM_1_22508596_FULL.pdf
│   │   │   ├── chr13_1_22508596_FULL.png
│   │   │   ├── chr13_1_22508596_TRI.pdf
│   │   │   ├── chr13_1_22508596_TRI.png
│   │   │   ├── chr14_1_17708411_FULL.pdf
│   │   │   ├── chr14_1_17708411_FULL.png
│   │   │   ├── chr14_1_17708411_TRI.pdf
│   │   │   ├── chr14_1_17708411_TRI.png
│   │   │   ├── chr15_1_22694466_FULL.pdf
│   │   │   ├── chr15_1_22694466_FULL.png
│   │   │   ├── chr15_1_22694466_TRI.pdf
│   │   │   ├── chr15_1_22694466_TRI.png
│   │   │   ├── chr21_1_16306378_FULL.pdf
│   │   │   ├── chr21_1_16306378_FULL.png
│   │   │   ├── chr21_1_16306378_TRI.pdf
│   │   │   ├── chr21_1_16306378_TRI.png
│   │   │   ├── chr22_1_20711065_FULL.pdf
│   │   │   ├── chr22_1_20711065_FULL.png
│   │   │   ├── chr22_1_20711065_TRI.pdf
│   │   │   └── chr22_1_20711065_TRI.png
│   │   └── unique_sample_name_hap2 (some can be missing)
│   │       ├── chr13_1_22508596_FULL.pdf
│   │       ├── chr13_1_22508596_FULL.png
│   │       ├── chr13_1_22508596_TRI.pdf
│   │       ├── chr13_1_22508596_TRI.png
│   │       ├── chr14_1_17708411_FULL.pdf
│   │       ├── chr14_1_17708411_FULL.png
│   │       ├── chr14_1_17708411_TRI.pdf
│   │       ├── chr14_1_17708411_TRI.png
│   │       ├── chr15_1_22694466_FULL.pdf
│   │       ├── chr15_1_22694466_FULL.png
│   │       ├── chr15_1_22694466_TRI.pdf
│   │       ├── chr15_1_22694466_TRI.png
│   │       ├── chr21_1_16306378_FULL.pdf
│   │       ├── chr21_1_16306378_FULL.png
│   │       ├── chr21_1_16306378_TRI.pdf
│   │       ├── chr21_1_16306378_TRI.png
│   │       ├── chr22_1_20711065_FULL.pdf
│   │       ├── chr22_1_20711065_FULL.png
│   │       ├── chr22_1_20711065_TRI.pdf
│   │       └── chr22_1_20711065_TRI.png
│   └── target_beds
│       ├── chr13_1_22508596.bed
│       ├── chr14_1_17708411.bed
│       ├── chr15_1_22694466.bed
│       ├── chr21_1_16306378.bed
│       └── chr22_1_20711065.bed
├── plots
│   ├── contigs
│   │   ├── unique_sample_name_hap1.dist.log.png
│   │   ├── unique_sample_name_hap1.scatter.log.png
│   │   ├── unique_sample_name_hap2.dist.log.png
│   │   └── unique_sample_name_hap2.scatter.log.png
│   ├── ideo
│   │   └── minimap2
│   │       ├── unique_sample_name_to_CHM13.ideoplot.pdf
│   │       └── unique_sample_name_to_HG38.ideoplot.pdf
│   └── ploidy (only for CHM13)
│       └── CHM13
│           └── minimap2
│               ├── unique_sample_name.ploidy.pdf
│               └── unique_sample_name.summary.txt
├── saffire
│   ├── CHM13
│   │   ├── reference
│   │   │   ├── CHM13.fa
│   │   │   ├── CHM13.fa.fai
│   │   │   └── CHM13.genome.txt
│   │   └── results
│   │       ├── unique_sample_name_hap1
│   │       │   ├── alignments
│   │       │   │   └── unique_sample_name_hap1.minimap2.paf
│   │       │   └── beds
│   │       │       └── unique_sample_name_hap1.minimap2.bed
│   │       └── unique_sample_name_hap2
│   │           ├── alignments
│   │           │   ├── unique_sample_name_hap2.minimap2.bam
│   │           │   └── unique_sample_name_hap2.minimap2.paf
│   │           └── beds
│   │               └── unique_sample_name_hap2.minimap2.bed
│   └── HG38
│       ├── reference
│       │   ├── HG38.fa
│       │   ├── HG38.fa.fai
│       │   └── HG38.genome.txt
│       └── results
│           ├── unique_sample_name_hap1
│           │   ├── alignments
│           │   │   └── unique_sample_name_hap1.minimap2.paf
│           │   └── beds
│           │       └── unique_sample_name_hap1.minimap2.bed
│           └── unique_sample_name_hap2
│               ├── alignments
│               │   ├── unique_sample_name_hap2.minimap2.bam
│               │   └── unique_sample_name_hap2.minimap2.paf
│               └── beds
│                   └── unique_sample_name_hap2.minimap2.bed
└── stats
    ├── seq_stats
    │   ├── unique_sample_name_hap1.contig.stats
    │   ├── unique_sample_name_hap1.scaftig.stats
    │   ├── unique_sample_name_hap2.contig.stats
    │   └── unique_sample_name_hap2.scaftig.stats
    └── telo
        ├── unique_sample_name_hap1.contig.telo.tbl
        ├── unique_sample_name_hap1.scaftig.telo.tbl
        ├── unique_sample_name_hap2.contig.telo.tbl
        └── unique_sample_name_hap2.scaftig.telo.tbl
```
## Partial Reulsts Options
 - get_saf : saffire results including asm to ref paf
 - get_compleasm : BUSCO results
 - get_cleaned_fasta : Cleaned fasta
 - get_stats : basic statistics such as N50
 - get_qv: QV score
 - get_plots : ploidy, projection, ModDot
```commandline
./runcluster 30 get_qv get_compleasm
```


## Overview
![pipeline vector](https://github.com/youngjun0827/assembly_qc/blob/main/rules/assembly_qc.workflow.png)
