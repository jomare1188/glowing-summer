# SNP Variant Calling Pipeline for Population Genomics

pipeline for SNP discovery and variant calling from raw Illumina paired-end reads, optimized for downstream population genetic analyses with **pixy**.

## Overview

This pipeline performs end-to-end variant calling from raw FASTQ files to filtered VCF files containing both **variant and invariant sites**, which are essential for accurate population genetic statistics with pixy (https://pixy.readthedocs.io/en/latest/).

### Key Features

- **All-sites VCF generation**: Produces VCF files with both variant and invariant sites required by pixy
- **Parallel processing**: Processes multiple samples simultaneously with configurable resource allocation
- **Smart checkpoint system**: Automatically resumes from interruptions and skips completed steps
- **Storage optimization**: Automatic cleanup of intermediate files to conserve disk space
- **Quality control**: Built-in metrics collection at every step

## Pipeline Workflow

```
Raw FASTQ reads
    ↓
1. QC analysis (fastqc, fastp)
    ↓
2. Read Alignment (BWA-MEM)
    ↓
3. Duplicate Marking (GATK MarkDuplicates)
    ↓
4. Variant Calling per sample (GATK HaplotypeCaller with -ERC GVCF)
    ↓
5. Joint Genotyping (GATK GenotypeGVCFs) - ALL sites
    ↓
5. Variant Filtering
    ├─→ Filtered variant sites only (SNPs)  -> for other populations genomic analyses
    ├─→ Filtered invariant sites only
    └─→ Merged all-sites VCF (for pixy)

```

## Requirements

### Software Used

| Tool | Version | Purpose |
|------|---------|---------|
| **fastp** | 1.0.1 | Quality control and read trimming |
| **BWA-MEM** | 2.3 | Read alignment to reference genome |
| **SAMtools** | 1.23 | BAM file manipulation and indexing |
| **GATK** | 4.6.2.0 | SNP calling and filtering |
| **BCFtools** | 1.23 | VCF file manipulation |
| **VCFtools** | 0.1.17 | Population genetics filtering |
| **PLINK** | 2.0.0a.6.9 | LD pruning and format conversion |
| **Pixy** | 2.0.0.beta14 | pi diversity estimation |
| **GNU parallel** | 20190922 | parallel processing |

### Input Data

- **Paired-end FASTQ files**: Illumina sequencing reads (gzip compressed)
  - Naming convention: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`
- **Reference genome**: https://doi.org/10.1093/gigascience/giaa093, FASTA format (will be indexed automatically)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/snp-variant-calling-pipeline.git
cd snp-variant-calling-pipeline
```

2. Ensure all dependencies are installed and in your PATH

3. Configure the pipeline by editing the configuration section in `snp_invariants.sh`:

```bash
# Input/Output paths
READS_DIR="/path/to/your/trimmed_reads"
REF_GENOME="/path/to/reference/genome.fa"
WORK_DIR="/path/to/output/directory"

# Computational resources
TOTAL_CORES=80
TOTAL_MEM=600  # GB
MAX_PARALLEL_SAMPLES=8
THREADS_PER_SAMPLE=10
```

## Usage

### Basic Usage

```bash
bash snp_invariants.sh
```

The pipeline will:
1. Automatically discover all FASTQ pairs in the input directory
2. Process samples in parallel
3. Generate joint-called VCF files
4. Apply quality filters
5. Produce population assignment templates for pixy

### Output Structure

```
results/
├── alignment/              # Sorted BAM files
├── markdup/                # Duplicate-marked BAM files
├── variants/
│   ├── raw/                # Raw GVCF files per sample
│   ├── all_sites/          # Joint-called all-sites VCF
│   │   ├── cohort.allsites.vcf.gz
│   │   ├── populations.txt
│   │   └── sample_list.txt
│   └── filtered/           # Filtered VCF files
│       ├── cohort.snps.final.vcf.gz
│       ├── cohort.invariant.final.vcf.gz
│       └── cohort.allsites.filtered.vcf.gz
├── metrics/                # QC statistics
└── logs_parallel/          # Detailed logs for each step


## Results

### QC Analaysis 
- Raw FastQC report: [raw_multiqc_report.html](qc/qc_results/raw_fastqc/raw_multiqc_report.html)
- Trimmed FastQC report (after fastp): [trimmed_multiqc_report.html](qc/qc_results/trimmed_fastqc/trimmed_multiqc_report.html)

> ⚠️   Note: GitHub does not render HTML reports. Download the file and open it locally to view the interactive MultiQC report.


|Raw Reads     | Trimmed reads|
|--------------|--------------|
|4.74 billion  | 4.73 billion |

99.71% retained

### Alignemnt to reference genome

| Sample | Mapped (%) |
|--------|------------|
| 547    | 68.82 |
| 548    | 73.30 |
| 549    | 16.65 |
| 550    | 68.66 |
| 551    | 70.73 |
| 552    | 77.13 |
| 553    | 10.71 |
| 554    | 63.19 |
| 555    | 23.35 |
| 556    | 51.35 |
| 557    | 6.02  |
| 558    | 6.90  |
| 559    | 75.34 |
| 560    | 17.26 |
| 561    | 36.01 |
| 562    | 29.70 |
| 563    | 16.51 |
| 564    | 11.39 |

Complete alignment and mark duplicate results:
[Metrics directory](snp_calling/snp_calling_pixy/results/metrics)

### Variant calling

- Total samples processed: 18 

 
|| all joined sites (variant and invariant) |  filtered variant sites (SNPs) | filtered invariant sites |
|-|-----------------------------------------|--------------------------------|--------------------------|
|File        | cohort.allsites.vcf.gz | cohort.snps.final.vcf.gz | cohort.invariant.final.vcf.gz (pixy) |
|ts/ tv      | 2.21                   | 2.62                     | 0                                    |
|SNPs        | 7,857,765              | 645,785                  | 0                                    |
|Invariants  | 211,002,973            | 0                        | 24,636,637                           |

- Hard filters

Variant sites

--QD (Quality by Depth): QD < 2.0
--FS (Fisher Strand) > 60.0
--MQ (Mapping Quality) < 40.0
--SOR (Strand Odds Ratio) > 3.0
--MQRankSum (Mapping Quality Rank Sum) < -12.5
--ReadPosRankSum (Read Position Rank Sum) < -8.0

Invariant sites

--DP (Depth of Coverage) < 20
--DP (Depth of Coverage) > 500
--MQ (Mapping Quality) < 40.0

Popfilters

Variant sites 

--remove-indels \
--hwe 0.01 \
--maf 0.01 \
--max-missing 0.8 \
--min-meanDP 20 \
--max-meanDP 500 \

## vcf files for downstream populations analysis

LD filter parameters
WINDOW_SIZE=50
LD_THRESHOLD=0.2
STEP_SIZE=10

ld_pruned_50K.snps.vcf.gz
ld_pruned_5K.snps.vcf.gz


### Population Genetic Statistics (Pixy Output)

*Add pixy results here:*

#### Nucleotide Diversity (π)
- Population A: 
- Population B: 

### Figures

*Add your plots/figures here*

<!-- END RESULTS SECTION -->

```

### Monitoring Progress

Monitor the main log in real-time:

```bash
tail -f logs_parallel/pipeline_main.log
```



- **Reference genome**: *Callicarpa* genome assembly: John P Hamilton, Grant T Godden, Emily Lanier, Wajid Waheed Bhat, Taliesin J Kinser, Brieanne Vaillancourt, Ha    iyan Wang, Joshua C Wood, Jiming Jiang, Pamela S Soltis, Douglas E Soltis, Bjoern Hamberger, C Robin Buell, Gene    ration of a chromosome-scale genome assembly of the insect-repellent terpenoid-producing Lamiaceae species, Callicarpa americana, GigaScience, Volume 9, Issue 9, September 2020, giaa093, https://doi.org/10.1093/gigascience/giaa093  (car_asm.fa)


| Sample | Mapped (%) |
|--------|------------|
| 547    | 68.82 |
| 548    | 73.30 |
| 549    | 16.65 |
| 550    | 68.66 |
| 551    | 70.73 |
| 552    | 77.13 |
| 553    | 10.71 |
| 554    | 63.19 |
| 555    | 23.35 |
| 556    | 51.35 |
| 557    | 6.02  |
| 558    | 6.90  |
| 559    | 75.34 |
| 560    | 17.26 |
| 561    | 36.01 |
| 562    | 29.70 |
| 563    | 16.51 |
| 564    | 11.39 |


---

| Filter | Threshold | Why? |
|--------|-----------|------|
| QD (Quality by Depth) | < 2.0 | Removes variants with low quality relative to depth |
| FS (Fisher Strand) | > 60.0 | Removes variants showing strand bias |
| MQ (Mapping Quality) | < 40.0 | Removes variants in hard-to-map regions |
| SOR (Strand Odds Ratio) | > 3.0 | Another strand bias metric |
| MQRankSum | < -12.5 | Removes variants where ref/alt reads have different mapping qualities |
| ReadPosRankSum | < -8.0 | Removes variants at ends of reads (often errors) |


#### Initial Quality Filtering
```bash
--remove-indels              # Keep only SNPs (no insertions/deletions)
--min-meanDP 20              # Minimum average depth of 20x coverage
--max-meanDP 100             # Maximum average depth of 100x (avoid repeats)
--minQ 30                    # Minimum quality score of 30
--max-missing 0.8            # Keep SNPs present in ≥80% of samples
--maf 0.05                   # Minimum allele frequency of 5%
--min-alleles 2              # Exactly 2 alleles (biallelic SNPs)
--max-alleles 2
```

**Why these filters?**
- **Depth filters (20-100x)**: Too low = unreliable calls; too high = repetitive regions
- **Max missing 0.8**: Ensures SNPs are present in most samples for accurate comparisons
- **MAF 0.05**: Rare variants can be unreliable and don't contribute much to population structure
- **Biallelic only**: Simplifies analysis and avoids complex multiallelic sites

---


### Output Files

#### vcf file with all sites (variant and invariant) for diversity calculation using pixy
- `/home/dmpachon/jorge/TATIANA/snp_calling/snp_calling_pixy/results/variants/all_sites/cohort.allsites.vcf.gz`




#### vcf filtered file (only variant sites and gen pop filters)
- /home/dmpachon/jorge/TATIANA/snp_calling/snp_calling_pixy/results/variants/filtered/cohort.snps.pass.vcf.gz
- 4,759,430 SNPs


## 💻 Software Used

| Tool | Version | Purpose |
|------|---------|---------|
| **fastp** | 1.0.1 | Quality control and read trimming |
| **BWA-MEM** | 2.3 | Read alignment to reference genome |
| **SAMtools** | 1.23 | BAM file manipulation and indexing |
| **GATK** | 4.6.2.0 | SNP calling and filtering |
| **BCFtools** | 1.23 | VCF file manipulation |
| **VCFtools** | 0.1.17 | Population genetics filtering |
| **PLINK** | 2.0.0a.6.9 | LD pruning and format conversion |
| **Pixy** | 2.0.0.beta14 | pi diversity estimation |

