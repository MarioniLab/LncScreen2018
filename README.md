# Screening for lncRNAs involved in mitosis

## Overview

This repository provides the code for the paper **TITLE YET TO BE DESCRIBED** by Stojic et al.

## RNA-seq

### Generating the results

To repeat the differential expression analysis, enter `rnaseq/` and follow these instructions:

1. Run `download.sh`, which will download the SDRF file from ArrayExpress using the accession code [E-MTAB-7432](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7432).
It will also download the processed data and place it in the relevant subdirectories.
2. Enter `analysis/` and run the following to compile the DE analysis report.

    ```r
    .use.local <- FALSE
    rmarkdown::render("de_analysis.Rmd")
    ``` 

3. Run `further_analysis.Rmd` to perform gene set tests.
4. Run `run_bg.R` (or `bedgraph_generator.R` directly) to create Bedgraph files.

### Realigning the sequence data 

If you want to re-align the sequence data and regenerate the count matrix:

1. Enter `genomes/builds` and follow the instructions in the README to create a `subread` index for hg38.
Also enter `genomes/annotation` and follow the instructions in the README to create a GTF file using Ensembl GRCm38 annotation.
2. Enter the `rnaseq/` subdirectory and download all FASTQ files from ArrayExpress (E-MTAB-7432) into `fastq/`.
3. Running  `mapme.sh` aligns Gzipped FASTQ files using the `subread` aligner.
4. Enter `analysis/` and run `count_me.sh` to generate the `count_me.R` script, which counts the number of reads mapped to each gene.
This yields `genic_counts.tsv` that can be used in place of the processed data for the DE analysis.

## CHART-seq

### Aligning the seuqencing data

Repeating the analysis requires re-aligning the sequencing data yourself.
Assuming you have already built a `subread` index for hg38:

1. Enter `chartseq/` and download all FASTQ files from ArrayExpress ([E-MTAB-7418](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7418)) into `fastq/`.
2. Rename each FASTQ file to its original name (see `Comment[SUBMITTED_FILE_NAME]` in the SDRF file).
3. Running `mapme.sh` aligns the Gzipped FASTQ files using the `subread` aligner.

### Repeating the analysis

We use a multi-resolution analysis with sliding windows to detect regions that are enriched in the antisense pulldown compared to the sense control:

- Enter `chartseq/analysis` and run `propagate.sh` to perform the DE analyses and to consolidate the results.
- Quality control reports can be generated with `run_qc.sh`.
- BedGraph files can be generated using `run_bg.sh`.

All scripts can be submitted via a SLURM job scheduler.
