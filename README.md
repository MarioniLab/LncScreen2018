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

