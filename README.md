# Screening for lncRNAs involved in mitosis

## Overview

This repository provides the code for the paper **A long noncoding RNA regulates microtubule behaviour during mitosis**
by Stojic et al.

## System requirements

Statistical analyses require [R version 3.5 or higher](https://www.r-project.org/) and several packages from the [Bioconductor project](https://bioconductor.org) (release 3.8 or higher):

```r
install.packages("BiocManager")
BiocManager::install(
    c(
        "csaw",  
        "BiocFileCache",
        "rtracklayer",
        "edgeR",
        "org.Hs.eg.db",
        "EnsDb.Hsapiens.v86",
        "BiocStyle"
    ),
    dependencies=TRUE)
```

Read alignment requires the [_subread_ aligner](http://subread.sourceforge.net/) version 1.6 or higher.
Installation of dependencies should take less than 10 minutes, depending on whether R package compilation is necessary.

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

This should take less than 10 minutes, except for step 4, which can take up to an hour.

### Realigning the sequence data 

If you want to re-align the sequence data and regenerate the count matrix:

1. Enter `genomes/builds` and follow the instructions in the README to create a `subread` index for hg38.
Also enter `genomes/annotation` and follow the instructions in the README to create a GTF file using Ensembl GRCm38 annotation.
2. Enter the `rnaseq/` subdirectory and download all FASTQ files from ArrayExpress (E-MTAB-7432) into `fastq/`.
3. Running  `mapme.sh` aligns Gzipped FASTQ files using the `subread` aligner.
4. Enter `analysis/` and run `count_me.sh` to generate the `count_me.R` script, which counts the number of reads mapped to each gene.
This yields `genic_counts.tsv` that can be used in place of the processed data for the DE analysis.

This will take several hours on a high-performance computing cluster.

## CHART-seq

### Aligning the sequencing data

Repeating the analysis requires re-aligning the sequencing data yourself.
Assuming you have already built a `subread` index for hg38:

1. Enter `chartseq/` and download all FASTQ files from ArrayExpress ([E-MTAB-7418](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7418)) into `fastq/`.
2. Rename each FASTQ file to its original name (see `Comment[SUBMITTED_FILE_NAME]` in the SDRF file).
3. Running `mapme.sh` aligns the Gzipped FASTQ files using the `subread` aligner.

This will take several hours on a high-performance computing cluster.

### Repeating the analysis

We use a multi-resolution analysis with sliding windows to detect regions that are enriched in the antisense pulldown compared to the sense control:

- Enter `chartseq/analysis` and run `propagate.sh` to perform the enrichment analyses and to consolidate the results.
- Quality control reports can be generated with `run_qc.sh`.
- BedGraph files can be generated using `run_bg.sh`.
- The sequence complementary analysis can be performed using `sequences/tppp_check.R`.

All scripts can be submitted via a SLURM job scheduler.
Each analysis will take 30 minutes to an hour on a single node.

## CUT&RUN

### Aligning the sequencing data

Repeating the analysis again requires re-aligning the sequencing data yourself.
Assuming you have already built a `subread` index for hg38:

1. Enter `cutnrun/` and download all FASTQ files from ArrayExpress ([E-MTAB-7419](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7419)) into `fastq/`.
2. Rename each FASTQ file to its original name (see `Comment[SUBMITTED_FILE_NAME]` in the SDRF file).
3. Running `mapme.sh` aligns the Gzipped FASTQ files using the `subread` aligner.

This will take several hours on a high-performance computing cluster.

### Repeating the analysis

Quality control reports can be generated with `run_qc.sh` and BedGraph files can be generated using `run_bg.sh`.
All scripts can be submitted via a SLURM job scheduler and will take 15 to 30 minutes on a single node.

For each mark, we use a multi-resolution analysis to deteect regions that are differentially bound upon lncRNA depletion.
Enter `cutnrun/analysis/<MARK>` for each mark of interest and call `propagate.sh` to perform the DB analyses and to consolidate the results.

## Proteomics

To reproduce the Poisson model analysis of the TMT data:

1. Enter `proteomics/` and create a `data` subdirectory.
2. Download raw proteomics data from [here](https://jmlab-gitlab.cruk.cam.ac.uk/publications/LncScreen2018-DataFiles) into `data` (preserve the subdirectory structure).
3. Run `analysis.Rmd` to test for significant differences upon lncRNA depletion.

This should take a few minutes.

## Integrating 'omics modalities

Once data from each modality has been analyzed, the results can be integrated across modalities using the various Rmarkdown files in `integration/`.
This focuses only on integration of RNA-seq data with each other modality.
For CUT&RUN, use `propagate.sh` to perform integration between RNA-seq and the DB results for each mark of interest.

## Other analyses

The `figures/` directory contains code required to reproduce some of the figures from available 'omics datasets.
Most of these are coverage tracks and assume that BedGraph files are available for each dataset.
The heatmap for RNA-seq data can also be regenerated.

The `miscellaneous/` directory contains scripts to generate other figures, generally for publicly available datasets (GTEx, TCGA).
