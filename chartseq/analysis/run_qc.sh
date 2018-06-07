#!/bin/bash
#SBATCH -o qc.out
#SBATCH -e qc.err
#SBATCH -n 1
#SBATCH --mem 16000

echo 'rmarkdown::render("data_qc.Rmd")' | R --no-save --slave
