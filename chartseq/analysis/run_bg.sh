#!/bin/bash
#SBATCH -o bg.out
#SBATCH -e bg.err
#SBATCH -n 1
#SBATCH --mem 16000

echo 'source("../../cutnrun/analysis/bedgraph_generator.R")' | R --no-save --slave
