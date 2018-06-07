#!/bin/bash
#SBATCH -o bg.log
#SBATCH -e bg.err
#SBATCH -n 1
#SBATCH --mem 16000

echo 'source("bedgraph_generator.R")' | R --no-save --slave
