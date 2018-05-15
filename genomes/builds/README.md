# Sequence sources

- hg38 sequence was obtained from UCSC (via `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa`)
- R64 sequence was obtained from Ensembl (via `/scratchb/bioinformatics/reference_data/reference_genomes/saccharomyces_cerevisiae/R64-1-1/fasta/sce.R64-1-1.fa`)

# Genome builds

This combines the hg38 build of the human genome with yeast spike-ins:

```sh
sbatch << EOT
#!/bin/bash
#SBATCH -o log.hg38_sc.out
#SBATCH -e log.hg38_sc.err
#SBATCH -n 1
#SBATCH --mem 16000
subread-buildindex -o hg38_sc /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa \
    /scratchb/bioinformatics/reference_data/reference_genomes/saccharomyces_cerevisiae/R64-1-1/fasta/sce.R64-1-1.fa
EOT
```
