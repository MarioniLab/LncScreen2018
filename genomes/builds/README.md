# Sequence sources

- hg38 sequence was obtained from UCSC (via `/scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa`)

# Genome builds

This combines the hg38 build of the human genome:

```sh
subread-buildindex -o hg38 /scratchb/bioinformatics/reference_data/reference_genomes/homo_sapiens/hg38/fasta/hsa.hg38.fa
```
