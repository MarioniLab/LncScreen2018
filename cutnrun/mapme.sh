ispet=1
extra="-g"
fastq=$(ls fastq/* | grep ".fq.gz")
genome=../genomes/builds/hg38_sc
source ../tools/multi_align.sh
