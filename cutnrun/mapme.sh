ispet=1
fastq=$(ls fastq/* | grep ".fq.gz")
genome=../genomes/builds/hg38
source ../tools/multi_align.sh
