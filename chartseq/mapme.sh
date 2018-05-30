ispet=1
extra="-g"
fastq=($(ls merged/*.fq.gz))
genome=../genomes/builds/hg38
source ${HOME}/Code/mapping/multi_align.sh
