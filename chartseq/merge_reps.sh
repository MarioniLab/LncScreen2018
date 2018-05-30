#!/bin/bash
#SBATCH -e merge.log
#SBATCH -o merge.out
#SBATCH -n 1
#SBATCH --mem 2000

cd fastq/
allprefix=$(ls | grep "fq.gz$" | sed "s/.s_[0-9].r_[12].fq.gz$//g" | sort | uniq)
for thing in ${allprefix}
do 
    curfiles1=$(ls | grep "_1\\.fq.gz$" | grep "${thing}")
    curfiles2=$(ls | grep "_2\\.fq.gz$" | grep "${thing}")
    cat ${curfiles1} > ${thing}_1.fq.gz
    rm ${curfiles1}
    cat ${curfiles2} > ${thing}_2.fq.gz
    rm ${curfiles2}
done

