cd fastq/
allprefix=$(ls | grep "fq.gz$" | sed "s/.s_[0-9].r_[12].fq.gz$//g" | sort | uniq)
for thing in ${allfiles}
do 
    curfiles1=$(ls | grep "_1\\.fq.gz$" | grep "${allprefix}")
    curfiles2=$(ls | grep "_2\\.fq.gz$" | grep "${allprefix}")
    cat ${curfiles1} > ${allprefix}_1.fq.gz
    rm ${curfiles1}
    cat ${curfiles2} > ${allprefix}_2.fq.gz
    rm ${curfiles2}
done

