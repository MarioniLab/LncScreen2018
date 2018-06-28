for mark in H3K4me3 H3K36me3 H3K27ac H3K27me3 antiRabbit
do 
    newfile=cutnrun_rnaseq_${mark}.Rmd
    cat cutnrun_rnaseq_template.Rmd | sed "s/#_MARK_#/${mark}/g" > ${newfile}
    echo "rmarkdown::render('${newfile}')" | R --slave --vanilla
done
