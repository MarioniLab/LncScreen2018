#!/bin/bash
#SBATCH -o log.out
#SBATCH -e log.err
#SBATCH -n 1
#SBATCH --mem 16000

set -e
set -u

# Creating and executing the width-specific DE reports.
for w in 500 1000 2000
do
    newfile=de_${w}.Rmd 
    cat ../template_de.Rmd | sed "s/#_WIDTH_#/${w}/g" \
        | sed "s/#_MARK_#/H3K36me3/g" \
        | sed 's/#_FILES_#/"SLX-15763.iPCRtagT004.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT009.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT032.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT043.CC6ULANXX.s_7.r.bam",\n    "SLX-15764.D708_D504.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D710_D501.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D711_D502.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D712_D503.CC6ULANXX.s_8.r.bam"/' \
        > ${newfile}
    echo "rmarkdown::render('${newfile}')" | Rdevel --no-save --slave
done

# Instantiating and executing the consolidation report.
newfile=cons.Rmd
cat ../template_cons.Rmd | sed "s/#_MARK_#/H3K36me3/g" > ${newfile}
echo "rmarkdown::render('${newfile}')" | Rdevel --no-save --slave

