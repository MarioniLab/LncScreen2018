set -e
set -u

# Creating and executing the width-specific DE reports.
for w in 150 500 1000
do
    newfile=de_${w}.Rmd 
    cat ../template_de.Rmd | sed "s/#_WIDTH_#/${w}/g" \
        | sed "s/#_MARK_#/antiRabbit/g" \
        | sed 's/#_FILES_#/"SLX-15763.iPCRtagT005.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT010.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT033.CC6ULANXX.s_7.r.bam",\n    "SLX-15763.iPCRtagT044.CC6ULANXX.s_7.r.bam",\n    "SLX-15764.D709_D501.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D710_D502.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D711_D503.CC6ULANXX.s_8.r.bam",\n    "SLX-15764.D712_D504.CC6ULANXX.s_8.r.bam"/' \
        > ${newfile}
    echo "rmarkdown::render('${newfile}')" | R --no-save --slave
done

# Instantiating and executing the consolidation report.
newfile=cons.Rmd
cat ../template_cons.Rmd | sed "s/#_MARK_#/antiRabbit/g" > ${newfile}
echo "rmarkdown::render('${newfile}')" | R --no-save --slave

