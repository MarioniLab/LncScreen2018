set -e
set -u

# Creating and executing the width-specific DE reports.
for w in 150 500 1000
do
    newfile=de_${w}.Rmd 
    cat ../template_de.Rmd | sed "s/#_WIDTH_#/${w}/g" \
        | sed "s/#_MARK_#/H3K4me3/g" \
        > ${newfile}
    echo "rmarkdown::render('${newfile}')" | R --no-save --slave
done

# Instantiating and executing the consolidation report.
newfile=cons.Rmd
cat ../template_cons.Rmd | sed "s/#_MARK_#/H3K4me3/g" > ${newfile}
echo "rmarkdown::render('${newfile}')" | R --no-save --slave

