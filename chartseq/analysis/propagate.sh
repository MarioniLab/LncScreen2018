set -e
set -u

# Creating and executing the width-specific DE reports.
for w in 150 500 1000
do
    newfile=de_${w}.Rmd 
    cat template_de.Rmd | sed "s/#_WIDTH_#/${w}/g" > ${newfile}
    echo "rmarkdown::render('${newfile}')" | Rdevel --no-save --slave
done

# Instantiating and executing the consolidation report.
echo "rmarkdown::render('cons.Rmd')" | Rdevel --no-save --slave

