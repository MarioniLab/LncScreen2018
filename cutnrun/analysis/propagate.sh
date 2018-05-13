if [[ ! -e h3k4me3 ]]
then
    mkdir h3k4me3
fi

for w in 150 500 1000
do 
    cat template_de.Rmd | sed "s/#_WIDTH_#/${w}/g" > h3k4me3/de_${w}.Rmd 
done

