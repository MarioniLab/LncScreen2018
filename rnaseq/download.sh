# Downloads processed data and SDRF files from ArrayExpress.

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7432/E-MTAB-7432.sdrf.txt
wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7432/E-MTAB-7432.processed.1.zip

unzip E-MTAB-7432.processed.1.zip
mv genic_counts.tsv analysis/
rm E-MTAB-7432.processed.1.zip
