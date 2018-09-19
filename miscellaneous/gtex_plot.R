# Loads in GTEx TPMs for the specified gene, and plots them by tissue.
# Data taken from https://gtexportal.org/home/datasets, and processed with:
# zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | head -3 | tail -1 > of_interest.tsv
# zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | egrep "^ENSG00000231711|ENSG00000265096|ENSG00000171368" >> of_interest.tsv

exprs <- read.delim("of_interest.tsv", check.names=FALSE)
meta <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

genes <- c(LINC00899="ENSG00000231711", C1qTNF1AS1="ENSG00000265096", TPPP="ENSG00000171368")
for (x in seq_along(genes)) {
    keep <- grep(genes[x], exprs$Name)
    lnc.exprs <- as.numeric(exprs[keep,-(1:2)])

    locale <- meta$SMTS[match(colnames(exprs)[-(1:2)], meta$SAMPID)]
    by.location <- split(lnc.exprs, locale)
    by.location <- by.location[names(by.location)!=""]
    by.location <- by.location[order(-sapply(by.location, median))]

    pdf(paste0(genes[x], "_GTEx.pdf"), height=6, width=6)
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    boxplot(by.location, ylab=paste(names(genes)[x], "expression (TPM)"), las=2, col=rev(viridis::viridis(length(by.location))), main=genes[x])
    dev.off()
}
