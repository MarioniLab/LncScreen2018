# Loads in GTEx TPMs for the specified gene, and plots them by tissue.
# Data taken from https://gtexportal.org/home/datasets, and processed with:
# zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | head -3 | tail -1 > of_interest.tsv
# zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | egrep "^ENSG00000231711|ENSG00000265096|ENSG00000171368" >> of_interest.tsv

exprs <- read.delim("of_interest.tsv", check.names=FALSE)
meta <- read.table("GTEx_v7_Annotations_SampleAttributesDS.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t", fill=TRUE, quote="")
locale <- meta$SMTS[match(colnames(exprs)[-(1:2)], meta$SAMPID)]

genes <- c(LINC00899="ENSG00000231711", C1qTNF1AS1="ENSG00000265096", TPPP="ENSG00000171368")
for (x in seq_along(genes)) {
    keep <- grep(genes[x], exprs$Name)
    cur.exprs <- as.numeric(exprs[keep,-(1:2)])

    by.location <- split(cur.exprs, locale)
    by.location <- by.location[names(by.location)!=""]
    by.location <- by.location[order(-sapply(by.location, median))]

    colors <- rev(viridis::viridis(length(by.location)))
    by.rgb <- col2rgb(colors) / sqrt(2)
    darker.colors <- rgb(by.rgb[1,], by.rgb[2,], by.rgb[3,], maxColorValue=255)

    pdf(paste0(genes[x], "_GTEx.pdf"), height=6, width=8)
    par(mar=c(10.1, 4.1, 4.1, 2.1))
    boxplot(by.location, ylab=paste(names(genes)[x], "expression (TPM)"), las=2, col=colors, main=genes[x], 
            border=darker.colors, pch=16, cex=0.5)
    dev.off()
}

# Plotting TPPP against the LINC in all GTEx samples.
lnc.keep <- grep(genes[1], exprs$Name)
lnc.exprs <- as.numeric(exprs[lnc.keep,-(1:2)])
tppp.keep <- grep(genes[3], exprs$Name)
tppp.exprs <- as.numeric(exprs[tppp.keep,-(1:2)])

pdf("TPPP_vs_linc_GTEx.pdf")
plot(lnc.exprs+1, tppp.exprs+1, log="xy", xlab="LINC0899 (TPM + 1)", ylab="TPPP (TPM + 1)",
    cex=0.8, col=ifelse(locale=="Brain", "red", "grey80"), pch=16)        
dev.off()
