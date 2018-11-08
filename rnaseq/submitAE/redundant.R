already_there <- read.table("https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5308/E-MTAB-5308.sdrf.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
current <- read.table("sdrf.tsv", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

m <- match(sub("p([12]).fq.gz", "_\\1.fq.gz", current[, "Array Data File"]), already_there[,"Comment[SUBMITTED_FILE_NAME]"])
current$`Comment[Existing ENA run]` <- already_there[m,"Comment[ENA_RUN]"]

write.table(current, file="sdrf2.tsv", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
