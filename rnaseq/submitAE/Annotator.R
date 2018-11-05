########################################################################################
# Coordinating the metadata obtainer with the file names.

meta.data <- read.table("../metadata.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
md5.data <- read.table("../fastq/md5.all", stringsAsFactors=FALSE)

do.numbers <- sub("_.*", "", md5.data[,2])
m <- match(do.numbers, meta.data$Library)
stopifnot(all(!is.na(m)))

final <- meta.data[m,]
final$FileName <- md5.data[,2]
final$MD5sum <- md5.data[,1]

########################################################################################
# Constructing the sdrf.tsv file.

output <- list()
output[["Source Name"]] <- final$Library
output[["Characteristics[organism]"]] <- "Homo sapiens"
output[["Characteristics[cell line]"]] <- "HeLa"
output[["Characteristics[loss of function method]"]] <- final$LOF
output[["Characteristics[genotype]"]] <- final$Genotype
output[["Material Type"]] <- "cell"
output[[paste0(rep(c("Protocol REF", "Performer"), 6), collapse="\t")]] <- paste0(c("P-MTAB-53234", "Lovorka Stojic",
                                                                                    "P-MTAB-53235", "Lovorka Stojic",
                                                                                    "P-MTAB-53236", "Lovorka Stojic",
                                                                                    "P-MTAB-53241", "Lovorka Stojic",
                                                                                    "P-MTAB-53237", "Lovorka Stojic",
                                                                                    "P-MTAB-53238","Lovorka Stojic"
                                                                                    ), collapse="\t")
output[["Extract Name"]] <- final$Library
output[["Comment[LIBRARY_LAYOUT]"]] <- "PAIRED"
output[["Comment[LIBRARY_SELECTION]"]] <- "Inverse rRNA"
output[["Comment[LIBRARY_SOURCE]"]] <- "TRANSCRIPTOMIC"
output[["Comment[LIBRARY_STRAND]"]] <- "first strand"
output[["Comment[LIBRARY_STRATEGY]"]] <- "ssRNA-seq"
output[["Comment[NOMINAL_LENGTH]"]] <- 295
output[["Comment[NOMINAL_SDEV]"]] <- 25
output[["Comment[ORIENTATION]"]] <- "5'-3'-3'-5'"
output[["Protocol REF\tPerformer"]] <- "P-MTAB-53239\tLovorka Stojic"
output[["Assay Name"]] <- sub("_.*_", "_", sub("\\.r", "", sub("[_p]?[12].fq.gz", "", final$FileName)))
output[["Technology Type"]] <- "sequencing assay"
output[["Comment[batch number]"]] <- final$Batch
output[["Comment[sequencing date]"]] <- paste0(substr(final$Date, 1, 4), "-", substr(final$Date, 5, 6), "-", substr(final$Date, 7, 8))
output[["Comment[experiment number]"]] <- final$Experiment
output[["Array Data File"]] <- final$FileName
output[["Protocol REF"]] <- "P-MTAB-53240"
output[["Derived Array Data File"]] <- "lncRNA_counts.tsv"
output[["Comment[MD5]"]] <- final$MD5sum
output[["Factor Value[loss of function method]"]] <- final$LOF
output[["Factor Value[genotype]"]] <- final$Genotype
output[["Factor Value[compound]"]] <- sub("271", "LINC000899", final$Compound)

output$check.names <- FALSE
sdrf <- do.call(data.frame, output)
write.table(file="sdrf.tsv", sdrf, row.names=FALSE, sep="\t", quote=FALSE)
