# This script pulls out the sequence of all TPPP-related binding sites.

res <- read.table("../db_regions.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
res <- res[res$FDR <= 0.3,]

library(GenomicRanges)
keep <- grepl("TPPP", res$overlap) | grepl("TPPP", res$left) | grepl("TPPP", res$right) 
chosen <- with(res[keep,], GRanges(seqnames, IRanges(start, end)))

library(BSgenome.Hsapiens.UCSC.hg38)
ref <- getSeq(BSgenome.Hsapiens.UCSC.hg38, chosen)
names(ref) <- as.character(chosen)
writeXStringSet(ref, file="TPPP.fa")

# We also pull out the sequence of the lncRNA.

library(EnsDb.Hsapiens.v86)
all.tx <- select(EnsDb.Hsapiens.v86, keytype="GENEID", keys="ENSG00000231711", columns="TXID")
trans <- exonsBy(EnsDb.Hsapiens.v86, by="tx")[all.tx$TXID]
seqlevels(trans) <- paste0("chr", seqlevels(trans))
genome(trans) <- "hg38"

transeq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, trans)
transeq <- lapply(transeq, paste, collapse="")
writeXStringSet(DNAStringSet(unlist(transeq)), file="LINC00899.fa")
