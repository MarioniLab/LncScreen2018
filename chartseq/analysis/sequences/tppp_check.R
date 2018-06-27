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
transeq <- DNAStringSet(unlist(transeq))
writeXStringSet(transeq, file="LINC00899.fa")

# Comparing the best pairwise alignments with and without shuffling.

pdf("alignments.pdf")
for (x in names(ref)) {
    for (y in names(transeq)) {
        new.refs <- character(100)
        seq.src <- strsplit(as.character(ref[[x]]), "")[[1]]
        for (i in seq_along(new.refs)) {
            new.refs[i] <- paste(seq.src[sample(length(seq.src))], collapse="")
        }
        aln.g.f <- pairwiseAlignment(DNAStringSet(new.refs), subject=transeq[y], type="local", scoreOnly=TRUE)
        aln.g.r <- pairwiseAlignment(reverseComplement(DNAStringSet(new.refs)), subject=transeq[y], type="local", scoreOnly=TRUE)
        hist(pmax(aln.g.f, aln.g.r), xlab="Score", col="grey80", breaks=20, main=paste(x, "on", y))

        aln.f <- pairwiseAlignment(ref[x], subject=transeq[y], type="local", scoreOnly=TRUE)
        aln.r <- pairwiseAlignment(reverseComplement(ref[x]), subject=transeq[y], type="local", scoreOnly=TRUE)
        abline(v=pmax(aln.f, aln.r), col="red", lwd=2, lty=2)
    }
}
dev.off()