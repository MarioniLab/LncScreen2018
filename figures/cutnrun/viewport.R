library(GenomicRanges)
all.regions <- list(TPPP=GRanges("chr5", IRanges(659862,695000)),
                    TPPP_zoom=GRanges("chr5", IRanges(691616, 694968)),
                    `271`=GRanges("chr22", IRanges(46041009, 46044853)))

# This file selection choice is _extremely_ fragile, and depends on the ordering of bedgraph file names.
# Fortunately, in this case, the ordering is consistent across experiments with different marks.
all.choices <- list(`RNAi 271`=c(1,5), `RNAi control`=c(2,6),
    `LNA 271`=c(3,7), `LNA control`=c(4, 8))

colors <- list(`RNAi 271`="lightblue", `RNAi control`="lightgrey",
    `LNA 271`="lightblue", `LNA control`="lightgrey")

bounds <- list(
    TPPP=c(H3K4me3=3, H3K27ac=5, H3K27me3=2, H3K36me3=2, antiRabbit=2),
    TPPP_zoom=c(H3K4me3=3, H3K27ac=5, H3K27me3=2, H3K36me3=2, antiRabbit=2),
    `271`=c(H3K4me3=5, H3K27ac=3, H3K27me3=2, H3K36me3=2, antiRabbit=2)
)

# Pre-work.
library(EnsDb.Hsapiens.v86)
all.tx <- exonsBy(EnsDb.Hsapiens.v86, "tx")
seqlevels(all.tx) <- paste0("chr", seqlevels(all.tx))

library(Gviz)
gtrack <- GenomeAxisTrack()
library(rtracklayer) 

# Iterating across all specified loci.    
for (loc in names(all.regions)) {
    REQUESTED <- all.regions[[loc]]
    CHR <- as.character(seqnames(REQUESTED))
    FORWARD <- as.logical(strand(REQUESTED)=="+")

    # Setting up the annotation:
    itrack <- IdeogramTrack(genome = "hg38", chromosome = CHR)
    current.tx <- all.tx[seqnames(all.tx)==CHR]
    current.tx <- current.tx[lengths(current.tx) > 0]

    nexons <- lengths(current.tx)
    txid <- names(current.tx)    
    anno <- select(EnsDb.Hsapiens.v86, keys=txid, keytype="TXID", columns=c("GENEID", "SYMBOL"))

    all.exons <- unlist(current.tx)
    all.exons$feature <- rep(anno$GENEID, nexons)
    all.exons$transcript <- rep(txid, nexons)
    all.exons$symbol <- rep(anno$SYMBOL, nexons)
    all.exons$exon <- all.exons$exon_id
    all.exons$exon_id <- NULL
    all.exons$exon_rank <- NULL
    
    grtrack <- GeneRegionTrack(all.exons, genome = "hg38", chromosome = CHR, name = "Gene Model")
    
    # Importing the bedgraph files on the current chromosome.
    for (mark in c("H3K4me3", "H3K36me3", "H3K27ac", "H3K27me3", "antiRabbit")) { 
        bg.files <- list.files(sprintf("../../cutnrun/analysis/%s/bedgraph", mark), full=TRUE, pattern=".bedgraph.gz$")
        bg.files <- sort(bg.files)
        stopifnot(length(bg.files)==8)

        all.tracks <- list()
        mean.tracks <- list()
        cur.bound <- bounds[[loc]][mark]

        for (comparison in names(all.choices)) {
            CUR.FILES <- bg.files[all.choices[[comparison]]]
            cur.tracks <- vector("list", length(CUR.FILES))
            cur.color <- colors[[comparison]]
                
            mean.cov <- 0
            for (i in seq_along(CUR.FILES)) {
                cmd <- sprintf("zcat %s | grep '^%s'", CUR.FILES[i], CHR)
                X <- import(pipe(cmd), format="bedgraph") # Restricting the input, otherwise memory usage explodes.
                seqlevels(X) <- CHR
    
                X0 <- restrict(X, start(REQUESTED), end(REQUESTED))
                cur.tracks[[i]] <- DataTrack(range = X0, genome = "hg38",  
                    type = "histogram", chromosome = CHR, name = paste(comparison, i),
                    col.histogram = cur.color, fill.histogram = cur.color, ylim=c(0, cur.bound))
                
                mean.cov <- mean.cov + Rle(values=X$score, lengths=width(X))
            }
    
            all.tracks[[comparison]] <- cur.tracks
            mean.cov <- mean.cov/length(CUR.FILES)
            ends <- cumsum(runLength(mean.cov))
            mean.range <- GRanges(CHR, IRanges(ends - runLength(mean.cov), ends))
            mean.range$score <- runValue(mean.cov)
    
            mean0 <- restrict(mean.range, start(REQUESTED), end(REQUESTED))
            mean.tracks[[comparison]] <- DataTrack(range = mean0, genome = "hg38",
                type = "histogram", chromosome = CHR, name = comparison,
                col.histogram = cur.color, fill.histogram = cur.color, ylim=c(0, cur.bound))
        }
    
        # Making the plot.
        pdf(paste0(loc, "_", mark, "_all.pdf"), height=12, width=8)
        plotTracks(c(list(itrack, gtrack), unlist(all.tracks), list(grtrack)), from=start(REQUESTED), to=end(REQUESTED), 
                transcriptAnnotation = "symbol")
        dev.off()

        pdf(paste0(loc, "_", mark, "_mean.pdf"), height=12, width=8)
        plotTracks(c(list(itrack, gtrack), mean.tracks, list(grtrack)), from=start(REQUESTED), to=end(REQUESTED), 
                transcriptAnnotation = "symbol")
        dev.off()
    }
}
