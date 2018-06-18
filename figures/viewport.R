library(GenomicRanges)
REQUESTED <- GRanges("chr5", IRanges(672501 - 2e4, 676750 + 2e4))

# Setting up the annotation:

CHR <- as.character(seqnames(REQUESTED))

library(Gviz)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg38", chromosome = CHR)

library(EnsDb.Hsapiens.v86)
current.tx <- exonsBy(EnsDb.Hsapiens.v86, "tx")
seqlevels(current.tx) <- paste0("chr", seqlevels(current.tx))
current.tx <- current.tx[seqnames(current.tx)==CHR]
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

# Importing the bedgraph files for CHART:

chart.files <- file.path("../chartseq/analysis/bedgraph", c(
    "SLX-14490.iPCRtagT001.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT003.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT005.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT007.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT009.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT002.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT004.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT006.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT008.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT010.HL73MBBXX.bedgraph.gz"
))

condition <- rep(c("CHART 271", "control"), each=5)
chart.label <- paste(condition, "replicate", rep(1:5, 2))
color <- rep(c("forestgreen", "goldenrod"), each=5)

library(rtracklayer) 
chart.tracks <- vector("list", length(chart.files))    
for (i in seq_along(chart.files)) { 
    cmd <- sprintf("zcat %s | grep '^%s'", chart.files[i], CHR)
    X <- import(pipe(cmd), format="bedgraph") # Restricting the input, otherwise memory usage explodes.
    seqlevels(X) <- CHR
    X$score <- pmin(X$score, 3)
    chart.tracks[[i]] <- DataTrack(range = X, genome = "hg38",  type = "histogram", chromosome = CHR, name = chart.label[i], 
            col.histogram = color[i], fill.histogram = color[i])
}

# Making the plot.

pdf("CHART_TPPP.pdf", height=12, width=8)
plotTracks(c(list(itrack, gtrack), chart.tracks, list(grtrack)), from=start(REQUESTED), to=end(REQUESTED), 
    transcriptAnnotation = "symbol")
dev.off()

 
