library(GenomicRanges)
all.regions <- list(TPPP=GRanges("chr5", IRanges(659862,695000)),
                    TPPP_zoom=GRanges("chr5", IRanges(672501,676750)))

# Defining the files, conditions and colors.
bg.files <- list(`271 antisense`=c(
    "SLX-14490.iPCRtagT001.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT003.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT005.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT007.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT009.HL73MBBXX.bedgraph.gz"
), `271 sense`=c(
    "SLX-14490.iPCRtagT002.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT004.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT006.HL73MBBXX.bedgraph.gz",
    "SLX-14481.iPCRtagT008.HL73MBBXX.bedgraph.gz",
    "SLX-14490.iPCRtagT010.HL73MBBXX.bedgraph.gz"
))

colors <- c(`271 antisense`="lightblue", `271 sense`="lightgrey")
bounds <- c(TPPP=3, TPPP_zoom=3)

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
    all.tracks <- list()
    mean.tracks <- list()
    cur.bound <- bounds[[loc]]

    for (comparison in names(bg.files)) {
        CUR.FILES <- file.path("../../chartseq/analysis/bedgraph", bg.files[[comparison]])
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
    pdf(paste0(loc, "_all.pdf"), height=12, width=8)
    plotTracks(c(list(itrack, gtrack), unlist(all.tracks), list(grtrack)), from=start(REQUESTED), to=end(REQUESTED),
            transcriptAnnotation = "symbol")
    dev.off()

    pdf(paste0(loc, "_mean.pdf"), height=12, width=8)
    plotTracks(c(list(itrack, gtrack), mean.tracks, list(grtrack)), from=start(REQUESTED), to=end(REQUESTED),
            transcriptAnnotation = "symbol")
    dev.off()
}
