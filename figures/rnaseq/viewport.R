library(GenomicRanges)
all.regions <- list(TPPP=GRanges("chr5", IRanges(659862,693395)),
                    `271`=GRanges("chr22", IRanges(46041009, 46044853)))

all.files <- list(RNAi=list(`271`=c(
            "do8280_RNA_interference.wild_type_genotype.271_siRNA.batch_1",
            "do8288_RNA_interference.wild_type_genotype.271_siRNA.batch_1",
            "do8308_RNA_interference.wild_type_genotype.271_siRNA.batch_1",
            "do8316_RNA_interference.wild_type_genotype.271_siRNA.batch_1"
        ), Ambion=c(
            "do8271_RNA_interference.wild_type_genotype.Ambion_control.batch_1",
            "do8303_RNA_interference.wild_type_genotype.Ambion_control.batch_1",
            "do8307_RNA_interference.wild_type_genotype.Ambion_control.batch_1",
            "do8310_RNA_interference.wild_type_genotype.Ambion_control.batch_1"
        ), Dharmacon=c(
            "do8289_RNA_interference.wild_type_genotype.Dharmacon_control.batch_1",
            "do8292_RNA_interference.wild_type_genotype.Dharmacon_control.batch_1",
            "do8297_RNA_interference.wild_type_genotype.Dharmacon_control.batch_1"
        )
    ), LNA=list(`271`=c(
            "do12583_LNA.wild_type_genotype.271_LNA.batch_3",
            "do12594_LNA.wild_type_genotype.271_LNA.batch_3",
            "do12646_LNA.wild_type_genotype.271_LNA.batch_3",
            "do12652_LNA.wild_type_genotype.271_LNA.batch_3"
        ), `Negative A`=c(
            "do12587_LNA.wild_type_genotype.Negative_control_A.batch_3",
            "do12615_LNA.wild_type_genotype.Negative_control_A.batch_3",
            "do12632_LNA.wild_type_genotype.Negative_control_A.batch_3",
            "do12643_LNA.wild_type_genotype.Negative_control_A.batch_3"
        ), `Negative B`=c(
            "do12575_LNA.wild_type_genotype.Negative_control_B.batch_3",
            "do12609_LNA.wild_type_genotype.Negative_control_B.batch_3",
            "do12619_LNA.wild_type_genotype.Negative_control_B.batch_3",
            "do12649_LNA.wild_type_genotype.Negative_control_B.batch_3"
        )
    )
)

colors <- list(RNAi=c(`271`="red", Ambion="salmon", Dharmacon="violet"),
    LNA=c(`271`="blue", `Negative A`="dodgerblue", `Negative B`="lightblue"))

bounds <- list(TPPP=c(RNAi=1, LNA=1),
    `271`=c(RNAi=1, LNA=1))

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
    for (comparison in names(all.files)) {
        FILES <- all.files[[comparison]]
        all.tracks <- list()
        mean.tracks <- list()
        cur.bound <- bounds[[loc]][comparison]

        for (f in names(FILES)) {
            CUR.FILES <- FILES[[f]]
            cur.tracks <- vector("list", length(CUR.FILES))
            cur.color <- colors[[comparison]][f]
            
            mean.cov <- 0
            for (i in seq_along(CUR.FILES)) {
                fpath <- file.path("../../rnaseq/analysis/bedgraph", paste0(CUR.FILES[i], ifelse(FORWARD, "_R", "_F"), ".bedgraph.gz"))
                cmd <- sprintf("zcat %s | grep '^%s'", fpath, CHR)
                X <- import(pipe(cmd), format="bedgraph") # Restricting the input, otherwise memory usage explodes.
                seqlevels(X) <- CHR

                X0 <- restrict(X, start(REQUESTED), end(REQUESTED))
                cur.tracks[[i]] <- DataTrack(range = X0, genome = "hg38",  
                    type = "histogram", chromosome = CHR, name = paste(f, i),
                    col.histogram = cur.color, fill.histogram = cur.color, ylim=c(0, cur.bound))
                
                mean.cov <- mean.cov + Rle(values=X$score, lengths=width(X))
            }

            all.tracks[[f]] <- cur.tracks
            mean.cov <- mean.cov/length(CUR.FILES)
            ends <- cumsum(runLength(mean.cov))
            mean.range <- GRanges(CHR, IRanges(ends - runLength(mean.cov), ends))
            mean.range$score <- runValue(mean.cov)

            mean0 <- restrict(mean.range, start(REQUESTED), end(REQUESTED))
            mean.tracks[[f]] <- DataTrack(range = mean0, genome = "hg38",
                type = "histogram", chromosome = CHR, name = f,
                col.histogram = cur.color, fill.histogram = cur.color, ylim=c(0, cur.bound))
        }

        # Making the plot.
        pdf(paste0(loc, "_", comparison, "_all.pdf"), height=12, width=8)
        plotTracks(c(list(itrack, gtrack), unlist(all.tracks), list(grtrack)), from=start(REQUESTED), to=end(REQUESTED), 
                transcriptAnnotation = "symbol")
        dev.off()

        pdf(paste0(loc, "_", comparison, "_mean.pdf"), height=12, width=8)
        plotTracks(c(list(itrack, gtrack), mean.tracks, list(grtrack)), from=start(REQUESTED), to=end(REQUESTED), 
                transcriptAnnotation = "symbol")
        dev.off()
    }
}
