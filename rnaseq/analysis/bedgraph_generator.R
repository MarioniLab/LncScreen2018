library(edgeR)
input <- readRDS("object.rds") 
eff.lib.size <- input$samples$lib.size * input$samples$norm.factors
headers <- colnames(input)

library(csaw)
library(Rsamtools)
library(rtracklayer)
dir.create("bedgraph", showWarnings=FALSE)

for (f in seq_along(headers)) {
    all.bams <- list.files("../bam", pattern=paste0(headers[f], "_.*\\.bam$"), full=TRUE)
    all.chrs <- scanBamHeader(all.bams[1])[[1]][[1]]
    all.chrs <- all.chrs[grepl("^chr[0-9XYM]+$", names(all.chrs))]

    collected.F <- collected.R <- vector("list", length(all.chrs))
    for (i in seq_along(all.chrs)) {
        region <- GRanges(names(all.chrs)[i], IRanges(1, all.chrs[i]))
        read.cov.F <- read.cov.R <- 0L

        for (bf in all.bams) {
            # Using only 'second' for strand specificity. TruSeq's strand-specific protocol
            # is such that the second read is the one on the original strand of the mRNA.
            out <- extractReads(bf, region=region, param=readParam(minq=10, pe="second")) 
            read.cov.F <- read.cov.F + coverage(out[strand(out)=="+"])
            read.cov.R <- read.cov.R + coverage(out[strand(out)=="-"])
        }

        collected.F[[i]] <- read.cov.F / eff.lib.size[f] * 1e6
        collected.R[[i]] <- read.cov.R / eff.lib.size[f] * 1e6
    }

    # Compiling forward/reverse strand coverage.
    collected.F <- do.call(c, collected.F)
    names(collected.F) <- names(all.chrs)
    collected.R <- do.call(c, collected.R)
    names(collected.R) <- names(all.chrs)

    # Saving to file.
    new.name <- paste0(headers[f], "_", input$samples$group[f])
    new.file <- file.path("bedgraph", paste0(new.name, "_F.bedgraph.gz"))
    new.handle <- gzfile(new.file, open="wb")
    writeLines(sprintf("track type=bedGraph name=%s", new.name), con=new.handle)
    export(collected.F, con=new.handle, format="bedGraph")
    close(new.handle)

    new.file <- file.path("bedgraph", paste0(new.name, "_R.bedgraph.gz"))
    new.handle <- gzfile(new.file, open="wb")
    writeLines(sprintf("track type=bedGraph name=%s", new.name), con=new.handle)
    export(collected.R, con=new.handle, format="bedGraph")
    close(new.handle)
}

    
