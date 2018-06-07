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

    collected <- vector("list", length(all.chrs))
    for (i in seq_along(all.chrs)) {
        region <- GRanges(names(all.chrs)[i], IRanges(1, all.chrs[i]))
        read.cov <- 0L
        for (bf in all.bams) {
            out <- extractReads(bf, region=region, param=readParam(minq=10))
            read.cov <- read.cov + coverage(out)
        }
        collected[[i]] <- read.cov / eff.lib.size[f] * 1e6
    }
    collected <- do.call(c, collected)
    names(collected) <- names(all.chrs)

    new.name <- paste0(headers[f], "_", input$samples$group[f])
    new.file <- file.path("bedgraph", paste0(new.name, ".bedgraph.gz"))
    new.handle <- gzfile(new.file, open="wb")
    writeLines(sprintf("track type=bedGraph name=%s", new.name), con=new.handle)
    export(collected, con=new.handle, format="bedGraph")
    close(new.handle)
}

    
