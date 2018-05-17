library(rtracklayer)
black <- import("http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz")

library(csaw)
param <- readParam(dedup=TRUE, minq=10, discard=black, pe="both", max.frag=1000)

input <- read.table("norm_1000.tsv", header=TRUE, stringsAsFactor=FALSE)
bam.files <- input$bam.files
eff.lib.size <- input$totals * input$norm.factors

library(Rsamtools)
all.chrs <- scanBamHeader(bam.files[1])[[1]][[1]]
all.chrs <- all.chrs[grepl("^chr[0-9XY]+$", names(all.chrs))]

dir.create("bedgraph", showWarnings=FALSE)
for (f in seq_along(bam.files)) {
    collected <- vector("list", length(all.chrs))
    for (i in seq_along(all.chrs)) {
        region <- GRanges(names(all.chrs)[i], IRanges(1, all.chrs[i]))
        out <- extractReads(bam.files[f], region=region, param=param)
        read.cov <- coverage(out)
        collected[[i]] <- read.cov / eff.lib.size[f] * 1e6
    }

    collected <- do.call(c, collected)
    new.file <- file.path("bedgraph", sub("bam", "bedgraph.gz", basename(bam.files[f])))
    new.handle <- gzfile(new.file, open="wb")
    writeLines(sprintf("track type=bedGraph name=%s", basename(bam.files[f])), con=new.handle)
    export(collected, con=new.handle, format="bedGraph")
    close(new.handle)
}

    
