echo '
anno.files <- "../../genomes/annotation/hg38.gtf"
bam.files <- list.files("../bam", full=TRUE, pattern="bam$")
ispet <- TRUE
strandspec <- 2
' | cat - ../../tools/counter.R > count_me.R
R CMD BATCH --no-save count_me.R 

