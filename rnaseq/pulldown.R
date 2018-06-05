# This script runs through the metadata and copies the relevant FASTQ files from tier II storage, 
# along with their MD5 sums.

dir.create("fastq")
metadata <- read.delim("metadata.tsv", header=TRUE, stringsAsFactors=FALSE)
path <- "/mnt/nas2-data/jmlab/group_folders/lun01/Odom/lncRNA_mitosis"
total.md5 <- list()

#################################################################################

for (b in c(20160208, 20160907)) {
    current <- metadata[b==metadata$Date,]
    curpath <- file.path(path, paste0("real_", b))
    fqsrc <- file.path(curpath, "fastq")
    all.fastqs <- list.files(fqsrc, pattern="(fastq|fq).gz") 

    # Figuring out which fastqs to keep.
    keep <- logical(length(all.fastqs))
    for (index in current$Library) {
        keep <- keep | grepl(tolower(index), all.fastqs)
    }
    retain.fastqs <- all.fastqs[keep]
    file.copy(file.path(fqsrc, retain.fastqs), "fastq")

    # Figuring out which MD5 sums to keep.
    md5 <- read.table(file.path(fqsrc, "md5.all"), stringsAsFactors=FALSE)
    md5 <- md5[md5[,2] %in% retain.fastqs,]
    total.md5[[as.character(b)]] <- md5
}

#################################################################################

# Repeating for the last batch, which has a different naming convention.
b <- "20161212"
current <- metadata[b==metadata$Date,]
curpath <- file.path(path, paste0("real_", b))
fqsrc <- file.path(curpath, "fastq")

all.fastqs <- list.files(fqsrc, pattern="(fastq|fq).gz") 
metasrc <- read.csv(file.path(curpath, "metadata.csv"), stringsAsFactors=FALSE)
md5 <- read.table(file.path(fqsrc, "md5.all"), stringsAsFactors=FALSE)

# Figuring out which fastqs to keep.
collected.md5 <- vector("list", nrow(current))
names(collected.md5) <- current$Library

for (index in current$Library) {
    tag <- metasrc$Indexes[metasrc$DO.numbers==index]
    tag <- sub("-", "_", tag)

    keep <- grepl(tag, all.fastqs)
    keep.fastqs <- all.fastqs[keep]
    renamed <- paste0(tolower(index), "_", keep.fastqs)
    file.copy(file.path(fqsrc, keep.fastqs), file.path("fastq", renamed))

    curmd5 <- md5[md5[,2] %in% keep.fastqs,]
    curmd5[,2] <- paste0(tolower(index), "_", curmd5[,2])
    collected.md5[[index]] <- curmd5
}

#################################################################################

# Combining all MD5 sums into a single file.
total.md5[[b]] <- do.call(rbind, collected.md5)
total.md5 <- do.call(rbind, total.md5)
write.table(file="fastq/md5.all", total.md5, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
