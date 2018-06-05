# This script collates files across all of the batches of data in this experiment.

path <- "../../lncRNA_mitosis"

########################################################################################

firstlot <- read.table(file.path(path, "real_20160208/analysis/sample_metadata.txt"), header=TRUE, sep ="\t")
firstlot <- firstlot[,c(1,7)]

lib.num <- firstlot$Library
exp.num <- sub(".*exp", "\\1", firstlot$Individual)
exp.num[exp.num=="4b"] <- "5"
condition <- sub("[ _-]*exp.*", "", firstlot$Individual)
cleaned.first <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Keeping only necessary libraries.
keep <- grepl("Ambion", condition) | grepl("Dharamaco", condition) | grepl("siRNA", condition)
keep <- keep & (!grepl("3417", condition) & !grepl("289", condition))
keep <- keep & lib.num != "do8276"

cleaned.first <- cleaned.first[keep,] 
rownames(cleaned.first) <- NULL
cleaned.first$Date <- 20160208
cleaned.first$Batch <- "20160208"

# Setting LOF mode:
LOF <- "RNA inteference"

# Setting genotype:
genotype <- "wild type genotype"

# Setting treatment compound:
compound <- character(nrow(cleaned.first))

compound[cleaned.first$Condition=="Con si Ambion"] <- "Ambion control"
compound[cleaned.first$Condition=="Control Dharamaco"] <- "Dharmacon control"
compound[cleaned.first$Condition=="271 siRNA"] <- "271 siRNA"
 
# Creating a new group.
cleaned.first$LOF <- LOF
cleaned.first$Genotype <- genotype 
cleaned.first$Compound <- compound

########################################################################################

secondlot <- read.csv(file.path(path, "real_20160907/analysis/sample_metadata.csv"), header=TRUE)

lib.num <- tolower(secondlot$Library)
exp.num <- sub(".*exp", "\\1", secondlot$Sample)
exp.num[exp.num=="4b"] <- "5"
condition <- sub("[ _-]*exp.*", "", secondlot$Sample)
cleaned.second <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Keeping only necessary libraries.
keep <- condition %in% c("NegA", "NegB", "Consi", "Condh") | grepl("C1", condition) 
keep <- keep & (!grepl("LNAold", exp.num))

cleaned.second <- cleaned.second[keep,] 
rownames(cleaned.second) <- NULL
cleaned.second$Date <- 20160907
cleaned.second$Batch <- "20160907"

# Setting LOF mode:
LOF <- character(nrow(cleaned.second))
LOF[cleaned.second$Condition %in% c("NegA", "NegB", "C1_LNA1")] <- "LNA"
LOF[cleaned.second$Condition %in% c("Consi", "Condh", "C1_si")] <- "RNA interference"

# Setting genotype:
genotype <- "wild type genotype"

# Setting treatment compound:
compound <- character(nrow(cleaned.second))

compound[cleaned.second$Condition=="NegA"] <- "Negative control A"
compound[cleaned.second$Condition=="NegB"] <- "Negative control B"
compound[cleaned.second$Condition=="C1_LNA1"] <- "C1 LNA"

compound[cleaned.second$Condition=="Consi"] <- "Ambion control"
compound[cleaned.second$Condition=="Condh"] <- "Dharmacon control"
compound[cleaned.second$Condition=="C1_si"] <- "C1 siRNA"
 
# Creating a new group.
cleaned.second$LOF <- LOF
cleaned.second$Genotype <- genotype 
cleaned.second$Compound <- compound

########################################################################################

thirdlot <- read.csv(file.path(path, "real_20161212/analysis/metadata.csv"), header=TRUE)
thirdlot <- thirdlot[,1:2]

lib.num <- tolower(thirdlot$DO.numbers)
exp.num <- sub(".*exp", "\\1", thirdlot$Description)
exp.num[exp.num=="4b"] <- "5"
condition <- sub("[ _-]*exp.*", "", thirdlot$Description)
cleaned.third <- data.frame(Library=lib.num, Condition=condition, Experiment=exp.num) 

# Keeping only necessary libraries.
keep <- grepl("^271", condition) | grepl("^Neg", condition)

cleaned.third <- cleaned.third[keep,] 
rownames(cleaned.third) <- NULL
cleaned.third$Date <- 20161212
cleaned.third$Batch <- "20161212"

# Setting LOF mode:
LOF <- "LNA"

# Setting genotype:
genotype <- "wild type genotype"

# Setting treatment compound:
compound <- character(nrow(cleaned.third))

compound[cleaned.third$Condition=="NegA"] <- "Negative control A"
compound[cleaned.third$Condition=="NegB"] <- "Negative control B"
compound[cleaned.third$Condition=="271 exon1_LNA2"] <- "271 LNA2"
compound[cleaned.third$Condition=="271 exon1_LNA3"] <- "271 LNA3"
 
# Creating a new group.
cleaned.third$LOF <- LOF
cleaned.third$Genotype <- genotype 
cleaned.third$Compound <- compound

########################################################################################

combined <- rbind(cleaned.first, cleaned.second, cleaned.third)
combined$Batch <- as.integer(factor(combined$Batch))
new.group <- paste0(combined$LOF, ".", combined$Genotype, ".", combined$Compound, ".batch_", combined$Batch)
new.group <- gsub("[- ]", "_", new.group)
combined$Condition <- new.group

write.table(combined, file="metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
