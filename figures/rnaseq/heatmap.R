# Creating the heatmap of DE genes for 271.

library(edgeR)
y <- readRDS("../rnaseq/analysis/object.rds")
sig <- read.table("../rnaseq/analysis/results_lfc/cons_271_all.txt", header=TRUE, sep="\t", row.names=1)
sig <- sig[sig$adj.P.Val <= 0.05,]

to.use <- c(
    "LNA.wild_type_genotype.271_LNA.batch_3",
    "LNA.wild_type_genotype.Negative_control_A.batch_3",
    "LNA.wild_type_genotype.Negative_control_B.batch_3", 
    "RNA_interference.wild_type_genotype.271_siRNA.batch_1",
    "RNA_interference.wild_type_genotype.Ambion_control.batch_1",
    "RNA_interference.wild_type_genotype.Dharmacon_control.batch_1"
)

# Using all significant genes, plus the LINC itself.
y <- y[c(rownames(sig), "ENSG00000231711"), y$samples$group %in% to.use]
adjc <- cpm(y, log=TRUE, prior.count=3)

colnames(adjc) <- y$samples$group
adjc <- adjc[,order(match(colnames(adjc), to.use))]
rownames(adjc) <- y$genes$Symbol

# Centering so that log-expression values are log-fold changes against controls.
library(pheatmap)
in.lna <- colnames(adjc) %in% to.use[1:3]
in.rnai <- colnames(adjc) %in% to.use[4:6]
adjc[,in.lna] <- adjc[,in.lna] - rowMeans(adjc[,colnames(adjc) %in% to.use[2:3]])
adjc[,in.rnai] <- adjc[,in.rnai] - rowMeans(adjc[,colnames(adjc) %in% to.use[5:6]])

# Capping intensity to preserve visibility around zero.
LIM <- 2
adjc[adjc > LIM] <- LIM
adjc[adjc < -LIM] <- -LIM
pdf("heat_271_rna.pdf")
pheatmap(adjc, cluster_cols=FALSE, breaks=seq(-LIM, LIM, length.out=101))
dev.off()
