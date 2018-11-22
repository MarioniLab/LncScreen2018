# This plots the average mitotic delay against the average knockdown efficiency for all depletion experiments.
# Raw data were obtained from Lovorka Stojic, where each value represents the average of 4-5 replicates.
# It also plots the average TPPP upregulation in each of these experiments.

lna.identity <- c("Cells", "OLD 271, ex1", "Exon1 4910", "Intron1 4916", "Intron2 4909", "Intron3 4912")
lna.delay <- c(1.006550701, 2.7607248557, 1.2661597163, 1.5986431166, 1.0295081167, 1.0568579716) 
lna.knockdown <- c(0.94337296, 0.2500654143, 0.61413758, 0.455920475, 0.68882896, 0.84169272)
lna.tppp <- c(1.36889805, 14.4084081111, 1.04975056, 3.0008166, 1.32256392, 0.99434702)

rnai.identity <- c("271 si", "cells")
rnai.delay <- c(2.0417657873, 1.1432664602)
rnai.knockdown <- c(0.022930295, 1.1187240875)
rnai.tppp <- c(14.3089043077, 1.1103203154)

identity <- c(lna.identity, rnai.identity)
delay <- c(lna.delay, rnai.delay)
knockdown <- c(lna.knockdown, rnai.knockdown)
tppp <- c(lna.tppp, rnai.tppp)    

pdf("combined.pdf")
plot(knockdown, tppp, pch=16, xlab="Knockdown (vs control)", ylab="TPP (vs control)", log="y")
text(knockdown, tppp, identity, xpd=TRUE, col="red", pos=4)

plot(tppp, delay, pch=16, xlab="TPPP (vs control)", ylab="Mitotic delay (vs control)", log="x")
text(tppp, delay, identity, xpd=TRUE, col="red", pos=4)
dev.off()
