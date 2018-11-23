##########################################

library(recount)

# download_study("TCGA")
load("TCGA/rse_gene.Rdata")

se <- rse_gene
libsize <- colSums(assay(se, withDimnames=FALSE))

##########################################
# Looking at all the cancers of the world and plotting relative change in expression.

library(edgeR)
library(survival)    

collected.diff <- collected.surv <- list() 

for (cancer in unique(se$cgc_file_disease_type)) {
    typekeep <- which(se$cgc_file_disease_type==cancer)
    disease <- (se$gdc_cases.samples.sample_type=="Primary Tumor")[typekeep]
    normal <- (se$gdc_cases.samples.sample_type=="Solid Tissue Normal")[typekeep]
    if (!any(disease) || !any(normal)) { 
        next
    }

    for (gene in c("ENSG00000231711.2", "ENSG00000265096.1", "ENSG00000171368.11")) {
        current.expr <- assay(se, withDimnames=FALSE)[gene,typekeep,drop=FALSE]
        current.expr <- cpm(current.expr, lib.size=libsize[typekeep], log=TRUE, prior.count=5)
        ylab <- paste("Log2-CPM for", gene)

        pdf(paste0(gsub(" ", "_", cancer), "-", gene, ".pdf"))
        all.exprs <- list(Normal=current.expr[normal], Disease=current.expr[disease])
        boxplot(all.exprs, ylab=ylab, main=cancer)
        if (all(lengths(all.exprs) >= 2L)) {
            tested <- t.test(all.exprs[[1]], all.exprs[[2]])
            collected.diff <- c(collected.diff, list(data.frame(Cancer=cancer, Gene=gene, P=tested$p.value)))
        }

        disease.expr <- current.expr[disease]
        status <- se$gdc_cases.diagnoses.vital_status[typekeep][disease]
        d2d <-se$gdc_cases.diagnoses.days_to_death[typekeep][disease]
        sex <- se$gdc_cases.demographic.gender[typekeep][disease]
        age <- se$gdc_cases.diagnoses.age_at_diagnosis[typekeep][disease]

        # Set the censoring time for living people.
        is.live <- which(status=="alive")
        d2d[is.live] <- se$gdc_cases.diagnoses.days_to_last_follow_up[typekeep][disease][is.live]

        # Specify expression boundaries.
        labels <- c("Low", "Intermediate", "High")
        thresholds <- quantile(disease.expr, c(0.25, 0.75))
        expr.level <- factor(labels[findInterval(disease.expr, thresholds)+1], levels=labels)

        km_fit <- survfit(Surv(d2d, status=="dead") ~ expr.level)
        colors <- c("red","salmon", "orange")
        plot(km_fit, xlab="Time", ylab="Survival", col=colors, lwd=2, main=cancer)
        legend("bottomleft", col=colors, lwd=2, legend=paste(labels,
            c(sprintf("[%.1f, %.1f]", min(disease.expr), thresholds[1]), 
                sprintf("(%.1f, %.1f]", thresholds[1], thresholds[2]), 
                sprintf("(%.1f, %.1f]", thresholds[2], max(disease.expr))),
            paste0("(n = ", table(expr.level), ")")
            )
        )

        Z.expr <- expr.level # qnorm((rank(disease.expr)-0.5)/length(disease.expr))
        if (length(unique(sex))>1L) {
            cx_fit <- coxph(Surv(d2d, status=="dead") ~ Z.expr + sex + age)
            cx_fit0 <- coxph(Surv(d2d, status=="dead") ~ sex + age)
        } else {
            cx_fit <- coxph(Surv(d2d, status=="dead") ~ Z.expr + age)
            cx_fit0 <- coxph(Surv(d2d, status=="dead") ~ age)
        }
        tested <- anova(cx_fit, cx_fit0)
        collected.surv <- c(collected.surv, list(data.frame(Cancer=cancer, Gene=gene, P=tested[["P(>|Chi|)"]][2])))
        dev.off()
    }
}

##########################################
# Examining the statistics.

tab.diff <- do.call(rbind, collected.diff)
tab.diff$adj.P.Val <- p.adjust(tab.diff$P, method="holm")
tab.diff <- tab.diff[order(tab.diff$P),]
write.table(tab.diff, file="tcga_diff.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

tab.surv <- do.call(rbind, collected.surv)
tab.surv$adj.P.Val <- p.adjust(tab.surv$P, method="holm")
tab.surv <- tab.surv[order(tab.surv$P),]
write.table(tab.surv, file="tcga_surv.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
