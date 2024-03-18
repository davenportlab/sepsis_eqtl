################################################################################
#
# 7. Plot individual interaction QTL against SRSq
#
################################################################################

library(lme4)
library(interactions)
library(jtools)
library(ggplot2)

load("eqtl_int_files.rda")
exp <- t(exp)

srs.int <- readRDS("srs_int_results_cond_incl.rds")
srs.int$fdr <- p.adjust(srs.int$Interaction_pval, method="fdr")
srs.sig <- subset(srs.int, fdr<0.05)
srs.sig <- srs.sig[order(srs.sig$fdr), ]

gene.info <- read.delim("gene_info_20412.txt")
srs.sig$Symbol <- gene.info$gene_name[match(srs.sig$Gene, gene.info$gene_id)]

sepsis.specific <- read.delim("sepsis_specificity_categories.txt")
sepsis.en <- subset(sepsis.specific, shared == "Bigger effect in GAinS")

# Subset to eGenes that are enhanced in sepsis and SRS1
srs.sig <- subset(srs.sig, Symbol %in% sepsis.en$gene_name)
srs.sig <- subset(srs.sig, sign(Interaction_beta)*sign(eQTL_beta) == 1)

srs.info <- read.delim("full-gains-SRS-predictions_mNN-RF.tsv")

pdf("Fig_2B_SRSq_interaction_plots.pdf", useDingbats = F, onefile = T, height=5, width=6)
for(g in 1:nrow(srs.sig)){
  gene <- as.character(srs.sig[g, "Gene"])
  gene.name <- as.character(srs.sig[g, "Symbol"])
  snp <- as.character(srs.sig[g, "SNP"])
  
  # reduce data to non-missing genotypes for plotting
  expr <- exp[complete.cases(geno.int[, snp]), gene]
  covs.red <- cbind(covs[complete.cases(geno.int[, snp]), 1:11], 
                    peer.factors[complete.cases(geno.int[, snp]), 1:20]) # 31 cols
  geno.snp <- geno.int[complete.cases(geno.int[, snp]), snp]
  srs <- covs[complete.cases(geno.int[, snp]), 12]
  gainsid.red <- GAinSID[complete.cases(geno.int[, snp])]
  
  df <- data.frame(expr, covs.red, geno.snp, srs, gainsid.red)
  
  lmm.model <- lmer(expr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10 + 
                      X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + 
                      X20 + X21 + X22 + X23 + X24 + X25 + X26 + X27 + X28 + 
                      X29 + X30 + X31 + geno.snp*srs + (1|gainsid.red),
                    REML=FALSE, data=df)
  gg <- interact_plot(lmm.model, pred = geno.snp, modx = srs, 
                      plot.points = TRUE, jitter = c(0.1, 0), 
                      colors=c("red", "dodgerblue"), centered = "none",
                      partial.residuals = T)
  
  # extract adjusted expression data
  tmp <- ggplot_build(gg)
  tmp <- tmp$data
  tmp <- tmp[[2]]
  tmp$geno <- as.factor(df$geno.snp)
  exp.data <- data.frame("Genotype"=tmp$geno,
                         "Exprn"=df$expr,
                         "AdjustedExprn"=tmp$y,
                         "SRS"=tmp$group)
  rownames(exp.data) <- rownames(df)
  # Add in SRSq
  exp.data$SRSq <- srs.info$SRSq[match(rownames(exp.data), srs.info$Sample_id)]
  
  # Fit lmm with SRSq
  lmm <- lmer(AdjustedExprn ~ Genotype*SRSq + (1|gainsid.red),
              REML=FALSE, data=exp.data)
  lmm <- summary(lmm)
  lmm <- lmm$coefficients
  # calculate lines from coefficients
  effect.lines <- data.frame("Genotype"=c(0,1,2),
                             "x1"=c(0, 0, 0),
                             "xend"=c(1, 1, 1),
                             "y1"=c(lmm[1, 1], lmm[1, 1] + lmm[2, 1], lmm[1, 1] + lmm[3, 1]),
                             "yend"=c(lmm[1, 1] + lmm[4, 1], 
                                      lmm[1, 1] + lmm[2, 1] + lmm[4, 1] + lmm[5, 1], 
                                      lmm[1, 1] + lmm[3, 1] + lmm[4, 1] + lmm[6, 1]))
  
  ggplot(exp.data, aes(SRSq, AdjustedExprn)) +
    geom_point(aes(colour=Genotype)) +
    theme_bw() + 
    ylab(label=paste(gene.name, "adjusted expression")) +
    geom_segment(data=effect.lines, 
                 aes(x=x1, y=y1, xend=xend, yend=yend, colour=as.factor(Genotype)))
}
