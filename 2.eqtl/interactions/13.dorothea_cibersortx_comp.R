################################################################################
#
# 13. Dorothea transcription factor activity inference and comparison
#
################################################################################

library(decoupleR)

# Change rownames of expression data to gene symbols
counts <- read.delim("logcpm_864_20412_hla.txt")
gene.info <- read.delim("gene_info_20412.txt")
counts$Gene <- gene.info$gene_name[match(rownames(counts), 
                                         gene.info$gene_id)]
counts <- counts[!(is.na(counts$Gene)), ]
counts <- counts[!duplicated(counts$Gene), ]
rownames(counts) <- counts$Gene
counts$Gene <- NULL
counts <- as.matrix(counts)

# read in SRS information
srs <- read.delim("full-gains-SRS-predictions_mNN-RF.tsv")
srs <- subset(srs, Assay == "RNA-seq")
# all(srs$Sample_id == colnames(counts))

net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))

sample_acts <- decouple(
   mat = counts,
   network = net,
   .source = "source",
   .target = "target",
   consensus_score = TRUE,
   minsize = 5)
# write.table(sample_acts, "sample_activities_consensus_ABC.txt", sep="\t")

# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'consensus') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs that differ between SRS groups
# lmm given repeat samples
srs.diff <- apply(sample_acts_mat, 2, function(x){
  srs <- srs$SRS[match(rownames(sample_acts_mat), srs$Sample_id)]
  srs <- as.factor(srs == "SRS1")
  gainsid <- substr(rownames(sample_acts_mat), 1, 
                    nchar(rownames(sample_acts_mat))-2)
  
  res <- lme4::lmer(x ~ srs + (1|gainsid), REML = F)
  null <- lme4::lmer(x ~ (1|gainsid), REML = F)
  
  return(c(summary(res)$coefficients[2, 1], anova(null, res)$'Pr(>Chisq)'[2]))
})

# permutation
perm.results <- data.frame("Perm"=1:1000, n.sig=NA)
gainsid <- substr(rownames(sample_acts_mat), 1, nchar(rownames(sample_acts_mat))-2)

for(p in 1:1000){
  print(p)
  srs.diff <- apply(sample_acts_mat, 2, function(x){
    srs.p <- sample(srs$SRS[match(rownames(sample_acts_mat), srs$Sample_id)],
                    size=nrow(srs), replace = F)
    srs.p <- as.factor(srs.p == "SRS1")
    res <- lme4::lmer(x ~ srs.p + (1|gainsid), REML = F)
    null <- lme4::lmer(x ~ (1|gainsid), REML = F)

    return(c(summary(res)$coefficients[2, 1], anova(null,
                                                  res)$'Pr(>Chisq)'[2]))
  })
  srs.diff <- data.frame(t(srs.diff))
  colnames(srs.diff) <- c("beta", "srs.diff")
  srs.diff$FDR <- p.adjust(srs.diff$srs.diff, method="fdr")
  perm.results[p, 2] <- nrow(subset(srs.diff, FDR<0.05))
}
# pdf("DoRothEA_permutation_histogram.pdf")
ggplot(perm.results, aes(n.sig)) +
  geom_histogram(bins=100, colour="black", fill="grey") +
  theme_bw() +
  coord_cartesian(expand=FALSE) +
  # scale_x_continuous(limits = c(0, 170)) +
  geom_vline(xintercept=253, col="red", lty=2) +
  xlab("Number of TF motifs significantly enriched in SRS interaction QTL") +
  ylab("Number of permutations")
# dev.off()
range(perm.results$n.sig) # 0-2

srs.diff <- data.frame(t(srs.diff))
colnames(srs.diff) <- c("beta", "srs.diff")
srs.diff$source <- rownames(srs.diff)
srs.diff$FDR <- p.adjust(srs.diff$srs.diff, method="fdr")

ggplot(srs.diff, aes(beta, -log10(FDR))) +
  geom_point(aes(colour=FDR<0.05)) +
  scale_colour_manual(values = c("lightgrey", "red")) +
  theme_bw() +
  geom_text_repel(aes(label=source))

write.table(srs.diff, "Dorothea_SRS_differential_activity_summ_stats.txt", sep="\t")


#######################
# compare to SNP2TFBS results
dorothea <- srs.diff
dorothea$TF <- toupper(dorothea$source)

# results from SNP2TFBS
srs.ld <- read.delim("SRS_ints_TFBS_enrichment.txt")
srs.ld$tf <- toupper(srs.ld$tf)
srs.ld$fullmotif <- srs.ld$tf
srs.ld <- separate_rows(srs.ld, tf, sep="_")
x <- list("SRS interactions enriched TFBS"=sort(as.character(srs.ld$tf[which(srs.ld$int_fdr < 0.05)])),
          "Dorothea SRS differential activity"=sort(dorothea$TF[which(dorothea$FDR < 0.05)]))
ggvenn(x, show_elements = F)
ggvenn(x, show_elements = T, label_sep = "\n")

dorothea.srs <- merge(dorothea, srs.ld, by.x="TF", by.y="tf")
dorothea.srs$Significance <- paste((dorothea.srs$FDR <0.05), 
                                   (dorothea.srs$int_fdr < 0.05))
dorothea.srs$Significance <- gsub("FALSE FALSE", "Neither",
                                  dorothea.srs$Significance)
dorothea.srs$Significance <- gsub("FALSE TRUE", "SNP2TFBS only",
                                  dorothea.srs$Significance)
dorothea.srs$Significance <- gsub("TRUE FALSE", "DoRothEA only",
                                  dorothea.srs$Significance)
dorothea.srs$Significance <- gsub("TRUE TRUE", "DoRothEA and SNP2TFBS",
                                  dorothea.srs$Significance)
dorothea.srs$label <- ifelse(dorothea.srs$TF == dorothea.srs$fullmotif,
                             dorothea.srs$TF, 
                             paste0(dorothea.srs$TF, "\n(", 
                                    dorothea.srs$fullmotif, ")"))
# pdf("TFBS_enrichment_vs_dorothea_activity_significance.pdf")
ggplot(dorothea.srs, aes(-log10(FDR), -log10(int_fdr))) +
  geom_point(aes(colour=Significance)) +
  scale_colour_manual(values=c("darkblue", "lightblue", "lightgrey",
                               "lightgrey",  "bisque2", "black")) +
  geom_hline(yintercept = -log10(0.05), lty=2) +
  geom_vline(xintercept = -log10(0.05), lty=2) +
  theme_bw() +
  xlab("Significance of Dorothea activity difference (FDR)") +
  ylab("Significance of SNP2TFBS motif enrichment (FDR)") +
  geom_text_repel(size = 2, aes(label=label, 
                                colour=(Significance == "DoRothEA and SNP2TFBS")), 
                  fontface="italic") +
  theme(legend.position = "none") +
  ggtitle("SRS interactions: putative driver transcription factors")
# dev.off()



###################
# compare to cell counts
cell.counts <- read.csv("CIBERSORTx_results_AK.csv", row.names=1)
sort(apply(cell.counts[, 1:23], 2, median))
cell.order <- names(sort(apply(cell.counts[, 1:23], 2, median)))
cibersort.m <- melt(as.matrix(cell.counts[order(cell.counts$S100A8.9_hi_neutrophils), 1:23]))
cibersort.m$Var2 <- factor(cibersort.m$Var2, levels = cell.order)

srs <- read.delim("full-gains-SRS-predictions_mNN-RF.tsv")
cibersort.m$SRS <- srs$SRS[match(cibersort.m$Var1, srs$Sample_id)]
cell.counts$SRS <- srs$SRS[match(rownames(cell.counts), srs$Sample_id)]
cell.counts$SRS <- ifelse(cell.counts$SRS == "SRS1", "SRS1", "Non-SRS1")
cell.counts$SRSq <- srs$SRSq[match(rownames(cell.counts), srs$Sample_id)]
cell.counts.m <- reshape2::melt(cell.counts[, c(1:23, 27), ])

tf.activity <- read.delim("gains_eqtl_tf_activities.txt")
tf.activity <- t(sample_acts_mat)
tf.activity.m <- reshape2::melt(t(tf.activity))

cell.counts <- cell.counts[match(colnames(tf.activity), rownames(cell.counts)), ]
colnames(tf.activity)[is.na(rownames(cell.counts))]

# restrict to enriched TFBS
diff.tfs <- tf.activity[match(unique(dorothea.srs$TF[which(dorothea.srs$FDR < 0.05 &
                                                             dorothea.srs$int_fdr < 0.05)]), 
                              rownames(tf.activity)), ]
diff.tfs <- diff.tfs[complete.cases(diff.tfs), ]

# Restrict to first sample only and test for correlation between cell props and TF activity
cell.counts$patient <- substr(rownames(cell.counts), 1,
                              nchar(rownames(cell.counts))-2)
patients <- cell.counts$patient
cell.counts$patients <- NULL
first.samples <- !(duplicated(patients))

cell.tf.cor <- list()
for(i in 1:23){
  cell.type <- cell.counts[first.samples, i]
  cell.type.res <- apply(diff.tfs[, first.samples], 1, function(x){
    res <- cor.test(x, cell.type, method = "spearman")
    return(c(res$estimate, res$p.value))
  })
  cell.type.res <- data.frame(t(cell.type.res))
  colnames(cell.type.res) <- c("Correlation", "PValue")
  cell.type.res$TF <- rownames(cell.type.res)
  cell.type.res$CellType <- colnames(cell.counts)[i]
  cell.tf.cor[[i]] <- cell.type.res
}
cell.tf.cor <- rbindlist(cell.tf.cor)
cell.tf.cor$padj <- p.adjust(cell.tf.cor$PValue, method="bonferroni")
# write.table(cell.tf.cor, "cell_proportion_estimates_correlations_TF_activities.txt", sep="\t")

corr.mat <- dcast(cell.tf.cor, formula = TF ~ CellType, value.var="Correlation")
rownames(corr.mat) <- corr.mat$TF
corr.mat$TF <- NULL

# get clusting order
hm <- pheatmap(corr.mat, angle_col = "45", na_col = "white")
tf.order <- hm$tree_col
cell.order <- hm$tree_row

cell.tf.cor$Correlation[which(cell.tf.cor$padj > 0.05)] <- NA
corr.mat <- dcast(cell.tf.cor, formula = TF ~ CellType, value.var="Correlation")
rownames(corr.mat) <- corr.mat$TF
corr.mat$TF <- NULL

# pdf("Fig_3e_cell_corr_tfs_heatmap_1sample.pdf", height=12, useDingbats = F)
pheatmap(corr.mat, angle_col = "45", na_col = "white", cluster_rows = cell.order,
         cluster_cols = tf.order)
# dev.off()
