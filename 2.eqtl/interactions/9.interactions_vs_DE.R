################################################################################
#
# 9. Do SRS interactions associate with change in gene expression?
#
################################################################################

de <- read.delim("SRS_DE.txt") # from SepstratifieR

srs.int <- readRDS("srs_int_results_cond_incl.rds")
srs.int <- srs.int[complete.cases(srs.int), ]
srs.int$fdr <- p.adjust(srs.int$Interaction_pval, method="fdr")
srs.de.int <- merge(de, srs.int, by.x="row.names", by.y="Gene")

fisher.test(table(srs.de.int$adj.P.Val < 0.05, srs.de.int$fdr < 0.05), alternative="greater")
# 0.0004587

# For genes significant in both, compare interaction type to direction of DE
table(ifelse(sign(srs.de.int$Interaction_beta*srs.de.int$eQTL_beta)==1, 
             "Magnifier", "Dampener")[which(srs.de.int$fdr < 0.05 & srs.de.int$adj.P.Val < 0.05)],
      sign(srs.de.int$logFC[which(srs.de.int$fdr < 0.05 & srs.de.int$adj.P.Val < 0.05)]))
#           -1   1
# Dampener  366 105
# Magnifier 283 628

# much more likely for eQTL effect sizes to increase as gene expression increases
fisher.test(table(ifelse(sign(srs.de.int$Interaction_beta*srs.de.int$eQTL_beta)==1, 
                         "Magnifier", "Dampener")[which(srs.de.int$fdr < 0.05 & srs.de.int$adj.P.Val < 0.05 )],
                  sign(srs.de.int$logFC[which(srs.de.int$fdr < 0.05 & srs.de.int$adj.P.Val < 0.05)])))$p.value
