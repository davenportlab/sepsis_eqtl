################################################################################
#
# 6. Summary of interaction effects
#
################################################################################

srs.int <- readRDS("srs_int_results_cond_incl.rds")
diag.int <- readRDS("diagnosis_int_results_cond_incl.rds")
cell.int <- readRDS("cellprops_int_results_cond_incl.rds")

# Remove the pairs where there weren't enough minor allele homozygotes to test
srs.int <- srs.int[complete.cases(srs.int), ]
# Correct p-values for multiple testing
srs.int$FDR <- p.adjust(srs.int$Interaction_pval, method="fdr")
# Categorise interactions as those where the main eQTL effect is magnified or dampened
srs.int$Type <- ifelse(srs.int$eQTL_beta*srs.int$Interaction_beta >0, "Magnifier", "Dampener")
# Get the significant results
srs.int <- srs.int[which(srs.int$FDR < 0.05), c(1:2, 9, 13, 14)]

diag.int <- diag.int[complete.cases(diag.int), ]
diag.int$FDR <- p.adjust(diag.int$Interaction_pval, method="fdr")
diag.int$Type <- ifelse(diag.int$eQTL_beta*diag.int$Interaction_beta >0, "Magnifier", "Dampener")
diag.int <- diag.int[which(diag.int$FDR < 0.05), c(1:2, 9, 13, 14)]

neut.int <- cell.int[!is.na(cell.int$Interaction_pval_neutrophil), ]
neut.int <- cell.int[complete.cases(cell.int), ]
neut.int <- neut.int[, c(1:12)]
neut.int$FDR <- p.adjust(neut.int$Interaction_pval_neutrophil, method="fdr")
neut.int$Type <- ifelse(neut.int$eQTL_beta*neut.int$Interaction_beta_neutrophil >0, "Magnifier", "Dampener")
neut.int <- neut.int[which(neut.int$FDR < 0.05), c(1:2, 9, 13, 14)]
colnames(neut.int)[3] <- "Interaction_beta"

lymph.int <- cell.int[!is.na(cell.int$Interaction_pval_lymphocytes), ]
lymph.int <- cell.int[complete.cases(cell.int), ]
lymph.int <- lymph.int[, c(1:2, 23:32)]
lymph.int$FDR <- p.adjust(lymph.int$Interaction_pval_lymphocytes, method="fdr")
lymph.int$Type <- ifelse(lymph.int$eQTL_beta*lymph.int$Interaction_beta_lymphocytes >0, "Magnifier", "Dampener")
lymph.int <- lymph.int[which(lymph.int$FDR < 0.05), c(1:2, 9, 13, 14)]
colnames(lymph.int)[3] <- "Interaction_beta"

mono.int <- cell.int[!is.na(cell.int$Interaction_beta_monocytes), ]
mono.int <- cell.int[complete.cases(cell.int), ]
mono.int <- mono.int[, c(1:2, 13:22)]
mono.int$FDR <- p.adjust(mono.int$Interaction_pval_monocytes, method="fdr")
mono.int$Type <- ifelse(mono.int$eQTL_beta*mono.int$Interaction_beta_monocytes >0, "Magnifier", "Dampener")
mono.int <- mono.int[which(mono.int$FDR < 0.05), c(1:2, 9, 13, 14)]
colnames(mono.int)[3] <- "Interaction_beta"

sig.int <- rbind(srs.int, diag.int, neut.int, lymph.int, mono.int)
sig.int$Variable <- c(rep("SRS", nrow(srs.int)),
                      rep("Diagnosis", nrow(diag.int)),
                      rep("Neutrophil", nrow(neut.int)),
                      rep("Lymphocyte", nrow(lymph.int)),
                      rep("Monocyte", nrow(mono.int)))
write.table(sig.int, "int_summary.txt", sep="\t", quote=F)

################################################################################
# Upset plot for interaction effects

library(tidyverse)
library(ggupset)
library(ggplot2)

# read in list of significant interaction effects
int.summary <- read.delim("int_summary.txt")
int.summary$pair <- paste0(int.summary$Gene, "_", int.summary$SNP)
# keep just pair, interaction variable and type
int.summary <- int.summary[, c(5:7)]
int.summary$Variable_Type <- paste0(int.summary$Variable, "_", int.summary$Type)

int.summary %>%
  as_tibble() %>%
  group_by(pair) %>%
  summarize(Variable_Type = list(paste0(Variable, "_", Type))) %>%
  # plot upset
  ggplot(aes(x=Variable_Type)) +
  ggtitle("Interaction effects") +
  geom_bar() +
  scale_x_mergelist(sep = "-") +
  scale_x_upset(order_by = "freq") +
  axis_combmatrix(sep = "-", override_plotting_function = function(df){
    df %>%
      mutate(magnifier = case_when(
        ! observed ~ "not observed",
        str_detect(single_label, "Magnifier") ~ "Magnifier",
        observed ~ "Dampener"
      )) %>%
      ggplot(aes(x= at, y= single_label)) +
      geom_rect(aes(fill= index %% 2 == 0), ymin=df$index-0.5, ymax=df$index+0.5, xmin=0, xmax=1) +
      geom_point(aes(color= magnifier), size = 3) +
      geom_line(data= function(dat) dat[dat$observed, ,drop=FALSE], aes(group = labels), size= 1.2) +
      ylab("") + xlab("") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_fill_manual(values= c(`TRUE` = "white", `FALSE` = "#F7F7F7")) +
      scale_color_manual(values= c("Magnifier" = "red", "Dampener" = "cornflowerblue",
                                   "not observed" = "lightgrey")) +
      guides(color="none", fill="none") +
      theme(
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0, "pt"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank()
      )
  }) +
  geom_text(stat='count', aes(label=..count..), hjust=-0.1, angle=90) +
  theme_bw() +
  theme_combmatrix(combmatrix.label.make_space = TRUE) +
  ylim(0, 502)

################################################################################

# Plot all interaction effects

library(qvalue)
library(RColorBrewer)
library(ggrepel)
library(XGR)
library(lme4)
library(reshape2)

RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
gene.info <- read.delim("gene_info_20412.txt")

# SRS
srs.int <- readRDS("srs_int_results_cond_incl.rds")
srs.int <- srs.int[complete.cases(srs.int), ]
srs.int$FDR <- p.adjust(srs.int$Interaction_pval, method="fdr")
table(srs.int$FDR < 0.05) # 1578
srs.int$Sig <- srs.int$FDR < 0.05
srs.int$Symbol <- gene.info$gene_name[match(srs.int$Gene, gene.info$gene_id)]

srs.int$Type <- (srs.int$eQTL_beta*srs.int$Interaction_beta) > 0
srs.int$Type <- gsub("TRUE", "Magnifier", srs.int$Type)
srs.int$Type <- gsub("FALSE", "Dampener", srs.int$Type)
table(srs.int$Type, srs.int$Sig)

srs.sig <- subset(srs.int, FDR < 0.05)
srs.sig <- srs.sig[order(srs.sig$FDR), ]

pdf("Fig_2A_srs_intns.pdf", useDingbats = F)
ggplot(srs.int, aes(eQTL_beta, Interaction_beta)) +
  geom_point(aes(colour=Sig)) +
  theme_bw() +
  geom_text_repel(data=srs.int[match(srs.sig$Gene[1:20], srs.int$Gene), ], 
                  aes(label=Symbol),
                  fontface="italic") +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  ggtitle(label = "SRS") + theme(legend.position = "none")
dev.off()

# Diagnosis
diag.int <- readRDS("diagnosis_int_results_cond_incl.rds")
diag.int <- diag.int[complete.cases(diag.int), ]
diag.int$FDR <- p.adjust(diag.int$Interaction_pval, method="fdr")
table(diag.int$FDR < 0.05) # 166
diag.int$Sig <- diag.int$FDR < 0.05
diag.int$Symbol <- gene.info$gene_name[match(diag.int$Gene, gene.info$gene_id)]

diagnosis.sig <- subset(diagnosis.int, FDR < 0.05)

pdf("Fig_1E_diagnosis_intns.pdf", useDingbats = F)
ggplot(diagnosis.int, aes(eQTL_beta, Interaction_beta)) +
  geom_point(aes(colour=Sig)) +
  theme_bw() +
  geom_text_repel(data=diagnosis.int[match(diagnosis.sig$Gene[1:20], diagnosis.int$Gene), ], 
                  aes(label=Symbol),
                  fontface="italic") +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  xlab("eQTL effect size") +
  ylab("Genotype x source interaction effect size") +
  ggtitle(label = "Diagnosis") + theme(legend.position = "none")
dev.off()

# Cell proportions
cell.int <- readRDS("cellprops_int_results_cond_incl.rds")
cell.int <- cell.int[complete.cases(cell.int), ]
cell.int$FDR_n <- p.adjust(cell.int$Interaction_pval_neutrophil, method="fdr")
cell.int$FDR_l <- p.adjust(cell.int$Interaction_pval_lymphocytes, method="fdr")
cell.int$FDR_m <- p.adjust(cell.int$Interaction_pval_monocytes, method="fdr")
cell.int$Sig_neutrophils <- cell.int$FDR_n < 0.05
cell.int$Sig_lymphocytes <- cell.int$FDR_l < 0.05
cell.int$Sig_monocytes <- cell.int$FDR_m < 0.05

cell.int$Symbol <- gene.info$gene_name[match(cell.int$Gene, gene.info$gene_id)]

table(cell.int$FDR_n < 0.05) # 1073
table(cell.int$FDR_l < 0.05) # 1013
table(cell.int$FDR_m < 0.05) # 608

write.table(cell.int, "cellprops_int_results_cond_incl.txt", sep="\t", row.names=FALSE)

pdf("FigS11_cell_proportion_intns.pdf", useDingbats = F, onefile=T)
ggplot(cell.int, aes(eQTL_beta_neutrophil, Interaction_beta_neutrophil)) +
  geom_point(aes(colour=FDR_n<0.05)) +
  theme_bw() +
  geom_text_repel(data=cell.int[order(cell.int$FDR_n), ][1:20, ], 
                  aes(label=Symbol),
                  fontface="italic") +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(label = "Neutrophil proportion") + theme(legend.position = "none") +
  xlab("eQTL effect size") + ylab("Interaction effect size")

ggplot(cell.int, aes(eQTL_beta_lymphocytes, Interaction_beta_lymphocytes)) +
  geom_point(aes(colour=FDR_l<0.05)) +
  theme_bw() +
  geom_text_repel(data=cell.int[order(cell.int$FDR_l), ][1:20, ], 
                  aes(label=Symbol),
                  fontface="italic") +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(label = "Lymphocyte proportion") + theme(legend.position = "none") +
  xlab("eQTL effect size") + ylab("Interaction effect size")

ggplot(cell.int, aes(eQTL_beta_monocytes, Interaction_beta_monocytes)) +
  geom_point(aes(colour=FDR_m<0.05)) +
  theme_bw() +
  geom_text_repel(data=cell.int[order(cell.int$FDR_m), ][1:20, ], 
                  aes(label=Symbol),
                  fontface="italic") +
  scale_color_manual(values=c("grey", "red")) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(label = "Monocyte proportion") + theme(legend.position = "none") +
  xlab("eQTL effect size") + ylab("Interaction effect size")
dev.off()

################################################################################
# Pathway enrichment

# SRS
eTerm <- xEnricherGenes(srs.sig$Symbol, background=srs.int$Symbol,
                        ontology = "MsigdbC2REACTOME")
eTerm.res <- xEnrichViewer(eTerm, top_num = "auto")

eTerm <- xEnricherGenes(srs.int$Symbol[srs.int$Sig == TRUE & srs.int$Type == "Magnifier"], 
                        background=srs.int$Symbol,
                        ontology = "MsigdbC2REACTOME")
eTerm.res <- xEnrichViewer(eTerm, top_num = "auto")

eTerm <- xEnricherGenes(srs.int$Symbol[srs.int$Sig == TRUE & srs.int$Type == "Dampener"], 
                        background=srs.int$Symbol,
                        ontology = "MsigdbC2REACTOME")
eTerm.res <- xEnrichViewer(eTerm, top_num = "auto")

# Diagnosis
eTerm <- xEnricherGenes(diagnosis.sig$Symbol, background=diagnosis.int$Symbol,
                        ontology = "MsigdbC2REACTOME")
eTerm.res <- xEnrichViewer(eTerm, top_num = "auto")

data <- diagnosis.sig[, c("Symbol", "Interaction_beta")]
subg_path <- xSubneterGenes(data=data, network="STRING_high",
                            subnet.size=50, RData.location=RData.location)

subg_path <- xSubneterGenes(data=data, network="PCommonsDN_medium",
                            subnet.size=50, RData.location=RData.location)
subg <- subg_path
pattern <- log2(as.numeric(V(subg)$significance))
xVisNet(g=subg, pattern=pattern, vertex.shape="circle", newpage=F)

################################################################################

# Permutation results

library(ggplot2)

srs.perm <- read.delim("srs_permutation_results.txt", header=F)
srs.perm$Permutation <- 1:nrow(srs.perm)
true.result <- read.delim("int_summary.txt")
true.result <- nrow(subset(true.result, Variable == "SRS"))

pdf("Fig_S10_srs_intns_perm.pdf", useDingbats = F)
ggplot(srs.perm, aes(V1)) +
  geom_histogram(bins = 100, colour="black", fill="grey") +
  theme_bw() +
  coord_cartesian(expand=FALSE) +
  geom_vline(xintercept=true.result, col="red", lty=2) +
  xlab("Number of eQTL with a significant interaction with SRS1 status") +
  ylab("Number of permutations")
dev.off()

diag.perm <- read.delim("diag_permutation_results.txt", header=F)
diag.perm$Permutation <- 1:nrow(diag.perm)
true.result <- read.delim("int_summary.txt")
true.result <- nrow(subset(true.result, Variable == "Diagnosis"))

pdf("Fig_S7_diagnosis_intns_perm.pdf", useDingbats = F)
ggplot(diag.perm, aes(V1)) +
  geom_histogram(bins=100, colour="black", fill="grey") +
  theme_bw() +
  coord_cartesian(expand=FALSE) +
  geom_vline(xintercept=true.result, col="red", lty=2) +
  xlab("Number of eQTL with a significant interaction with source of sepsis") +
  ylab("Number of permutations")
dev.off()
