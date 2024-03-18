################################################################################
#
# Investigation of sepsis-specific eQTL effects
#
################################################################################

options(stringsAsFactors = FALSE)
library(XGR)
library(ggplot2)
library(gridExtra)
library(data.table)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

#################################################################################

# We compare the GAinS results (lead SNP for each significant eGene) to those from
# GTEx whole blood EUR individuals where the same SNP-gene pair was tested 
# (n=8137). Only a subset were significant in the GTEx analysis but they are well
# correlated in general.

mashr.input <- read.delim("gtex_gains_comp_lead_for_mashr.txt")
ggplot(mashr.input, aes(beta/se, slope/slope_se)) +
  geom_point() +
  theme_bw() +
  xlab("GAinS Z-score") +
  ylab("GTEx Z-score")

ggplot(subset(mashr.input, Sig==TRUE), aes(beta, slope)) +
  geom_point(aes(colour=gtex.sig)) +
  scale_colour_manual(values=c("gold", "darkgreen")) +
  theme_bw() +
  ggtitle("Input for mashr: initial effect sizes for GAinS lead eSNPs") +
  xlab("GAinS beta") +
  ylab("GTEx slope") +
  geom_abline(intercept = 0, slope=1, lty=2)

# Mashr returns posterior effect sizes and lfsr (local false sign rate) for each 
# SNP-gene pair in GAinS and GTEx, which can be used to assess significance and 
# sharing of eQTL across conditions. There are 15 GAinS eQTL that were deemed not
# significant in either case in the mashr outputs. These borderline cases are 
# not considered for sepsis-specificity.

mashr.results <- read.delim("gains_gtex_mashr_results.txt")
table(mashr.results$GAinS.lfsr < 0.05)
mashr.results$shared <- NA
mashr.results$shared[which(mashr.results$GAinS.lfsr > 0.05)] <- "not significant"

# "Shared effects" are defined by default as those where the
# direction of effect is the same and within a factor of 0.5 of each other; i.e. 
# the ratio of GAinS effect/GTEx effect is >0.5 but < 1/0.5=2. 

mashr.results$shared[is.na(mashr.results$shared)] <- ifelse(
  mashr.results$diff[is.na(mashr.results$shared)] == "FALSE", "shared", "not shared")

# We then want to think about the eQTL that don't meet this definition of sharing. 
# We term all "context-specific" but want to identify the subset that have bigger 
# effect sizes in sepsis/HV and those with a convincing opposite direction of effect 
# for the same SNP. 

# The eQTL with an opposite direction of effect include some that aren't significant
# in GTEx i.e. are a sepsis-specific effect. We therefore define any eQTL 
# significant in GTEx but with opposite direction of
# effect to GAinS as a specific subset and those not significant in GTEx as sepsis-magnified

ggplot(subset(mashr.results, shared == "not shared"), aes(GAinS, GTEx)) +
  geom_point(aes(colour=(ratio < 0))) +
  scale_colour_manual(values=c("lightgrey", "red")) +
  theme_bw() +
  ggtitle("Not shared: opposite effects") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

ggplot(subset(mashr.results, ratio < 0), aes(GAinS, GTEx)) +
  geom_point(aes(colour=GTEx.lfsr < 0.05)) +
  scale_colour_manual(values=c("lightgrey", "red")) +
  theme_bw() +
  ggtitle("Not shared: opposite effects") +
    xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

mashr.results$shared[which(mashr.results$ratio < 0 &
                                mashr.results$GTEx.lfsr < 0.05 &
                                mashr.results$GAinS.lfsr < 0.05)] <- "Opposite direction of effect"
ggplot(subset(mashr.results, shared != "shared"), aes(GAinS, GTEx)) +
  geom_point(aes(colour=(shared=="Opposite direction of effect"))) +
  scale_colour_manual(values=c("lightgrey", "red")) +
  theme_bw() +
  ggtitle("Not shared: convincing opposite effects") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

# Then we look at the remaining cases.

ggplot(subset(mashr.results, shared == "not shared"), 
       aes(GAinS, GTEx)) +
  geom_point(colour="darkblue") +
  theme_bw() +
  ggtitle("Not shared") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_abline(intercept = 0, slope=1, lty=2) +
    geom_abline(intercept = 0, slope=0.5, lty=2) +
    geom_abline(intercept = 0, slope=1/0.5, lty=2)

# We are particularly interested in the "sepsis-magnified" eQTL, so those where the
# effect size is bigger in GAinS than GTEx. Some of these may not be significant 
# in GTEx but it makes sense to consider these sepsis-magnified too, as there is a
# significant effect in GAinS only (regardless of the sign of the GTEx effect size).

ggplot(subset(mashr.results, shared == "not shared"), 
       aes(GAinS, GTEx)) +
  geom_point(aes(colour=GTEx.lfsr < 0.05)) +
  scale_color_manual(values=c("lightgrey", "blue")) +
  theme_bw() +
  ggtitle("Not shared") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_abline(intercept = 0, slope=1, lty=2) +
    geom_abline(intercept = 0, slope=0.5, lty=2) +
    geom_abline(intercept = 0, slope=1/0.5, lty=2)

mashr.results$shared[which(mashr.results$ratio < 0 & 
                             mashr.results$GTEx.lfsr >= 0.05 &
                             mashr.results$GAinS.lfsr < 0.05)] <- "Only significant in GAinS"
mashr.results$shared[which(mashr.results$shared == "not shared")] <- ifelse(
  mashr.results$ratio[which(mashr.results$shared == "not shared")] >2, 
                                  "Bigger effect in GAinS", 
                                  "Bigger effect in GTEx")

ggplot(subset(mashr.results, shared=="Bigger effect in GAinS"), aes(GAinS, GTEx)) +
  geom_point(aes(colour=GTEx.lfsr < 0.05)) +
  scale_colour_manual(values=c("gold", "darkgreen")) +
  theme_bw() +
  ggtitle("Not shared: sepsis-magnified") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_abline(intercept = 0, slope=1, lty=2) +
    geom_abline(intercept = 0, slope=0.5, lty=2) +
    geom_abline(intercept = 0, slope=1/0.5, lty=2)

# Those with a smaller effect in sepsis are all significant in GTEx

ggplot(subset(mashr.results, shared=="Bigger effect in GTEx"), aes(GAinS, GTEx)) +
  geom_point(aes(colour=GTEx.lfsr < 0.05)) +
  scale_colour_manual(values=c("blue")) +
  theme_bw() +
  ggtitle("Not shared: sepsis-reduced") +
    xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_abline(intercept = 0, slope=1, lty=2) +
    geom_abline(intercept = 0, slope=0.5, lty=2) +
    geom_abline(intercept = 0, slope=1/0.5, lty=2)

# Final categories:
# eQTL that are significant in sepsis are classed as "shared" (effect size same
# direction and within a factor of 0.5 of the GTEx effect size) or "context specific"
# (not fulfilling these criteria). Of the context-specific eQTL, those significant 
# in both GAinS and GTex but with opposite directions of effects are classed as 
# "opposite effects", and the remainder are divided into those with bigger effects/
# only significant in GTEx or in sepsis ("sepsis-dampened/magnified"). Those non
# significant in GTEx are labelled separately

pdf("Fig_1F_mashr_scatterplot.pdf", width=7, height=5)
ggplot(mashr.results, aes(GAinS, GTEx)) +
  geom_point(aes(colour=shared)) +
  scale_colour_manual(values=c("darkgreen", "darkblue", "lightgrey", "darkseagreen3", "gold", "darkgrey")) +
  theme_bw() +
  ggtitle("Comparison of sepsis and GTEx effect sizes with mashr") +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline(intercept = 0, slope=1, lty=2) +
  geom_abline(intercept = 0, slope=0.5, lty=2) +
  geom_abline(intercept = 0, slope=1/0.5, lty=2)
dev.off()

gene.info <- read.delim("../../../Data/gene_info_20412.txt")
mashr.results$gene <- mashr.input$gene[match(rownames(mashr.results), mashr.input$Pairs)]
mashr.results$gene_name <- gene.info$gene_name[match(mashr.results$gene, gene.info$gene_id)]
mashr.results$snps <- mashr.input$snps[match(rownames(mashr.results), mashr.input$Pairs)]
write.table(mashr.results, "sepsis_specificity_categories.txt", sep="\t", quote=F)

############################ 
# Enrichment analysis
mashr <- subset(mashr.results, shared != "shared" & shared != "not significant")
opp <- subset(mashr.results, shared == "Opposite direction of effect")
background <- mashr.results

# SNP enrichment
snp.list <- as.character(unique(mashr$snps))
background.snps <- as.character(unique(background$snps))
eTerm_EF <- xEnricherSNPs(snp.list,
                          background=background.snps,
                          ontology="EF_disease",
                          LD.r2=0.8, true.path.rule=T, ontology.algorithm="lea",
                          RData.location=RData.location)
xEnrichViewer(eTerm_EF)

# Gene enrichment
ss.genes <- unique(mashr$gene_name)
ss.genes <- unique(opp$gene_name)
background.genes <- gene.info$gene_name[match(background$gene, gene.info$gene_id)]

eTerm <- xEnricherGenes(ss.genes, background=background.genes,
                        ontology = "MsigdbC2REACTOME")
xEnrichViewer(eTerm)

# network analysis - needs pvalue rather than fold change
data <- mashr[, c("gene_name", "BF.FDR")]
subg_path <- xSubneterGenes(data=data, RData.location=RData.location)
subg <- subg_path
pattern <- -log10(as.numeric(V(subg)$significance))
xVisNet(g=subg, pattern=pattern, vertex.shape="circle", newpage=F, legend=F, 
        colorbar = F,vertex.size = 7)

#################### magnified in sepsis
# restrict to those with larger effect sizes in sepsis
ggplot(mashr.results, 
       aes(GAinS, GTEx, colour=(shared == "Bigger effect in GAinS"))) +
  geom_point() +
  theme_bw() +
  xlab("GAinS posterior effect size") +
  ylab("GTEx posterior effect size") +
  scale_color_manual(values=c("lightgrey", "darkblue"))

up.sepsis <- mashr[which(mashr$shared == "Bigger effect in GAinS"), ]
ssup.genes <- unique(up.sepsis$gene_name)
eTerm <- xEnricherGenes(ssup.genes, background=background.genes,
                        ontology = "MsigdbC2REACTOME")
xEnrichViewer(eTerm)
eTerm <- xEnricherGenes(ssup.genes, background=background.genes,
                        ontology = "GOBP")
xEnrichViewer(eTerm)

opp.genes <- mashr$gene_name[which(mashr$shared == "Opposite direction of effect")]
eTerm <- xEnricherGenes(opp.genes, background=background.genes,
                        ontology = "MsigdbC2REACTOME")
xEnrichViewer(eTerm)

# are there networks within the sepsis magnified egenes?
x <- signif(as.numeric(E(sem.sim)$weight), digits=2)
edge.width <- 1 + (x-min(x))/(max(x)-min(x))*3
xVisNet(g=sem.sim, vertex.shape="circle", edge.width=edge.width,
        edge.label=x, edge.label.cex=0.7, vertex.label.color="black")
output <- igraph::get.data.frame(sem.sim, what="edges")

data <- mashr[which(mashr$ratio >1), c("gene_name", "BF.FDR")]
subg_path <- xSubneterGenes(data=data, RData.location=RData.location)
subg <- subg_path
pattern <- -log10(as.numeric(V(subg)$significance))
xVisNet(g=subg, pattern=pattern, vertex.shape="circle", newpage=F, legend=T, 
        vertex.size = 7)

#############################################################################
# overlap with eQTL interactions - which have been run on the conditional
# signals so should only consider 1st signal/on the gene level
ints <- read.delim("../interactions/int_summary.txt")
# restrict to primary signal
cond <- read.delim("../conditional/conditional_eQTL_results_final.txt")
ints$pair <- paste0(ints$Gene, "_", ints$SNP)
cond$pair <- paste0(cond$Gene, "_", cond$SNP)
ints <- subset(ints, pair %in% cond$pair[which(cond$Number == 1)])

# reduce to those tested vs gtex
tested.genes <- unique(mashr.input$gene[which(mashr.input$Sig == TRUE)])
int.enrich <- data.frame("Gene"=tested.genes,
                          "Int"=tested.genes %in% ints$Gene)

# how many were also sepsis specific
int.enrich$ss <- int.enrich$Gene %in% mashr$gene[mashr$shared != "shared"]
table(int.enrich$Int, int.enrich$ss)
fisher.test(table(int.enrich$Int, int.enrich$ss))

# how many were magnified in sepsis
int.enrich$ssup <- int.enrich$Gene %in% up.sepsis$gene
table(int.enrich$Int, int.enrich$ssup)
fisher.test(table(int.enrich$Int, int.enrich$ssup))

# specific interaction types
ints$TypeVariable <- paste0(ints$Type, ints$Variable)
int.types <- unique(ints$TypeVariable)
ints$ss <- ints$Gene %in% mashr$gene
ints$ssup <- ints$Gene %in% up.sepsis$gene
table(ints$ss, ints$TypeVariable)
table(ints$ssup, ints$TypeVariable)
mashr.results$ss <- rownames(mashr.results) %in% rownames(mashr)[mashr$shared != "shared"]
mashr.results$ssup <- rownames(mashr.results) %in% rownames(up.sepsis)
mashr.results$opp <- mashr.results$shared == "Opposite direction of effect"

mashr.results$gene <- mashr.input$gene[match(rownames(mashr.results), mashr.input$Pairs)]

results.ints <- data.frame("inttype"=int.types,
                           "pvalss"=NA,
                           "pvalsssup"=NA)
for(i in 1:length(int.types)){
  print(int.types[i])
  ints.red <- subset(ints, TypeVariable == int.types[i])
  mashr.results$inti <- mashr.results$gene %in% ints.red$Gene
  print(table(mashr.results$inti, mashr.results$ss))
  results.ints$pvalss[i] <- (fisher.test(table(mashr.results$inti, mashr.results$ss)))$p.value
  print(table(mashr.results$inti, mashr.results$ssup))
  results.ints$pvalsssup[i] <- (fisher.test(table(mashr.results$inti, mashr.results$ssup)))$p.value
}
results.ints$fdrss <- p.adjust(results.ints$pvalss, method="fdr")
results.ints$fdrsup <- p.adjust(results.ints$pvalsssup, method="fdr")
results.ints

################################################################################################

# Some further comparisons of the SNPs/eQTLs with shared and specific effects

primary.peak.results <- read.delim("ciseqtl_eigenMT_corrected.txt")
primary.peak.results <- subset(primary.peak.results, Sig == TRUE)

ss <- read.delim("sepsis_specificity_categories.txt")
ss <- subset(ss, shared %in% c("Bigger effect in GAinS",
                               "Bigger effect in GTEx",
                               "shared"))
primary.peak.results$pair <- paste0("chr", primary.peak.results$chr, "_", 
                                    primary.peak.results$SNPpos, "_", primary.peak.results$gene)
primary.peak.results <- merge(primary.peak.results, ss, 
                              by.x = "pair",
                              by.y="row.names")

gex <- read.delim("gex_rnaseq_for_eqtl_hla.txt")
gene.info <- read.delim("gene_info_20412.txt")
gene.info$TSS <- gene.info$start
gene.info$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]

freq <- read.table("Genotyping_for_rna-seq_eQTL.frq", header=T)

# Compare MAF
primary.peak.results$MAF <- freq$MAF[match(primary.peak.results$snps.x,
                                           freq$SNP)]
pdf("FigS9_sepsisdependency_MAF.pdf", width=7, height=5)
ggplot(subset(primary.peak.results, shared %in% c("Bigger effect in GAinS",
                                                  "shared")),
              aes(MAF*100, group=shared, colour=shared)) +
  geom_density() +
  theme_bw() +
  xlab("Minor allele frequency (%)")
dev.off()

aggregate(primary.peak.results$MAF, by=list(primary.peak.results$shared),
          median)

# Mann-Whitney test - independent samples
wilcox.test(MAF ~ shared, data=subset(primary.peak.results, 
                                      shared %in% c("Bigger effect in GAinS",
                                                  "shared")))

# sepsis specific eSNPs are lower in frequency than the eSNPs for eQTL shared with GTEx.

# TSS
primary.peak.results$start <- gene.info$start[match(primary.peak.results$gene.x,
                                                    gene.info$gene_id)]
primary.peak.results$end <- gene.info$end[match(primary.peak.results$gene.x,
                                                    gene.info$gene_id)]
primary.peak.results$TSS <- gene.info$TSS[match(primary.peak.results$gene.x,
                                                    gene.info$gene_id)]
primary.peak.results$DistanceFromTSS <- primary.peak.results$SNPpos - primary.peak.results$TSS
ggplot(primary.peak.results, aes(DistanceFromTSS, abs(beta))) +
  geom_point() +
  theme_bw()

pdf("FigS9_sepsisdependency_TSS.pdf", width=7, height=5)
ggplot(subset(primary.peak.results, shared %in% c("Bigger effect in GAinS",
                                                  "shared")),
       aes(DistanceFromTSS, group=shared, 
                                 colour=shared)) +
  geom_density(alpha=0.3) + 
  scale_x_continuous(labels = scales::comma) +
  xlab("Distance from eGene transcriptional start site") +
  theme_bw()
dev.off()

pdf("FigS9_sepsisdependency_promoter.pdf", width=7, height=5)
ggplot(subset(primary.peak.results, shared %in% c("Bigger effect in GAinS",
                                                  "shared")), 
       aes(abs(DistanceFromTSS), group=shared, 
                                 colour=shared)) +
  geom_density() + 
  scale_x_log10(labels = scales::comma) +
  theme_bw()
dev.off()

wilcox.test(abs(DistanceFromTSS) ~ shared, 
             data=subset(primary.peak.results, shared %in% c("Bigger effect in GAinS", "shared")))
