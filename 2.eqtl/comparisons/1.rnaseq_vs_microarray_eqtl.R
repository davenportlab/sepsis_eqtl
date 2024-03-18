################################################################################
#
# 1. Replication of microarray eQTL results
#
################################################################################

# NB: 134 samples used in the published microarray eQTL mapping were repeated in
# the RNA-seq cohort. To explore replication, we re-ran the initial RNA-seq eQTL
# mapping without these samples and used the same thresholds to determine significance

# reduce to non-overlapping samples
rnaseq.unique <- read.delim("../replication_nodups/rnaseq_unique_samples.txt")
 
for(chr in 1:22){
  # read in the data
  load(paste("../data/eqtl_files_", chr, ".rda", sep=""))

  to.use <- match(rnaseq.unique$x, colnames(exp))
  to.use <- to.use[!is.na(to.use)]
  # non-genotyped samples already excluded ==> 688 samples :)

  exp <- exp[, to.use]
  covs <- covs[to.use, ]
  peer.factors <- peer.factors[to.use, ]
  GAinSID <- GAinSID[to.use]
  geno <- geno[to.use, ]

  eQTL_files <- paste("../data/eqtl_files_nodups_", chr, ".rda", sep="")
  save(list=c("exp", "geno", "covs", "GAinSID", "peer.factors"),
       file = eQTL_files)
}

################################################################################

# use these files as input to cis-eqtl mapping (script 2.ciseqtl_mapping.R) and collate results as before

################################################################################

# compare to microarray cis eqtl
library(ggplot2)
library(ggrastr)

# Results from Davenport et al 2016 - all pairs tested with p value < 0.01
microarray <- read.delim("../replication_nodups/cis eQTL 30 PCs 240 samples.txt")
microarray$Sig <- microarray$FDR < 0.05
microarray$Pair <- paste0(microarray$SNP, "_", microarray$gene)

# get allele tested
bim <- read.table("../replication_nodups/genotyping_for_eQTL_published.bim")
microarray$MinorAllele <- bim$V5[match(microarray$SNP, bim$V2)]
microarray$MajorAllele <- bim$V6[match(microarray$SNP, bim$V2)]

# read in reduced RNA-seq results (non-conditioned as not done in microarray data)
gains.eqtl <- readRDS("../replication_nodups/ciseqtl_all.rds")
gains.eqtl$Sig <- gains.eqtl$pvalue < gains.eqtl$threshold
table(gains.eqtl$Sig)

gene.info <- read.delim("../data//gene_info_864_20412_hla.txt",
                        stringsAsFactors = F)
gains.eqtl$Symbol <- gene.info$gene_name[match(gains.eqtl$gene,
                                               gene.info$gene_id)]
gains.eqtl$Pair <- paste0(gains.eqtl$snps, "_", gains.eqtl$Symbol)

# get alelle information
bim <- read.table("../data/genotyping_for_rna-seq_eQTL.bim")
gains.eqtl$MinorAllele <- bim$V5[match(gains.eqtl$snps, bim$V2)]
gains.eqtl$MajorAllele <- bim$V6[match(gains.eqtl$snps, bim$V2)]

# merge
comp <- merge(gains.eqtl, microarray, by="Pair")
comp$Sig <- paste(comp$Sig.x, comp$Sig.y)

# Check that major and minor alleles are the same
comp$action <- "NA"
comp[which(comp$MinorAllele.x == comp$MinorAllele.y), "action"] <- "same"
# Some alleles can't be distinguished because of strand
comp[which(comp$MinorAllele.x=="A" & comp$MajorAllele.x=="T" |
             comp$MinorAllele.x=="T" & comp$MajorAllele.x=="A" |
             comp$MinorAllele.x=="C" & comp$MajorAllele.x=="G" |
             comp$MinorAllele.x=="G" & comp$MajorAllele.x=="C"), 
     "action"] <-"strand"
# Same alleles but different strand
comp[which(comp$MinorAllele.x=="A" & comp$MinorAllele.y=="T" & comp$action != "strand" |
             comp$MinorAllele.x=="T" & comp$MinorAllele.y=="A" & comp$action != "strand" |
             comp$MinorAllele.x=="C" & comp$MinorAllele.y=="G" & comp$action != "strand" |
             comp$MinorAllele.x=="G" & comp$MinorAllele.y=="C" & comp$action != "strand"),
     "action"] <-"same"
# Other allele tested: flip direction of effects
comp[which(comp$action == "NA"), "action"] <- "flip"
table(comp$action, useNA = "ifany")
comp$beta.x.2 <- comp$beta.x
comp[which(comp$action == "flip"), "beta.x.2"] <- 
  comp[which(comp$action == "flip"), "beta.x.2"] * (-1)

# Remove SNPs which can't be distinguished because of strand
comp$beta.x.2[which(comp$action == "strand")] <- NA
comp <- comp[complete.cases(comp), ]

cor.test(comp$beta.x.2, comp$beta.y) # 0.70

pdf("FigS6_microarray_replication_sig_nodups.pdf", useDingbats = F, width=7, height=5)
ggplot(comp, aes(x=beta.y, y=beta.x.2) ) +
  rasterize(geom_point(aes(colour=Sig %in% c("TRUE FALSE", "TRUE TRUE"))), dpi=1000) +
  scale_colour_manual(values=c("lightgrey", "darkblue"), name="Significant in RNA-seq") +
  xlab("Microarray beta value for SNP-gene pair") +
  ylab("RNA-seq beta value for SNP-gene pair") +
  xlim(c(-2, 4)) + ylim(c(-2, 4)) +
  geom_abline(intercept = 0, slope=1, lty=2) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  annotate("text", label="85,678 SNP-gene pairs tested in both\n Pearson's r=0.70", x=1, y=-1.5)
dev.off()

comp <- comp[, c("SNP", "MinorAllele.x", "MajorAllele.x", "gene.x", "Symbol", "beta.x.2", "pvalue", "Sig.x",
                 "beta.y", "p.value", "FDR", "Sig.y")]
colnames(comp) <- c("SNP", "Minor Allele", "Major Allele", "Gene", "Symbol", "RNAseq_beta", "RNAseq_pvalue",
                    "RNAseq_significant", "Microarray_beta", "Microarray_pvalue", "Microarray_FDR", 
                    "Microarray_significant")
write.table(comp, "microarray_comparison_noreplicates.txt", sep="\t", row.names = FALSE)
