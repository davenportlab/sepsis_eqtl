################################################################################
#
# 2. Comparison to GTEx whole blood (EUR)
#
################################################################################

library(data.table)
library(ggplot2)
options(stringsAsFactors = FALSE)
# GTEx files downloaded as a shared resource: resources/GTEx/AnalysisV8

# Sepsis lead SNP-gene pairs
gains <- read.delim("ciseqtl_eigenMT_corrected.txt")
gains.bim <- data.frame(fread("genotyping_for_rna-seq_eQTL.bim"))
m1 <- match(gains$snps, gains.bim$V2)
gains$minor <- gains.bim$V5[m1]
gains$major <- gains.bim$V6[m1]

# SNP in same format as GTEx
gains$SNP <- paste0("chr", gains$chr, "_", gains$SNPpos)
gains$Pairs <- paste0(gains$SNP, "_", gains$gene)

# Extract GTEx results if the SNP-gene pair was a lead SNP-gene pair in GAinS
gtex.files <- list.files(path="resources/GTEx/AnalysisV8/GTEx_Analysis_v8_EUR_eQTL_all_associations_csv/",
                          pattern="Whole_Blood.v8.EUR.allpairs.chr*", full.names = TRUE)
gtex.results <- lapply(gtex.files, function(x){
  gtex <- fread(x)
  gtex.variants <- unlist(strsplit(gtex$variant_id, "\\_"))
  gtex$SNP <- paste0(gtex.variants[seq(1, nrow(gtex)*5, 5)], "_",
                     gtex.variants[seq(2, nrow(gtex)*5, 5)])
  gtex$REF <- gtex.variants[seq(3, nrow(gtex)*5, 5)]
  gtex$ALT <- gtex.variants[seq(4, nrow(gtex)*5, 5)]
  gtex$phenotype_id <- unlist(strsplit(gtex$phenotype_id, "\\."))[seq(1, nrow(gtex)*2, 2)]
  gtex$Pairs <- paste0(gtex$SNP, "_", gtex$phenotype_id)
  gtex <- subset(gtex, Pairs %in% gains$Pairs)
  print(dim(gtex))
  gtex <- as.data.frame(gtex)
  return(gtex)
})

gtex.results.df <- rbindlist(gtex.results)
write.table(gtex.results.df, "../gtex/gtex_snp_gene_pairs_in_gains_lead.txt", sep="\t", quote=F)
gtex.results.df <- read.delim("Results/RNASeq/sepsis_specificity/gtex_snp_gene_pairs_in_gains_lead.txt")

comp <- merge(gtex.results.df, gains, by="Pairs")

# check alleles/direction of effect
comp$action <- "NA"
comp[which(comp$ALT == comp$minor), "action"] <- "same"

# Alleles can't distinguish because of strand
comp[which(comp$ALT=="A" & comp$REF=="T" |
             comp$ALT=="T" & comp$REF=="A" |
             comp$ALT=="C" & comp$REF=="G" |
             comp$ALT=="G" & comp$REF=="C"), "action"] <-"strand"

# Same but different strand
comp[which(comp$ALT=="A" & comp$minor=="T" & comp$action != "strand" |
             comp$ALT=="T" & comp$minor=="A" & comp$action != "strand" |
             comp$ALT=="C" & comp$minor=="G" & comp$action != "strand" |
             comp$ALT=="G" & comp$minor=="C" & comp$action != "strand"),
     "action"] <-"same"

# Flip effects for those that are opposite
comp[which(comp$action == "NA"), "action"] <- "flip"
comp[which(comp$action == "flip"), "beta"] <- comp[which(comp$action == "flip"), "beta"] * (-1)

# Remove SNPs which can't be distinguished because of strand
comp$beta[which(comp$action == "strand")] <- NA
comp <- comp[complete.cases(comp), ]

# use gtex eGene qval to decide if gene is an eGene, then use nominal threshold 
# to determine if a significant eSNP.
gtex.egenes <- read.delim("Whole_Blood.v8.EUR.egenes.txt.gz")
gtex.egenes <- subset(gtex.egenes, qval < 0.05)

gtex.egenes$phenotype_id <- unlist(strsplit(gtex.egenes$phenotype_id, "\\."))[seq(1, nrow(gtex.egenes)*2, 2)]
comp$gtex.egene <- comp$phenotype_id %in% gtex.egenes$phenotype_id
table(comp$gtex.egene) # 6519 significant

comp$nominal_threshold <- gtex.egenes$pval_nominal_threshold[match(comp$phenotype_id,
                                                                   gtex.egenes$phenotype_id)]
comp$gtex.sig <- comp$pval_nominal <= comp$nominal_threshold & comp$gtex.egene == TRUE
comp$gtex.sig[is.na(comp$gtex.sig)] <- "FALSE"
table(comp$gtex.sig, comp$Sig)

# FALSE TRUE
# FALSE  6004 3147
# TRUE    208 4990
write.table(comp, "../gtex/gtex_gains_comp_lead_for_mashr.txt", sep="\t", quote=F)

# reduce to significant results from GAinS
comp <- subset(comp, Sig == TRUE)
cor.test(comp$beta/comp$se, comp$slope/comp$slope_se)

pdf("S8_GTEx.pdf")
ggplot(comp, aes(beta/se, slope/slope_se)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("GAinS Z score") + ylab("GTEx Z score") +
  annotate("text", label="8,137 lead SNP-gene pairs tested in GTEx\n4,990 (61.3%) significant\nPearson's r=0.923", 
           x=-50, y=40) +
  xlim(c(-80, 80)) + ylim(c(-64, 64))
dev.off()