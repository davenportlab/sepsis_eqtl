################################################################################
#
# 3. Comparison to eQTLGen
#
################################################################################

options(stringsAsFactors = FALSE)

library(data.table)
library(ggplot2)

# get sepsis lead snp-gene pairs also tested in eQTLGen and see if significant in each dataset

# Sepsis data
gains.eqtl <- read.delim("ciseqtl_eigenMT_corrected.txt")

# eQTLGen data
eqtlgen <- fread("eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")

# Reduce both to shared genes and SNPs
gains.eqtl$Pair <- paste0(gains.eqtl$gene, "_", gains.eqtl$snps)
eqtlgen$Pair <- paste0(eqtlgen$Gene, "_", eqtlgen$SNP)
gains.eqtl <- gains.eqtl[which(gains.eqtl$Pair %in% eqtlgen$Pair), ]
eqtlgen <- eqtlgen[match(gains.eqtl$Pair, eqtlgen$Pair), ]

# Add GAinS allele information to confirm same SNPs/minor allele assessed
gains.bim <- data.frame(fread("genotyping_for_rna-seq_eQTL.bim"))
m1 <- match(gains.eqtl$snps, gains.bim$V2)
gains.eqtl$Allele1 <- gains.bim$V5[m1] # minor
gains.eqtl$Allele2 <- gains.bim$V6[m1] # major

# Merge (already in same order from above match)
eqtlgen.results <- cbind(gains.eqtl, eqtlgen)

# Check that major and minor alleles are the same
eqtlgen.results$action <- "NA"
eqtlgen.results[which(eqtlgen.results$AssessedAllele == eqtlgen.results$Allele1),
        "action"] <- "same"

# Alleles where strand can't be distinguished
eqtlgen.results[which(eqtlgen.results$AssessedAllele=="A" & 
                eqtlgen.results$OtherAllele=="T" |
                eqtlgen.results$AssessedAllele=="T" & 
                eqtlgen.results$OtherAllele=="A" |
                eqtlgen.results$AssessedAllele=="C" & 
                eqtlgen.results$OtherAllele=="G" |
                eqtlgen.results$AssessedAllele=="G" & 
                eqtlgen.results$OtherAllele=="C"), 
        "action"] <-"strand"

# Same alleles but different strand
eqtlgen.results[which(eqtlgen.results$AssessedAllele=="A" & 
                eqtlgen.results$Allele1=="T" & 
                eqtlgen.results$action != "strand" |
                eqtlgen.results$AssessedAllele=="T" & 
                eqtlgen.results$Allele1=="A" & 
                eqtlgen.results$action != "strand" |
                eqtlgen.results$AssessedAllele=="C" & 
                eqtlgen.results$Allele1=="G" & 
                eqtlgen.results$action != "strand" |
                eqtlgen.results$AssessedAllele=="G" & 
                eqtlgen.results$Allele1=="C" & 
                eqtlgen.results$action != "strand"),
        "action"] <-"same"

# Flip effects where opposite
eqtlgen.results[which(eqtlgen.results$action == "NA"), "action"] <- "flip"
eqtlgen.results$eQTL_beta.2 <- eqtlgen.results$eQTL_beta
eqtlgen.results[which(eqtlgen.results$action == "flip"), "eQTL_beta.2"] <- eqtlgen.results[which(eqtlgen.results$action == "flip"), "eQTL_beta.2"] * (-1)

# Remove SNPs which can't be distinguished because of strand
eqtlgen.results$eQTL_beta.2[which(eqtlgen.results$action == "strand")] <- NA
eqtlgen.results <- eqtlgen.results[complete.cases(eqtlgen.results), ]

# Mark if lead GAinS association is significant
eqtlgen.results$Sig <- eqtlgen.results$pvalue < eqtlgen.results$threshold
table(eqtlgen.results$Sig, eqtlgen.results$BonferroniP < 0.05)
write.table(eqtlgen.results, "gains_lead_vs_eqtlgen.results_stats.txt", sep="\t", quote=F)

# Reduce to significant sepsis results
eqtlgen.results <- subset(eqtlgen.results, Sig == TRUE)
cor.test(eqtlgen.results$eQTL_beta.2/eqtlgen.results$se, eqtlgen.results$Zscore)

pdf("S8_eQTLgen.pdf")
ggplot(eqtlgen.results, aes(eQTL_beta.2/se, Zscore)) +
  geom_point() +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("GAinS Z score") + ylab("eQTLgen Z score") +
  annotate("text", label="7,751 lead SNP-gene pairs tested in eQTLgen\n6,464 (83%) significant\nPearson's r=0.822", 
           x=-50, y=170) +
  xlim(c(-80, 80)) + ylim(c(-195, 195))
dev.off()
