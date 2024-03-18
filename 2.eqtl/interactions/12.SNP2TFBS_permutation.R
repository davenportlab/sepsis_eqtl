# SNP2TFBS SRS permutation analysis

library(dplyr)
library(tidyr)
library(TFBSTools)
library(JASPAR2014)
library(ggplot2)

setwd("Results/RNASeq/interactions/")

# Get info on interactions: background should be all tested for that variable
srs.int <- readRDS("srs_int_results_cond_incl.rds")
srs.int <- srs.int[complete.cases(srs.int), ]
srs.int$FDR <- p.adjust(srs.int$Interaction_pval, method="fdr")
srs.int$Type <- ifelse(srs.int$eQTL_beta*srs.int$Interaction_beta > 0,
                       "Magnifier", "Dampener")
srs.int$Pair <- paste0(srs.int$Gene, "_", srs.int$SNP)

# sig subset
srs.sig <- subset(srs.int, FDR < 0.05)

# LD info
ld.info <- read.delim("conditional_snps.80r2.tags.tsv")
ld.info$TAGS[which(ld.info$TAGS == "NONE")] <- ld.info$SNP[which(ld.info$TAGS == "NONE")]

ld.key <- separate_rows(ld.info[, c("SNP", "TAGS")], TAGS, sep="\\|")
ld.key$TAGS[which(ld.key$TAGS == "NONE")] <- ld.key$SNP[which(ld.key$TAGS == "NONE")]
# add the eSNP as a tag SNP for itself
ld.key <- rbind(ld.key, data.frame("SNP"=ld.info$SNP,
                                   "TAGS"=ld.info$SNP))

# TFBS info
tf.matches <- read.delim("all_cond_snps_ld_tf_matches.txt", header=F)
colnames(tf.matches) <- c("snp", "chr", "start", "end", "ref", "alt", "nmatches",
                          "tf", "scorediff")
tf.matches$nmatches <- gsub("MATCH=", "", tf.matches$nmatches)
tf.matches$tf <- gsub("TF=", "", tf.matches$tf)
tf.matches <- separate_rows(tf.matches, tf, scorediff, sep=",")
tf.matches$scorediff <- as.numeric(gsub("ScoreDiff=", "", tf.matches$scorediff))
tf.matches$tf <- toupper(tf.matches$tf)

# Get list of Jaspar 2014 motifs and restrict to human
opts <- list()
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- F
opts[["tax_group"]] <- "vertebrates"
PFMatrixList <- getMatrixSet(JASPAR2014, opts)
# 205 vertebrates i.e. matches SNP2TFBS

human.jaspar2014 <- data.frame("ID"=names(PFMatrixList),
                               "Name"=NA,
                               "Species"=NA)
for(i in 1:length(PFMatrixList)){
  human.jaspar2014$ID[i] <- PFMatrixList@listData[[i]]@ID
  human.jaspar2014$Name[i] <- PFMatrixList@listData[[i]]@name
  spp <- names(PFMatrixList@listData[[i]]@tags$species)
  human.jaspar2014$Species[i] <- ifelse(is.null(spp), "NA", spp)
}
human.jaspar2014$Name <- gsub("::", "_", human.jaspar2014$Name)
human.jaspar2014$Name <- toupper(human.jaspar2014$Name)
nonhuman.jaspar2014 <- human.jaspar2014[!grepl("9606", human.jaspar2014$Species), ]
human.jaspar2014 <- human.jaspar2014[grepl("9606", human.jaspar2014$Species), ]

tf.matches <- subset(tf.matches, tf %in% human.jaspar2014$Name)

# JASPAR 2024 filter removes 3 more:
to.rm.24 <- c("MZF1_5-13", "BRCA1", "SMAD2_SMAD3_SMAD4")
tf.matches <- subset(tf.matches, !(tf %in% to.rm.24))

# restrict to cases where there is a change in the motif score?
# tf.matches <- subset(tf.matches, abs(scorediff)>0)

ints <- srs.int
int.snps <- as.character(ints$SNP[which(ints$FDR < 0.05)])
sig.egenes <- as.character(ints$Pair[which(ints$FDR < 0.05)])  
intstested <- as.character(ints$SNP)
  
ld.info.snps <- ld.info[match(int.snps, ld.info$SNP), ]
int.snps.ld <- unique(c(int.snps, unlist(strsplit(ld.info.snps$TAGS, split="\\|"))))
  
all.snps <- ld.info[match(intstested, ld.info$SNP), ]
all.snps.ld <- unique(c(intstested, unlist(strsplit(all.snps$TAGS, split = "\\|"))))
  
int.snps <- int.snps.ld
intstested <- all.snps.ld

tf.matches <- subset(tf.matches, snp %in% intstested)

# get lead sepsis SNP to link LD SNPs to eGenes
snp.gene.pairs <- ints[, c(1:2, 13:15)]
ld.key <- merge(ld.key, snp.gene.pairs, by.x="SNP", by.y="SNP")
tf.matches <- merge(tf.matches, ld.key, by.x="snp", by.y="TAGS")
colnames(tf.matches)[10] <- "lead_snp"

# Enrichment on eQTL signal level
egenes <- snp.gene.pairs$Pair
tf.cont.table <- data.frame("eGene"=egenes,
                            "lead_snp"=snp.gene.pairs$SNP,
                            "Interaction"=egenes %in% sig.egenes)
  
# For all TFs, test for enrichment for interactions
tfs <- unique(tf.matches$tf)
results.tf <- data.frame("tf"=tfs,
                           "int"=NA)
perm.results <- data.frame("perm"=1:1001,
                           "n.sig"=NA)

set.seed("23112022")
for(p in 1:1001){
  print(p)
  if(p>1){
    # permute interaction label in contingency table
    tf.cont.table$Interaction <- sample(tf.cont.table$Interaction, 
                                        size=nrow(tf.cont.table), replace = F)
  }
  # test for TFBS enrichment in interaction eQTL
  for(i in 1:length(tfs)){
    tf.i <- subset(tf.matches, tf == tfs[i])
    tf.cont.table$tfi <- tf.cont.table$lead_snp %in% tf.i$lead_snp
  if(sum(as.numeric(tf.cont.table$tfi)) >= 1){
    results.tf$int[i] <- (fisher.test(table(tf.cont.table$tfi, 
                                            tf.cont.table$Interaction, 
                                            useNA = "ifany"), alternative = "greater"))$p.value
    }
  }
  # report number of significant enrichments
  results.tf$int_fdr <- p.adjust(results.tf$int, method="fdr")
  results.sig <- subset(results.tf, int_fdr < 0.05)
  perm.results[p, 2] <- nrow(results.sig)
}

# pdf("TFBS_permutation_histogram_updated.pdf")
ggplot(perm.results, aes(n.sig)) +
  geom_histogram(bins=60, colour="black", fill="grey") +
  theme_bw() +
  coord_cartesian(expand=FALSE) +
  # scale_x_continuous(limits = c(0, 170)) +
  geom_vline(xintercept=56, col="red", lty=2) +
  xlab("Number of TF motifs significantly enriched in SRS interaction QTL") +
  ylab("Number of permutations")
# dev.off()
