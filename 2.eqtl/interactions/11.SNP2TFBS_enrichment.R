################################################################################
#
# 11. Transcription factor identification for SRS interaction QTL
#
################################################################################

library(data.table)
library(tidyr)
library(dplyr)
library(tibble)

library(ggplot2)
library(ggrepel)
library(ggvenn)
library(pheatmap)
library(lineup)
library(RColorBrewer)
library(gridExtra)

library(JASPAR2014)
library(JASPAR2024)
library(TFBSTools)

options(stringsAsFactors = FALSE)

# Get info on interactions: background should be all tested for that variable
srs.int <- readRDS("srs_int_results_cond_incl.rds")
srs.int <- srs.int[complete.cases(srs.int), ]
srs.int$FDR <- p.adjust(srs.int$Interaction_pval, method="fdr")
srs.int$Type <- ifelse(srs.int$eQTL_beta*srs.int$Interaction_beta > 0,
                       "Magnifier", "Dampener")
srs.int$Pair <- paste0(srs.int$Gene, "_", srs.int$SNP)

# significant subset
srs.sig <- subset(srs.int, FDR < 0.05)

# LD info
ld.info <- read.delim("conditional_snps.80r2.tags.tsv")
ld.info$TAGS[which(ld.info$TAGS == "NONE")] <- ld.info$SNP[which(ld.info$TAGS == "NONE")]

ld.key <- separate_rows(ld.info[, c("SNP", "TAGS")], TAGS, sep="\\|")
ld.key$TAGS[which(ld.key$TAGS == "NONE")] <- ld.key$SNP[which(ld.key$TAGS == "NONE")]
# add the eSNP as a tag SNP for itself
ld.key <- rbind(ld.key, data.frame("SNP"=ld.info$SNP,
                                   "TAGS"=ld.info$SNP))

# TFBS info - input all eSNPS and LD partners to SNP2TFBS
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
opts[["all_versions"]] <- FALSE
opts[["tax_group"]] <- "vertebrates"
PFMatrixList <- getMatrixSet(JASPAR2014, opts)
# 205 vertebrates i.e. matches SNP2TFBS when all_versions set to FALSE

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

human.jaspar2014 <- human.jaspar2014[grepl("9606", human.jaspar2014$Species), ]
# write.table(human.jaspar2014, "JASPAR_human_TFBS.txt", sep="\t")
tf.matches <- subset(tf.matches, tf %in% human.jaspar2014$Name)

# JASPAR 2024 filter removes 3 more:
to.rm.24 <- c("MZF1_5-13", "BRCA1", "SMAD2_SMAD3_SMAD4")
tf.matches <- subset(tf.matches, !(tf %in% to.rm.24))

################################
# Test for enrichment of interaction eSNPs for each TF
ints <- srs.int
int.snps <- as.character(ints$SNP[which(ints$FDR < 0.05)])
sig.egenes <- as.character(ints$Pair[which(ints$FDR < 0.05)])  
intstested <- as.character(ints$SNP) # background

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

# Total number of overlaps
num.tf.hits <- tf.matches %>% group_by(Pair) %>%
  summarise(no_rows = length(Pair))
num.tf.hits$FDR <- ints$FDR[match(num.tf.hits$Pair, ints$Pair)]
wilcox.test(num.tf.hits$no_rows ~ num.tf.hits$FDR<0.05)
aggregate(num.tf.hits$no_rows, by=list(num.tf.hits$FDR<0.05), median)


# Enrichment on eQTL signal level
egenes <- snp.gene.pairs$Pair
tf.cont.table <- data.frame("eGene"=egenes,
                            "lead_snp"=snp.gene.pairs$SNP,
                            "Interaction"=egenes %in% sig.egenes)

# For all TFs, test for enrichment for interactions
tfs <- unique(tf.matches$tf)
results.tf <- data.frame("tf"=tfs,
                         "int"=NA,
                         "int_or"=NA,
                         "ntarget"=NA,
                         "nbackground"=NA)

for(i in 1:length(tfs)){
  tf.i <- subset(tf.matches, tf == tfs[i])
  tf.cont.table$tfi <- tf.cont.table$lead_snp %in% tf.i$lead_snp
  # If there is at least 1 overlap do the test
  if(sum(as.numeric(tf.cont.table$tfi)) >= 1){
      results.tf$int[i] <- (fisher.test(table(tf.cont.table$tfi, 
                                              tf.cont.table$Interaction, 
                                              useNA = "ifany"), alternative = "greater"))$p.value
      results.tf$int_or[i] <- (fisher.test(table(tf.cont.table$tfi, 
                                                 tf.cont.table$Interaction, 
                                                 useNA = "ifany"), alternative = "greater"))$estimate
      results.tf$ntarget[i] <- nrow(subset(tf.cont.table, tfi == TRUE & Interaction == TRUE))
      results.tf$nbackground[i] <- nrow(subset(tf.cont.table, tfi == TRUE))
    }
}

results.tf <- results.tf[complete.cases(results.tf), ]
results.tf$int_fdr <- p.adjust(results.tf$int, method="fdr")
  
results.sig <- subset(results.tf, int_fdr < 0.05)
results.sig <- results.sig[order(results.sig$int_or, decreasing = T), ]
results.sig$tf <- factor(results.sig$tf, levels=unique(results.sig$tf))
    
gg <- ggplot(results.sig, aes(tf, int_or, colour=-log10(int_fdr))) +
    geom_point() +
    theme_bw() + scale_color_gradient(low="gold", high="darkblue",
                                      limits=c(0, max(-log10(results.sig$int_fdr)))) +
    ylab("Fold enrichment") +
    xlab("Transcription factor") +
    scale_x_discrete(guide = guide_axis(angle = 45))

write.table(results.tf, "SRS_ints_TFBS_enrichment.txt", sep="\t")
