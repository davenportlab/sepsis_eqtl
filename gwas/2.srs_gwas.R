##############################################################################

# GWAS for SRS1 ever vs never

##############################################################################

# AIM: to test genome-wide for genetic variants associated with a patient's SRS status
# Because SRS can change over time, and we have a different number of samples for different patients,
# we encode SRS as SRS1 ever vs never

##############################################################################

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(qqman)
library(reshape2)
library(tidyverse)

options(stringsAsFactors = FALSE)

gwsig <- 5e-8
gwsugg <- 5e-5

##############################################################################
# Match IDs between genotyping and samples with SRS assignments
fam <- read.table("All_genotyping_merged_filtered_b38_refiltered_rsID.fam")
key <- read.delim("GEx_geno_data_availability.txt")
srs <- read.delim("full-gains-SRS-predictions_mNN-RF.tsv")
srs$Patient <- key$GenotypingPatientID[match(srs$Patient, key$SRSPatientID)]

table(srs$Assay)
# Microarray       qPCR    RNA-seq 
# 676        115        864

by(srs, srs$Assay, function(x) length(unique(x$Patient))) # number of patients with genotyping in each
# srs$Assay: Microarray
# [1] 514
#   srs$Assay: qPCR
# [1] 107
#   srs$Assay: RNA-seq
# [1] 667

# Make SRS status SRS1-centric i.e. SRS1=case, SRS2 or 3=control
srs$SRS1 <- ifelse(srs$SRS == "SRS1", "1", "0")

# Are there any repeat samples discrepant between assays?
# samples with SRS assignments from different assays
samples.multi.assays <- names((rowSums(table(srs$Sample_id, 
                                             srs$Assay)>0)))[rowSums(table(srs$Sample_id, 
                                                                           srs$Assay)>0) > 1]
srs.multi <- srs[which(srs$Sample_id %in% samples.multi.assays), ]
# Look for cases where the same sample has different assignments
samples.srs.diff <- names(which(table(srs.multi$Sample_id, srs.multi$SRS1)[,1] > 0 &
                                  table(srs.multi$Sample_id, srs.multi$SRS1)[,2] > 0))
# only 10 samples with discrepant assignments across platforms, most involve qPCR replicate

# Use RNA-seq then microarray then qPCR assignments
srs <- rbind(srs[which(srs$Assay == "RNA-seq"), ],
             srs[which(srs$Assay == "Microarray"), ],
             srs[which(srs$Assay == "qPCR"), ])
srs <- srs[!duplicated(srs$Sample_id), ] # 1412 unique samples from 1044 patients

table(srs$Assay)
# Microarray       qPCR    RNA-seq 
# 541          7        864

by(srs, srs$Assay, function(x) length(unique(x$Patient))) # number of patients in each (patients can be duplicated too)
# 382, 7, 667

# how many patients have multiple assays used
multiassay <- data.frame(table(srs$Patient, srs$Assay))
multiassay <- subset(multiassay, Freq != 0)
table(duplicated(multiassay$Var1)) # only 7

# Number of patients with SRS and genotyping data
srs$Genotyping <- srs$Patient %in% fam$V1 # 997
table(srs$Assay[which(srs$Genotyping == TRUE)])
# Microarray       qPCR    RNA-seq 
# 506          5        823 

###############################
# Phenotype: SRS1ever vs never
###############################

srs1ever <- data.frame(table(srs$Patient, srs$SRS1)) 
# for each patient, how many timepoints SRS1 (1) or not (0)
srs1ever <- subset(srs1ever, Var2 == 1 & Freq > 0)
srs$SRS1ever <- srs$Patient %in% srs1ever$Var1

# add phenotype to fam file
fam$V6 <- srs$SRS1ever[match(fam$V1, srs$Patient)]
table(fam$V6, useNA = "ifany")
# FALSE  TRUE  <NA> 
#   557   440   171

pheno <- fam[, c(1:2, 6)] # GCTA case/control format 1/0
pheno$V6 <- as.integer(fam$V6)
# https://yanglab.westlake.edu.cn/software/gcta/#Tutorial
 write.table(pheno, "All_genotyping_SRS1ever_gcta.phen",
            row.names = F, col.names = F, quote=F, sep="\t")
pheno$V6 <- pheno$V6+1
pheno$V6[which(is.na(fam$V6))] <- "-9"
table(pheno$V6) # case=2, control=1, NA=-9 (Plink format)
# https://www.cog-genomics.org/plink/1.9/input#pheno
write.table(pheno, "All_genotyping_SRS1ever_plink.txt",
            row.names = F, col.names = F, quote=F, sep="\t")

##############################################################################
# Plink: logistic/linear regression
##############################################################################

# Phenotype: SRS1ever vs never
# Covariates: age, age2, sex, genotyping PCs
# plink --bfile All_genotyping_merged_filtered_b38_refiltered_rsID --logistic hide-covar --allow-extra-chr --maf 0.01 --geno 0.02 --chr 1-22 --out srs_gwas_logistic_sex_age_pcs_srs1ever --allow-no-sex --covar All_genotyping_covs_for_gwas.cov --adjust --pheno  All_genotyping_SRS1ever_plink.txt --ci 0.95

# Genomic inflation est. lambda (based on median chisq) = 1.01975

srs1ever.plink <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic",
                                    header=T)
table(srs1ever.plink$P < gwsig) # 0
table(srs1ever.plink$P < gwsugg) # 155

# png("Fig_S1B_SRS1ever_qqplot.png")
qq(srs1ever.plink$P)
# dev.off()

# ~/plink --bfile All_genotyping_merged_filtered_b38_refiltered_rsID --clump srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic --clump-r2 0.5 --clump-p1 0.00005 --clump-p2 0.00005 --clump-kb 250 --out srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic --allow-extra-chr --chr 1-22

srs1ever.plink.clumped <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic.clumped", header=T)
srs1ever.plink.clumped$SP2 <- gsub("\\(1)", "", srs1ever.plink.clumped$SP2)
srs1ever.plink.clumped$SP2 <- gsub("NONE", "", srs1ever.plink.clumped$SP2)

df <- srs1ever.plink %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=as.numeric(max(BP, na.rm=T))) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(srs1ever.plink, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf <- df %>% group_by(CHR) %>% 
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

p1 <- ggplot(df, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  rasterize(geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3), dpi=100) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # Add highlighted points
  geom_point(data=subset(df, SNP %in% srs1ever.plink.clumped$SNP), 
             color="orange", size=2) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  # scale_y_continuous(expand = c(0, 2) ) +
  ylim(c(0, -log10(1e-8))) +
  xlab("SNP position") +
  ylab("-log10(P value)") +
  # Custom the theme:
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = -log10(5e-5), lty=2) +
  geom_hline(yintercept = -log10(5e-8), lty=2)

pdf("Fig_1B_SRS1ever_plink_manhattan.pdf", width = 14, height=5)
print(p1)
dev.off()
