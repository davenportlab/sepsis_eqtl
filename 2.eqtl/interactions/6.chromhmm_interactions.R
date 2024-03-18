################################################################################
#
# 8. Interaction QTL - ChromHMM overlap
#
################################################################################

# Aim: to investigate whether interaction eSNPs are enriched in different types 
# of chromatin across cell types using ChromHMM annotations from Roadmap for 
# primary peripheral blood cells. 

library(XGR)
library(tidyverse)
library(ggupset)
library(ggplot2)
library(ggvenn)

RData.location <- "http://galahad.well.ox.ac.uk/bigdata"
options(stringsAsFactors = FALSE)

# Download and extract ChromHMM info for primary blood cells
info <- xRDataLoader('EpigenomeAtlas_15Segments_info',
                     RData.location=RData.location)
GR.annotations <- paste0('EpigenomeAtlas_15Segments_', names(info))
names(GR.annotations) <- info
GR.annotations <- GR.annotations[c(29:30, 32, 34:35, 37:40, 43:48)]

# Get info on interacions
int.summary <- read.delim("int_summary.txt")
int.summary$Variable_Type <- paste0(int.summary$Variable, "_", int.summary$Type)
int.summary$pair <- paste0(int.summary$Gene, "_", int.summary$SNP)

# Get background set (all eQTL tested for interactions)
int.tested <- readRDS("cellprops_int_results_cond_incl.rds")
ints <- as.character(unique(int.summary$SNP))
intstested <- as.character(unique(int.tested$SNP))

# Get b19 positions for those SNPs
b19pos.ints <- xSNPlocations(unique(ints), GR.SNP = "dbSNP_Common")
b19pos.intstested <- xSNPlocations(unique(intstested), GR.SNP = "dbSNP_Common")

# All interaction QTL vs no interaction QTL
ls_df <- lapply(1:length(GR.annotations), function(i){
  GR.annotation <- GR.annotations[i]
  message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                  as.character(Sys.time())), appendLF=T)
  df <- xGRviaGenomicAnno(data.file=b19pos.ints, 
                          background.file = b19pos.intstested,
                          format.file="GRanges",
                          p.tail="one-tail",
                          GR.annotation=GR.annotation, 
                          RData.location=RData.location, verbose=F)
  df$group <- names(GR.annotation)
  return(df)
})
df <- do.call(rbind, ls_df)
df$fc[which(df$adjp > 0.05)] <- NA
for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
df$name <- gsub(" from peripheral blood", "", df$name)
ggplot(df, aes(name, group)) + 
  geom_tile(aes(fill=log2(fc))) + 
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("Interactions vs non-interactions")

# Cell proportion interactions vs no cell proportion interaction

# Neutrophil magnifiers
neut.mag <- int.summary$SNP[int.summary$Variable_Type == "Neutrophil_Magnifier"]
b19pos.neutmag <- xSNPlocations(unique(neut.mag), GR.SNP = "dbSNP_Common")

ls_df <- lapply(1:length(GR.annotations), function(i){
  GR.annotation <- GR.annotations[i]
  message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                  as.character(Sys.time())), appendLF=T)
  df <- xGRviaGenomicAnno(data.file=b19pos.neutmag, 
                          background.file = b19pos.intstested,
                          format.file="GRanges",
                          p.tail="one-tail",
                          GR.annotation=GR.annotation, 
                          RData.location=RData.location, verbose=F)
  df$group <- names(GR.annotation)
  return(df)
})
df <- do.call(rbind, ls_df)
df$fc[which(df$adjp > 0.05)] <- NA
for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
df$name <- gsub(" from peripheral blood", "", df$name)
ggplot(df, aes(name, group)) + 
  geom_tile(aes(fill=log2(fc))) + 
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("Neutrophil magnifiers vs all tested")

# lymphocytes
lymph.mag <- int.summary$SNP[int.summary$Variable_Type == "Lymphocyte_Magnifier"]
b19pos.lymphmag <- xSNPlocations(unique(lymph.mag), GR.SNP = "dbSNP_Common")

ls_df <- lapply(1:length(GR.annotations), function(i){
  GR.annotation <- GR.annotations[i]
  message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                  as.character(Sys.time())), appendLF=T)
  df <- xGRviaGenomicAnno(data.file=b19pos.lymphmag, 
                          background.file = b19pos.intstested,
                          format.file="GRanges",
                          p.tail="one-tail",
                          GR.annotation=GR.annotation, 
                          RData.location=RData.location, verbose=F)
  df$group <- names(GR.annotation)
  return(df)
})
df <- do.call(rbind, ls_df)
df$fc[which(df$adjp > 0.05)] <- NA
for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
df$name <- gsub(" from peripheral blood", "", df$name)
ggplot(df, aes(name, group)) + 
  geom_tile(aes(fill=log2(fc))) + 
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("Lymphocyte magnifiers vs all tested")

# monocytes
mono.mag <- int.summary$SNP[int.summary$Variable_Type == "Monocyte_Magnifier"]
b19pos.monomag <- xSNPlocations(unique(mono.mag), GR.SNP = "dbSNP_Common")

ls_df <- lapply(1:length(GR.annotations), function(i){
  GR.annotation <- GR.annotations[i]
  message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                  as.character(Sys.time())), appendLF=T)
  df <- xGRviaGenomicAnno(data.file=b19pos.monomag, 
                          background.file = b19pos.intstested,
                          format.file="GRanges",
                          p.tail="one-tail",
                          GR.annotation=GR.annotation, 
                          RData.location=RData.location, verbose=F)
  df$group <- names(GR.annotation)
  return(df)
})
df <- do.call(rbind, ls_df)
df$fc[which(df$adjp > 0.05)] <- NA
for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
df$name <- gsub(" from peripheral blood", "", df$name)
ggplot(df, aes(name, group)) + 
  geom_tile(aes(fill=log2(fc))) + 
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("Monocyte magnifiers vs all tested")

# SRS interactions
srs <- int.summary$SNP[int.summary$Variable == "SRS"]
b19pos.srs <- xSNPlocations(unique(srs), GR.SNP = "dbSNP_Common")

ls_df <- lapply(1:length(GR.annotations), function(i){
  GR.annotation <- GR.annotations[i]
  message(sprintf("Analysing '%s' (%s) ...", names(GR.annotation),
                  as.character(Sys.time())), appendLF=T)
  df <- xGRviaGenomicAnno(data.file=b19pos.srs, 
                          background.file = b19pos.intstested,
                          format.file="GRanges",
                          p.tail="one-tail",
                          GR.annotation=GR.annotation, 
                          RData.location=RData.location, verbose=F)
  df$group <- names(GR.annotation)
  return(df)
})
df <- do.call(rbind, ls_df)
df$fc[which(df$adjp > 0.05)] <- NA
for(i in 1:9){df$name <- gsub(paste0("E", i, "_"), paste0("E0", i, "_"), df$name)}
df$name <- gsub(" from peripheral blood", "", df$name)

pdf("Fig_2D_heatmap_chromhmm_srs_ints.pdf", useDingbats = F, height=5, width=10)
ggplot(df, aes(name, group)) + 
  geom_tile(aes(fill=log2(fc))) + 
  scale_fill_gradient2(low="red", mid="white", high="blue") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("SRS interactions vs all tested")
dev.off()
