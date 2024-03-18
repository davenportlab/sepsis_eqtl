################################################################################
#
# 8. Module QTL-GWAS colocalisation
#
################################################################################# 

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(coloc)
library(susieR)

setwd("../modqtl_rerun/")

# Load ModQTL data
geno.bim <- fread("../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim")
colnames(geno.bim) <- c("chr", "snp", "cM", "pos", "minor_allele", "major_allele")

module.qtl.snps <- read.delim("mqtl_full_summary_statistics_snps.txt")
module.qtl.snps$pair <- paste0(module.qtl.snps$ME, "_", module.qtl.snps$SNP)
lead.mqtl <- read.csv("updated_lead_mqtl.csv")

module.qtl.sum <- do.call(rbind, lapply(list.files("./", pattern="*_1.tsv"), function(file.name) {
  fread(file.name) %>%
    as.data.frame() %>%
    dplyr::select(snp=1, beta=2, se=3, t=4, p=5) %>%
    dplyr::mutate(module.qtl=gsub(".tsv", "", file.name)) %>%
    dplyr::mutate(module=gsub("_1", "", module.qtl))
}))

module.qtl.sum <- subset(module.qtl.sum, module %in% lead.mqtl$me)
module.qtl.sum$pair <- paste0(module.qtl.sum$module, "_", module.qtl.sum$snp)
module.qtl.sum <- merge(module.qtl.sum, module.qtl.snps, by="pair")
module.qtl.sum$pair <- NULL
module.qtl.sum$module.qtl <- paste0(module.qtl.sum$module.qtl, "-", module.qtl.sum$QTL.ID)
module.qtl.sum$SNP <- NULL
module.qtl.sum$QTL.ID <- NULL

iupac.map <- matrix(
  c(
    "A", "M", "R", "W",
    "M", "C", "S", "Y",
    "R", "S", "G", "K",
    "W", "Y", "K", "T"
  ), 
  nrow=4
)
rownames(iupac.map) <- diag(iupac.map)
colnames(iupac.map) <- diag(iupac.map)

alleles.iupac <- Vectorize(function(x, y) {
  
  if (x %in% rownames(iupac.map) && y %in% colnames(iupac.map)) {
    return(iupac.map[x, y])
  } else {
    return(NA)
  }
})

module.qtl.sum <- module.qtl.sum %>% merge(., geno.bim, by="snp") %>%
  dplyr::mutate(IUPAC = alleles.iupac(minor_allele, major_allele))

mqtl.geno <- fread("../../nikhil/data/genotypes/eigengene_sva_ss_genotypes.raw", sep=" ", drop=2:6) %>%
  as.data.frame()

rownames(mqtl.geno) <- mqtl.geno$FID
mqtl.geno$FID <- NULL
colnames(mqtl.geno) <- gsub("_.*$", "", colnames(mqtl.geno))

# Limit to the ModQTL of interest
module.qtl.sum <- subset(module.qtl.sum, module == "ME_47")
module.qtl.sum$SNPID <- paste0(module.qtl.sum$chr, ":", module.qtl.sum$pos, "_", module.qtl.sum$major_allele, "_", module.qtl.sum$minor_allele)

# Load module eigengenes
eigengenes <- read.csv("../../nikhil/expression/gene_expression/eigengenes.multiple.csv", row.names=1)

# Load EBI SNPs that are Module QTL
ibd <- fread("../../nikhil/data/EBI_GWAS_Catalog/ibd_build37_59957_20161107_corrected.txt", sep="\t")
# ebi.studies <- fread("/lustre/scratch125/humgen/resources/NHGRI-EBI_GWAS_catalogue/harmonised_files_download_2022-07-19/28067908-GCST004131-EFO_0003767.h.tsv.gz", header=TRUE, quote="") %>%
#   as.data.frame()

ebi.me <- merge(module.qtl.sum, ebi.studies, by.x="snp", by.y="hm_rsid")
ibd.me <- merge(module.qtl.sum, ibd, by="snp")

ibd.snps <- unique(ibd.me$snp)

module.qtl.set <- module.qtl.sum %>%
  dplyr::filter(snp %in% ibd.snps)

module.qtl.set <- unique(module.qtl.set$module.qtl)

all.snps <- module.qtl.sum %>%
  dplyr::filter(module.qtl %in% module.qtl.set)

all.snps <- all.snps$snp

ibd <- ibd %>%
  dplyr::filter(snp %in% all.snps) %>%
  dplyr::filter(StdErr > 0) %>%
  dplyr::mutate(Allele1 = str_to_upper(Allele1), Allele2 = str_to_upper(Allele2)) %>%
  dplyr::mutate(IUPAC = alleles.iupac(Allele1, Allele2))

mqtl.locus.info = module.qtl.sum %>%
    dplyr::select(m.snp = snp, m.beta=beta, m.se=se, m.pos=pos, major_allele, minor_allele, IUPAC)
  
all.info = ibd %>%
    dplyr::filter(snp %in% mqtl.locus.info$m.snp) %>%
    merge(
      ., mqtl.locus.info, 
      by.x=c("snp", "IUPAC"),
      by.y=c("m.snp", "IUPAC")
    ) %>%
    dplyr::select(snp=snp, m.beta, m.se, pos=m.pos, g.beta=Effect, g.se=StdErr)

module.qtl = list()
module.qtl$beta = all.info$m.beta
module.qtl$varbeta = all.info$m.se^2
module.qtl$snp = all.info$snp
module.qtl$position = all.info$pos
module.qtl$type = "quant"
module.qtl$sdY = sd(eigengenes[, "ME_47_1"], na.rm=TRUE)
 
gwas.assoc = list()
gwas.assoc$beta = all.info$g.beta
gwas.assoc$varbeta = all.info$g.se^2
gwas.assoc$snp = all.info$snp
gwas.assoc$position = all.info$pos
gwas.assoc$type = "cc"

abf.res = coloc.abf(gwas.assoc, module.qtl)

data.frame(t(abf.res$summary)) %>%
  dplyr::mutate(Module.QTL = module.qtl.id) %>%
  do.call(rbind, .) %>%
  dplyr::mutate(PP3plusPP4 = PP.H3.abf + PP.H4.abf) %>%
  dplyr::mutate(COLOC.Factor = PP.H4.abf / PP3plusPP4) %>%
  dplyr::mutate(Colocalise = (PP3plusPP4 > 0.25) & (COLOC.Factor > 0.7))  
