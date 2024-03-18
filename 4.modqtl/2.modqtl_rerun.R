################################################################################
#
# 2. Module QTL mapping
#
################################################################################

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

library(data.table)
library(lme4)
library(parallel)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
genotypes.file <- "eigengene_sva_genotypes.raw"
design.matrix.file <- "mapping_data.csv"
ME <- args[1] # 1:106
ME <- paste0("ME_", ME, "_1")
output_dir <- "../modqtl_rerun/"
output <- paste0(output_dir, ME, ".tsv")

#----------------------------------------------------------
# Load mQTL Data
#----------------------------------------------------------
# Design Matrix
design.matrix <- read.csv(design.matrix.file)

# Genotype Matrix
genotypes <- fread(genotypes.file, sep=" ", drop=2:6)

# Clean Genotype Matrix
patient.sample.match <- match(design.matrix$GAinS.ID, genotypes$FID)
genotypes <- genotypes[patient.sample.match, ]

colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL

# Filter to just the conditional cis-eQTL signal SNPs
mqtl.snp.table <- read.csv("../../nikhil/expression/eigengene_sva/mqtl_snp_table.csv")
cond.snps <- unique(mqtl.snp.table$snps[which(mqtl.snp.table$source == "Conditional cis-eQTL SNP")])
cond.snps <- cond.snps[cond.snps %in% colnames(genotypes)]

genotypes <- genotypes %>% select(cond.snps)
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- design.matrix$Sample.ID

genotype.ids <- colnames(genotypes)

#----------------------------------------------------------
# Linear Mixed Model
#----------------------------------------------------------

# Create a full design matrix with all genotypes
genotypes <- cbind(design.matrix, genotypes)

all.vars <- colnames(design.matrix)
eigens <- colnames(design.matrix)[grepl("^ME", colnames(design.matrix))]
covs <- setdiff(setdiff(all.vars, eigens), c("Sample.ID", "GAinS.ID"))

results <- rbindlist(mclapply(genotype.ids, function(snp) {
  # for ME == 1
  
  variant.design <- genotypes[,c(ME, snp, covs, "GAinS.ID")]
  
  f.null <- as.formula(paste0(ME, "~", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
  model.null <- lmer(f.null, data=variant.design, REML=FALSE)
  
  f.alt <- as.formula(paste0(ME, "~`", snp , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
  model.test <- lmer(f.alt, data=variant.design, REML=FALSE)
  
  if (!all(complete.cases(variant.design[, snp]))) {
    model.null <- update(model.null, subset=complete.cases(variant.design[, snp]))
    model.test <- update(model.test, subset=complete.cases(variant.design[, snp]))
  }
  
  data.frame(matrix(
    data=c(
      summary(model.test)$coefficients[2, ],
      anova(model.null, model.test)["model.test", "Pr(>Chisq)"]
    ),
    nrow=1, ncol=4
  ))
}))

#----------------------------------------------------------
# Write Output
#----------------------------------------------------------

# Results has the following columns:
#   SNP
#   Beta
#   Standard Error
#   t Value
#   P-Value from ANOVA

results <- cbind(genotype.ids, results)

write.table(results, output, quote=F, row.names=F, col.names=F, sep="\t")
