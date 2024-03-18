#----------------------------------------------------------
# Author: Nikhil Milind
# Created: 25 November 2021
#----------------------------------------------------------

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

if (!require("lme4")) {
    install.packages("lme4")
}

library(tidyverse)
library(data.table)
library(lme4)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
genotypes.file <- args[1]
design.matrix.file <- args[2] # covariates
locus <- args[3]
snps.file <- args[4]
conditional.snps.file <- args[5]

#----------------------------------------------------------
# Load eQTL Data
#----------------------------------------------------------

# Design Matrix
design.matrix <- read.csv(design.matrix.file)

# Genotype Matrix
genotypes <- fread(genotypes.file, sep=" ", drop=2:6)

# QTL SNPs
snps <- read.table(snps.file, header=T)
conditional.snps <- read.table(conditional.snps.file)[, 1]

# Clean Genotype Matrix
patient.sample.match <- match(design.matrix$GAinS.ID, genotypes$FID)
genotypes <- genotypes[patient.sample.match,]
colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- design.matrix$Sample.ID

# Gene Expression Matrix
gene.exp <- fread("logcpm_864_20412_hla.txt") %>% as.data.frame()
rownames(gene.exp) <- gene.exp[, 1]
gene.exp <- as.data.frame(t(gene.exp[, -1]))
gene.exp <- gene.exp[design.matrix$Sample.ID,]
design.matrix <- cbind(gene.exp[, locus], design.matrix)
colnames(design.matrix)[1] <- locus

#----------------------------------------------------------
# Linear Mixed Model
#----------------------------------------------------------

# Create a full design matrix with all genotypes
data.mtx <- cbind(design.matrix, genotypes[, intersect(snps$snps, colnames(genotypes))])

for (signal in 1:length(conditional.snps)) {

    conditional.snp.set <- setdiff(conditional.snps, conditional.snps[signal])

    snp.set <- setdiff(intersect(snps$snps, colnames(genotypes)), conditional.snp.set)

    all.vars <- colnames(design.matrix)
    covs <- setdiff(setdiff(all.vars, locus), c("Sample.ID", "GAinS.ID"))

    results <- rbindlist(mclapply(snp.set, function(snp) {

        variant.design <- data.mtx[,c(locus, snp, conditional.snp.set, covs, "GAinS.ID")]

        conditional.snps.str <- paste0(paste0("`", conditional.snp.set, "`"), collapse="+")

        f.null <- as.formula(paste0(locus, "~", conditional.snps.str, "+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
        model.null <- lmer(f.null, data=variant.design, REML=FALSE)

        f.alt <- as.formula(paste0(locus, "~`", snp , "`+", conditional.snps.str, "+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
        model.test <- lmer(f.alt, data=variant.design, REML=FALSE)

        if (!all(complete.cases(variant.design[, snp]))) {
            model.null <- update(model.null, subset=complete.cases(variant.design[, snp]))
            model.test <- update(model.test, subset=complete.cases(variant.design[, snp]))
        }

        data.frame(matrix(
            data=c(
                summary(model.test)$coefficients[snp, ],
                anova(model.null, model.test)["model.test", "Pr(>Chisq)"]
            ),
            nrow=1, ncol=4
        ))
    }, mc.cores=16))

    results <- cbind(snp.set, results) %>%
        as.data.frame() %>%
        dplyr::select(snps=1, beta=2, se=3, t=4, pval=5) %>%
        merge(., snps, by="snps") %>%
        dplyr::mutate(gene=locus, signal=signal) %>%
        dplyr::select(gene, signal, chr, snps, SNPpos, beta, se, pval) %>%
        dplyr::filter(!is.na(pval))

    write.table(results, paste0(locus, "-", signal, ".tsv"), quote=F, row.names=F, col.names=F, sep="\t")
}
