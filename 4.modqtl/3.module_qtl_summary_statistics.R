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

library(data.table)
library(lme4)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
genotypes.file <- args[1]
design.matrix.file <- args[2]
qtl.snps.file <- args[3]
ME <- args[4]

#----------------------------------------------------------
# Load mQTL Data
#----------------------------------------------------------

# Design Matrix
design.matrix <- read.csv(design.matrix.file)

# Genotype Matrix
genotypes <- fread(genotypes.file, sep=" ", drop=2:6)

# QTL SNPs
qtl.snps <- read.table(qtl.snps.file, sep="\t", header=TRUE)
qtl.snps <- qtl.snps[qtl.snps$ME == ME, ]
qtl.snps$SNP <- gsub("\\:", "_", qtl.snps$SNP)

qtl.ids <- unique(qtl.snps$QTL.ID)

# Clean Genotype Matrix
patient.sample.match <- match(design.matrix$GAinS.ID, genotypes$FID)
genotypes <- genotypes[patient.sample.match,]
colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- design.matrix$Sample.ID

#----------------------------------------------------------
# Linear Mixed Model
#----------------------------------------------------------

for (qtl.id in qtl.ids) {

    qtl.snp.set <- qtl.snps$SNP[qtl.snps$QTL.ID == qtl.id]
    qtl.snp.set <- intersect(qtl.snp.set, colnames(genotypes))

    # Create a full design matrix with all genotypes
    data.mtx <- cbind(design.matrix, genotypes[, qtl.snp.set])

    all.vars <- colnames(design.matrix)
    eigens <- colnames(design.matrix)[grepl("^ME", colnames(design.matrix))]
    covs <- setdiff(setdiff(all.vars, eigens), c("Sample.ID", "GAinS.ID"))

    results <- rbindlist(mclapply(qtl.snp.set, function(snp) {

        ME.pc.1 = paste0(ME, "_1")

        variant.design <- data.mtx[,c(ME.pc.1, snp, covs, "GAinS.ID")]

        f.null <- as.formula(paste0(ME.pc.1, "~", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
        model.null <- lmer(f.null, data=variant.design, REML=FALSE)

        f.alt <- as.formula(paste0(ME.pc.1, "~`", snp , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
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
    }))

    results <- cbind(qtl.snp.set, results)

    write.table(results, paste0(ME, "-", qtl.id, ".tsv"), quote=F, row.names=F, col.names=F, sep="\t")
}
