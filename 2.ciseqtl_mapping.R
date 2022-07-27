################################################################################
#
# 2. cis-eQTL mapping: for one gene identified by job index
#
################################################################################

.libPaths("~/R/x86_64-pc-linux-gnu-library/3.6/")
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
# first argument is the gene for which to map eQTL
jobindex <- as.numeric(args[1])
# second argument is the number of PEER factors to include in the model
n.peer <- as.numeric(args[2])

library(lme4)
library(data.table)

################################################################################

# read in the list of SNP-gene pairs
pairs <- fread("../data/gene_snp_pairs_cis_rnaseq.txt", stringsAsFactors = F)

# use the job index to identify the gene for which to map eqtl
genes <- unique(pairs$Gene)
gene <- genes[jobindex]

# extract the list of SNPs to test for this gene
pairs.eqtl <- pairs[which(pairs$Gene == gene), ]
pairs.eqtl$SNP <- gsub(":", ".", pairs.eqtl$SNP)
irange <- 1:nrow(pairs.eqtl)
chr <- pairs.eqtl[1, 3]

# check if results file exists
done <- file.exists(paste0("../cisresults/chr", chr, "/", gene, ".rds"))
print(done)
if(done == FALSE){
  
# read in the data
load(paste("../data/eqtl_files_", chr, ".rda", sep=""))

results <- rbindlist(lapply(irange, function(i){
  # for each SNP, fit a LMM for gene expression including genotyping and compare
  # to a null model including SRS, diagnosis, cell proportions, genotyping PCs and PEER factors
  model.null <- lmer(as.numeric(exp[gene, ]) ~
                       covs +
                       peer.factors[, 1:n.peer] +
                       (1|GAinSID),
                     REML=FALSE) 
  
  model.test <- lmer(as.numeric(exp[gene, ]) ~
                       as.numeric(geno[, as.character(pairs.eqtl[i, 2])]) +
                       covs + 
                       peer.factors[, 1:n.peer] +
                       (1|GAinSID),
                     REML=FALSE)
  
  #If there are no missing genotypes
  if (all(complete.cases(geno[, as.character(pairs.eqtl[i, 2])]))){
    data.frame(matrix(data=c(summary(model.test)$coefficients[2, ],
                             anova(model.null, 
                                   model.test)$'Pr(>Chisq)'[2]), nrow = 1, ncol=4))
    
    #Update the models if there are missing genotypes
  } else {
    model.null.subset <- update(model.null, 
                                subset=complete.cases(geno[, as.character(pairs.eqtl[i, 2])]))
    model.test.subset <- update(model.test, 
                                subset=complete.cases(geno[, as.character(pairs.eqtl[i, 2])]))
    data.frame(matrix(data=c(summary(model.test.subset)$coefficients[2,],
                             anova(model.null.subset, 
                                   model.test.subset)$'Pr(>Chisq)'[2]), nrow = 1, ncol=4))
  }
}))

# Make results table for this gene
colnames(results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "eQTL_pval")

results <- data.frame(Gene = gene, SNP = pairs.eqtl[irange, 2], results)

output_dir <- paste0("../cisresults/chr", chr)
if (!dir.exists(output_dir)){
  dir.create(output_dir)
}
saveRDS(results, paste0(output_dir, "/", gene, ".rds"))
}