################################################################################
#
# 2b. Source of sepsis interactions permutation analysis
#
################################################################################

library(lme4)
set.seed(31052022)

# Load data
load("eqtl_int_files.rda")
n.peer <- 20
irange <- 1:nrow(pairs.int) # 16049
args <- commandArgs(TRUE)
# first argument is the permutation number
jobindex <- as.numeric(args[1])

# permute diagnosis, this should still be the same across serial samples from the same patient
diag <- data.frame("GAinSID"=GAinSID, "diagnosis"=covs[, 11])
diag <- diag[!duplicated(diag$GAinSID), ]
diag.perm <- replicate(1000, sample(diag$diagnosis, nrow(diag), replace=FALSE))
covs[, 11] <- diag.perm[match(GAinSID, diag$GAinSID), jobindex]

results <- do.call(rbind, lapply(irange, function(i){
  pairs.int <- data.frame(pairs.int)
  # is there more than one signal for this gene?
  gene <- pairs.int[i, 1]
  signals <- pairs.int[pairs.int$Gene == gene, ]
  
  if(nrow(signals) == 1){
    model.null <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + 
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)
    
    model.test <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + 
                         peer.factors[, 1:n.peer] +
                         geno.int[, pairs.int[i, 2]]*covs[, 11] +
                         (1|GAinSID),
                       REML=FALSE)
  } else {
    other.snps <- signals$SNP[signals$SNP != pairs.int[i, 2]]
    model.null <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + 
                         geno.int[, other.snps] +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)
    
    model.test <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + 
                         geno.int[, other.snps] +
                         peer.factors[, 1:n.peer] +
                         geno.int[, pairs.int[i, 2]]*covs[, 11] +
                         (1|GAinSID),
                       REML=FALSE)
  }
  
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 11] == unique(covs[, 11])[1])])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 11] == unique(covs[, 11])[2])])) > 1){
    
    # Just store the interaction p-value
    anova(model.null, model.test)$'Pr(>Chisq)'[2]
    
  } else {
    NA
  }
}))

# And store how many are significant for each permutation
results <- results[!is.na(results)]
fdr <- p.adjust(results, method="fdr")
n.sig <- length(which(fdr < 0.05))
cat(n.sig, file="diag_permutation_results.txt", append=TRUE, fill=TRUE)
