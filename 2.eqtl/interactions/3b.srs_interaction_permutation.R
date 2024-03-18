################################################################################
#
# 3b. SRS interactions permutation analysis
#
################################################################################

library(lme4)
set.seed(31052022)

# Load data
load("eqtl_int_files.rda")
n.peer <- 20
irange <- 1:nrow(pairs.int) # 16049
args <- commandArgs(TRUE)
# Run as job array; first argument is the permutation number
jobindex <- as.numeric(args[1])

srs.perm <- replicate(1000, sample(covs[, 12], 823, replace=FALSE))
covs[, 12] <- srs.perm[, jobindex]

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
                         geno.int[, pairs.int[i, 2]]*covs[, 12] +
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
                         geno.int[, pairs.int[i, 2]]*covs[, 12] +
                         (1|GAinSID),
                       REML=FALSE)
  }
  
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 12] == unique(covs[, 12])[1])])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 12] == unique(covs[, 12])[2])])) > 1){
    
    # Just store the p-value
    anova(model.null, model.test)$'Pr(>Chisq)'[2]
    
  } else {
    NA
  }
}))

# And store how many are significant for each permutation
results <- results[!is.na(results)]
fdr <- p.adjust(results, method="fdr")
n.sig <- length(which(fdr < 0.05))
cat(n.sig, file="srs_permutation_results.txt", append=TRUE, fill=TRUE)
