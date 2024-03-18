################################################################################
#
# 3. SRS interactions
#
################################################################################

library(lme4)

# Load data
load("eqtl_int_files.rda")

# set number of PEER factors (same as original eQTL mapping)
n.peer <- 20

# Number of SNP-gene pairs that will be tested
irange <- 1:nrow(pairs.int) # 16049

# Test each SNP-gene pair in turn for an interaction
results <- do.call(rbind, lapply(irange, function(i){
  # is there more than one signal for the gene involved in this pair?
  gene <- pairs.int[i, 1]
  signals <- pairs.int[pairs.int$Gene == gene, ]
  
  # If only 1 signal
  if(nrow(signals) == 1){
    # Compare eQTL model to model with interaction term (col 12 of covariates is SRS1 status)
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
    # If there are multiple signals for the gene include other SNPs in models
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
  
  # Check that there are enough patients in each group that are minor allele homs
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                       & covs[, 12] == unique(covs[, 12])[1])])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                       & covs[, 12] == unique(covs[, 12])[2])])) > 1){
    
    c(summary(model.test)$coefficients[2,],
      summary(model.test)$coefficients["covs12",],
      summary(model.test)$coefficients[dim(summary(model.test)$coefficients)[1],],
      anova(model.null, model.test)$'Pr(>Chisq)'[2])
    
  } else {
    rep(NA, 10)
  }
}))

colnames(results) <- c("eQTL_beta",
                       "eQTL_SE",
                       "eQTL_t",
                       "DE_beta",
                       "DE_SE",
                       "DE_t",
                       "Interaction_beta",
                       "Interaction_SE",
                       "Interaction_t",
                       "Interaction_pval")

results <- data.frame(Gene = pairs.int[irange, 1],
                      SNP = pairs.int[irange, 2],
                      results)

saveRDS(results, "srs_int_results_cond_incl.rds")
