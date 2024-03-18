################################################################################
#
# 4. Cell proportion interactions
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
  signals <- pairs.int[which(pairs.int$Gene == gene), ]
  
  # If only 1 signal
  if(nrow(signals) == 1){
    # Compare eQTL model to model with interaction term (col 1/2/3 of covariates is neut/lymph/mono proportion)
    model.null <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + 
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)
    
    model.n <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 1] +
                         (1|GAinSID),
                       REML=FALSE)
    
    model.m <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 3] +
                      (1|GAinSID),
                    REML=FALSE)
    
    model.l <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 2] +
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
    
    model.n <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      geno.int[, other.snps] +
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 1] +
                      (1|GAinSID),
                    REML=FALSE)
    
    model.m <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      geno.int[, other.snps] +
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 3] +
                      (1|GAinSID),
                    REML=FALSE)
    
    model.l <- lmer(exp[pairs.int[i, 1], ] ~
                      geno.int[, pairs.int[i, 2]] +
                      covs + 
                      geno.int[, other.snps] +
                      peer.factors[, 1:n.peer] +
                      geno.int[, pairs.int[i, 2]]*covs[, 2] +
                      (1|GAinSID),
                    REML=FALSE)
  }
  
  # Check that there are enough patients in each half of neutrophil distribution that are minor allele homs  
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 1] > 0)])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 1] < 0)])) > 1){
    
    n <- c(summary(model.n)$coefficients[2, ],
           summary(model.n)$coefficients[3,],
           summary(model.n)$coefficients[dim(summary(model.n)$coefficients)[1],],
           anova(model.null, model.n)$'Pr(>Chisq)'[2])
    
  } else {
    n <- rep(NA, 10)
  }
  # Check that there are enough patients in each half of monocyte distribution that are minor allele homs  
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 3] > 0)])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 3] < 0)])) > 1){
    
    m <- c(summary(model.m)$coefficients[2, ],
           summary(model.m)$coefficients[5,],
           summary(model.m)$coefficients[dim(summary(model.n)$coefficients)[1],],
           anova(model.null, model.m)$'Pr(>Chisq)'[2])
    
  } else {
    m <- rep(NA, 10)
  }
  # Check that there are enough patients in each half of lymphocyte distribution that are minor allele homs  
  if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 2] > 0)])) > 1 &
      length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                  & covs[, 2] < 0)])) > 1){
    
    l <- c(summary(model.l)$coefficients[2, ],
           summary(model.l)$coefficients[4,],
           summary(model.l)$coefficients[dim(summary(model.n)$coefficients)[1],],
           anova(model.null, model.l)$'Pr(>Chisq)'[2])
    
  } else {
    l <- rep(NA, 10)
  }
  c(n, m, l)
}))

colnames(results) <- c("eQTL_beta_neutrophil",
                       "eQTL_SE_neutrophil",
                       "eQTL_t_neutrophil",
                       "DE_beta_neutrophil",
                       "DE_SE_neutrophil",
                       "DE_t_neutrophil",
                       "Interaction_beta_neutrophil",
                       "Interaction_SE_neutrophil",
                       "Interaction_t_neutrophil",
                       "Interaction_pval_neutrophil",
                       "eQTL_beta_monocytes",
                       "eQTL_SE_monocytes",
                       "eQTL_t_monocytes",
                       "DE_beta_monocytes",
                       "DE_SE_monocytes",
                       "DE_t_monocytes",
                       "Interaction_beta_monocytes",
                       "Interaction_SE_monocytes",
                       "Interaction_t_monocytes",
                       "Interaction_pval_monocytes",
                       "eQTL_beta_lymphocytes",
                       "eQTL_SE_lymphocytes",
                       "eQTL_t_lymphocytes",
                       "DE_beta_lymphocytes",
                       "DE_SE_lymphocytes",
                       "DE_t_lymphocytes",
                       "Interaction_beta_lymphocytes",
                       "Interaction_SE_lymphocytes",
                       "Interaction_t_lymphocytes",
                       "Interaction_pval_lymphocytes")

results <- data.frame(
  Gene = pairs.int[irange, 1],
  SNP = pairs.int[irange, 2],
  results
)

saveRDS(results, "cellprops_int_results_cond_incl.rds")
