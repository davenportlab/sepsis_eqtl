################################################################################
#
# 5. Sex interactions (additional analysis for review and comparison)
#
################################################################################library(lme4)

library(lme4)

# Load data
load("eqtl_int_files.rda")

# set number of PEER factors (same as original eQTL mapping)
n.peer <- 20

# Number of SNP-gene pairs that will be tested
irange <- 1:nrow(pairs.int) # 16049

# Add sex information
sex.info <- read.table("All_genotyping_discretecovs_for_gwas_nodiagorpcs.cov",
                       header=F)
sex <- sex.info$V3[match(GAinSID, sex.info$V1)]
sex <- sex-1

# get threshold for significance (check baseline significance with sex added to model)
thresholds <- read.delim("nominal_pval_thresholds.txt")

# Test each SNP-gene pair in turn for an interaction
results <- do.call(rbind, lapply(irange, function(i){
  # is there more than one signal for this gene?
  gene <- pairs.int[i, 1]
  signals <- pairs.int[pairs.int$Gene == gene, ]
  # get significance threshold
  threshold <- thresholds$threshold[match(gene, thresholds$gene)]
  
  if(nrow(signals) == 1){
    
    # check if SNP still significant with sex in the model
    model.null <- lmer(exp[pairs.int[i, 1], ] ~
                         covs + sex +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)

    model.test <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, pairs.int[i, 2]] +
                         covs + sex +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)

      if (all(complete.cases(geno.int[, pairs.int[i, 2]]))){
        check.sig <- anova(model.null, model.test)$'Pr(>Chisq)'[2] < threshold
        #Update the models if there are missing genotypes
      } else {
        model.null.subset <- update(model.null,
                                    subset=complete.cases(geno.int[, pairs.int[i, 2]]
                                    ))
        model.test.subset <- update(model.test,
                                    subset=complete.cases(geno.int[, pairs.int[i, 2]]
                                    ))
        check.sig <- anova(model.null.subset, model.test.subset)$'Pr(>Chisq)'[2] < threshold
      }
    # if still significant
    # check.sig <- TRUE
    if(check.sig){
      model.null <- lmer(exp[pairs.int[i, 1], ] ~
                           geno.int[, pairs.int[i, 2]] +
                           covs + sex +
                           peer.factors[, 1:n.peer] +
                           (1|GAinSID),
                         REML=FALSE)
      
      model.test <- lmer(exp[pairs.int[i, 1], ] ~
                           geno.int[, pairs.int[i, 2]] +
                           covs + sex +
                           peer.factors[, 1:n.peer] +
                           geno.int[, pairs.int[i, 2]]*sex +
                           (1|GAinSID),
                         REML=FALSE)
      
      # Check that there are enough patients in each group that are minor allele homs
      if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                      & sex == unique(sex)[1])])) > 1 &
          length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                      & sex == unique(sex)[2])])) > 1){
        
        c(summary(model.test)$coefficients[2,],
          summary(model.test)$coefficients["sex",],
          summary(model.test)$coefficients[dim(summary(model.test)$coefficients)[1],],
          anova(model.null, model.test)$'Pr(>Chisq)'[2])
        
      } else {
        rep(NA, 10)
      }
    # if baseline no longer significant  
    } else {
      rep(NA, 10)
    }
    
  } else {
    # If there are multiple signals for the gene include other SNPs in models
    other.snps <- signals$SNP[signals$SNP != pairs.int[i, 2]]
    
    check if SNP still significant with sex in the model
    model.null <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, other.snps] +
                         covs + sex +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)

    model.test <- lmer(exp[pairs.int[i, 1], ] ~
                         geno.int[, other.snps] +
                         geno.int[, pairs.int[i, 2]] +
                         covs + sex +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)

    if (all(complete.cases(geno.int[, pairs.int[i, 2]]))){
      check.sig <- anova(model.null, model.test)$'Pr(>Chisq)'[2] < threshold
      #Update the models if there are missing genotypes
    } else {
      model.null.subset <- update(model.null,
                                  subset=complete.cases(geno.int[, pairs.int[i, 2]]
                                  ))
      model.test.subset <- update(model.test,
                                  subset=complete.cases(geno.int[, pairs.int[i, 2]]
                                  ))
      check.sig <- anova(model.null.subset, model.test.subset)$'Pr(>Chisq)'[2] < threshold
    }
    check.sig <- TRUE
    if(check.sig == TRUE){
      
      model.null <- lmer(exp[pairs.int[i, 1], ] ~
                           geno.int[, pairs.int[i, 2]] +
                           covs + sex +
                           geno.int[, other.snps] +
                           peer.factors[, 1:n.peer] +
                           (1|GAinSID),
                         REML=FALSE)
      
      model.test <- lmer(exp[pairs.int[i, 1], ] ~
                           geno.int[, pairs.int[i, 2]] +
                           covs + sex +
                           geno.int[, other.snps] +
                           peer.factors[, 1:n.peer] +
                           geno.int[, pairs.int[i, 2]]*sex +
                           (1|GAinSID),
                         REML=FALSE)
      
      if (length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                      & sex == unique(sex)[1])])) > 1 &
          length(unique(GAinSID[which(geno.int[, pairs.int[i, 2]] == 2
                                      & sex == unique(sex)[2])])) > 1){
        
        c(summary(model.test)$coefficients[2,],
          summary(model.test)$coefficients["sex",],
          summary(model.test)$coefficients[dim(summary(model.test)$coefficients)[1],],
          anova(model.null, model.test)$'Pr(>Chisq)'[2])
        
      } else {
        rep(NA, 10)
      }
      
    } else {
      rep(NA, 10)
    }
    
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

saveRDS(results, "sex_int_results_cond_incl.rds")

################################################################################

sex.int.cond <- readRDS("sex_int_results_cond_incl.rds")
sex.int.cond <- sex.int.cond[complete.cases(sex.int.cond), ]
sex.int.cond$FDR <- p.adjust(sex.int.cond$Interaction_pval, method="fdr")
table(sex.int.cond$FDR < 0.05) # 9
sex.int.cond$Sig <- sex.int.cond$FDR < 0.05

gene.info <- read.delim("gene_info_864_20412_hla.txt")
sex.int.cond$Symbol <- gene.info$gene_name[match(sex.int.cond$Gene, gene.info$gene_id)]
write.table(sex.int.cond, "sex_int_results_cond_incl.txt", sep="\t")
