################################################################################
#
# 4. Conditional analysis: identification of independent eQTL signals
#
################################################################################

.libPaths("/nfs/users/nfs_k/kb21/R/x86_64-pc-linux-gnu-library/3.6")
options(stringsAsFactors = FALSE)
args <- commandArgs(TRUE)
# first argument is the gene for which to map eQTL
jobindex <- as.numeric(args[1])
# second argument is the number of PEER factors to include in the model
n.peer <- as.numeric(args[2])

library(lme4)
library(data.table)

################################################################################

# read in results from the first pass analysis and use the job index to select gene
primary.peak.results <- read.delim("../cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt",
                                   stringsAsFactors = F)
this.gene <- primary.peak.results[jobindex, ]
gene <- this.gene$gene

# get threshold for significance
thresholds <- read.delim("../cisresults/nominal_pval_thresholds.txt")
threshold <- thresholds$threshold[match(gene, thresholds$gene)]

# set this output (the peak association) as the most recent result
last.result <- this.gene
# is there a significant association in the most recent results?
last.was.sig <- min(last.result$pvalue) <= threshold

# if it was determined to be significant, add it to a list of snps for the forward pass
if(this.gene$Sig == TRUE){
  sig.snps <- this.gene[, c("snps", "gene", "beta", "se", "pvalue")]
  colnames(sig.snps) <- c("SNP", "Gene", "eQTL_beta", "eQTL_SE", "pvalue")
  
  # extract list of snp-gene pairs to test and read in data for this gene
  pairs <- fread("../data/gene_snp_pairs_cis_rnaseq.txt", stringsAsFactors = F)
  pairs.eqtl <- pairs[pairs$Gene == this.gene$gene, ]
  pairs.eqtl$SNP <- gsub(":", ".", pairs.eqtl$SNP)
  print(dim(pairs.eqtl))
  chr <- pairs.eqtl[1, 3]
  # read in eqtl files
  load(paste("../data/eqtl_files_", chr, ".rda", sep=""))
}

# forward step to identify additional signals for eGenes
while(last.was.sig == TRUE){
  # remove all peak snps from the list of pairs to test
  pairs.eqtl.it <- pairs.eqtl[!(pairs.eqtl$SNP %in% sig.snps$SNP), ]
  dim(pairs.eqtl.it)
  irange <- 1:nrow(pairs.eqtl.it)
  
  # rerun eqtl mapping including peak SNPs in base model
  sig.snps.geno <- geno[, sig.snps[, 1], drop=FALSE]
  results <- rbindlist(lapply(irange, function(i){
    
    model.null <- lmer(as.numeric(exp[gene, ]) ~
                         sig.snps.geno +
                         covs +
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE) 
    model.test <- lmer(as.numeric(exp[gene, ]) ~
                         as.numeric(geno[, as.character(pairs.eqtl.it[i, 2])]) +
                         sig.snps.geno +
                         covs + 
                         peer.factors[, 1:n.peer] +
                         (1|GAinSID),
                       REML=FALSE)
    
    # if the test SNP is identical to a previously significant SNP the models are the same
    if (any(is.na(fixef(model.test, add.dropped=T)[3:(ncol(sig.snps.geno)+2)]))){
      data.frame(matrix(data=c(summary(model.test)$coefficients[2, ], 1), nrow=1, ncol=4))
      
      # If there are no missing genotypes
    } else if (all(complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]))){
      data.frame(matrix(data=c(summary(model.test)$coefficients[2, ],
                               anova(model.null, model.test)$'Pr(>Chisq)'[2]), nrow=1, ncol=4))
      
      #Update the models if there are missing genotypes
    } else {
      model.null.subset <- update(model.null, 
                                  subset=complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]
                                  ))
      model.test.subset <- update(model.test, 
                                  subset=complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]
                                  ))
      data.frame(matrix(data=c(summary(model.test.subset)$coefficients[2,],
                               anova(model.null.subset, model.test.subset)$'Pr(>Chisq)'[2]), nrow=1, ncol=4))
    }
  }))
  
  colnames(results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "pvalue")
  last.result <- data.frame(
    Gene = gene,
    SNP = pairs.eqtl.it[irange, 2],
    results
  )
  
  # update: is there a significant assn in the most recent results
  last.was.sig <- min(last.result$pvalue) <= threshold
  
  # if yes, add peak sig SNP to list of sig SNPs
  if(last.was.sig == TRUE){
    sig.snps <- rbind(sig.snps,
                      last.result[which.min(last.result$pvalue), 
                                  c("SNP", "Gene", "eQTL_beta", "eQTL_SE", "pvalue")])
  }
}

# when there are no significant SNPs in the last iteration, stop

# If there are significant SNPs, save the forward selection results
if(exists("sig.snps")){
  output_dir <- paste0("../cisresults/conditionalanalysis/chr", chr)
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
  }
  saveRDS(sig.snps, paste0(output_dir, "/forward_", gene, ".rds"))
  
  n.snps <- nrow(sig.snps)
  # If there is only 1 significant SNP the backward step isn't necessary
  if(n.snps == 1){
    sig.snps.backward <- sig.snps
  } else {
    # backward selection if there are multiple significant SNPs
    sig.snps.backward <- rbindlist(lapply(1:n.snps, function(s){
      test.snp <- sig.snps[s, ]
      other.snps <- sig.snps[-s, ]
      other.snps.geno <- geno[, other.snps[, 1], drop=F]
      
      # scan for cis-eQTLs controlling for all but one SNP from the forward pass
      pairs.eqtl.it <- pairs.eqtl[!(pairs.eqtl$SNP %in% other.snps[, 1]), ]
      dim(pairs.eqtl.it)
      irange <- 1:nrow(pairs.eqtl.it)
      backward.results <- rbindlist(lapply(irange, function(i){
        
        model.null <- lmer(as.numeric(exp[gene, ]) ~
                             other.snps.geno +
                             covs +
                             peer.factors[, 1:n.peer] +
                             (1|GAinSID),
                           REML=FALSE) 
        
        model.test <- lmer(as.numeric(exp[gene, ]) ~
                             as.numeric(geno[, as.character(pairs.eqtl.it[i, 2])]) +
                             other.snps.geno +
                             covs +
                             peer.factors[, 1:n.peer] +
                             (1|GAinSID),
                           REML=FALSE)
        
        # if the test SNP is identical to a previously significant SNP the models are the same
        if (any(is.na(fixef(model.test, add.dropped=T)[3:(ncol(other.snps.geno)+2)]))){
          data.frame(matrix(data=c(summary(model.test)$coefficients[2, ], 1), nrow=1, ncol=4))
          
          # If there are no missing genotypes
        } else if (all(complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]))){
          data.frame(matrix(data=c(summary(model.test)$coefficients[2, ],
                                   anova(model.null, model.test)$'Pr(>Chisq)'[2]), nrow=1, ncol=4))
          
          # Update the models if there are missing genotypes
        } else {
          model.null.subset <- update(model.null, 
                                      subset=complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]
                                      ))
          model.test.subset <- update(model.test, 
                                      subset=complete.cases(geno[, as.character(pairs.eqtl.it[i, 2])]
                                      ))
          data.frame(matrix(data=c(summary(model.test.subset)$coefficients[2,],
                                   anova(model.null.subset, model.test.subset)$'Pr(>Chisq)'[2]), nrow=1, ncol=4))
        }
      }))
      
      colnames(backward.results) <- c("eQTL_beta", "eQTL_SE", "eQTL_t", "pvalue")
      backward.results <- data.frame(Gene = gene, SNP = pairs.eqtl.it[irange, 2],
                                     backward.results)
      backward.results$Sig <- backward.results$pvalue <= threshold
      
      # the lead variant from this scan, which controls for all other signals 
      # found in the forward stage, was chosen as the
      # variant that represents the signal best in the full model.
      if(any(backward.results$Sig == TRUE)){
        sig.snps.backward <- backward.results[which.min(backward.results$pvalue), 
                                              c("SNP", "Gene", "eQTL_beta", "eQTL_SE", "pvalue")]
        sig.snps.backward
      } else {
        # If no variant was significant at the gene-level threshold the variant in 
        # question was dropped
        sig.snps.backward <- data.frame(matrix(data=rep("NA", 5), nrow=1, ncol=5))
        colnames(sig.snps.backward) <- c("SNP", "Gene", "eQTL_beta", "eQTL_SE", "pvalue")
        sig.snps.backward
      }
    }))
    
    # This gives the list of signal SNPs
    # Nikhil generated conditional summary stats for each signal SNP finally
    # adjusting for the other signal SNPs identified in the backwards pass
    # (rather than those from the forward pass as here)
    
  }
  saveRDS(sig.snps.backward, paste0(output_dir, "/backward_", gene, ".rds"))
} else {print("Not significant")}
