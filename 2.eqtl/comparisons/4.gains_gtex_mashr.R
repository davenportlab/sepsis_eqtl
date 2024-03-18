################################################################################
#
# 4. Compare GAinS and GTEx EUR results using mashr
#
################################################################################

# Following the mashr eQTL vignette

library(data.table)
library(ashr)
library(mashr)
library(ggplot2)
set.seed(14022022)

# strong results - pairs significant in GAinS and also tested in GTEx
gtex.lead <- read.delim("gtex_gains_comp_lead_for_mashr.txt") # from gtex_comparison.R
gtex.lead <- subset(gtex.lead, Sig == TRUE)
gains.strong <- data.frame("beta"=gtex.lead$beta, "SE"=gtex.lead$se)
gtex.strong <- data.frame("beta"=gtex.lead$slope, "SE"=gtex.lead$slope_se)

# random results - compare all results (initial pass) and sub-sample
gains <- readRDS("../cisresults/ciseqtl_all.rds")
gains$SNP <- paste0("chr", gains$chr, "_", gains$SNPpos)
gains$Pairs <- paste0(gains$SNP, "_", gains$gene)
gains.bim <- data.frame(fread("../data/genotyping_for_rna-seq_eQTL.bim"))
m1 <- match(gains$snps, gains.bim$V2)
gains$minor <- gains.bim$V5[m1] # minor
gains$major <- gains.bim$V6[m1] # major

gtex.files <- list.files(path="resources/GTEx/AnalysisV8/GTEx_Analysis_v8_EUR_eQTL_all_associations_csv/",
                        pattern="Whole_Blood.v8.EUR.allpairs.chr*", full.names = TRUE)
gtex.results <- lapply(gtex.files, function(x){
  gtex <- fread(x)
  gtex.variants <- unlist(strsplit(gtex$variant_id, "\\_"))
  gtex$SNP <- paste0(gtex.variants[seq(1, nrow(gtex)*5, 5)], "_",
                     gtex.variants[seq(2, nrow(gtex)*5, 5)])
  gtex$REF <- gtex.variants[seq(3, nrow(gtex)*5, 5)]
  gtex$ALT <- gtex.variants[seq(4, nrow(gtex)*5, 5)]
  gtex$phenotype_id <- unlist(strsplit(gtex$phenotype_id, "\\."))[seq(1, nrow(gtex)*2, 2)]
  gtex$Pairs <- paste0(gtex$SNP, "_", gtex$phenotype_id)
  gtex <- subset(gtex, Pairs %in% gains$Pairs)
  print(dim(gtex))
  gtex <- as.data.frame(gtex)
  return(gtex)
})

gtex.results <- rbindlist(gtex.results)
gtex.results <- merge(gtex.results, gains, by="Pairs")

# check alleles/direction of effect
gtex.results$action <- "NA"
gtex.results[which(gtex.results$ALT == gtex.results$minor), "action"] <- "same"

# Alleles can't distinguish because of strand
gtex.results[which(gtex.results$ALT=="A" & gtex.results$REF=="T" |
             gtex.results$ALT=="T" & gtex.results$REF=="A" |
             gtex.results$ALT=="C" & gtex.results$REF=="G" |
             gtex.results$ALT=="G" & gtex.results$REF=="C"), "action"] <-"strand"

# Same but different strand
gtex.results[which(gtex.results$ALT=="A" & gtex.results$minor=="T" & gtex.results$action != "strand" |
             gtex.results$ALT=="T" & gtex.results$minor=="A" & gtex.results$action != "strand" |
             gtex.results$ALT=="C" & gtex.results$minor=="G" & gtex.results$action != "strand" |
             gtex.results$ALT=="G" & gtex.results$minor=="C" & gtex.results$action != "strand"),
     "action"] <-"same"

gtex.results[which(gtex.results$action == "NA"), "action"] <- "flip"
gtex.results$beta.2 <- gtex.results$beta
gtex.results[which(gtex.results$action == "flip"), "beta.2"] <- gtex.results[which(gtex.results$action == "flip"), "beta.2"] * (-1)

# Remove SNPs which can't be distinguished because of strand
gtex.results <- gtex.results[which(gtex.results$action != "strand"), ]
gtex <- gtex.results[, c(1,7:9, 11:18, 22, 25:27)]

# write.table(gtex, "../gtex/gtex_snp_gene_pairs_vs_gains_all_red.txt", sep="\t",
# quote=F, row.names=F)

# sub-sample random results from gtex-gains comparison
gtex <- gtex[sample(1:nrow(gtex), size=nrow(gtex)*0.01), ]
gains.random <- data.frame("beta"=gtex$beta, "SE"=gtex$se)
gtex.random <- data.frame("beta"=gtex$slope, "SE"=gtex$slope_se)

# set up data for mash
effects.strong <- as.matrix(cbind(gains.strong$beta, gtex.strong$beta))
ses.strong <- as.matrix(cbind(gains.strong$SE, gtex.strong$SE))
colnames(effects.strong) <- c("GAinS", "GTEx")
colnames(ses.strong) <- c("GAinS", "GTEx")
data.strong <- mash_set_data(effects.strong, ses.strong)

effects.random <- as.matrix(cbind(gains.random$beta, gtex.random$beta))
ses.random <- as.matrix(cbind(gains.random$SE, gtex.random$SE))
colnames(effects.random) <- c("GAinS", "GTEx")
colnames(ses.random) <- c("GAinS", "GTEx")
data.random <- mash_set_data(effects.random, ses.random)

# correlation structure from random subset
Vhat = estimate_null_correlation_simple(data.random)

# Now we can set up our main data objects with this correlation structure in place:
data.random = mash_set_data(effects.random, ses.random, V=Vhat)
data.strong = mash_set_data(effects.strong, ses.strong, V=Vhat)

# Now we use the strong tests to set up data-driven covariances.
U.pca = cov_pca(data.strong, 2)
U.ed = cov_ed(data.strong, U.pca)

# Now we fit mash to the random tests using both data-driven and canonical covariances. 
# (Remember the Crucial Rule! We have to fit using a random set of tests, and not a 
# dataset that is enriched for strong tests.) The outputlevel=1 option means that it 
# will not compute posterior summaries for these tests (which saves time).

U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
# - Computing 384052 x 301 likelihood matrix.
# - Likelihood calculations took 95.62 seconds.
# - Fitting model with 301 mixture components.
# - Model fitting took 759.76 seconds.

# Now we can compute posterior summaries etc for any subset of tests using the 
# above mash fit. Here we do this for the strong tests. We do this using the same 
# mash function as above, but we specify to use the fit from the previous run of 
# mash by specifying g=get_fitted_g(m), fixg=TRUE. (In mash the parameter g is 
# used to denote the mixture model which we learned above.)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
# - Computing 8137 x 301 likelihood matrix.
# - Likelihood calculations took 1.12 seconds.
# - Computing posterior matrices.
# - Computation allocated took 0.11 seconds.

print(get_pairwise_sharing(m2)) # the same sign and within a factor 0.5 of each other
#           GAinS      GTEx
# GAinS 1.0000000 0.7301587
# GTEx  0.7301587 1.0000000
# list of different ones:
post.mean <- data.frame(get_pm(m2))
rownames(post.mean) <- gtex.lead$Pairs

post.mean$ratio <- post.mean$GAinS/post.mean$GTEx
post.mean$diff <- post.mean$ratio < 0.5 | post.mean$ratio > (1/0.5)
table(post.mean$diff)

post.sig <- get_lfsr(m2)
colnames(post.sig) <- c("GAinS.lfsr", "GTEx.lfsr")
post.mean <- data.frame(post.mean, post.sig)

ggplot(post.mean, aes(GAinS, GTEx)) +
  geom_point(aes(colour=diff)) +
  scale_color_manual(values=c("lightgrey", "darkblue")) +
  theme_bw()

write.table(post.mean, "gains_gtex_mashr_results.txt", sep="\t", quote=F)

gtex.lead$diff.mashr <- post.mean$diff
gtex.lead.diff.mashr <- subset(gtex.lead, diff.mashr == TRUE)
write.table(gtex.lead.diff.mashr, "sepsis_specific_gains_lead_mashr.txt", 
            sep="\t", quote=F)

print(get_pairwise_sharing(m2, factor=0)) # same sign
# GAinS      GTEx
# GAinS 1.0000000 0.9716993
# GTEx  0.9716993 1.0000000
