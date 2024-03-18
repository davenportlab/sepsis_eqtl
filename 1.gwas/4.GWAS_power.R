##############################################################################
# GWAS replication power calculation
##############################################################################

# Aim: to estimate the sample size required to replicate the SRS1ever vs never
# GWAS results

# https://cran.r-project.org/web/packages/genpwr/vignettes/vignette.html

library(devtools)
# install_github("camillemmoore/Power_Genetics", subdir="genpwr")
library(genpwr)
library(ggplot2)

# Extract full results for the signal SNPs for each locus
gwas.res <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic.clumped", header=T)
gwas.full <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic", header=T)
gwas.full <- subset(gwas.full, SNP %in% gwas.res$SNP)
rm(gwas.res)

# Get MAF information
frq <- read.table("gains_b38.afreq")
gwas.full$maf <- frq$V5[match(gwas.full$SNP, frq$V2)]

# For each signal SNP, use the observed MAF and OR in the power calculation
ss <- genpwr.calc(calc = "ss", 
                  model = "logistic", 
                  ge.interaction = NULL,
                  OR=gwas.full$OR, 
                  Case.Rate=0.4, # estimated occurence of SRS1ever phenotype
                  k=NULL,
                  MAF=gwas.full$maf, 
                  Power=0.8, Alpha=0.05,
                  True.Model="Additive", 
                  Test.Model="Additive")

# function calculated sample size for all pairs of MAF and OR; reduce to the 
# specific pairs for the GWAS results
gwas.full$pair <- paste0(gwas.full$OR, "_", gwas.full$maf)
ss$pair <- paste0(ss$OR, "_", ss$MAF)
ss <- subset(ss, pair %in% gwas.full$pair)

ggplot(ss, aes(MAF, abs(log2(OR)), colour=N_total_at_Alpha_0.05)) +
  geom_point() +
  theme_bw() +
  ylab("log2(Odds Ratio)") +
  xlab("MAF") +
  scale_colour_continuous(type="viridis")
