##############################################################################
# GCTA
##############################################################################

# Aim: to estimate the heritability of SRS and the variance explained by other
# covariates using imputed genotyping data
# https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata

# Input files:
# Genotyping data: All_genotyping_merged_filtered_b38_refiltered_rsID_autosomes
# Phenotype: All_genotyping_SRS1ever_gcta.phen
# Discrete covariates: sex All_genotyping_discretecovs_for_gwas.cov
# Quantitative covariates: age, age2, genotyping PCs 1-7 All_genotyping_continuouscovs_for_gwas.cov

# --reml-est-fix gives estimates of covariate effects

# LD prune the data
plink --bfile All_genotyping_merged_filtered_b38_refiltered_rsID_autosomes --indep-pairwise 50 5 0.2 --out All_genotyping_merged_filtered_b38_refiltered_rsID_autosomes

# Make GRM
gcta64 --bfile All_genotyping_merged_filtered_b38_refiltered_rsID_autosomes --autosome --extract All_genotyping_merged_filtered_b38_refiltered_rsID_autosomes.prune.in --make-grm --out All_genotyping_merged_filtered_b38_refiltered_rsID_pruned --thread-num 10 --maf 0.01

# Estimate SNP heritability
gcta64 --grm All_genotyping_merged_filtered_b38_refiltered_rsID_pruned --pheno All_genotyping_SRS1ever_gcta.phen --reml --out h2_srs1ever_covspruned --covar All_genotyping_discretecovs_for_gwas.cov --qcovar All_genotyping_continuouscovs_for_gwas.cov --reml-est-fix
# https://gcta.freeforums.net/thread/211/estimating-fixed-effects-gcta-greml

####################################################################################
# FINAL RESULT - plot in R
library(ggplot2)
hsq <- read.table("h2_srs1ever_covspruned.hsq", header=TRUE,
                  fill=TRUE, na.strings = "NA")
hsq.p <- as.numeric(hsq[9, 2])
hsq.prop <- as.numeric(hsq[4, 2])
hsq.ci <- hsq[4, 3]
hsq <- hsq[1:3, ]
hsq$Source <- factor(hsq$Source, levels=c("V(G)", "V(e)", "Vp"))
hsq$Variance <- as.numeric(hsq$Variance)
hsq$SE <- as.numeric(hsq$SE)

pdf("Fig_S1_heritability_SRS1ever_agesexPCs_pruned.pdf")
ggplot(hsq, aes(Source, Variance, ymin=Variance-SE, ymax=Variance+SE)) +
  geom_pointrange() +
  coord_flip() +
  theme_bw() +
  geom_text(aes(y=Variance[1], x=1.2,
                label=paste0("V(G)/Vp=", signif(hsq.prop, 2), " (",
                             signif(hsq.prop-hsq.ci, 2),
                             "-", signif(hsq.prop+hsq.ci, 2), ")\n p=",
                             signif(hsq.p, 2))))
dev.off()
