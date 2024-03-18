##############################################################################
# PEER
##############################################################################

# Aim: to calculate PEER factors for the RNA-seq dataset to use in eQTL mapping
# Hold out specific covariates that we plan to test for interactions so these
# are not captured by PEER factors

options(stringsAsFactors = FALSE)
library(RNOmni)
library(peer)
set.seed("27082021")

# read in expression data
expr <- read.delim("../data/gex_rnaseq_for_eqtl_hla.txt") # only samples with genotyping data
expr <- t(expr)

# read in covariates file (sampleID, patient ID, diagnosis, cell proportions)
covs <- read.delim("../data/covs_for_eqtl_823_updated_srs.txt", row.names = 1)
covs <- covs[match(rownames(expr), rownames(covs)), ]

# rank inverse normal transform the cell proportions
covs$Neutrophils[!(is.na(covs$Neutrophils))] <- rankNorm(covs$Neutrophils[!(is.na(covs$Neutrophils))], k = 3/8)
covs$Lymphocytes[!(is.na(covs$Lymphocytes))] <- rankNorm(covs$Lymphocytes[!(is.na(covs$Lymphocytes))], k = 3/8)
covs$Monocytes[!(is.na(covs$Monocytes))] <- rankNorm(covs$Monocytes[!(is.na(covs$Monocytes))], k = 3/8)
# Set missing values to the mean (0)
covs$Neutrophils[is.na(covs$Neutrophils)] <- 0
covs$Lymphocytes[is.na(covs$Lymphocytes)] <- 0
covs$Monocytes[is.na(covs$Monocytes)] <- 0

# write.table(covs, "../data/covs_for_eqtl_823_updated_srs_int_cells.txt", sep="\t")

covs <- covs[, c("Neutrophils", "Lymphocytes", "Monocytes", 
                 "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                 "Diagnosis", "SRS1")]
# Diagnosis: 0=CAP, 1=FP. SRS1: 0=non-SRS1, 1=SRS1

model <- PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model, 30)

# PEER can automatically include an additional factor (covariate) to account for 
# the mean expression. For most use cases, including the mean effect is likely 
# to be a good choice.
PEER_setAdd_mean(model, TRUE)
PEER_setCovariates(model, as.matrix(covs))
# This sets the first C factors to be fixed to the observed covariates, and 
# extends the hidden factor matrix X to have additional C columns. 

PEER_setNmax_iterations(model, 500)
PEER_update(model)
# Converged (var(residuals)) after 138 iterations

factors <- PEER_getX(model)
weights <- PEER_getW(model)
precision <- PEER_getAlpha(model)
residuals <- PEER_getResiduals(model)

save(model, factors, weights, precision, residuals, 
     file="../peer/peer_factors_covs_int_cells_823_hla.Rda")

var.names <- c("Neutrophils", "Lymphocytes", "Monocytes", 
               "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
               "Diagnosis", "SRS1", "Intercept", "PEER_1", 
               "PEER_2", "PEER_3", "PEER_4", "PEER_5", "PEER_6", 
               "PEER_7", "PEER_8", "PEER_9", "PEER_10", "PEER_11", 
               "PEER_12", "PEER_13", "PEER_14", "PEER_15", "PEER_16",
               "PEER_17", "PEER_18", "PEER_19", "PEER_20", "PEER_21", 
               "PEER_22", "PEER_23", "PEER_24", "PEER_25", "PEER_26", 
               "PEER_27", "PEER_28", "PEER_29", "PEER_30")
precision <- data.frame("precision"=precision, 
                        "variable"=var.names)
precision$variable <- factor(precision$variable, 
                             levels = var.names)
# convert from inverse variance to variance
ggplot(precision, aes(variable, 1/precision)) +
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
