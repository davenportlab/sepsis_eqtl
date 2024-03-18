################################################################################
#
# 1. Make files for interaction analysis
#
################################################################################

# Take the lead SNP-gene pairs for all conditional eQTL signals to test for interactions
res <- read.delim("conditional_eQTL_results_final.txt")
snps.sig <- unique(res$SNP)

# Subset the genotyping data to just those SNPs
geno.int <- do.call(cbind, lapply(1:22, function(i){
  load(paste0("../data/eqtl_files_", i, ".rda"))
  geno[, which(colnames(geno) %in% snps.sig)]
}))

# Load the rest of the required data for interaction mapping
load("../data/eqtl_files_22.rda")

# Make a list of SNP-gene pairs to test
pairs.int <- res[, c("Gene", "SNP")]
pairs.int$Gene <- as.character(pairs.int$Gene)
pairs.int$SNP <- as.character(pairs.int$SNP)
dim(pairs.int) #  16049     2
rownames(pairs.int) <- 1:nrow(pairs.int)

# Make the R data file containing everything needed for interaction testing
int_file <- "eqtl_int_files.rda"
save(list=c("exp", "geno.int", "covs", "GAinSID", "peer.factors", "pairs.int"),
     file = int_file)
