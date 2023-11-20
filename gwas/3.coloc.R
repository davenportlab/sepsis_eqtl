#################################################################################
#
# Colocalisation of SRS GWAS and Sepsis RNA-seq eQTL hits
#
#################################################################################

# AIM: to find loci that are significant in both the SRS GWAS and eQTL and test 
# if there is likely to be a shared causal effect (colocalisation)

# STRATEGY:

# Start with the GWAS summary stats
# Find all that pass the genome-wide suggestive (1e-05) threshold
# Look for overlap with conditioned eQTL: same SNP set tested in both so don't expand to LD SNPs

# Taking the GWAS hits with an overlap
# Reduce to independent loci from clumping

# To test for shared effect with coloc:
# Reduce to SNPs tested in both
# Look up eQTL results for the same SNPs 
# Pull out MAF etc from the full cohort to calculate posterior probabliities
# Run coloc on each locus: eQTL results are conditioned on all other signals so fits the 
# assumption of 1 causal SNP

#################################################################################

options(stringsAsFactors = FALSE)

library(coloc)
library(rlist)
library(gdata)
library(data.table)
library(ggplot2)
library(gridExtra)

# read in SRS GWAS results
srs.gwas <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic", header=T)

# find all that pass the genome-wide suggestive (1e-05) threshold
gwas.hits <- subset(srs.gwas, P < 5e-5)

# read in conditioned eQTL results and subset to significant results
eqtl <- readRDS("../cisresults/conditionalanalysis/all_associations_conditional.rds")
sig.eqtl <- subset(eqtl, P_Value < threshold)

# look for overap between GWAS and eQTL
table(gwas.hits$SNP %in% sig.eqtl$SNP)
combine.a <- merge(gwas.hits, sig.eqtl, by.x="SNP", by.y="SNP")

# get gene names
gene.info <- read.delim("gene_info_20412.txt")
combine.a$Symbol <- gene.info$gene_name[match(combine.a$Gene, gene.info$gene_id)]
unique(combine.a$Symbol)

# remove duplicates (same GWAS signal and eGene and rank)
# assign gwas SNPs to independent loci and test each once
clumped <- read.table("srs_gwas_logistic_sex_age_pcs_srs1ever.assoc.logistic.clumped", header=T)
for(i in 1:nrow(combine.a)){
  if(combine.a$SNP[i] %in% clumped$SNP){
    combine.a$GWASLoci[i] <- clumped$SNP[match(combine.a$SNP[i], clumped$SNP)]
  } else {
    combine.a$GWASLoci[i] <- clumped$SNP[grepl(combine.a$SNP[i], clumped$SP2)]
  }
}

combine.a$gwas_eqtl_pair <- paste0(combine.a$GWASLoci, "_", combine.a$Gene, "_", combine.a$Signal)
combine <- combine.a[!duplicated(combine.a$gwas_eqtl_pair), ]
unique(c(combine$GWASLoci, combine$SNP))
# 5 loci

# get r2 information for local association plots
# plink --allow-extra-chr --bfile genotyping_for_rna-seq_eQTL --out rs60143665 --ld-snp rs60143665 --ld-window 99999 --ld-window-kb 1100 --ld-window-r2 0 --r2
# plink --allow-extra-chr --bfile genotyping_for_rna-seq_eQTL --out rs34100 --ld-snp rs34100 --ld-window 99999 --ld-window-kb 1100 --ld-window-r2 0 --r2
# plink --allow-extra-chr --bfile genotyping_for_rna-seq_eQTL --out rs10405668 --ld-snp rs10405668 --ld-window 99999 --ld-window-kb 1100 --ld-window-r2 0 --r2
# plink --allow-extra-chr --bfile genotyping_for_rna-seq_eQTL --out rs34099 --ld-snp rs34099 --ld-window 99999 --ld-window-kb 1100 --ld-window-r2 0 --r2
# plink --allow-extra-chr --bfile genotyping_for_rna-seq_eQTL --out rs9524848 --ld-snp rs9524848 --ld-window 99999 --ld-window-kb 1100 --ld-window-r2 0 --r2

# and MAF for coloc
maf <- read.delim("gains_b38.afreq", sep="")

# For each signal, test for colocalisation
coloc.results <- list()

# using all snps in 1Mb window
merging.distance <- 1000000

pdf("srs_gwas_coloc_results.pdf", onefile = T, useDingbats=F, width=10)
for(i in 1:nrow(combine)){
  tryCatch({
    # get SNPs in window from GWAS result
    gwas.locus <- data.frame(srs.gwas[(srs.gwas$CHR == combine$CHR[i] &
                                         srs.gwas$BP %in% (combine$BP[i] - merging.distance):(combine$BP[i] + merging.distance)), ])
    
    # For these SNPs, add the eQTL stats for that eGene signal
    eGene <- combine$Gene[i]
    rank <- combine$Signal[i]
    gwas.locus <- merge(gwas.locus, eqtl[which(eqtl$Gene == eGene & Signal == rank), ], 
                        by="SNP")
    
    # And add allele frequency information
    gwas.locus$MAF <- maf$ALT_FREQS[match(gwas.locus$SNP, maf$ID)]
    
    # GWAS data in coloc format
    gwas <- list(pvalues=gwas.locus$P,
                 N=997,
                 MAF=gwas.locus$MAF, 
                 type="cc",
                 beta=log(gwas.locus$OR),
                 s=0.44)
    # eQTL data in coloc format
    eqtl.loc <- list(pvalues=gwas.locus$P_Value,
                     N=638,
                     MAF=gwas.locus$MAF,
                     beta=gwas.locus$Beta,
                     varbeta=gwas.locus$SE.y^2,
                     type="quant")
    
    # Run coloc
    coloc.results[[i]] <- coloc.abf(gwas, eqtl.loc)$summary
    names(coloc.results)[i] <- eGene
    
    # Plots
    # Look up LD info from the lead SNP
    snp <- gwas.locus[which.min(gwas.locus$P), "SNP"]
    r2 <- read.table(paste0(snp, ".ld"), header=T)
    gwas.locus$r2 <- r2$R2[match(gwas.locus$SNP, r2$SNP_B)]
    gwas.locus$LD <- cut(gwas.locus$r2, breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))
    
    gg1 <- ggplot(gwas.locus, aes(-log10(P), -log10(P_Value))) +
      geom_point(aes(colour=LD)) +
      scale_colour_manual(values=c("darkblue", "skyblue", "gold", "orange", "red"), drop=F) +
      theme_bw() +
      xlab("SRS GWAS pvalue") +
      ylab(paste0(combine$Symbol[i], " eQTL pvalue"))
    gg2 <- ggplot(gwas.locus, aes(BP, -log10(P))) +
      geom_point(aes(colour=LD)) +
      scale_colour_manual(values=c("darkblue", "skyblue", "gold", "orange", "red"), drop=F) +
      theme_bw() +
      ggtitle("SRS GWAS")
    gg3 <- ggplot(gwas.locus, aes(BP, -log10(P_Value))) +
      geom_point(aes(colour=LD)) +
      scale_colour_manual(values=c("darkblue", "skyblue", "gold", "orange", "red"), drop=F) +
      theme_bw() +
      ggtitle(paste0(combine$Symbol[i], " eQTL"))
    lay <- rbind(c(1,1,2,2),
                 c(1,1,3,3))
    grid.arrange(grobs = list(gg1, gg2, gg3), layout_matrix = lay)
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
dev.off()

coloc.results.red <- coloc.results[lapply(coloc.results, length) > 0]
coloc.results.df <- lapply(coloc.results.red, function(x){
  x <- data.frame(t(data.frame(x)))
})

coloc.results.df <- rbindlist(coloc.results.df)
apply(coloc.results.df[, 2:6], 2, function(x){
  length(which(x > 0.7))
})

# PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
# 1         0         1         0         2 

write.table(coloc.results.df, "../SRS_GWAS_coloc_results.txt", sep="\t")
