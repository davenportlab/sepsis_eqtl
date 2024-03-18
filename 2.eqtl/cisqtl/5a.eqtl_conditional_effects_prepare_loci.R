#----------------------------------------------------------
# Prepare Loci for Conditional Effects
# Created: 20 March 2022
#----------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
chr = args[1]

#----------------------------------------------------------
# Environment
#----------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(data.table)
library(parallel)
library(foreach)

# SNP information
geno.bim <- fread("All_genotyping_merged_filtered_b38_refiltered_rsID.bim", sep="\t") %>%
    as.data.frame()
colnames(geno.bim) <- c("chr", "snp", "cM", "Position", "minor_allele", "major_allele")
geno.bim$snp <- gsub(":", ".", geno.bim$snp)

# Conditional Analysis Results
gene.info <- read.table("gene_info_864_20412_hla.txt")
cis.eqtl.conditional <- readRDS("conditional_eQTL_results_final.rds") %>%
    merge(., gene.info, by.x="Gene", by.y=0) %>%
    dplyr::filter(seqnames == chr) %>%
    dplyr::select(SNP, Gene, eQTL_beta, eQTL_SE, pvalue, Number)

cis.eqtl.conditional.loci <- split(cis.eqtl.conditional, cis.eqtl.conditional$Gene)

# Summary statistics from initial pass
cis.eqtl.summary <- readRDS("ciseqtl_all.rds")
cis.eqtl.summary <- cis.eqtl.summary[cis.eqtl.summary$gene %in% cis.eqtl.conditional$Gene,]
cis.eqtl.summary <- merge(cis.eqtl.summary, geno.bim, by.x="snps", by.y="snp")
cis.eqtl.summary$chr.y <- NULL
setnames(cis.eqtl.summary, "chr.x", "chr")
cis.eqtl.loci <- split(cis.eqtl.summary, cis.eqtl.summary$gene)

# Number of signals at each locus
cis.eqtl.conditional.n.signals <- cis.eqtl.conditional %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarize(Signals=n())

#----------------------------------------------------------
# Single Signal Summary Statistics
#----------------------------------------------------------

single.signal.stats <- do.call(rbind, mclapply(cis.eqtl.conditional.n.signals$Gene[cis.eqtl.conditional.n.signals$Signals == 1], function(locus) {

    cis.eqtl.loci[[locus]] %>%
        dplyr::mutate(Gene=locus, Signal=1) %>%
        dplyr::select(Gene, Signal, chr, snps, SNPpos, beta, se, pvalue)
}, mc.cores=16))

write.table(single.signal.stats, "single_signal_stats.txt", sep="\t", row.names=F, col.names=F, quote=F)

#----------------------------------------------------------
# Multiple Signal Summary Statistics
#----------------------------------------------------------

# Only identify conditional summary statistics for loci with more than one signal

doParallel::registerDoParallel(cores=16)
foreach(locus=cis.eqtl.conditional.n.signals$Gene[cis.eqtl.conditional.n.signals$Signals > 1]) %dopar% {

    # Save SNPs at eQTL locus
    cis.eqtl.loci[[locus]] %>%
        dplyr::select(chr, snps, SNPpos) %>%
        write.table(., paste0(locus, ".snps.txt"), row.names=F, col.names=T, quote=F)

    # Save conditional SNP list
    write.table(cis.eqtl.conditional.loci[[locus]]$SNP, paste0(locus, ".conditional_snps.txt"), row.names=F, col.names=F, quote=F)
}

#----------------------------------------------------------
# Design Matrix
#----------------------------------------------------------

held.out <- c("Neutrophils", "Lymphocytes", "Monocytes", paste0("PC", 1:7), "Diagnosis", "SRS1")
peer <- paste0("PEER_", 1:20)

covs <- read.table("covs_and_peer_factors.txt") %>%
    dplyr::mutate(Sample.ID=gsub("^GA", "", rownames(.))) %>%
    dplyr::mutate(GAinS.ID=gsub("\\_.", "", Sample.ID)) %>%
    dplyr::select(Sample.ID, GAinS.ID, any_of(held.out), any_of(peer))

write.csv(covs, "design_matrix.csv", row.names=F)
