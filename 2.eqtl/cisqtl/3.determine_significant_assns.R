################################################################################
#
# 3. Hierarchical multiple testing control: identify significant cis-eQTL from 1st pass
#
################################################################################

options(stringsAsFactors = FALSE)
library(data.table)

################################################################################

# combine results from first pass mapping by chromosome
cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(length(temp_eqtl))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  eqtl <- data.frame("snps" = eqtl$SNP, "gene" = eqtl$Gene,
                     "statistic"=eqtl$eQTL_t, "pvalue" = eqtl$eQTL_pval,
                     "beta"=eqtl$eQTL_beta, 'se'=eqtl$eQTL_SE)
  eqtl$snps <- gsub(":", ".", eqtl$snps)
  print(dim(eqtl))
  write.table(eqtl, output_filename, sep="\t", row.names=F, quote=F)
}

# make files for eigenMT
bim <- data.frame(fread("../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim",
                        header=FALSE))

gene.info <- read.delim("../data/gene_info_20412.txt")
genes <- gene.info[, 1:3]
colnames(genes) <- c("chrom", "TSS", "TES")
genes$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]
genes$TES[which(gene.info$strand == "-")] <- gene.info$start[which(gene.info$strand == "-")]
# write.table(genes, "../cisresults/eigenMT/gene_pos_for_eigenMT.txt")

for(i in 1:22){
  # eqtl results
  cat.rds(paste0("../cisresults/chr", i), "*rds",
          paste0("../cisresults/eigenMT/eqtl_results_for_eigenMT_chr", i, ".txt"))
  
  res <- data.frame(fread(paste0("../cisresults/eigenMT/eqtl_results_for_eigenMT_chr", i, ".txt")))
  
  # genotyping data
  load(paste0("../data/eqtl_files_", i, ".rda"))
  geno <- as.data.table(geno)
  geno <- geno[which(!duplicated(GAinSID))]
  geno <- t(geno)
  colnames(geno) <- GAinSID[which(!duplicated(GAinSID))]
  write.table(geno, paste0("../cisresults/eigenMT/geno_for_eigenMT_chr", i, ".txt"),
              sep="\t", quote=F)
  
  # snp positions
  chr.bim <- subset(bim, V1 == i)
  chr.bim$V2 <- gsub(":", ".", chr.bim$V2)
  geno.pos <- chr.bim[match(res$snps, chr.bim$V2), c(2, 1, 4)]
  write.table(geno.pos, paste0("../cisresults/eigenMT/geno_pos_for_eigenMT_chr", i, ".txt"), sep="\t",
              row.names=F, quote=F)
  
  # gene positions
  genes.pos <- subset(genes, chrom == i)
  genes.pos <- genes.pos[match(unique(res$gene), rownames(genes.pos)), ]
  write.table(genes.pos, paste0("../cisresults/eigenMT/gene_pos_for_eigenMT_chr", i, ".txt"), 
              sep="\t", quote=F)
  
  # add to res file
  res$chr <- i
  res$SNPpos <- geno.pos[match(res$snps, geno.pos$V2), 3]
  res$TSS <- genes.pos$TSS[match(res$gene, rownames(genes.pos))]
  write.table(res, paste0("../cisresults/eigenMT/eqtl_results_for_eigenMT_chr", i, ".txt"),
              sep="\t", quote=F, row.names=F)
}

################################################################################

# Run eigenMT to determine local significance

# for CHR in {1..22}
# do
# python eigenMT.py --CHROM $CHR --QTL eQTL_results_for_eigenMT_chr$CHR.txt --GEN geno_for_eigenMT_chr$CHR.txt --GENPOS geno_pos_for_eigenMT_chr$CHR.txt --OUT eigenMT_chr$CHR --window 200 --PHEPOS gene_pos_for_eigenMT_chr$CHR.txt --cis_dist 1000000
# done

################################################################################

# Hierarchical correction to identify eGenes

# read in lead results and correct for number of genes
eigen.res <- list.files(path = ".", pattern = "^eigenMT_chr*", full.names = T)
eigen.res <- rbindlist(lapply(eigen.res, read.delim))
eigen.res$BF.FDR <- p.adjust(eigen.res$BF, method="fdr")
eigen.res$Sig <- eigen.res$BF.FDR < 0.05
write.table(eigen.res, "cis-eQTL_eigenMT_corrected.txt", sep="\t")

################################################################################

# Identify all significant associations for eGenes

# calculate the global significance threshold: 
# the locally corrected p value corresponding to q-value=0.05
ub <- min(eigen.res$BF[eigen.res$BF.FDR > 0.05])  # smallest p-value above FDR
lb <- max(eigen.res$BF[eigen.res$BF.FDR <= 0.05])  # largest p-value below FDR
fdr.threshold <- (lb+ub)/2

# determine the nominal pvalue threshold for each gene
thresholds <- rbindlist(lapply(1:nrow(eigen.res), function(g){
  gene <- as.character(eigen.res[g, 2])
  n.tests <- eigen.res$TESTS[g]
  threshold <- fdr.threshold/n.tests
  return(data.frame(gene, n.tests, threshold))
}))
write.table(thresholds, "../cisresults/nominal_pval_thresholds.txt", sep="\t",
            row.names=F, quote=F)
eigen.res$threshold <- thresholds$threshold
write.table(eigen.res, "../cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt", sep="\t")

full.res <- list.files(path = "../cisresults/eigenMT/", 
                       pattern = "^eqtl_results_for_eigenMT_chr", 
                       full.names = T)

full.res <- rbindlist(lapply(full.res, read.delim, stringsAsFactors=F))
full.res$threshold <- eigen.res$threshold[match(full.res$gene, eigen.res$gene)]
saveRDS(full.res, "../cisresults/ciseqtl_all.rds")

all.sig.eqtl <- full.res[which(full.res$pvalue <= full.res$threshold), ]
saveRDS(all.sig.eqtl, "../cisresults/cisqtl_all_significant.rds")
