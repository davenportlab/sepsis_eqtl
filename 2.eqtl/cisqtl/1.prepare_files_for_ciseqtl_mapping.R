################################################################################
#
# 1. Prepare files for eQTL mapping
#
################################################################################

options(stringsAsFactors = FALSE)
set.seed("27082021")

library(data.table)
library(peer)
library(RNOmni)
library(edgeR)
library(genefilter)
library(lineup)
library(ggplot2)
library(reshape2)

################################################################################

# read in the RNAseq data with updated HLA counts (all samples)
gex <- read.delim("../../hla/alt.4.3.mhc_mapping.counts/3.new_gene_counts.nM8.60681.txt", row.names=1)
gex.info <- read.delim("../data/sample_key.txt")

# restrict to samples passing QC and update column names with GAinS IDs
gex <- gex[, match(gex.info$SangerSampleID, colnames(gex))]
colnames(gex) <- gex.info$GAinSID

# gene filtering (at least 10 reads in 5% of samples)
f1 <- kOverA(43, 10)
flist <- filterfun(f1)
exprs <- genefilter(gex, flist)
counts.red <- gex[exprs, ]
dim(counts.red) # 20412 864

# normalisation
y <- DGEList(counts=counts.red, samples = gex.info)
y <- calcNormFactors(y, method="TMM")
cpms <- cpm(y, log = F)
logcpm <- log(cpms + 1, 2)
# write.table(logcpm, "../data/logcpm_864_20412_hla.txt", sep="\t")
rm(gex, counts.red, cpms, y, exprs)

# prepare gene position information file
gtf <- as.data.frame(rtracklayer::import("../data/Homo_sapiens.GRCh38.99.gtf"))
gene.info.filt <- gtf[match(rownames(logcpm), gtf$gene_id), ]
rownames(gene.info.filt) <- rownames(logcpm)
gene.info.filt$gene_id <- rownames(logcpm)
gene.info.filt <- gene.info.filt[, c("seqnames", "start", "end", "strand", 
                                     "gene_id", "gene_name")]
gene.info.filt$seqnames <- as.character(gene.info.filt$seqnames)
gene.info.filt$seqnames[gene.info.filt$seqnames == "X"] <- 23
# write.table(gene.info.filt, "../data/gene_info_20412.txt", sep = "\t")
rm(gtf, gene.info.filt)

# Get data on patients with both gene expression and genotyping data available
fam <- read.table("../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.fam", stringsAsFactors = F)
gex.info$PatientID <- substr(gex.info$GAinSID, start = 1, stop = nchar(gex.info$GAinSID)-2)
gex.info$Genotyped <- gex.info$PatientID %in% fam$V1
gex.info.eqtl <- gex.info[gex.info$PatientID %in% fam$V1, ]
gex.eqtl <- logcpm[, gex.info$PatientID %in% fam$V1]
write.table(gex.eqtl, "../data/gex_rnaseq_for_eqtl_hla.txt", sep="\t")
write.table(gex.info.eqtl, "../data/gex_info_rnaseq_for_eqtl", sep="\t")
rm(fam, gex.info, logcpm)

pts <- unique(gex.info.eqtl$PatientID)
pts <- data.frame("FID"=pts, "IID"=pts) # 638
write.table(pts, "../data/rnaseq_patients_with_genotyping.txt", sep="\t", quote=F,
            row.names = F)

system('plink --bfile ../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID --keep ../data/rnaseq_patients_with_genotyping.txt --make-bed --out ../data/genotyping_for_rna-seq_eQTL --allow-extra-chr --maf 0.01')

################################################################################

# Make a list of SNP-gene pairs to test (SNPs 1Mb either side of TSS for each gene)

################################################################################

# read in SNP and gene positional information
snp.info <- data.frame(fread("../data/genotyping_for_rna-seq_eQTL.bim", 
                             header=F, stringsAsFactors = F))
genes <- read.delim("../data/gene_info_20412.txt", stringsAsFactors = F)

# flip start and end for genes on - strand
genes$start[which(genes$strand == "-")] <- genes$end[which(genes$strand == "-")]
genes$end[which(genes$strand == "-")] <- genes$start[which(genes$strand == "-")]

# divide up genes and SNPs by chromosome
genes.list <- vector("list", 22)
chr.bim.list <- vector("list", 22)

for(i in 1:22){
  genes.list[[i]] <- subset(genes, seqnames == i)
  rownames(genes.list[[i]]) <- genes.list[[i]][, "gene_id"]
  
  chr.bim.list[[i]] <- snp.info[which(snp.info$V1 == i), ]
  chr.bim.list[[i]] <- chr.bim.list[[i]][, c(1, 2, 4)]
  rownames(chr.bim.list[[i]]) <- chr.bim.list[[i]][, 2]
}

# make a list of SNP-gene pairs to test
results <- vector("list", 22)
# for each chromosome
for(c in 1:22){
  results.c <- vector("list", nrow(genes.list[[c]]))
  # for each gene get SNPs in a 1Mb window
  for (i in 1:nrow(genes.list[[c]])) {
    snp.names <- rownames(chr.bim.list[[c]])[which(as.numeric(genes.list[[c]][i, 1]) == chr.bim.list[[c]][, 1]
                                                   & (chr.bim.list[[c]][, 3] > (genes.list[[c]][i, 2] - 1000000)
                                                      &  chr.bim.list[[c]][, 3] < (genes.list[[c]][i, 2] + 1000000)))]
    if(length(snp.names) > 1){
      results.c[[i]] <- data.frame(rownames(genes.list[[c]])[i], snp.names, 
                                   rep(as.numeric(genes.list[[c]][i, 1]), length(snp.names)))
    }
  }
  results[[c]] <- rbindlist(results.c)
  colnames(results[[c]]) <- c("Gene", "SNP", "Chr")
}

results.2 <- rbindlist(results)
write.table(results.2, "../data/gene_snp_pairs_cis_rnaseq.txt", sep="\t", 
            quote=FALSE, row.names=FALSE)

snps <- unique(as.character(results.2$SNP))
write.table(snps, "../data/snps_for_cis_eqtl_rnaseq.txt", sep="\t", quote=F, 
            row.names = F, col.names = F)
rm(chr.bim.list, genes, genes.list, results, results.2, results.c, snp.info, snp.names, snps)

# make genotyping files
system("for CHR in {1..22} ; do plink --bfile ../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID --keep ../data/rnaseq_patients_with_genotyping.txt --extract ../data/snps_for_cis_eqtl_rnaseq.txt --recode A --out ../data/genotyping_for_rna-seq_eqtl_${CHR} --allow-extra-chr --chr ${CHR} --maf 0.01 ; done")

################################################################################
# assemble all files

exp <- read.delim("../data/gex_rnaseq_for_eqtl_hla.txt")
exp <- as.matrix(exp)
info <- read.delim("../data/covs_for_eqtl_823_updated_srs_int_cells.txt")
all(colnames(exp) == rownames(info))
GAinSID <- info$GAinSID
load("../peer/peer_factors_covs_int_cells_823_hla.Rda")
covs <- as.matrix(factors[, 1:12])
peer.factors <- as.matrix(factors[, 14:43])

for(i in 1:22){
  geno <- data.frame(fread(paste("../data/genotyping_for_rna-seq_eqtl_", i, ".raw", sep=""),
                           sep=" ", drop = 2:6))
  m1 <- match(GAinSID, geno$FID)
  geno <- geno[m1, ]
  rownames(geno) <- rownames(info)
  colnames(geno) <- gsub("X", "", colnames(geno))
  colnames(geno) <- substr(colnames(geno), 1, nchar(colnames(geno))-2)
  geno[, 1] <- NULL
  geno <- as.matrix(geno)
  
  eQTL_files <- paste("../data/eqtl_files_", i, ".rda", sep="")
  save(list=c("exp", "geno", "covs", "GAinSID", "peer.factors"),
       file = eQTL_files)
  
}
