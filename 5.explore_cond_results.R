################################################################################
#
# 5. Combine conditional results
#
################################################################################

options(stringsAsFactors = FALSE)
library(data.table)
library(ggplot2)

################################################################################

res.path <- "../cisresults/conditionalanalysis/"
file.pattern <- "backward*"
temp_eqtl <- list.files(path=res.path, pattern = file.pattern, 
                        full.names = F, recursive = T)
temp_eqtl <- unlist(temp_eqtl)
for(i in 1:22){
  temp_eqtl <- gsub(paste0("chr", i, "/backward_"), "", temp_eqtl)
}
temp_eqtl <- gsub(".rds", "", temp_eqtl)

################################################################################

# Should include all significant eGenes (10618)

# Check if any jobs failed
# primary.peak.results <- read.delim("../cisresults/eigenMT/ciseqtl_eigenMT_corrected.txt")
# results.avail <- primary.peak.results$gene %in% temp_eqtl
# fails <- which(results.avail == FALSE & primary.peak.results$Sig == TRUE)
# length(fails)
# write.table(fails, "Fails.txt", sep="\t", row.names=F, quote=F, col.names=F)

# # If failures are found, resubmit failed jobs and redo combining results
# for i in `cat Fails.txt`
# do
#   bsub -q "long" -o logs/conditionalciseqtl.${i}.stdout -e logs/conditionalciseqtl.${i}.stderr -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" -M4000 /software/R-3.6.1/bin/Rscript --vanilla 4.conditionalanalysis.R ${i} 20
# done

################################################################################

# combine results
cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, 
                          recursive=T, full.names = T)
  print(length(temp_eqtl))
  eqtl <- list()
  eqtl <- lapply(temp_eqtl, function(x){
    x <- readRDS(x)
    x <- x[!(x$SNP == "NA"), ]
    x$Number <- 1:nrow(x)
    return(x)
  })
  eqtl <- rbindlist(eqtl)
  eqtl <- data.frame(eqtl)
  print(dim(eqtl))
  saveRDS(eqtl, output_filename)
  return(eqtl)
}

forward <- cat.rds("../cisresults/conditionalanalysis/", "forward*", 
                   "../cisresults/conditionalanalysis/conditional_eQTL_results_forward.rds")
write.table(forward, "../cisresults/conditionalanalysis/conditional_eQTL_results_forward.txt", sep="\t")

backward <- cat.rds("../cisresults/conditionalanalysis/", "backward*", 
                    "../cisresults/conditionalanalysis/conditional_eQTL_results_final.rds")
write.table(backward, "../cisresults/conditionalanalysis/conditional_eQTL_results_final.txt", sep="\t")

################################################################################

# Nikhil then generated summary statistics for each signal, conditioning on the 
# other conditional signal SNPs for that gene
# These are the most accurate beta estimates: update final output file with these
# Also select those passing significance threshold?

cond.inpt <- read.delim("../cisresults/conditionalanalysis/conditional_eQTL_results_final.txt")
cond.inpt$pair <- paste0(cond.inpt$SNP, "_", cond.inpt$Gene)
thresholds <- read.delim("../cisresults/nominal_pval_thresholds.txt")

results <- cond.inpt[!duplicated(cond.inpt$pair), ]
# write.table(results, "../cisresults/conditionalanalysis/conditional_eQTL_results_final.txt",
#             sep="\t")

results.chr <- list()
all.sig <- list()
for(i in 1:22){
  summ.stats <- read.delim(paste0("../../nikhil/colocalization/cis_eqtl/conditional_effects/LMM/chr", 
                                  i, "_conditional_cis_eQTL_summary_statistics.tsv"))
  # summ.stats$SNP <- gsub("_A_C", "", summ.stats$SNP)
  summ.stats$pair <- paste0(summ.stats$SNP, "_", summ.stats$Gene)
  results.chr.i <- merge(results, summ.stats, by="pair")
  results.chr.i <- results.chr.i[order(results.chr.i$Number), ]
  results.chr[[i]] <- results.chr.i[!duplicated(results.chr.i$pair), ]
  summ.stats$threshold <- thresholds$threshold[match(summ.stats$Gene, thresholds$gene)]
  all.sig[[i]] <- summ.stats[which(summ.stats$P_Value < summ.stats$threshold), ]
}
final.summ.stats <- rbindlist(results.chr)
all.sig <- rbindlist(all.sig)
all.sig$pair <- NULL
saveRDS(all.sig, "all_significant_associations_conditional.rds")
all.sig$signal <- paste0(all.sig$Gene, "_", all.sig$Signal)
n.per.signal <- data.frame(table(all.sig$signal))
n.per.signal$Var1 <- as.character(n.per.signal$Var1)
n.per.signal$rank <- as.factor(all.sig$Signal[match(n.per.signal$Var1, all.sig$signal)])
ggplot(n.per.signal, aes(rank, Freq)) +
  geom_boxplot() +
  theme_bw()

ggplot(final.summ.stats, aes(eQTL_beta, Beta)) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1, intercept = 0)

cond.inpt <- final.summ.stats[, c("SNP.x", "Gene.x", "Beta", "SE", "P_Value", "Number")]
colnames(cond.inpt) <- c("SNP", "Gene", "Beta", "SE", "P_Value", "Rank")

####################################################################################

# summary of results
table(table(cond.inpt$Gene))
table(table(cond.inpt$Gene)>1)

pdf("../cisresults/conditionalanalysis/conditional_analysis_plots.pdf", onefile = T, useDingbats = F)
library(ggplot2)
n.sig <- data.frame(table(cond.inpt$Gene))
colnames(n.sig) <- c("Gene", "NumberAssns")
ggplot(n.sig, aes(NumberAssns)) +
  geom_histogram(bins=10, fill="lightgrey", colour="black") +
  theme_bw() +
  ggtitle("Number of independent signals per gene")
median(n.sig$NumberAssns)

ggplot(cond.inpt, aes(Rank, Beta)) +
  geom_point(position=position_jitter(width=0.2)) +
  theme_bw() +
  ggtitle("Effect size across conditionally independent signals")

ggplot(cond.inpt, aes(Beta, -log10(P_Value))) +
  geom_point(aes(colour=as.factor(Rank))) + theme_bw() +
  ggtitle("Effect size across conditionally independent signals")

cond.forward <- read.delim("../cisresults/conditionalanalysis//conditional_eQTL_results_forward.txt")
ggplot(cond.forward, aes(Number, eQTL_beta)) +
  geom_point(position=position_jitter(width=0.2)) +
  theme_bw() +
  ggtitle("Forward step of conditional analysis")

ggplot(cond.forward, aes(eQTL_beta, -log10(pvalue))) +
  geom_point() +
  theme_bw() +
  ggtitle("Forward step of conditional analysis")

cond.forward <- cond.forward[match(paste0(cond.inpt$SNP, cond.inpt$Gene),
                                   paste0(cond.forward$SNP, cond.forward$Gene)), ]

plot(cond.forward$eQTL_beta, cond.inpt$Beta)
plot(cond.forward$Number, cond.inpt$Rank)

# calculate distance between SNP and TSS: positive value if downstream
bim <- fread("../data/genotyping_for_rna-seq_eQTL.bim")
bim$V2 <- gsub(":", ".", bim$V2)
cond.inpt$SNPpos <- bim$V4[match(cond.inpt$SNP, bim$V2)]
gene.info <- read.delim("../data/gene_info_20412.txt")
# table(gene.info$start < gene.info$end) # Yes therefore not flipped for strand
gene.info$TSS <- gene.info$start
gene.info$TSS[which(gene.info$strand == "-")] <- gene.info$end[which(gene.info$strand == "-")]
cond.inpt$TSS <- gene.info$TSS[match(cond.inpt$Gene, gene.info$gene_id)]
gene.info$TES <- gene.info$end
gene.info$TES[which(gene.info$strand == "-")] <- gene.info$start[which(gene.info$strand == "-")]
cond.inpt$TES <- gene.info$TES[match(cond.inpt$Gene, gene.info$gene_id)]

cond.inpt$dist <- cond.inpt$SNPpos - cond.inpt$TSS

# for genes on the - strand the sign of the difference needs to be flipped
cond.inpt$strand <- gene.info$strand[match(cond.inpt$Gene, gene.info$gene_id)]
cond.inpt$dist[which(cond.inpt$strand == "-")] <- cond.inpt$dist[which(cond.inpt$strand == "-")]*-1
range(cond.inpt$dist)
cond.inpt$Rank <- as.factor(cond.inpt$Rank)
ggplot(cond.inpt, aes(dist, group=Rank, fill=Rank)) +
  geom_density(alpha=0.2) +
  theme_bw() +
  ggtitle("Distance from TSS to SNP") +
  geom_vline(xintercept = 0, lty=2)

ggplot(cond.inpt, aes(Rank, abs(dist))) +
  geom_boxplot() +
  # geom_point() +
  theme_bw() + 
  ggtitle("Distance from TSS to SNP")

ggplot(cond.inpt, aes(Rank, abs(dist))) +
  geom_violin(trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75)) +
  # geom_point() +
  theme_bw() + 
  ggtitle("Distance from TSS to SNP")

ggplot(cond.inpt, aes(dist, abs(Beta))) +
  geom_point() +
  theme_bw() +
  facet_grid(Rank ~ .) + 
  ggtitle("Distance from TSS to SNP")

eigenMT <- read.delim("../cisresults/nominal_pval_thresholds.txt")
ggplot(eigenMT, aes(1, -log10(threshold))) +
  geom_boxplot() +
  theme_bw()
cond.inpt$ntests <- eigenMT$n.tests[match(cond.inpt$Gene, eigenMT$gene)]
cond.inpt.max <- cond.inpt[order(cond.inpt$Rank, decreasing = T), ]
cond.inpt.max <- cond.inpt.max[!duplicated(cond.inpt.max$Gene), ]
ggplot(cond.inpt.max, aes(Rank, ntests)) +
  theme_bw() +
  geom_boxplot() +
  ggtitle("Estimated number of tests performed per gene vs number of signals detected")
# geom_point(position=position_jitter(width=0.2))
ggplot(cond.inpt.max, aes(Rank, ntests)) +
  theme_bw() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  ggtitle("Estimated number of tests performed per gene vs number of signals detected")

snp.frq <- read.delim("../../Genotyping/gains_b38.afreq")
cond.inpt$MAF <- snp.frq$ALT_FREQS[match(cond.inpt$SNP, snp.frq$ID)]
ggplot(cond.inpt, aes(Rank, MAF)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Allele frequency against order of signal detection")

ggplot(cond.inpt, aes(Beta, MAF)) +
  geom_point(aes(colour=as.numeric(Rank)), alpha=0.5) +
  theme_bw() +
  scale_color_gradient2(low = "blue", mid = "gold", high = "red", midpoint = 5)
dev.off()

write.table(cond.inpt, "../cisresults/conditionalanalysis/conditional_eQTL_results_summary_stats.txt",
            sep="\t", quote=F)
