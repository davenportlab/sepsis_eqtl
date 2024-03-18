################################################################################
#
# 6. Combine conditional results and explore
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

cond.inpt <- read.delim("../cisresults/conditionalanalysis/conditional_eQTL_results_final.txt")
cond.inpt$pair <- paste0(cond.inpt$SNP, "_", cond.inpt$Gene)
thresholds <- read.delim("../cisresults/nominal_pval_thresholds.txt")

results <- cond.inpt[!duplicated(cond.inpt$pair), ]
# write.table(results, "../cisresults/conditionalanalysis/conditional_eQTL_results_final.txt",
#             sep="\t")

################################################################################

# We then generated summary statistics for each signal, conditioning on the 
# other conditional signal SNPs for that gene (see scripts 5a and 5b)
# These are the most accurate beta estimates and the final output files

################################################################################

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

# Summary of results

# How many independent associations do we have for each gene?
table(table(cond.inpt$Gene))
table(table(cond.inpt$Gene)>1)
length(unique(cond.inpt$SNP))

n.sig <- data.frame(table(cond.inpt$Gene))
colnames(n.sig) <- c("Gene", "NumberAssns")
median(n.sig$NumberAssns)
n.sig$NumberAssns <- as.factor(n.sig$NumberAssns)

pdf("Fig_1D_conditionalsignals_histogram.pdf")
ggplot(n.sig, aes(NumberAssns)) +
  geom_bar(fill="lightgrey", colour="black") +
  theme_bw() +
  coord_cartesian(expand = FALSE) +
  xlab("Number of independent associations\n per gene") +
  ylab("Number of genes")
dev.off()

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
# table(gene.info$start < gene.info$end) # Always therefore not flipped for strand
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

cond.inpt$Signal <- cond.inpt$Rank
cond.inpt$Signal[cond.inpt$Signal >2] <- "3+"

pdf("Fig_S4_TSS_dist_density.pdf", width=7, height=5)
ggplot(cond.inpt, aes(dist, group=Signal, colour=Signal)) +
  geom_density(alpha=0.2) +
  theme_bw() +
  scale_x_continuous(labels = scales::comma) +
  xlab("Distance from TSS to lead eSNP")
dev.off()

write.table(cond.inpt, "../cisresults/conditionalanalysis/conditional_eQTL_results_summary_stats.txt",
            sep="\t", quote=F)
