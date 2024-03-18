################################################################################
#
# 2a Compare eQTL results after including increasing number of PEER factors in the model
#
################################################################################

library(data.table)
library(ggplot2)
options(stringsAsFactors = FALSE)

# Function to read in set of eQTL results
cat.rds <- function(my.path, input_file_pattern, output_filename){
  temp_eqtl <- list.files(path=my.path, pattern = input_file_pattern, full.names = T)
  print(paste("Number of files:", length(temp_eqtl)))
  eqtl <- lapply(temp_eqtl, readRDS)
  eqtl <- rbindlist(eqtl)
  print(c("eQTL results:", dim(eqtl)))
  saveRDS(eqtl, output_filename)
}

n.peer <- seq(0, 30, 5)

for(i in n.peer){
  cat.rds(paste0("results/peer", i), "*.rds",
          paste0("results/chr1_results_", i, "peer_factors.rds"))
}

peer.res <- data.frame("NPeerFactors"=n.peer,
                       "NeGenes"=NA,
                       "NSigPairs"=NA)

for(i in n.peer){
  eqtl <- readRDS(paste0("results/chr1_results_", i, "peer_factors.rds"))
  print(length(unique(eqtl$Gene)))
  eqtl.o <- eqtl[order(eqtl$eQTL_pval), ]
  eqtl.o$padj <- p.adjust(eqtl.o$eQTL_pval, method="bonferroni")
  print(table(eqtl.o$padj < 0.05))
  peer.res[p, 3] <- length(which(eqtl.o$padj < 0.05))
  eqtl.u <- subset(eqtl.o, !duplicated(Gene))
  print(table(eqtl.u$padj < 0.05))
  peer.res[p, 2] <- length(which(eqtl.u$padj < 0.05))
}

peer.res$DiffeGenes <- c(peer.res$NeGenes[1], diff(peer.res$NeGenes, lag=1))
peer.res$DiffSigPairs <- c(peer.res$NSigPairs[1], diff(peer.res$NSigPairs, lag=1))
write.table(peer.res, "chr1_eqtl_number_by_peer.txt", sep="\t")

p1 <- ggplot(peer.res, aes(NPeerFactors, NeGenes)) +
  geom_point() +
  geom_text(aes(y=NeGenes + 10, label=NeGenes)) +
  geom_line() +
  theme_bw() +
  xlab("Number of PEER factors included") +
  ylab("Number of chr1 eGenes identified")

p2 <- ggplot(peer.res, aes(NPeerFactors, DiffeGenes)) +
  geom_point() +
  geom_text(aes(y=DiffeGenes + 5, label=DiffeGenes)) +
  geom_line() +
  theme_bw() +
  ylim(c(0, 330)) +
  xlab("Number of PEER factors included") +
  ylab("Increase in number of chr1 eGenes identified")

pdf("FigS3_PEER_factor_test.pdf", width=10, height=5)
gridExtra::grid.arrange(p1, p2, nrow=1)
dev.off()
