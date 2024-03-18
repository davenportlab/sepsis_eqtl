################################################################################
#
# 3. Module regulon annotation
#
################################################################################

library(data.table)
library(pheatmap)
library(decoupleR)
options(stringsAsFactors = FALSE)

modules <- read.csv("modules.csv")

# Dorothea regulons
net <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))

# all genes for background
gene.info <- read.delim("gene_info_20412.txt")
# limit regulons to measured genes
net <- subset(net, target %in% gene.info$gene_name)
gene.info$regulon <- gene.info$gene_name %in% net$target

regulon.enrichment <- list()

# Test each module for enrichment of each regulon
for(i in 1:106){
  module <- paste0("Module_", i)
  genes <- subset(modules, Module == module)
  
  regulon.enrichment.module <- data.frame(Module=module,
                                             TF=unique(net$source),
                                             OR=NA,
                                             Pval=NA)
  
  for(r in 1:length(unique(net$source))){
    tf <- unique(net$source)[r]
    regulon <- subset(net, source == tf)
    # only test regulons with at least 5 genes assessed
    if(nrow(regulon)>4 & nrow(genes)>4){
      cont.table <- table(unique(gene.info$gene_name) %in% genes$Symbol,
                          unique(gene.info$gene_name) %in% regulon$target)
      # Add a filter for number of genes
      if(dim(cont.table)[2] == 2){
        if(cont.table[2, 2] >2){
          results <- fisher.test(cont.table, alternative = "greater")
          regulon.enrichment.module$OR[r] <- results$estimate
          regulon.enrichment.module$Pval[r] <- results$p.value
        }
      }
    }
  }
  regulon.enrichment.module <- regulon.enrichment.module[complete.cases(regulon.enrichment.module),]
  regulon.enrichment.module$FDR <- p.adjust(regulon.enrichment.module$Pval, method="fdr")
  regulon.enrichment[[module]] <- regulon.enrichment.module
}

enrichment.results <- rbindlist(regulon.enrichment)
sig.results <- subset(enrichment.results, FDR<0.05)

enrichment.results$FDRall <- p.adjust(enrichment.results$Pval, method="fdr")
subset(enrichment.results, FDRall<0.05)
sig.results <- subset(enrichment.results, FDRall<0.05)
