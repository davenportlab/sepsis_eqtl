################################################################################
#
# 2. Module Pathway annotation
#
################################################################################

library(tidyverse)
library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata"

options(stringsAsFactors = FALSE)

# Load data
gene.info <- read.table("gene_info_20412.txt")
modules <- read.csv("modules.csv")
eigengenes <- read.csv("eigengenes.csv", row.names=1)
variance.explained <- read.csv("variance.explained.csv")

# Setup
module.names <- paste0("Module_", 1:dim(eigengenes)[2])
module.list <- lapply(module.names, function(module.name) {
  module.info = modules %>%
    dplyr::filter(Module==module.name) %>%
    merge(., gene.info, by.x="Gene", by.y="gene_id") %>%
    dplyr::select(Ensembl.ID=Gene, Gene.Name=gene_name)
})
names(module.list) <- module.names

background = modules$Symbol
full.gene.id.map <- do.call(rbind, module.list)

# Reactome term enrichment
reactome.all <- lapply(names(module.list), function(module) {
  
  eTerm <- xEnricherGenes(module.list[[module]]$Gene.Name, background=background,
                          ontology = "MsigdbC2REACTOME")
  xEnrichViewer(eTerm, top_num = "auto", details = T) %>%
    as.data.frame() %>%
    dplyr::mutate(Module = module) %>%
    dplyr::select(Module, everything())
}) %>%
  do.call(rbind, .)

reactome.all.filtered <- reactome.all %>%
  dplyr::select(
    Module, Pathway=name,
    Number_of_member_genes=nAnno, Number_in_module=nOverlap,
    Fold_Change=fc, P_Value=pvalue, Adjusted_P_Value=adjp,
    Module_genes=members_Overlap, Member_genes=members_Anno) %>% 
  dplyr::filter(Adjusted_P_Value < 0.05)

write.csv(reactome.all.filtered, "WGCNA_Reactome_XGR.csv", row.names=F)
