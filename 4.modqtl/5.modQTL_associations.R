################################################################################
#
# 5. Module QTL Phenotype Associations etc
#
################################################################################

# Module QTL only phenotype heatmap

all.modqtl <- read.csv("updated_all_mqtl.csv")

modules <- read.csv("modules.csv")
modules$Module <- gsub("Module", "ME", modules$Module)
gene.info <- read.delim("gene_info_20412.txt")
modules$Symbol <- gene.info$gene_name[match(modules$Gene, gene.info$gene_id)]

for(i in 1:nrow(all.modqtl)){
  module <- all.modqtl$me[i]
  symbols <- modules$Symbol[modules$Module == module]
  all.modqtl$Genes[i] <- paste0(symbols, collapse=",")
}

modqtl.anno <- read.delim("conditional_eQTL_results_final.txt")
modqtl.anno$symbol <- gene.info$gene_name[match(modqtl.anno$Gene, gene.info$gene_id)]

for(i in 1:nrow(all.modqtl)){
  snp <- all.modqtl$snp[i]
  egenes <- modqtl.anno$symbol[which(modqtl.anno$SNP == snp)]
  egenes <- unique(egenes)
  all.modqtl$cisegene[i] <- ifelse(all(is.na(egenes)), "NA", paste0(egenes, collapse=","))
  all.modqtl$egeneinmodule[i] <- any(egenes %in% unlist(strsplit(all.modqtl$Genes[i], ",")))
}

lead.modqtl <- all.modqtl[order(all.modqtl$p), ]
lead.modqtl <- lead.modqtl[!duplicated(paste(lead.modqtl$me, lead.modqtl$chr)), ]

for(i in 1:nrow(lead.modqtl)){
  cisegenes <- all.modqtl$cisegene[which(all.modqtl$me == lead.modqtl$me[i])]
  lead.modqtl$cisegenes[i] <- paste0(cisegenes, collapse=",")
}
table(all.modqtl$me, all.modqtl$egeneinmodule)
table(lead.modqtl$me, lead.modqtl$egeneinmodule)

mod.assns <- subset(association.tbl, Adjusted.P.Value < 0.05)
table(mod.assns$Association.Variable.Type)
mod.assns <- subset(mod.assns, Eigengene %in% all.modqtl$me)

for(i in 1:nrow(all.modqtl)){
  me <- all.modqtl$me[i]
  assns <- mod.assns$Association.Variable[which(mod.assns$Eigengene == me)]
  assns <- unique(assns)
  all.modqtl$assns[i] <- ifelse(all(is.na(assns)), "NA", paste0(assns, collapse=","))
}

# ModQTL endophenotype heatmap
association.tbl <- read.delim("module_endophenotype_associations.txt")

cell.prop.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Cell Proportion") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Association.Variable, Statistic) %>%
  tidyr::spread(Association.Variable, Statistic)
cell.prop.data <- subset(cell.prop.data, Eigengene %in% all.modqtl$me)
cell.prop <- as.matrix(cell.prop.data[, c("Neutrophils", "Lymphocytes", "Monocytes")])
rownames(cell.prop) <- cell.prop.data$Eigengene

srs.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "SRS") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, SRS=Statistic) %>%
  dplyr::arrange(Eigengene)
srs.data <- subset(srs.data, Eigengene %in% all.modqtl$me)
srsq <- as.matrix(srs.data[, "SRS", drop=F])
rownames(srsq) <- srs.data$Eigengene

diagnosis.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Diagnosis") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Diagnosis=Statistic) %>%
  dplyr::arrange(Eigengene)
diagnosis.data <- subset(diagnosis.data, Eigengene %in% all.modqtl$me)
diagnosis <- as.matrix(diagnosis.data[, "Diagnosis", drop=F])
rownames(diagnosis) <- diagnosis.data$Eigengene

outcome.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Outcome") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Outcome=Statistic) %>%
  dplyr::arrange(Eigengene)
outcome.data <- subset(outcome.data, Eigengene %in% all.modqtl$me)
outcome <- as.matrix(outcome.data[, "Outcome", drop=F])
rownames(outcome) <- outcome.data$Eigengene

time.point.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Day") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Time.Point=Statistic) %>%
  dplyr::arrange(Eigengene)
time.point.data <- subset(time.point.data, Eigengene %in% all.modqtl$me)
time.point <- as.matrix(time.point.data[, "Time.Point", drop=F])
rownames(time.point) <- time.point.data$Eigengene

col.fun.beta <- colorRamp2(seq(-0.06, 0.06, length.out=11), 
                           rev(brewer.pal("RdBu", n=11)))
col.fun.coxph <- colorRamp2(seq(-20, 20, length.out=11), 
                            rev(brewer.pal("BrBG", n=11)))

h1 <- Heatmap(
  cell.prop, name="Beta",
  col=col.fun.beta, na_col="white",
  rect_gp=gpar(col="white", lwd=2),
  cluster_rows=F, cluster_columns=F
)

h2 <- Heatmap(
  srsq, name="Beta",
  col=col.fun.beta, na_col="white",
  rect_gp=gpar(col="white", lwd=2)
)

h3 <- Heatmap(
  diagnosis, name="Beta",
  col=col.fun.beta, na_col="white",
  rect_gp=gpar(col="white", lwd=2)
)

h4 <- Heatmap(
  time.point, name="Beta",
  col=col.fun.beta, na_col="white",
  rect_gp=gpar(col="white", lwd=2)
)

h5 <- Heatmap(
  outcome, name="LogOdds",
  col=col.fun.coxph, na_col="white",
  rect_gp=gpar(col="white", lwd=2)
)
# pdf("Fig4B_module_association_heatmap_modqtl.pdf",
# width=4, height=6, useDingbats = F)
h1 + h2 + h3 + h4 + h5
# dev.off()

# Combine all information
ebi <- read.delim("modqtl_gwas_overlap.txt")
lead.modqtl.ebi <- lead.modqtl %>%
  merge(., ebi, by="me") %>%
  group_by(me, snp) %>%
  summarise(gwas = toString(DISEASE.TRAIT)) %>%
  ungroup()

lead.modqtl.all <- lead.modqtl %>%
  merge(., lead.modqtl.ebi, by="snp",all.x = T)

# add in pathways
lead.modqtl.all$Reactome <- reactome$Reactome[match(lead.modqtl.all$me.x, reactome$Module)]
lead.modqtl.all$assns <- all.modqtl$assns[match(lead.modqtl.all$me.x, all.modqtl$me)]

lead.modqtl.all <- lead.modqtl.all[, c("me.x",	"snp",	"beta",	"p",	"chr", "Genes",
                                       "cisegenes",	"egeneinmodule",	"assns", "gwas", "Reactome")]

write.table(lead.modqtl.all, "ModQTL_all_information_update.txt", sep="\t")
