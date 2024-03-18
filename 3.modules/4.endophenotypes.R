################################################################################
#
# 4. Module Phenotype Associations
#
################################################################################

if (!requireNamespace("patchwork")) {
  install.packages("devtools")
  devtools::install_github("thomasp85/patchwork")
}

if (!requireNamespace("rstatix")) {
  install.packages("rstatix")
}

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(rstatix)
library(UpSetR)
library(ComplexHeatmap)
library(circlize)

modules <- read.csv("modules.csv")
eigengenes <- read.csv("eigengenes.csv", row.names=1)
variance.explained <- read.csv("variance.explained.csv")

logcpm <- read.delim("logcpm_864_20412_hla.txt")
gene.info <- read.delim("gene_info_20412.txt")

# Phenotypes
srs.info <- read.delim("RNASeq/Covs_for_PEER.txt")
cell.props <- read.delim("Cell_props_823_updated.txt")

srs.info$Neutrophils[is.na(srs.info$Neutrophils)] <- cell.props$Neutrophils[match(srs.info$Sample[is.na(srs.info$Neutrophils)], cell.props$Sample)]
srs.info$Lymphocytes[is.na(srs.info$Lymphocytes)] <- cell.props$Lymphocytes[match(srs.info$Sample[is.na(srs.info$Lymphocytes)], cell.props$Sample)]
srs.info$Monocytes[is.na(srs.info$Monocytes)] <- cell.props$Monocytes[match(srs.info$Sample[is.na(srs.info$Monocytes)], cell.props$Sample)]

srs.info$Sample <- gsub("GA", "", srs.info$Sample)
srs.info <- srs.info[match(rownames(eigengenes), srs.info$Sample), ]
srs.info$supplier_name <- rownames(eigengenes)


##### SRS
# update SRSq to SRS1 vs non-SRS1 to match rest of paper and run all tests as
# LMM for consistency
srs.info$SRS1 <- srs.info$SRS == 1

eigen.srs.data <- srs.info %>%
  dplyr::select(GAinSID, SRS1, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.srs.data[, i+3] ~ SRS1 + (1|GAinSID), data=eigen.srs.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.srs.data[, i+3] ~ (1|GAinSID), data=eigen.srs.data,
                     REML=FALSE)
  test.res <- anova(null, test)
  
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

# write.table(test.results, "Module_SRS_associations.txt", sep="\t")

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="SRS1", Association.Variable.Type="SRS", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()

srs.assns <- associated.eigengenes

#### Diagnosis
eigen.diag.data <- srs.info %>%
  dplyr::select(GAinSID, Diagnosis, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)
eigen.diag.data$Diagnosis <- (eigen.diag.data$Diagnosis -1)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.diag.data[, i+3] ~ Diagnosis + (1|GAinSID), data=eigen.diag.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.diag.data[, i+3] ~ (1|GAinSID), data=eigen.diag.data,
                     REML=FALSE)
  test.res <- anova(null, test)
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="Diagnosis", Association.Variable.Type="Diagnosis", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()

diag.assns <- associated.eigengenes


###### Cell proportions
# Inverse normal transformed cell proportions for three broad leukocyte lineages
# rank inverse normal transform the cell proportions
library(RNOmni)
srs.info$Neutrophils[!(is.na(srs.info$Neutrophils))] <- RankNorm(srs.info$Neutrophils[!(is.na(srs.info$Neutrophils))])
srs.info$Lymphocytes[!(is.na(srs.info$Lymphocytes))] <- RankNorm(srs.info$Lymphocytes[!(is.na(srs.info$Lymphocytes))])
srs.info$Monocytes[!(is.na(srs.info$Monocytes))] <- RankNorm(srs.info$Monocytes[!(is.na(srs.info$Monocytes))])

# replace missing values with mean
srs.info$Neutrophils[is.na(srs.info$Neutrophils)] <- 0
srs.info$Lymphocytes[is.na(srs.info$Lymphocytes)] <- 0
srs.info$Monocytes[is.na(srs.info$Monocytes)] <- 0

# Neutrophils
eigen.neut.data <- srs.info %>%
  dplyr::select(GAinSID, Neutrophils, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.neut.data[, i+3] ~ Neutrophils + (1|GAinSID), data=eigen.neut.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.neut.data[, i+3] ~ (1|GAinSID), data=eigen.neut.data,
                     REML=FALSE)
    test.res <- anova(null, test)
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="Neutrophils", Association.Variable.Type="Cell Proportion", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()
neut.assns <- associated.eigengenes

# Lymphocytes
eigen.lymph.data <- srs.info %>%
  dplyr::select(GAinSID, Lymphocytes, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.lymph.data[, i+3] ~ Lymphocytes + (1|GAinSID), data=eigen.lymph.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.lymph.data[, i+3] ~ (1|GAinSID), data=eigen.lymph.data,
                     REML=FALSE)
  test.res <- anova(null, test)
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="Lymphocytes", Association.Variable.Type="Cell Proportion", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()
lymph.assns <- associated.eigengenes

# Monocytes
eigen.mono.data <- srs.info %>%
  dplyr::select(GAinSID, Monocytes, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.mono.data[, i+3] ~ Monocytes + (1|GAinSID), data=eigen.mono.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.mono.data[, i+3] ~ (1|GAinSID), data=eigen.mono.data,
                     REML=FALSE)
  test.res <- anova(null, test)
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="Monocytes", Association.Variable.Type="Cell Proportion", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()
mono.assns <- associated.eigengenes


#########
# Timepoints

eigen.time.data <- srs.info %>%
  dplyr::select(GAinSID, Day, supplier_name) %>%
  merge(., eigengenes, by.x="supplier_name", by.y=0)

test.results <- list()
test.results[["P.Value"]] <- list()
test.results[["Beta"]] <- list()

for (i in 1:ncol(eigengenes)) {
  eigengene <- colnames(eigengenes)[i]
  test = lme4::lmer(eigen.time.data[, i+3] ~ Day + (1|GAinSID), data=eigen.time.data,
                    REML = FALSE)
  null =  lme4::lmer(eigen.time.data[, i+3] ~ (1|GAinSID), data=eigen.time.data,
                     REML=FALSE)
  test.res <- anova(null, test)
  test.results[["Beta"]][[eigengene]] <- summary(test)$coefficients[2, 1]
  test.results[["P.Value"]][[eigengene]] <- test.res$`Pr(>Chisq)`[2]
}

test.results[["P.Value"]] <- unlist(test.results[["P.Value"]])
test.results[["Beta"]] <- unlist(test.results[["Beta"]])

associated.eigengenes <- as.data.frame(test.results) %>%
  dplyr::mutate(Eigengene=rownames(.)) %>%
  dplyr::mutate(Association.Variable="Day", Association.Variable.Type="Day", Statistic.Type="Beta") %>%
  dplyr::mutate(Adjusted.P.Value=p.adjust(P.Value, method="bonferroni")) %>%
  dplyr::select(Eigengene, Association.Variable, Association.Variable.Type, Statistic=Beta, Statistic.Type, P.Value, Adjusted.P.Value) %>%
  dplyr::arrange(P.Value)

head(associated.eigengenes)
ggplot(associated.eigengenes, aes(Statistic, -log10(Adjusted.P.Value))) + geom_point() + theme_bw()
time.assns <- associated.eigengenes

########
# Survival: The 28-day survival of patients was measured in the GAinS cohort. 
# Association between this patient outcome and each eigengene was tested using a 
# Cox proportional hazards model as implemented in the survival R package using 
# the coxph function. For each patient, the value of the eigengene at the last 
# time point recorded was used as a predictor for the survival function.

# Use Nikhil's results
mod.assns <- read.delim("estimates.all.csv")
surv.assns <- subset(mod.assns, Association.Variable == "Outcome")
surv.assns$Statistic.Type <- "CoxPH"

# Heatmap
association.tbl <- dplyr::bind_rows(
  neut.assns,
  lymph.assns,
  mono.assns,
  srs.assns,
  time.assns,
  diag.assns,
  surv.assns
)
rownames(association.tbl) <- NULL
association.tbl$Adjusted.P.Value <- p.adjust(association.tbl$P.Value, method="bonferroni")
write.table(association.tbl, "module_endophenotype_associations.txt", sep="\t")

cell.prop.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Cell Proportion") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Association.Variable, Statistic) %>%
  tidyr::spread(Association.Variable, Statistic)

cell.prop <- as.matrix(cell.prop.data[, c("Neutrophils", "Lymphocytes", "Monocytes")])
rownames(cell.prop) <- cell.prop.data$Eigengene

srs.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "SRS") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, SRS=Statistic) %>%
  dplyr::arrange(Eigengene)

srsq <- as.matrix(srs.data[, "SRS", drop=F])
rownames(srsq) <- srs.data$Eigengene

diagnosis.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Diagnosis") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Diagnosis=Statistic) %>%
  dplyr::arrange(Eigengene)

diagnosis <- as.matrix(diagnosis.data[, "Diagnosis", drop=F])
rownames(diagnosis) <- diagnosis.data$Eigengene

outcome.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Outcome") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Outcome=Statistic) %>%
  dplyr::arrange(Eigengene)

outcome <- as.matrix(outcome.data[, "Outcome", drop=F])
rownames(outcome) <- outcome.data$Eigengene

time.point.data <- association.tbl %>%
  dplyr::mutate(Eigengene = factor(Eigengene, levels=paste0("ME_", 1:length(unique(Eigengene))))) %>%
  dplyr::filter(Association.Variable.Type == "Day") %>%
  dplyr::mutate(Statistic = ifelse(Adjusted.P.Value < 5e-2, Statistic, NA)) %>%
  dplyr::select(Eigengene, Time.Point=Statistic) %>%
  dplyr::arrange(Eigengene)

time.point <- as.matrix(time.point.data[, "Time.Point", drop=F])
rownames(time.point) <- time.point.data$Eigengene

col.fun.beta <- colorRamp2(breaks=seq(-0.06, 0.06, length.out=11), 
                           colors=rev(brewer.pal("RdBu", n=11)))
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
options(repr.plot.width=5, repr.plot.height=20)

# pdf("FigS15_module_association_heatmap.pdf", width=6, height=20)
h1 + h2 + h3 + h4 + h5
# dev.off()

########
# xCell heatmap
xcell.info <- read.csv("xCell_Aran_et_al_Additional_File_1.csv")
filtered.cell.type.enrichment <- read.csv("module.xcell.signature.enrichment.csv")

plot.data <- filtered.cell.type.enrichment %>%
  dplyr::mutate(Median.Log.2.Enrichment=log2(Median.Enrichment))  %>%
  merge(., xcell.info, by.x="Cell.Type", by.y="Cell.types") %>%
  dplyr::filter(Subgroup %in% c("HSC", "Lymphoid", "Myeloid"))

h.data <- plot.data %>%
  dplyr::select(Module, Full.name, Median.Enrichment) %>%
  tidyr::spread(Module, Median.Enrichment, fill=0)
rownames(h.data) <- h.data$Full.name
h.data$Full.name <- NULL
h.data <- as.matrix(h.data)

h <- hclust(dist(t(h.data)))
h.o <- hclust(dist(h.data))

# pdf("Fig_S15_Xcell_marker_enrichment.pdf", 
#     width=15.5, height=12)

plot.data %>%
    dplyr::mutate(Module=factor(Module, levels=colnames(h.data)[h$order])) %>%
    dplyr::mutate(Full.Name=factor(Full.name, levels=rownames(h.data)[h.o$order])) %>%
    ggplot() +
    geom_tile(aes(x=Module, y=Full.name, fill=Median.Log.2.Enrichment)) +
    scale_fill_distiller(palette="Blues", direction=1, limits=c(0, max(plot.data$Median.Log.2.Enrichment))) +
    ylab("xCell Signature") + xlab("Module") +
    labs(fill=bquote("Median Log"[2]*"(Enrichment)")) +
    facet_grid(Subgroup~., scales="free_y", space="free_y") +
    theme_bw(base_size=18) +
    theme(panel.border=element_blank(),
          panel.grid=element_blank(),
          axis.line.x.bottom=element_line(color="black", size=0.25),
          axis.line.y.left=element_line(color="black", size=0.25),
          legend.position="bottom",
          strip.background=element_rect(fill="#EEEEEE", color="white", size=0.25)) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), 
    legend.position="right", panel.grid.major=element_line(size=0.05))
# dev.off()

# And sepsis 
cibersort.info <- read.csv("kwok_2022_cell_type.csv")
cibersort.markers <- read.delim("kwok_2022_sepsis_cell_markers.txt")

cibersort.markers.enrichment <- read.csv("sepsis_cell_markers.csv")

plot.data <- cibersort.markers.enrichment %>%
  dplyr::mutate(Log.2.Enrichment=log2(Enrichment)) %>%
  merge(., cibersort.info, by.x="Cell.Type", by.y="ID")

h.data <- plot.data %>%
  dplyr::select(Module, Name, Enrichment) %>%
  tidyr::spread(Module, Enrichment, fill=0)
rownames(h.data) <- h.data$Name
h.data$Name <- NULL
h.data <- as.matrix(h.data)

h <- hclust(dist(t(h.data)))
h.o <- hclust(dist(h.data))

# pdf("~/Documents/GAinS_eQTL/paper/figures/Fig_4A_sepsis_cell_marker_enrichment.pdf",
#     width=16, height=8)

plot.data %>%
  dplyr::mutate(Module=factor(Module, levels=colnames(h.data)[h$order])) %>%
  dplyr::mutate(Name=factor(Name, levels=rownames(h.data)[h.o$order])) %>%
  # dplyr::filter(gsub("Module", "ME", Module) %in% lead.modqtl$me) %>%
  ggplot() +
  geom_tile(aes(x=Module, y=Name, fill=Log.2.Enrichment)) +
  scale_fill_distiller(palette="Blues", direction=1, limits=c(0, max(plot.data$Log.2.Enrichment))) +
  ylab("Cell Type") + xlab("Module") +
  labs(fill=bquote("Log"[2]*"(Enrichment)")) +
  facet_grid(Lineage~., scales="free_y", space="free_y") +
  theme_bw(base_size=18) +
  theme(panel.border=element_blank(),
        panel.grid=element_blank(),
        axis.line.x.bottom=element_line(color="black", size=0.25),
        axis.line.y.left=element_line(color="black", size=0.25),
        legend.position="bottom",
        strip.background=element_rect(fill="#EEEEEE", color="white", size=0.25)) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="right", panel.grid.major=element_line(size=0.05))

# dev.off()
