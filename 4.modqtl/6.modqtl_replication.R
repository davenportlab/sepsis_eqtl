################################################################################
#
# 6. Module QTL replication in microarray data
#
################################################################################# 

# 676 samples from 514 patients, of which 135 duplicates and 541 (382 patients) unique
# 641 samples from 493 patients with genotyping data. 
# 135 (all have genotyping) should be dropped for the replication to give 506 samples, 361 patients
options(stringsAsFactors = FALSE)

library(tidyverse)
library(SummarizedExperiment)
library(data.table)
library(lme4)
library(gridExtra)

exp.micro <- readRDS("../../nikhil/data/GAinS_Microarray/gains_full_microarray_dedup_norm_combat_average-per-gene.rds")
modules <- read.csv("../../nikhil/expression/gene_expression/modules.csv")
eigengenes <- read.csv("../../nikhil/expression/gene_expression/eigengenes.multiple.csv", row.names=1)

micro.info <- data.frame(exp.micro@colData) # 676 samples
length(unique(micro.info$GenotypingID)) # 514 patients, 676 samples
table(micro.info$SampleID %in% rownames(eigengenes))
# 135 duplicated samples in RNAseq, 541 unique

micro.info.unique <- micro.info[(micro.info$SampleID %in% rownames(eigengenes)), ]
length(unique(micro.info.unique$GenotypingID))
# 382 individuals after removing the duplicated samples

load("../peer/microarray/PEER_factors_covs_norm_cells_641_microarray.Rda") # all samples with genotyping
colnames(factors) <- c(colnames(covs), "Dummy", paste0("PEER", 1:30))

# Extract replication genotypes from all samples
# awk -F ',' 'NR > 1 { print $1; }' updated_all_mqtl.csv | sort | uniq > rep.snps.txt
# plink --bfile ../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID \
# --extract rep.snps.txt \
# --recode A \
# --out ./eigengene_sva_rep_genotypes \
# --allow-extra-chr \
# --maf 0.01
# 241 variants (all present)

geno <- fread("../modqtl_rerun/eigengene_sva_rep_genotypes.raw", drop=2:6) %>%
  as.data.frame()

rownames(geno) <- geno$FID
geno$FID <- NULL
snp.names <- colnames(geno) # to check minor allele
colnames(geno) <- gsub("_.$", "", colnames(geno))

mqtl <- read.csv("../modqtl_rerun/updated_all_mqtl.csv")
lead.mqtl <- read.csv("../modqtl_rerun/updated_lead_mqtl.csv")

micro.eigens <- lapply(unique(modules$Module), function(module) {
  
  module.genes = modules$Gene[modules$Module == module]
  if (sum(rowData(exp.micro)$Ensembl_ID %in% module.genes) < 5) {
    return(NULL)
  }

  mod.exp = exp.micro[rowData(exp.micro)$Ensembl_ID %in% module.genes, ]
  module.svd = svd(scale(t(assays(mod.exp)[[1]])))
  module.top.eigens = module.svd$u[, 1:5]
  rownames(module.top.eigens) = colData(mod.exp)$SampleID
  colnames(module.top.eigens) = paste0(gsub("Module_", "Micro_ME_", module), "_", 1:5)

  return(module.top.eigens)
}) %>%
  do.call(cbind, .)

rep.mes <- unique(gsub("_.$", "", gsub("Micro_ME_", "", colnames(micro.eigens))))
length(rep.mes) # 95

all.eigens = merge(micro.eigens, eigengenes, by=0)[, -1]
# 135 samples with both microarray and RNA-seq

rhos <- list()
gg <- list()

for(i in 1:length(rep.mes)){
  me.n <- rep.mes[i]
  module.eigengene = all.eigens[, paste0("ME_", me.n, "_1")]
  module.micro.eigengene = all.eigens[, paste0("Micro_ME_", me.n, "_1")]
  rhos[[i]] <- cor(module.eigengene, module.micro.eigengene, method="spearman")
  df <- data.frame("RNASeq_Module"=module.eigengene,
                   "Microarray_Module"=module.micro.eigengene)
  
  gg[[i]] <- ggplot(df, aes(RNASeq_Module, Microarray_Module)) +
    geom_point() +
    geom_smooth(formula='y~x', method="lm") +
    theme_bw() + ggtitle(label=paste0("ME_", me.n))
}

# pdf("ME_replication_correlation.pdf", onefile=T, useDingbats=F)
do.call(grid.arrange, gg[1:9])
do.call(grid.arrange, gg[10:18])
do.call(grid.arrange, gg[19:27])
do.call(grid.arrange, gg[28:36])
do.call(grid.arrange, gg[37:45])
do.call(grid.arrange, gg[46:54])
do.call(grid.arrange, gg[55:63])
do.call(grid.arrange, gg[64:72])
do.call(grid.arrange, gg[73:81])
do.call(grid.arrange, gg[82:90])
do.call(grid.arrange, gg[91:95])

rhos <- unlist(rhos)
names(rhos) <- paste0("ME_", rep.mes)
hist(rhos, breaks = 100)
# dev.off()

###########
# ModQTL replication
design.mtx <- colData(exp.micro)[, c("GAinSID", "SampleID")] # 676 samples

design.mtx <- merge(design.mtx, micro.eigens, by.x="SampleID", by.y=0) %>%
  merge(., geno, by.x="GAinSID", by.y=0) %>%
  merge(., factors, by.x="SampleID", by.y=0) %>%
  dplyr::filter(!(SampleID %in% rownames(eigengenes)))
# 506 samples not in RNA-seq with genotyping
length(unique(design.mtx$GAinSID)) # 361

covs = c(
  c("SRS", "Diagnosis", "Neutrophils", "Lymphocytes", "Monocytes"),
  paste0("PEER", 1:20), paste0("PC", 1:7)
)

mqtl.reps <- do.call(rbind, lapply(1:nrow(lead.mqtl), function(i) {
  
  snp = as.character(lead.mqtl[i, "snp"])
  
  if (!(snp %in% colnames(geno))) {
    print(as.character(lead.mqtl[i, "me"]))
    print(snp)
    print("SNP not present")
    return(NULL)
  }
  
  module = gsub("ME_", "Module_", lead.mqtl[i, "me"])
  me = paste0(gsub("Module_", "Micro_ME_", module), "_1")
  
  if (!(me %in% colnames(design.mtx))) {
    print(module)
    print("Module not replicable")
    return(NULL)
  }
  
  variant.design <- design.mtx[,c(me, snp, covs, "GAinSID")]
  
  f.null <- as.formula(paste0(me, "~", paste0(covs, collapse="+"), "+(1|GAinSID)"))
  model.null <- lmer(f.null, data=variant.design, REML=FALSE)
  
  f.alt <- as.formula(paste0(me, "~`", snp , "`+", paste0(covs, collapse="+"), "+(1|GAinSID)"))
  model.test <- lmer(f.alt, data=variant.design, REML=FALSE)
  
  if (!all(complete.cases(variant.design[, snp]))) {
    model.null <- update(model.null, subset=complete.cases(variant.design[, snp]))
    model.test <- update(model.test, subset=complete.cases(variant.design[, snp]))
  }
  
  data.frame(matrix(
    data=c(
      me, snp,
      summary(model.test)$coefficients[snp, ],
      anova(model.null, model.test)["model.test", "Pr(>Chisq)"]
    ),
    nrow=1, ncol=6
  ))
})) %>%
  as.data.frame() %>%
  dplyr::select(me=1, snp=2, beta=3, se=4, t=5, p=6) %>%
  dplyr::mutate(beta=as.numeric(beta), se=as.numeric(se), t=as.numeric(t), p=as.numeric(p))

# Of the 32 lead mQTL, only 29 had enough genes to generate a module.
# Module 97, 88, 103 not replicable

table(mqtl.reps$p < 0.05) # 19 replicate
rep.res <- lead.mqtl %>%
  dplyr::ungroup() %>%
  dplyr::select(me, snp, beta.original=beta, se.original=se, me.original=me) %>%
  dplyr::mutate(me=paste0("Micro_", me.original, "_1")) %>%
  merge(., mqtl.reps, by=c("me", "snp")) %>%
  dplyr::mutate(rho=rhos[me.original]) %>%
  dplyr::mutate(p=p.adjust(p, method="BH")) %>%
  dplyr::mutate(significant=p < 0.05) %>%
  dplyr::mutate(matching=sign(beta.original) * sign(rho) == sign(beta))

wilcox.test(abs(rep.res$rho[which(rep.res$significant)]), abs(rep.res$rho[-which(rep.res$significant)]), exact=F)
table(rep.res$significant, rep.res$matching)
# 19 are significant but only 16 are matching direction

# Is minor allele the same?
lead.mqtl$minor.rep <- substr(snp.names, nchar(snp.names), nchar(snp.names))[match(lead.mqtl$snp,
                                                                                   colnames(geno))]
# All same minor allele.
sum(rep.res$matching[rep.res$significant])

# pdf("../modqtl_rerun/FigS18_modQTL_replication_rho.pdf", useDingbats = F)
rep.res %>%
  dplyr::mutate(significant=ifelse(significant & matching, "Replicated", "Not Replicated")) %>%
  ggplot(aes(x=significant, y=abs(rho))) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitter(width=0.2)) +
  scale_fill_manual(values=c("royalblue1", "royalblue4")) +
  guides(fill="none") +
  ylab("|Spearman's Rho|") +
  theme_bw() +
  theme(axis.title.x=element_blank())
# dev.off()

# pdf("../modqtl_rerun/Fig_S17_modQTL_replication_forest_plot.pdf", useDingbats = F, height = 8)
rep.res %>%
  # dplyr::filter(significant, matching) %>%
  dplyr::mutate(me.snp=paste0(me.original, "-", snp)) %>%
  dplyr::mutate(beta=beta * sign(rho)) %>%
  dplyr::select(me.snp, beta.original, se.original, beta, se) %>%
  dplyr::arrange(desc(beta.original)) %>%
  dplyr::mutate(me.snp=factor(me.snp, levels=me.snp)) %>%
  tidyr::pivot_longer(!me.snp, names_to="col.name") %>%
  dplyr::mutate(statistic=gsub(".original", "", col.name)) %>%
  dplyr::mutate(dataset=ifelse(grepl("original", col.name), "RNA-seq", "Microarray")) %>%
  dplyr::select(me.snp, statistic, value, dataset) %>%
  tidyr::pivot_wider(names_from=statistic, values_from=value) %>%
  dplyr::mutate(dataset=factor(dataset, levels=c("RNA-seq", "Microarray"))) %>%
  dplyr::mutate(beta.low=beta-(1.96 * se), beta.high=beta+(1.96*se)) %>%
  ggplot() +
  geom_hline(aes(yintercept=0)) +
  geom_point(aes(x=me.snp, y=beta, color=dataset), position=position_dodge(0.9), size=I(2)) +
  geom_linerange(aes(x=me.snp, ymin=beta.low, ymax=beta.high, color=dataset), position=position_dodge(0.9)) +
  scale_color_manual("Dataset", values=c("darkblue", "cornflowerblue")) +
  facet_grid(me.snp ~ ., scales="free_y") +
  coord_flip() +
  ylab("ModQTL effect size") +
  xlab("Module-SNP pair") +
  theme_bw() +
  theme(
    strip.background.y=element_blank(), strip.text.y=element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
# dev.off()
