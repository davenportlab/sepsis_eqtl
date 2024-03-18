################################################################################
#
# 4. Module QTL summary
#
################################################################################

# Module QTL
library(data.table)
library(vctrs)
library(tidyverse)
library(parallel)
library(GenomicRanges)
library(circlize)
library(lme4)
library(jtools)
library(ggrepel)
library(mediation)
library(gridExtra)

options(stringsAsFactors = FALSE)

#####
# Preparing files for mapping
covs <- read.table("../data/covs_and_peer_factors.txt")
held.out <- c("Neutrophils", "Lymphocytes", "Monocytes", paste0("PC", 1:7), "Diagnosis", "SRS1")
peer <- paste0("PEER_", 1:20)

covs <- covs %>%
  dplyr::select(Sample.ID, any_of(held.out), any_of(peer)) %>%
  as.data.frame()

rownames(covs) <- covs$Sample.ID
covs <- covs %>%
  dplyr::select(-Sample.ID)

eigengenes <- read.csv("../../nikhil/expression/gene_expression/eigengenes.multiple.csv", row.names=1)
eigengene.patients <- sapply(strsplit(rownames(eigengenes), "_"), function(x) { x[1] })

geno.fam <- fread("../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.fam") %>%
  dplyr::select(Family.ID=1, Individual.ID=2) %>%
  dplyr::mutate(GAinS.ID=gsub("^GA", "", Individual.ID)) %>%
  dplyr::filter(GAinS.ID %in% eigengene.patients) %>%
  unique()

samples.with.genotypes <- rownames(eigengenes)[
  sapply(
    strsplit(rownames(eigengenes), "_"),
    function(x) { x[1] %in% geno.fam$GAinS.ID }
  )
]

peer.included <- peer[1:20]
peer.included

cov.names <- c("Diagnosis", "SRS1", paste0("PC", 1:7), "Neutrophils", "Lymphocytes", "Monocytes")
cov.names

rand.effect.names <- c("GAinS.ID")
rand.effect.names

mapping.data <- merge(eigengenes, covs, by=0) %>%
  dplyr::filter(Row.names %in% samples.with.genotypes) %>%
  dplyr::select(Sample.ID=Row.names, everything()) %>%
  dplyr::mutate(GAinS.ID.NonPrefix=sapply(strsplit(Sample.ID, "_"), function(x) { x[1] })) %>%
  merge(., geno.fam, by.x="GAinS.ID.NonPrefix", by.y="GAinS.ID") %>%
  dplyr::select(GAinS.ID=Individual.ID, everything()) %>%
  dplyr::select(Sample.ID, any_of(rand.effect.names), any_of(colnames(eigengenes)), any_of(cov.names), any_of(peer.included))

colnames(mapping.data)

write.csv(mapping.data, "../modqtl_rerun/mapping_data.csv", row.names=F)

#####
# read in results
mqtl.snp.table <- read.csv("../../nikhil/expression/eigengene_sva/mqtl_snp_table.csv")
cond.snps <- unique(mqtl.snp.table$snps[which(mqtl.snp.table$source == "Conditional cis-eQTL SNP")])

geno.bim <- fread("../../Genotyping/All_genotyping_merged_filtered_b38_refiltered_rsID.bim") %>%
  as.data.frame()
colnames(geno.bim) <- c("chr", "snps", "cM", "pos", "minor", "major")
geno.bim$snps <- gsub(":", ".", geno.bim$snps)
geno <- geno.bim
colnames(geno) <- c("chr", "snp", "cM", "pos", "minor", "major")

me.assocs <- do.call(rbind, lapply(list.files("../modqtl_rerun/", 
                                              pattern="*_1.tsv"), function(file) {
                                                file.full = paste0("../modqtl_rerun/", file)
                                                me.assoc = fread(file.full, sep="\t", fill=TRUE) %>%
                                                  as.data.frame() %>%
                                                  dplyr::mutate(me = gsub("_1.tsv", "", file)) %>%
                                                  dplyr::select(snp=1, beta=2, se=3, t=4, p=5, me)
                                              }))

####### Plots for paper
# We decided that we need to add an additional filter for the SNPs i.e. same as
# the interactions (at least 4 minor allele homozygotes in total)

# Genotype Matrix
genotypes.file <- "../../nikhil/data/genotypes/eigengene_sva_genotypes.raw"
genotypes <- fread(genotypes.file, sep=" ", drop=2:6)

# Clean Genotype Matrix
genotypes$FID <- gsub("GA", "", genotypes$FID)
colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL
colnames(genotypes) <- gsub(":", ".", colnames(genotypes))

frqs <- apply(genotypes, 2, function(x){
  as.data.frame(table(x))
})

freqs <- lapply(1:length(frqs), function(x){
  frqs[[x]]$snp <- names(frqs[x])
  return(frqs[[x]])
})

frqs <- rbindlist(freqs)
frqsmin <- subset(frqs, x == 2)

mqtl.snp.table$majorhoms <- frqs$Freq[match(mqtl.snp.table$snps, frqs$snp)]
mqtl.snp.table$minorhoms <- frqsmin$Freq[match(mqtl.snp.table$snps, frqsmin$snp)]

me.assocs$majorhoms <- frqs$Freq[match(me.assocs$snp, frqs$snp)]
me.assocs$minorhoms <- frqsmin$Freq[match(me.assocs$snp, frqsmin$snp)]

me.assocs <- subset(me.assocs, minorhoms > 3)
mqtl.snp.table <- subset(mqtl.snp.table, minorhoms > 3)
mqtl.snp.table <- subset(mqtl.snp.table, source == "Conditional cis-eQTL SNP")

# Significance: 0.05. We performed tests for 106 individual eigengenes for 12335 SNPs. 
sig.thresh <- 0.05/(106*12335)

lambdas <- me.assocs %>%
  dplyr::group_by(me) %>%
  dplyr::summarize(lambda = median(qchisq(1 - p, 1)) / qchisq(0.5, 1))
range(lambdas$lambda)
# 0.9940493 1.1399296

num.modules <- length(unique(me.assocs$me)) # 106
num.snps <- length(unique(me.assocs$snp)) # 12335

mqtl.lead <- me.assocs %>%
  dplyr::filter(p < 0.05 / (num.modules * num.snps)) %>%
  merge(., geno, by="snp") %>%
  dplyr::group_by(me, chr) %>%
  dplyr::slice_min(p, n=1) %>%
  dplyr::arrange(p)
# 32

head(mqtl.lead) # 30 modules have a significant association.
# write.csv(mqtl.lead, "../modqtl_rerun/updated_lead_mqtl.csv", row.names=F, quote=F)

mqtl.all <- me.assocs %>%
  dplyr::filter(p < 0.05 / (num.modules * num.snps)) %>%
  merge(., geno, by="snp") %>%
  dplyr::group_by(me, chr) %>%
  dplyr::arrange(me, chr, pos)

head(mqtl.all) # 241 individual associations (fewer SNPs included in analysis)
# write.csv(mqtl.all, "../modqtl_rerun/updated_all_mqtl.csv", row.names=F, quote=F)

plot.data <- me.assocs %>%
  dplyr::filter(p < 0.05) %>%
  merge(., geno, by="snp") %>%
  dplyr::filter(chr %in% 1:22) %>%
  dplyr::mutate(chr=factor(chr, levels=1:22)) %>%
  dplyr::mutate(chr.type=ifelse(as.numeric(chr) %% 2 == 0, "Even", "Odd")) %>%
  dplyr::mutate(Log.10.p = -log10(p))

# Manhattan plot
plot.data %>%
  ggplot() +
  geom_point(aes(x=pos, y=Log.10.p, color=chr.type), size=I(0.5), alpha=0.5) +
  geom_hline(yintercept=-log10(0.05 / (num.modules * num.snps)), lty=3, linewidth=0.5) +
  scale_color_manual(values=c("Even" = "royalblue4", "Odd" = "royalblue1")) +
  facet_grid(~chr, scale="free_x", space="free_x") +
  guides(color="none") +
  xlab("Position") + ylab(expression('-log'[10]~'(p)')) +
  theme_minimal() +
  theme(
    strip.background=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.spacing=unit(0.1, "lines"),
    strip.text=element_blank(),
    axis.title.x=element_blank()
  )

# Module windows
windows <- me.assocs %>%
  dplyr::filter(p < 0.05 / (num.modules * num.snps)) %>%
  merge(., geno, by="snp") %>%
  dplyr::mutate(start = sapply(pos, function(x) { max(1, x - 10^6) }), end = pos + 10^6 - 1)

windows.merged <- lapply(split(windows, windows$me), function(x) {
  ranges <- makeGRangesFromDataFrame(x, seqnames.field="chr", start.field="start", end.field="end") %>%
    reduce() %>%
    as.data.frame()
  ranges$me <- x$me[1]
  return(ranges)
}) %>%
  do.call(rbind, .) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

geno.bim <- subset(geno.bim, snps %in% cond.snps)

snp.ranges <- makeGRangesFromDataFrame(
  as.data.frame(geno.bim), 
  seqnames.field="chr", start.field="pos", end.field="pos", 
  keep.extra.columns=TRUE
)

overlaps <- findOverlaps(windows.merged, snp.ranges)

module.qtl.snps <- cbind(
  as.data.frame(windows.merged)[overlaps@from,],
  as.data.frame(snp.ranges)[overlaps@to,]
) %>%
  as.data.frame() %>%
  dplyr::select(me, snps, chr=1, start=2, end=3) %>%
  dplyr::mutate(QTL.ID=paste0(chr, ":", start, "-", end)) %>%
  dplyr::select(ME=me, SNP=snps, QTL.ID)

head(module.qtl.snps)

# write.table(module.qtl.snps, "mqtl_full_summary_statistics_snps.txt", 
  # sep="\t", row.names=F, quote=F
# )

chr.lengths <- read.table("rna_seq_genomes/star2.7.8a_index_Homo_sapiens.GRCh38.99_100bp/chrNameLength.txt") %>%
  dplyr::mutate(Start=1) %>%
  dplyr::select(Chr=1, Start, End=2) %>%
  dplyr::filter(Chr %in% as.character(1:22)) %>%
  dplyr::mutate(Chr = factor(Chr, levels=as.character(1:22))) %>%
  dplyr::arrange(Chr)

modules <- read.csv("../../nikhil/expression/gene_expression/modules.csv") %>%
  dplyr::mutate(Eigengene=paste0("ME_", gsub(".*_", "", Module)))

plot.data <- as.data.frame(table(modules$Module)) %>%
  dplyr::select(Module=Var1, Frequency=Freq) %>%
  dplyr::filter(Module != "Unassigned")
plot.data$Module <- as.numeric(gsub("Module_", "", plot.data$Module))

# pdf("../modqtl_rerun/Module_size.pdf", useDingbats = F)
ggplot(plot.data) +
  geom_col(aes(x=Module, y=Frequency), fill="royalblue4", color="royalblue4", width=1, position=position_dodge(1)) +
  xlab("Module") + ylab("Number of Genes") +
  theme_bw()
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# dev.off()

gene.info.circos <- read.table("gene_info_864_20412_hla.txt") %>%
  dplyr::select(gene_id, seqnames, start, end, gene_name) %>%
  dplyr::filter(seqnames %in% as.character(1:22)) %>%
  dplyr::mutate(seqnames=factor(seqnames, levels=as.character(1:22)))


# Module stats
# Connectivity is a measure assigned to each gene that describes how well the 
# gene fits into a certain group of genes. Intra-module connectivity is a measure 
# of how well the gene belongs to the module, Conversely, inter-module connectivity 
# is a measure of how well the gene belongs to all other genes outside the module. 
# Genes that are more "central" in the module will tend to have a high intra-module 
# connectivity and a low inter-module connectivity.
connectivity <- read.csv("../../nikhil/expression/gene_expression/connectivity.csv")
ggplot(connectivity, aes(kWithin)) +
  geom_histogram() +
  facet_wrap(~ Module, scales = "free") +
  theme_bw()

# Plot correlation of module genes with module eigengenes, and eQTL boxplots (adjusted ME)
logcpm <- read.delim("../data/logcpm_864_20412_hla.txt")
eigengenes <- read.csv("../../nikhil/expression/gene_expression/eigengenes.csv")
# all(colnames(logcpm) == eigengenes$X)

design.matrix <- read.csv("../modqtl_rerun/mapping_data.csv")
genotypes <- fread("../../nikhil/data/genotypes/eigengene_sva_genotypes.raw", 
                   sep=" ", drop=2:6)

gene.info <- read.table("../data/gene_info_20412.txt") %>%
  dplyr::select(gene_id, seqnames, start, end, gene_name)

# Clean Genotype Matrix
genotypes$FID <- gsub("GA", "", genotypes$FID)
patient.sample.match <- match(design.matrix$GAinS.ID, genotypes$FID)
genotypes <- genotypes[patient.sample.match, ]
colnames(genotypes) <- gsub("X", "", colnames(genotypes))
colnames(genotypes) <- sapply(strsplit(colnames(genotypes), "_"), function(x) x[1])
genotypes[, 1] <- NULL
genotypes <- as.matrix(genotypes)
rownames(genotypes) <- design.matrix$Sample.ID

genotype.ids <- colnames(genotypes)

genotypes <- cbind(design.matrix, genotypes)
all.vars <- colnames(design.matrix)
eigens <- colnames(design.matrix)[grepl("^ME", colnames(design.matrix))]
covs <- setdiff(setdiff(all.vars, eigens), c("Sample.ID", "GAinS.ID"))

##### Circos plots
for (module in unique(mqtl.lead$me)[1:30]){
  pdf(paste0("../modqtl_rerun/", module, "_circos.pdf"), 
      useDingbats = F, onefile = T, width=10, height=10)

  # Circos plot  
  this.me <- module
  
  snp.df <- me.assocs %>%
    dplyr::filter(me == this.me) %>%
    merge(., geno, by="snp") %>%
    dplyr::filter(chr %in% 1:22) %>%
    dplyr::mutate(chr=factor(chr, levels=1:22)) %>%
    dplyr::mutate(chr.type=ifelse(as.numeric(chr) %% 2 == 0, "royalblue4", "royalblue1")) %>%
    dplyr::mutate(Log.10.p = -log10(p))
  
  qtl.df <- module.qtl.snps %>%
    # dplyr::select(ME, QTL.ID) %>%
    dplyr::filter(ME == this.me) %>%
    # unique() %>%
    dplyr::mutate(Chr=factor(gsub(":.*", "", QTL.ID), levels=as.character(1:22))) %>%
    dplyr::mutate(Range=gsub(".*:", "", QTL.ID)) %>%
    dplyr::mutate(Start=as.numeric(gsub("-.*", "", Range))) %>%
    dplyr::mutate(End=as.numeric(gsub(".*-", "", Range))) %>%
    dplyr::mutate(Midpoint=(Start + End) / 2)
  
  mod.genes <- merge(modules, gene.info.circos, by.x="Gene", by.y="gene_id") %>%
    merge(., qtl.df, by.x="Eigengene", by.y="ME") %>%
    dplyr::filter(Eigengene == this.me)
  
  circos.initialize(sectors=chr.lengths$Chr, xlim=chr.lengths[, c("Start", "End")])
  circos.par(points.overflow.warning=F)
  
  genes.sector = mod.genes %>%
    dplyr::select(seqnames, start, end, gene_name)  %>%
    unique()
  
  # If there is a cis-eGene, colour it differently
  # Take all associated SNPs for this module and look up their lead eGenes
  cis.eqtl <- mqtl.all %>%
    dplyr::select(me, snp, chr) %>%
    dplyr::filter(me == this.me) %>%
    merge(., mqtl.snp.table, by.x="snp", by.y="snps") %>%
    dplyr::filter(source == "Conditional cis-eQTL SNP") %>%
    merge(., gene.info.circos, by.x="egene", by.y="gene_id")
  
  label.color <- (genes.sector$gene_name %in% cis.eqtl$gene_name)
  label.color <- ifelse(label.color, "orange", "black")
  circos.genomicLabels(genes.sector, labels.column = 4, side = "outside",
                       col = label.color,
                       font=3)
  
  # circos.track(ylim=c(0, max(snp.df$Log.10.p)), track.height=0.2, bg.border=F, panel.fun=function(x, y) {
  circos.track(ylim=c(0, 1), track.height=0.05, bg.border=F, bg.col = "lightgrey", panel.fun=function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ycenter, CELL_META$sector.index, cex=0.8)
    
    # snp.sector = snp.df %>%
    #   dplyr::filter(chr == CELL_META$sector.index)
    # circos.points(snp.sector$pos, snp.sector$Log.10.p, sector.index=CELL_META$sector.index,
    #               col=rgb(col2rgb(snp.sector$chr.type)[1], col2rgb(snp.sector$chr.type)[2], col2rgb(snp.sector$chr.type)[3], max=255, alpha=125),
    #               cex=1, pch=20)
    
  })
  
  circos.track(ylim=c(0, 1), track.height=0.05, bg.border=F, panel.fun=function(x, y) {
    
    qtl.sector = qtl.df %>%
      dplyr::filter(Chr == CELL_META$sector.index)
    if (nrow(qtl.sector) > 0) {
      circos.segments(
        qtl.sector$Midpoint, 0, qtl.sector$Midpoint, 1, sector.index=CELL_META$sector.index,
        col="orange", lwd=2
      )
    }
  })
  
  red1 <- rgb(col2rgb("red")[1], col2rgb("red")[2], col2rgb("red")[3], max=255, alpha=125)
  rb1 <- rgb(col2rgb("royalblue1")[1], col2rgb("royalblue1")[2], col2rgb("royalblue1")[3], max=255, alpha=125)
  
  # If the beta is positive, colour blue, if negative colour red?
  # Use ME and SNP chr to match to lead mqtl table
  
  for (g in 1:nrow(mod.genes)) {
    this.beta <- mqtl.lead$beta[which(mqtl.lead$me == this.me &
                                        mqtl.lead$chr == mod.genes$Chr[g])]
    # if (mod.genes$Chr[g] != mod.genes$seqnames[g]) {
    circos.link(
      mod.genes$Chr[g], mod.genes$Midpoint[g],
      mod.genes$seqnames[g], mod.genes[g, c("start", "end")],
      col=ifelse(this.beta > 0, rb1, red1) 
      # lwd = abs(this.beta*100)
    )
    # }
  }
  
  # ME correlations
  module.genes <- modules$Gene[modules$Eigengene == module]
  module.gex <- t(logcpm[module.genes, ])
  me <- eigengenes[, module]
  df <- data.frame(me, module.gex)
  df.m <- reshape2::melt(df, id.vars="me")
  if(ncol(df) > 25){
    df.m <- df.m[which(df.m$variable %in% colnames(df)[2:26]), ]
  }
  df.m$gene.name <- gene.info$gene_name[match(df.m$variable, gene.info$gene_id)]
  gg <- ggplot(df.m, aes(me, value)) + 
    geom_point() +
    geom_smooth(formula = y ~ x, method='lm', se=F) +
    facet_wrap(~ gene.name, scales = "free") +
    ggtitle(module)
  print(gg)
  
  # eQTL boxplot
  snp <- as.character(mqtl.lead$snp[match(module, mqtl.lead$me)])
  
  variant.design <- genotypes[, c(paste0(module, "_1"), snp, covs, "GAinS.ID")]
  variant.design <- variant.design[!is.na(variant.design[, snp]), ]
  
  colnames(variant.design)[2] <- "SNP"
  
  f.alt <- as.formula(paste0(paste0(module, "_1"), "~`", "SNP" , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
  model.test <- lmer(f.alt, data=variant.design, REML=FALSE)
  
  summary(model.test)
  print(effect_plot(model.test, 
                    pred="SNP", x.label=snp,
                    interval = TRUE, 
                    plot.points = TRUE,
                    main.title = module,
                    jitter=c(0.1, 0),
                    partial.residuals=TRUE) + 
          theme_bw() +
            scale_x_continuous(breaks=c(0, 1, 2)))
  
  # genes within modules
  gg <- list()
  results <- list()
  par(mfrow=c(5, 5))
  
  for (gene in unique(colnames(df[, -1]))){
    variant.design <- genotypes[, c(paste0(module, "_1"), snp, covs, "GAinS.ID")]
    variant.design[, 1] <- df[match(rownames(genotypes), rownames(df)), gene]
    # variant.design <- variant.design[!is.na(variant.design[, snp]), ]
    
    colnames(variant.design)[1] <- "Module_gene"
    colnames(variant.design)[2] <- "SNP"
    
    f.alt <- as.formula(paste0("Module_gene", "~`", "SNP" , "`+", paste0(covs, collapse="+"), "+(1|GAinS.ID)"))
    set.seed("310323")
    model.test <- lmer(f.alt, data=variant.design, REML=FALSE)
    
    results[[gene]] <- data.frame(matrix(data=
        c(summary(model.test)$coefficients["SNP", ], "NA"),
      nrow=1, ncol=4), row.names = gene)
    
    gg[[gene]] <- effect_plot(model.test, 
                      pred="SNP", x.label=snp,
                      plot.points = TRUE,
                      main.title = gene.info$gene_name[match(gene, gene.info$gene_id)],
                      jitter=c(0.1, 0),
                      partial.residuals=TRUE) + 
            theme_bw() +
            scale_x_continuous(breaks=c(0, 1, 2))
  }

  par(mfrow=c(1, 1))
  do.call(grid.arrange, gg[which(names(gg) %in% unique(df.m$variable))])

  results <- do.call(rbind, results)
  results$kIn <- connectivity$kWithin[match(rownames(results), connectivity$X)]
  results$genename <- gene.info$gene_name[match(rownames(results), gene.info$gene_id)]
  results$genechr <- gene.info$seqnames[match(rownames(results), gene.info$gene_id)]
  results$local <- results$genechr == mqtl.lead$chr[match(module, mqtl.lead$me)]

  # plot kin against eQTL betas for each gene
  gg <- ggplot(results, aes(kIn, abs(as.numeric(X1)))) +
    geom_point(aes(colour=local)) +
    scale_shape_manual(values=c(1, 16))+
    theme_bw() +
    xlab("Within-module connectivity") +
    ylab("|eQTL Beta|") +
    geom_smooth(formula = 'y ~ x', method="lm") +
    geom_text_repel(aes(label=genename), size=2)
  print(gg)

  dev.off()
}

# Look up any GWAS variants
ebi <- fread("../modqtl_rerun/gwas_catalog_v1.0-associations_e105_r2022-04-07.tsv", header=TRUE, quote="") %>%
  as.data.frame()
ebi <- subset(ebi, SNPS %in% mqtl.all$snp)
gwas.ol <- merge(ebi, mqtl.all, by.x="SNPS", by.y="snp")
write.table(gwas.ol, "../modqtl_rerun/modqtl_gwas_overlap.txt", sep="\t")
