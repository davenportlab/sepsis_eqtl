################################################################################
#
# 10. Get LD tagging SNPs for eSNPs
#
################################################################################# 

# cis.eqtl.conditional <- readRDS("conditional_eQTL_results_final.rds")
# write.table(gsub("\\.", ":", unique(cis.eqtl.conditional$SNP)), "conditional_snps.txt", quote=F, col.names=F, row.names=F)

plink \\
--bfile All_genotyping_merged_filtered_b38_refiltered_rsID \\
--show-tags conditional_snps.txt \\
--allow-extra-chr \\
--tag-kb 1000 \\
--tag-r2 1 \\
--list-all

cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.100r2.tags.tsv

plink \\
--bfile All_genotyping_merged_filtered_b38_refiltered_rsID \\
--show-tags conditional_snps.txt \\
--allow-extra-chr \\
--tag-kb 1000 \\
--tag-r2 0.9 \\
--list-all

cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.90r2.tags.tsv

plink \\
--bfile All_genotyping_merged_filtered_b38_refiltered_rsID \\
--show-tags conditional_snps.txt \\
--allow-extra-chr \\
--tag-kb 1000 \\
--tag-r2 0.8 \\
--list-all

cat plink.tags.list | sed -e "s/[[:space:]]\\+/\\t/g" | sed -e "s/^\\t//g" > ${snps_file.getSimpleName()}.80r2.tags.tsv