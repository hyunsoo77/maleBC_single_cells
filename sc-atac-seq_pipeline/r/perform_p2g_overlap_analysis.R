#
# perform_p2g_overlap_analysis.R
# check cluster markers from scRNA-seq
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
#
#
# usage:
# called by analyze_cancer_specific_p2g.R
#
# reference:
#
# 0. turorials
# https://www.archrproject.com/bookdown
# 1. papers
# https://www.nature.com/articles/s41588-021-00790-6
# 2. codes
# https://github.com/RegnerM2015/scENDO_scOVAR_2020
#






#wilcox_DEGs <- readRDS("wilcox_DEGs.rds")
if (grepl("-bc|+bc|tnbc|-breast", args$cancer_type)) {

	wilcox_DEGs <- readRDS(sprintf("%s/wilcox_degs/%s_%s_wilcox_degs.rds", dir_seurat_obj, cancer_type, col_cluster_types))

} else if (grepl("-oc", args$cancer_type)) {

	dir_cancer_type <- "/home/hkim77/francolab.w/sc-rna-seq/ovarian_cancer/run-20210916"
	wilcox_DEGs <- readRDS(sprintf("%s/rds/ovar_wilcox_degs.rds", dir_cancer_type))

} else {

	wilcox_DEGs <- readRDS(sprintf("%s/wilcox_degs/%s_%s_wilcox_degs.rds", dir_seurat_obj, cancer_type, col_cluster_types))
	#stop(sprintf("no reference data defined for %s", cancer_type))

} # if




# epithealial cells
idx <- grep("pithel", levels(factor(wilcox_DEGs$cluster)))
idx.new <- grep("-Ciliated", levels(factor(wilcox_DEGs$cluster)))
idx <- c(idx, idx.new)
labels <- levels(factor(wilcox_DEGs$cluster))[idx]

# cluster markers from scRNA-seq
f_select_degs <- wilcox_DEGs$avg_logFC >= 1.0 & 
		 wilcox_DEGs$p_val_adj <= 0.01 &
		 wilcox_DEGs$cluster %in% labels

wilcox_DEGs <- wilcox_DEGs[f_select_degs,]

p2g.df.sub.plot.cancer.degs <- p2g.df.sub.plot[p2g.df.sub.plot$geneName %in% wilcox_DEGs$gene,]



# save
saveRDS(p2g.df.sub, sprintf("%s/cancer_specific_p2g_table_degs.rds", dir_rds))




