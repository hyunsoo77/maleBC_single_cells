#
# plot_sc_clusters.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Mar.
#
# requirement:
#   ArchR v1.0.1
#
# content:
#   get_df_cell_type()
#   sort_cluster_members()
#   get_color_for_cell_type()
#   update_nv_color_cell_type()
#   print_umap_cluster_types_with_seurat_obj()
#   print_umap_samples()
#   print_umap_clusters()
#   print_umap_cluster_types()
#   print_umap_meta.data()
#   print_umap_cell_cycle()
#
#   print_featureplot_markers()
#   featureplot_enrichment_analysis()
#
#   print_boxplot_total_cnv()
#
#   print_barplot_cell_type_percent()
#   print_barplot_cell_numbers_for_cell_types()
#
#   plot_browser_track()
#   heatmap_enrichment_analysis()
#
# usage:
# source("plot_sc_clusters.R")
#
# refernce:
# http://sape.inf.usi.ch/quick-reference/ggplot2/colour
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# http://applied-r.com/rcolorbrewer-palettes/
# https://www.google.com/search?q=rgb+color+picker
#


# for heatmap
source("./r/seurat/seurat_heatmap.R")

# for pattern_epi
source("./r/utilities_for_sc_analyses.R")


# font size
font_size_sc_clusters <- 18

# umap range
min_umap_range <- c(10, 10)

# raster
f_rasterize <- TRUE
options("ggrastr.default.dpi" = 300)



#nv_ep <- c("LEp_prog"="LEp_prog", "BEp"="BEp", "LEp"="LEp")
#nv_ep <- c("Luminal Epi. Prog."="LEp_prog", "Basal Epi."="BEp", "Luminal Epi."="LEp")

nv_ep <- c("Epithelial cells"="Epithelial cells", "Epi. Non-tumor"="Epi. Non-tumor", "Normal"="Normal", "LEp_prog"="LEp_prog", "LEp_secretory"="LEp_secretory", "LEp"="LEp", "LEp_hormone"="LEp_hormone", "BEp"="BEp", "BEp_MaSCs"="BEp_MaSCs", "BEp_myo"="BEp_myo", "Epi. Unassigned"="Epi. Unassigned", "Epi. Tumor"="Epi. Tumor", "Normal-like"="Normal-like", "NBL"="NBL", "Basal"="Basal", "CLow"="CLow", "Her2E"="Her2E", "LumA"="LumA" ,"LumB"="LumB")


#library(RColorBrewer)
epithelial.cols <- colorRampPalette(c("#a0e989", "#245719"))
epithelial.cols <- epithelial.cols(length(nv_ep)+1)


# nv_color
nv_color_epi_cancer <- c()
if (length(nv_ep) > 3) {
  nv_color_epi_cancer <- c("Epi. Unassigned"=epithelial.cols[8], "Epi. Tumor"=epithelial.cols[9], "Normal-like"=epithelial.cols[10], "NBL"=epithelial.cols[10], "Basal"=epithelial.cols[11], "CLow"=epithelial.cols[12], "Her2E"=epithelial.cols[13], "LumA"=epithelial.cols[14], "LumB"=epithelial.cols[15])
}

# define others
orange.cols <- colorRampPalette(c("#ff9d5c", "#ff6600"))

# nv_color
if (!exists("nv_color")) {
	nv_color <- c("Epithelial cells"=epithelial.cols[1], "Epi. Non-tumor"=epithelial.cols[1], "Normal"=epithelial.cols[1], "LEp_prog"=epithelial.cols[2], "LEp_secretory"=epithelial.cols[3], "LEp"=epithelial.cols[3], "LEp_hormone"=epithelial.cols[4], "BEp"=epithelial.cols[5], "BEp_MaSCs"=epithelial.cols[6], "BEp_myo"=epithelial.cols[7], nv_color_epi_cancer, "Keratinocytes"="khaki2", "Stromal cells"="#fabfd2", "Fibroblasts"="#fabfd2", "Smooth muscle"="#b47fe5", "Endothelial cells"="#93ceff", "Pericytes"="lightblue1", "Adipocytes"="yellow4", "Immune cells"="gray70", "T-cells"="gray60", "CD4+ T-cells"="gray60", "CD8+ T-cells"="gray50", "NK cells"="gray40", "Myeloid cells"="gold2", "Monocytes"="#c0c0c0", "DCs"="#d9d9d9", "Macrophages"="gold2", "Mast cells"="violetred1", "B-cells"="#b22222", "Others"="gray10", "Unknown"="black")

	cell_type_others <- c("Smooth muscle", "Pericytes")
	nv_color <- nv_color[!names(nv_color) %in% cell_type_others]
} # if



# get_df_cell_type
#
# input:
#    args$method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", "singler_blueprint_encode", ...}
#
# comment:
# nv_color is defined in sort_cluster_members.R
# this function is called for geom_bar()
#
# usage:
# df_cell_type <- get_df_cell_type(rna, args)
#
get_df_cell_type <- function(obj, args, col_cell_types=NULL, f_merge_stromal_cells=FALSE, f_merge_immune_cells=FALSE, n_log=0, log_cmd="display") {

  switch(class(obj),
 	"data.frame"={
		md <- obj %>% as.data.table
	},
	"Seurat"={
		meta.data <- obj@meta.data
		# extract meta data
		md <- rna@meta.data %>% as.data.table
		# the resulting md object has one "row" per cell
	},
	{}
  ) # switch

  if (is.null(col_cell_types)) {
    col_cell_types <- get_column_name_for_cell_types(obj, args)
  }

  if (f_merge_stromal_cells) {
	md[, (col_cell_types) := lapply(.SD, function(x) gsub(pattern_stromal_cells_for_replacement, "Stromal cells", x)), .SDcols=col_cell_types]
  }

  if (f_merge_immune_cells) {
	md[, (col_cell_types) := lapply(.SD, function(x) gsub(pattern_immune_cells_for_replacement, "Immune cells", x)), .SDcols=col_cell_types]
  }

  # count the number of cells per unique combinations of "Sample" and "seurat_clusters"
  #dt.rna1 <- md[, .N, by = c("Sample", col_cell_types)]

  # with additional casting after the counting
  df <- md[, .N, by = c("Sample", col_cell_types)] %>% 
    	#data.table::dcast(., Sample ~ SingleR.BED, value.var = "N") %>%
    	data.table::dcast(., sprintf("Sample ~ %s", col_cell_types), value.var = "N") %>%
	as.data.frame

  if (n_log > 0) {
    do.call(log_cmd, list(df) )
  }

  idx_col <- 2:ncol(df)
  cols <- colnames(df)[idx_col]
  f <- cols %in% names(nv_color)

  df_others <- df[, c("Sample", cols[!f]), drop=F]

  if (n_log > 0) {
    do.call(log_cmd, list(df_others) )
  }

  # df_cell_type
  df_cell_type <- df[, c("Sample", cols[f]), drop=F]

  if (ncol(df_others) > 1) {
    # update df$others
    cols <- 2:ncol(df_others)
    df_cell_type$Others <- apply(df_others[,cols,drop=F], 1, sum, na.rm=T)
  }

  if (n_log > 0) {
    do.call(log_cmd, list(df_cell_type) )
  }

  df_cell_type

} # get_df_cell_type






# sort_cluster_members
# input:
#   obj: 
#     (1) ArchRProject:
#     (2) data.frame, cluster.type information is located at col_cluster_types 
#     (3) factor: cell_type
#   args:
#     args$method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", "singler_blueprint_encode", ...}
#   col_cell_types: (default=NULL) column of cell types
#   col_cluster_types: {"cluster.type", ...} column of the most common cell type for each cluster 
#   str_umap_reduction: {"umap", ...}
#
# output:
#   list_sort$df: data.frame
#       UMAP_1, UMAP_2, Sample, cluster, cluster.type, cell_type
#	CNV.value, CNV.Pos
#	cluster_type: factor, gsub("^[0-9]+-", "", df$cluster.type)
#	cluster.new: factor
#   list_sort$df_cell_type
#       data.frame for count summary with col_cell_types (e.g. Epithelial cells, Fibroblass, ...)
#   list_sort$vec_cluster_num: vector, ordered by cell types (e.g. 7, 8, 0, ...)
#   list_sort$vec_cluster_name: vector, ordered by cell types (e.g. 7-Epithelial cells, 8-Epithelial cells, 0-LumA, ...)
#   list_sort$nv_color_cell_type: named vector (e.g. Epithelial cells: '#d9ffcf, LumA: #E48F6B, LumB: #FFB69A, ...)
#   list_sort$nv_color_cluster_type_ordered: named vector
#   list_sort$vec_cluster_name_epithelial
#   
#
# usage:
# list_sort <- sort_cluster_members(rna, col_cluster_types=col_cluster_types, str_umap_reduction=str_umap_reduction, f_merge_immune_cell=TRUE)
# rna_levels <- unique(rna.info_sort$vec_cluster_num)
# list_sort <- sort_cluster_members(factor(atac.sub$predictedGroup_ArchR))
sort_cluster_members <- function(obj, args, col_cell_types=NULL, col_cluster_types="cluster.type", str_umap_reduction="umap", f_merge_stromal_cells=FALSE, f_merge_immune_cells=FALSE, n_log_df_cell_type=1, n_log=0) {

  list_out <- NULL
  if (is.null(col_cell_types)) {
    # col_cell_types="SingleR.BED" when args$method_to_identify_cell_types="singler_blueprint_encode"
    col_cell_types <- get_column_name_for_cell_types(obj, args)
  }

  switch(class(obj),
  	"data.frame"={
		# The most common SingleR label in each cluster becomes the cluster label (see make_sc-rna-seq_merged_seurat_obj.R)
		df <- obj
		fac_cluster.type <- naturalsort::naturalfactor(df[, col_cluster_types])
	},

	"Seurat"={

		rna <- obj
  		list_out$df_cell_type <- get_df_cell_type(rna, args, col_cell_types, f_merge_stromal_cells, f_merge_immune_cells, n_log=n_log_df_cell_type)
		#df <- as.data.frame(rna@reductions$umap@cell.embeddings)
		#df <- Embeddings(object = rna[["str_umap_reduction"]]) %>% as.data.frame
		df <- as.data.frame(Embeddings(rna, reduction = str_umap_reduction))
		colnames(df) <- c("UMAP_1", "UMAP_2")
		#length(which(rownames(df)==rownames(rna@meta.data)))
		df$Sample <- rna$Sample
		df$cluster <- factor(rna@meta.data[,str_column_of_meta_data_cluster])
		df$cluster.type <- rna@meta.data[,col_cluster_types]
		df$cell_type <- rna@meta.data[,col_cell_types]
		df$phase <- rna@meta.data[,"Phase"]
		df$CNV.value <- rna$CNV.value
		df$CNV.Pos <- rna$CNV.Pos

		# manually annotate 23-cluster as smooth muscle
		#df$cluster.type <- str_replace_all(rna.df$cluster.type, "23-Stromal fibroblast", "23-Smooth muscle cells")

		df$cluster <- as.factor(df$cluster)
		df$cluster.type <- naturalsort::naturalfactor(df$cluster.type)
		#levels(df$cluster)
		#levels(df$cluster.type)

		fac_cluster.type <- df$cluster.type

	},

	"ArchRProject"={
		atac <- obj
		f_verbose <- getArchRVerbose()
		suppressMessages( addArchRVerbose(verbose = FALSE) )
		atac.df <- plotEmbedding(atac, colorBy = "cellColData", name = "Sample", embedding = "UMAP", logFile = createLogFile("plotEmbedding", logDir="log"))
		atac.df <- as.data.frame(atac.df$data) # columns: {"x", "y", "color"}

		df <- plotEmbedding(atac, colorBy = "cellColData", name = "predictedGroup_ArchR", embedding = "UMAP", logFile = createLogFile("plotEmbedding", logDir="log"))
		suppressMessages( addArchRVerbose(verbose = f_verbose) )

		df <- as.data.frame(df$data)
		colnames(df) <- c("UMAP_1", "UMAP_2", "color")
		df$Sample <- gsub(".*-", "", atac.df$color) # 1-4CC61L --> 4CC61L
		df$cluster.type <- sub(".*?-", "", df$color)
		df$cluster.type <- str_replace_all(df$cluster.type, "Epithelial$", "Epithelial cells")
		df$cluster <-  gsub("-.*", "", df$cluster.type)

		df$cluster <- as.factor(df$cluster)
		df$cluster.type <- naturalsort::naturalfactor(df$cluster.type)
		#levels(df$cluster.type) %in% levels(factor(atac$predictedGroup_ArchR)) --> all TRUE

		df$cell_type <- gsub("^[0-9]+-", "", df$cluster.type)
		# predictedGroup_ArchR has been done with relationship between scRNA-seq clusters and scATAC-seq clusters.
		list_out$df_cell_type <- get_df_cell_type(df, args, col_cell_types="cell_type", f_merge_stromal_cells, f_merge_immune_cells, n_log=n_log_df_cell_type)

		fac_cluster.type <- df$cluster.type

	},

	{
		# string vector
		fac_cluster.type <- naturalsort::naturalfactor(obj)
		df <- data.frame(cluster.type=fac_cluster.type)
  	}
  ) # switch
 
  # levels_cluster_type
  levels_cluster_type <- levels(fac_cluster.type)

  list_idx <- list()

  ### epithelial cells
  # HTAPP: Epithelial cell, HPCA: Epithelial_cells, BED: Epithelial cells
  epi <- grep("[Ee]pitheli", levels_cluster_type)
  epi.nontumor <- grep(sprintf("%s$", nv_ep[["Epi. Non-tumor"]]), levels_cluster_type)
  epi.normal <- grep(sprintf("%s$", nv_ep[["Normal"]]), levels_cluster_type)
  epi <- c(epi, epi.nontumor, epi.normal)

  epi.lep_prog <- grep(sprintf("%s$", nv_ep[["LEp_prog"]]), levels_cluster_type)
  epi.lep_secretory <- grep(sprintf("%s$", nv_ep[["LEp_secretory"]]), levels_cluster_type)
  epi.lep <- grep(sprintf("%s$", nv_ep[["LEp"]]), levels_cluster_type)
  epi.lep_hormone <- grep(sprintf("%s$", nv_ep[["LEp_hormone"]]), levels_cluster_type)
  epi.bep <- grep(sprintf("%s$", nv_ep[["BEp"]]), levels_cluster_type)
  epi.bep_mascs <- grep(sprintf("%s$", nv_ep[["BEp_MaSCs"]]), levels_cluster_type)
  epi.bep_myo <- grep(sprintf("%s$", nv_ep[["BEp_myo"]]), levels_cluster_type)

  epi.tumor <- c(); epi.unassigned <- c(); epi.normal_like <- c(); epi.basal <- c(); epi.clow <- c(); epi.her2e <- c(); epi.luma <- c(); epi.lumb <- c()
  if (length(nv_ep) > 3) {
    epi.unassigned <- grep(sprintf("%s$", nv_ep[["Epi. Unassigned"]]), levels_cluster_type)
    epi.tumor <- grep(sprintf("%s$", nv_ep[["Epi. Tumor"]]), levels_cluster_type)
    epi.normal_like <- grep(sprintf("%s$", nv_ep[["Normal-like"]]), levels_cluster_type)
    epi.nbl <- grep(sprintf("%s$", nv_ep[["NBL"]]), levels_cluster_type)
    epi.basal <- grep(sprintf("%s$", nv_ep[["Basal"]]), levels_cluster_type)
    epi.clow <- grep(sprintf("%s$", nv_ep[["CLow"]]), levels_cluster_type)
    epi.her2e <- grep(sprintf("%s$", nv_ep[["Her2E"]]), levels_cluster_type)
    epi.luma <- grep(sprintf("%s$", nv_ep[["LumA"]]), levels_cluster_type)
    epi.lumb <- grep(sprintf("%s$", nv_ep[["LumB"]]), levels_cluster_type)
  } # if

  epi.keratinocytes <- grep("Keratinocytes", levels_cluster_type)
  epi.new <- grep("-Ciliated", levels_cluster_type)

  #list_idx[["Epithelial cells"]] <- c(epi, epi.lep_prog, epi.bep, epi.lep, epi.normal_like, epi.basal, epi.her2, epi.luma, epi.lumb, epi.new)
  #list_idx[["Epithelial cells"]] <- c(epi, epi.lep_prog, epi.bep, epi.lep, epi.normal_like, epi.basal, epi.her2, epi.luma, epi.lumb, epi.keratinocytes, epi.new)
  #list_idx[["Epithelial cells"]] <- c(epi, epi.lep_prog, epi.lep_secretory, epi.lep, epi.lep_hormone, epi.bep, epi.bep_mascs, epi.bep_myo, epi.normal_like, epi.basal, epi.her2, epi.luma, epi.lumb, epi.new)
  list_idx[["Epithelial cells"]] <- c(epi, epi.lep_prog, epi.lep_secretory, epi.lep, epi.lep_hormone, epi.bep, epi.bep_mascs, epi.bep_myo, epi.unassigned, epi.tumor, epi.normal_like, epi.nbl, epi.basal, epi.clow, epi.her2e, epi.luma, epi.lumb, epi.new)

  list_normal_epi_idx <- list()
  list_normal_epi_idx[["Epi. Non-tumor"]] <- epi.nontumor
  list_normal_epi_idx[["Normal"]] <- epi.normal
  list_normal_epi_idx[["LEp_prog"]] <- epi.lep_prog
  list_normal_epi_idx[["LEp_secretory"]] <- epi.lep_secretory
  list_normal_epi_idx[["LEp"]] <- epi.lep
  list_normal_epi_idx[["LEp_hormone"]] <- epi.lep_hormone
  list_normal_epi_idx[["BEp"]] <- epi.bep
  list_normal_epi_idx[["BEp_MaSCs"]] <- epi.bep_mascs
  list_normal_epi_idx[["BEp_myo"]] <- epi.bep_myo
  list_normal_epi_idx[["Ep_others"]] <- c(epi.keratinocytes, epi.new)
  list_tumor_epi_idx <- list()
  if (length(nv_ep) > 3) {
    list_tumor_epi_idx[["Epi. Unassigned"]] <- epi.unassigned
    list_tumor_epi_idx[["Epi. Tumor"]] <- epi.tumor
    list_tumor_epi_idx[["Normal-like"]] <- epi.normal_like
    list_tumor_epi_idx[["NBL"]] <- epi.nbl
    list_tumor_epi_idx[["Basal"]] <- epi.basal
    list_tumor_epi_idx[["CLow"]] <- epi.clow
    list_tumor_epi_idx[["Her2E"]] <- epi.her2e
    list_tumor_epi_idx[["LumA"]] <- epi.luma
    list_tumor_epi_idx[["LumB"]] <- epi.lumb
  }

  # Keratinocytes
  list_idx[["Keratinocytes"]] <- epi.keratinocytes





  ### stromal cells
  if (f_merge_stromal_cells) {

    # stromal cells
    stromal <- grep(pattern_stromal_cells, levels_cluster_type)
    list_idx[["Stromal cells"]] <- stromal
    df$cluster.type <- gsub(pattern_stromal_cells_for_replacement, "Stromal cells", df$cluster.type)

  } else {

    list_idx[["Stromal cells"]] <- grep("[Ss]tromal", levels_cluster_type)

    # HTAPP: Fibroblast, HPCA: Fibroblasts, BED: Fibroblasts
    list_idx[["Fibroblasts"]] <- grep("[Ff]ibro", levels_cluster_type)

    # HTAPP: ??, HPCA: Smooth_muscle_cells, BED: Smooth muscle
    list_idx[["Smooth muscle"]] <- grep("[Ss]mooth", levels_cluster_type)

    # HTAPP: Endothelial cell, HPCA: Endothelial_cells, BED: Endothelial cells
    list_idx[["Endothelial cells"]] <- grep("[Ee]ndothel", levels_cluster_type)

    # Pericytes
    list_idx[["Pericytes"]] <- grep("[Pp]ericytes", levels_cluster_type)

    # Adipocytes
    list_idx[["Adipocytes"]] <- grep("[Aa]dipocytes", levels_cluster_type)

  } # if





  ### immune cells
  if (f_merge_immune_cells) {

    # immune cells
    immune <- grep(pattern_immune_cells, levels_cluster_type)
    list_idx[["Immune cells"]] <- immune
    df$cluster.type <- gsub(pattern_immune_cells_for_replacement, "Immune cells", df$cluster.type)

  } else {

    # HTAPP: T cell, HPCA: T_cells, BED: CD8+ T-cells
    #t <- grep("T cell", levels_cluster_type)
    t <- grep("T cell|T_cell|T-cell", levels_cluster_type)
    t.cd4 <- grep("CD4\\+ T", levels_cluster_type)
    t.cd8 <- grep("CD8\\+ T", levels_cluster_type)
    # HTAPP: NK cell, HPCA: NK_cell, BED: NK cells
    nk <- grep("NK cell|NK_cell", levels_cluster_type)
    # ??
    t.nk.new <- grep("Lym", levels_cluster_type)

    list_idx[["CD4+ T-cells"]] <- t.cd4
    list_idx[["CD8+ T-cells"]] <- t.cd8
    list_idx[["NK cells"]] <- nk
    list_idx[["T-cells"]] <- c(t)
    list_idx[["T/NK-cells"]] <- unique(c(t.cd4, t.cd8, nk, t, t.nk.new))

    # monocytes
    mono <- grep("[Mm]onocyte", levels_cluster_type)
    list_idx[["Monocytes"]] <- mono

    # dendritic cells
    dc <- grep("[Dd]endritic cell|DC", levels_cluster_type)
    list_idx[["DCs"]] <- dc

    # HTAPP: Macrophage, HPCA: Macrophage, BED: Macrophages
    list_idx[["Macrophages"]] <- grep("[Mm]acrophage", levels_cluster_type)

    # HTAPP: ???, HPCA: ???, BED: ???, make_sc-rna-seq_merged_seurat_obj.R renamed to Mast cells
    list_idx[["Mast cells"]] <- grep("[Mm]ast", levels_cluster_type)

    # HTAPP: B cell, HPCA: ??, BED: B-cells, make_sc-rna-seq_merged_seurat_obj.R renamed to B-cells
    #b <- grep("B cell", levels_cluster_type)
    list_idx[["B-cells"]] <- grep("B cell|B_cell|B-cell", levels_cluster_type)
  } # if


  # add here when any(is.na(list_sort$df$cluster_type))
  list_idx[["Others"]] <- grep("[Oo]thers", levels_cluster_type)
  list_idx[["Unknown"]] <- grep("[Uu]nknown", levels_cluster_type)


  cluster.types.idx <- unlist(list_idx)

  vec_cluster_num <- numeric(0)
  vec_cluster_name <- c()

  for(i in 1:length(cluster.types.idx)){

    # 0-Epithelial cell
    name <- levels_cluster_type[cluster.types.idx[i]]

    if (f_merge_stromal_cells) {
	name <- gsub(pattern_stromal_cells_for_replacement, "Stromal cells", name)
    }

    if (f_merge_immune_cells) {
	name <- gsub(pattern_immune_cells_for_replacement, "Immune cells", name)
    }

    #print(name)
    #print(gsub("-.*","",name))

    # 0-Epithelial cell --> 0
    new.name <- gsub("-.*", "", name)
    # 1_0-Epithelial cell
    if (grepl("_", new.name)) {
	new.num <- NA
    } else {
    	new.num <- as.numeric(new.name)
    }

    vec_cluster_num[i] <- new.num
    vec_cluster_name[i] <- name 
   
    str_type <- unique(gsub("[0-9]+-", "", name))  
    if (n_log > 0) {
      cat(sprintf("%s: %s\n", str_type, paste(new.num, collapse=", ")))  
    }
  } # for

  if (n_log > 0) {
    print(paste(vec_cluster_num, collapse=", "))
  }

  list_out$vec_cluster_num <- unique(vec_cluster_num)
  list_out$vec_cluster_name <- unique(vec_cluster_name)

  # cluster
  df$cluster <- gsub("-.*", "", df$cluster.type)

  # cluster_type
  df$cluster_type <- gsub("^[0-9]+-", "", df$cluster.type)
  levels.label <- unique(gsub("^[0-9]+-", "", list_out$vec_cluster_name))
  

  # substitution
  # LEp_prog --> Luminal Epi. Prog.
  idx <- match(df$cluster_type, names(nv_ep))
  f <- !is.na(idx); df$cluster_type[f] <- as.character(nv_ep[idx[f]])
  idx <- match(levels.label, names(nv_ep))
  f <- !is.na(idx); levels.label[f] <- as.character(nv_ep[idx[f]])

  df$cluster_type <- factor(df$cluster_type, levels=levels.label)
  list_out$nv_color_cell_type <- get_color_for_cell_type(list_idx, list_normal_epi_idx, list_tumor_epi_idx)

  # levels_tmp for ordering
  levels_tmp <- vec_cluster_num

  # sort
  switch(class(obj),
	"ArchRProject"={
  		# update predictedgroup_archr
		for ( cluster.type1 in levels(factor(atac$predictedGroup_ArchR))){
 			num <-  gsub("-.*", "", cluster.type1)
			idx <- match(num, levels_tmp)
			atac$predictedGroup_ArchR <- str_replace(atac$predictedGroup_ArchR, pattern = cluster.type1, replacement = paste0(idx, "_", atac$predictedGroup_ArchR))
  			#print("iter complete")
		} # for
		list_out$atac <- atac
		df$cluster.new <- factor(x = df$cluster, levels = unique(levels_tmp))
		list_out$df <- df
	},
	"data.frame"={
		# update df
    		df$cluster.new <- factor(x = df$cluster, levels = unique(levels_tmp))
    		list_out$df <- df
  	},
	{
		# update df
    		df$cluster.new <- factor(x = df$cluster, levels = unique(levels_tmp))
    		list_out$df <- df
	}
  ) # switch







  ##### colors


  ### epithelial
  # epi: light green - dark green
  epi.cols <- colorRampPalette(c("#a0e989", "#245719"))
  epi.cols <- epi.cols(length(list_idx[["Epithelial cells"]]))

  # keratinocytes: khaki2
  keratinocytes.cols <- c()
  len1 <- length(list_idx[["Keratinocytes"]])
  if (len1 > 0) {
    keratinocytes.cols <- rep("khaki2", len1)
  }


  ### stromal
  # stromal: light purple - dark purple
  stromal.cols <- colorRampPalette(c("#fabfd2", "#bf46ac"))
  stromal.cols <- stromal.cols(length(list_idx[["Stromal cells"]]))

  # fibroblasts: light purple - dark purple
  fibro.cols <- colorRampPalette(c("#fabfd2", "#bf46ac"))
  fibro.cols <- fibro.cols(length(list_idx[["Fibroblasts"]]))

  # smooth muscle
  smooth.cols <- c()
  len1 <- length(list_idx[["Smooth muscle"]])
  if (len1 > 0) {
    smooth.cols <- c("#b47fe5","#d8b7e8")
    smooth.cols <- smooth.cols[1:len1]
  }

  # endothelial: light blue - dark blue
  endo.cols <- colorRampPalette(c("#93ceff","#4a99ff"))
  endo.cols <- endo.cols(length(list_idx[["Endothelial cells"]]))

  # pericytes: light blue
  pericytes.cols <- c()
  len1 <- length(list_idx[["Pericytes"]])
  if (len1 > 0) {
    pericytes.cols <- rep("lightblue1", len1)
  }


  # adipocytes: light yellow
  adipocytes.cols <- c()
  len1 <- length(list_idx[["Adipocytes"]])
  if (len1 > 0) {
    adipocytes.cols <- rep("yellow4", len1)
  }


  ### immune cells
  immune.cols <- colorRampPalette(c("gray50", "gray30"))
  immune.cols <- immune.cols(length(list_idx[["Immune cells"]]))

  # T-cells: light gray - dark gray
  t.cols <- colorRampPalette(c("gray50", "gray30"))
  t.cols <- t.cols(length(list_idx[["T/NK-cells"]]))

  # monocytes: #c0c0c0
  mono.cols <- c()
  len1 <- length(list_idx[["Monocytes"]])
  if (len1 > 0) {
    mono.cols <- rep("#c0c0c0", len1)
  }

  # dendritic cells: light silver
  dc.cols <- c()
  len1 <- length(list_idx[["DCs"]])
  if (len1 > 0) {
    dc.cols <- rep("#d9d9d9", len1)
  }

  # macrophages: light orange - dark organge
  macro.cols <- colorRampPalette(c("gold2", "gold3"))
  macro.cols <- macro.cols(length(list_idx[["Macrophages"]]))
    
  # mast cells: violetred1
  mast.cols <- c()
  len1 <- length(list_idx[["Mast cells"]])
  if (len1 > 0) {
    #mast.cols <- "violetred1"
    #mast.cols <- mast.cols[1:length(mast)]  
    mast.cols <- rep("violetred1", len1)
  }

  # B-cells: dark red
  b.cols <- c()
  len1 <- length(list_idx[["B-cells"]])
  if (len1 > 0) {
    #b.cols <- c("#b22222","#cd5c5c")
    #b.cols <- b.cols[1:len1]
    b.cols <- rep("#b22222", len1)
  }



  ### others
  others.cols <- c()
  len1 <- length(list_idx[["Others"]])
  if (len1 > 0) {
    others.cols <- rep("black", len1)
  }

  # unknown
  unknown.cols <- c()
  len1 <- length(list_idx[["Unknown"]])
  if (len1 > 0) {
    unknown.cols <- rep("black", len1)
  }

  cols <- c(

	# epithelial
	epi.cols,
	keratinocytes.cols,

	# stromal
	stromal.cols,
	fibro.cols,
	smooth.cols,
	endo.cols,
	pericytes.cols,
	adipocytes.cols,

	# immune
	immune.cols,
	t.cols,
	mono.cols,
	dc.cols,
	macro.cols,
	mast.cols,
	b.cols,

	others.cols,
	unknown.cols 

    ) # cols


  list_out$nv_color_cluster_type_ordered <- cols
  names(list_out$nv_color_cluster_type_ordered) <- list_out$vec_cluster_name

  # select
  f <- grepl(pattern_epi, vec_cluster_name, ignore.case = F)
  list_out$vec_cluster_name_epithelial <- vec_cluster_name[f]


  # warning
  idx <- which(is.na(list_out$df$cluster_type))
  cluster.type.skipped <- unique(list_out$df[idx, "cluster.type"])
  if (length(cluster.type.skipped) > 0) {
  	cat(sprintf("\twarning: cluster.type.skipped: %s\n", paste(cluster.type.skipped, collapse=", ")))
  }

  list_out

} # sort_cluster_members








# get_color_for_cell_type
#
# comment:
# epithelial.cols was defined as a global variable, so it could be changed before calling this function when you used different epithelial.cols.
get_color_for_cell_type <- function(list_idx, list_normal_epi_idx, list_tumor_epi_idx) {

  idx_normal_epi <- which(sapply(list_normal_epi_idx, length) > 0)
  idx_tumor_epi <- which(sapply(list_tumor_epi_idx, length) > 0)
  if ((length(idx_normal_epi) > 1) && (length(idx_tumor_epi) == 0)) {
    # diversify epithelial cell types
    list_idx[["Epithelial cells"]] <- NULL
    list_idx <- c(list_normal_epi_idx, list_idx)
  }

  if (("Epi. Non-tumor" %in% names(list_normal_epi_idx)) && ("Epi. Tumor" %in% names(list_tumor_epi_idx))) {
    # diversify epithelial cell types
    list_idx[["Epithelial cells"]] <- NULL
    list_idx <- c(list_normal_epi_idx, list_idx)
  }

  if (length(idx_tumor_epi) > 0) {
    list_idx <- c(list_normal_epi_idx, list_tumor_epi_idx, list_idx)
  }

  # information of cell types: levels(rna.df$cluster_type)=c('Epithelial cells', 'Fibroblasts', 'Endothelial cells', 'CD8+ T-cells', 'NK cells', 'Macrophages', 'Mast cells', 'B-cells', 'Keratinocytes'}

  nv_color_epi_normal <- c(
		"Epithelial cells"=epithelial.cols[1],
		"Epi. Non-tumor"=epithelial.cols[1],
		"Normal"=epithelial.cols[1],
		"LEp_prog"=epithelial.cols[2],
		"LEp_secretory"=epithelial.cols[3],
		"LEp"=epithelial.cols[4],
		"LEp_hormone"=epithelial.cols[4],
		"BEp"=epithelial.cols[5],
		"BEp_MaSCs"=epithelial.cols[6],
		"BEp_myo"=epithelial.cols[7],
		"Ep_others"=epithelial.cols[1]
	)

  nv_color_epi_cancer <- c(
		"Epi. Unassigned"=epithelial.cols[8],
		"Epi. Tumor"=epithelial.cols[9],
		"Normal-like"=epithelial.cols[10],
		"NBL"=epithelial.cols[10],
		"Basal"=epithelial.cols[11],
		"CLow"=epithelial.cols[12],
		"Her2E"=epithelial.cols[13],
		"LumA"=epithelial.cols[14],
		"LumB"=epithelial.cols[15]
	)

  nv_color <- c(

		nv_color_epi_normal,
		nv_color_epi_cancer,
		"Keratinocytes"="khaki2",

		# stromal
		"Stromal cells"="#fabfd2",
		"Fibroblasts"="#fabfd2", 
		"Smooth muscle"="#b47fe5",
		"Endothelial cells"="#93ceff",
		"Pericytes"="lightblue1",
		"Adipocytes"="yellow4",

		# immune
		"Immune cells"="gray70",
		"T-cells"="gray60",
		"CD4+ T-cells"="gray60",
		"CD8+ T-cells"="gray50",
		"NK cells"="gray40",
		"Monocytes"="#c0c0c0",
		"DCs"="#d9d9d9",
		"Macrophages"="gold2",
		"Mast cells"="violetred1",
		"B-cells"="#b22222",

		# others
		"Others"="gray10",
		"Unknown"="black"

	)

  cell_types <- names(nv_color)

  # modify nv_color
  f <- grepl(" T[a-z]*-cell", cell_types)
  f <- sapply(list_idx[cell_types[f]], length) > 0
  if (any(f)) {
	nv_color <- nv_color[!names(nv_color) %in% c("T-cells")]
	cell_types <- names(nv_color)
  }

  f <- sapply(list_idx[cell_types], length) > 0
  #vec_color <- as.vector(nv_color[f])
  nv_color <- nv_color[f]

  nv_color

} # get_color_for_cell_type







# update_nv_color_cell_type
# usage:
# list_sort <- update_nv_color_cell_type(list_sort)
update_nv_color_cell_type <- function(list_sort) {

  idx <- match(names(list_sort$nv_color_cell_type), names(nv_color))
  f <- !is.na(idx)
  list_sort$nv_color_cell_type[f] <- nv_color[idx[f]]

  list_sort

} # nv_color_cell_type












### umap plots



# print_umap_cluster_types_with_seurat_obj
#
# input:
#   rna: seurat obj
#   str_umap_reduction: {["umap"], "umapharmony"}
#
# usage:
# gg <- print_umap_cluster_types_with_seurat_obj(rna, width=10, height=7)
# ggsave(sprintf("png/umap_%s_cell_types.png", str_condition), width = 10, height = 7, plot=gg)
print_umap_cluster_types_with_seurat_obj <- function(rna, str_umap_reduction = "umap", col_cluster_types="cluster.type", font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-1, xlim=NULL, ylim=NULL, width=10, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  #rna.df <- as.data.frame(rna@reductions$umap@cell.embeddings)
  #rna.df <- Embeddings(object = rna[[str_umap_reduction]]) %>% as.data.frame
  rna.df <- as.data.frame(Embeddings(rna, reduction = str_umap_reduction))

  rna.df$Sample <- rna$Sample
  rna.df$cluster.type <- as.factor(rna@meta.data[, col_cluster_types])
  rna.info_sort <- sort_cluster_members(rna.df)
  # update rna.df
  rna.df <- rna.info_sort$df

  # all(rownames(rna.df) == Cells(rna))
  rna$cluster_type <- rna.df$cluster_type

  # convert_cell_type_names
  rna$cluster_type <- convert_cell_type_names(rna$cluster_type)
  names(list_sort$nv_color_cell_type) <- convert_cell_type_names(names(list_sort$nv_color_cell_type))

  gg <- DimPlot(rna, 
	     reduction = str_umap_reduction,
	     group.by = "cluster_type", 
	     cols=rna.info_sort$nv_color_cell_type,
	     pt.size=0.1,
	     label=TRUE,
	     label.size = 6, label.color = "black")

  gg <- gg + ggtitle("") +
    guides(colour = guide_legend(title="cluster type",
			override.aes = list(size=6),
			ncol=1, bycol=TRUE)) +
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_text(size=font_size),
      legend.text=element_text(size=font_size.legend.text),
      #legend.key.size = unit(0.1, "cm")
     )

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_cluster_types%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width = width, repr.plot.height = height)
    print(gg)
  }

  gg

} # print_umap_cluster_types_with_seurat_obj













# print_umap_samples
# 
# input:
#   list_sort:
#   colors_samples:
#   pattern_cluster.type_removal:
#
# usage:
# gg <- print_umap_samples(list_sort, colors_samples=sampleColors, legend.position="bottom", ncol=2, width=7, height=9, str_condition=str_condition)
print_umap_samples <- function(list_sort, colors_samples, pattern_cluster.type_removal=NULL, title=NULL, legend.position="bottom", ncol=2, font_size=font_size_sc_clusters, font_size.legend.text=NULL, xlim=NULL, ylim=NULL, width=NULL, height=NULL, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  samples <- unique(as.character(list_sort$df$Sample))

  width <- 7
  if (is.null(height)) {
	height <- 7
	switch(legend.position,
		"bottom"={
			if (is.null(font_size.legend.text)) {
				if ((ncol == 1) || (length(samples) < 3)) {
					font_size.legend.text <- font_size_sc_clusters-4
				} else {
					font_size.legend.text <- font_size_sc_clusters-6
				}
			}
			height <- height + length(samples)*0.2 / ncol
		},
		"right"={
			if (is.null(font_size.legend.text)) {
				font_size.legend.text <- font_size_sc_clusters-1
			}
			width <- width + max(nchar(samples))*0.2
		},
		{}
	) # switch
	#cat(sprintf("\twidth=%g, height=%g\n", width, height))
  } # if


  df <- list_sort$df
  if (!is.null(pattern_cluster.type_removal)) {
	f <- grepl(pattern_cluster.type_removal, df$cluster.type)
	df <- df[!f,]
  } # if


  gg <- ggplot(df,
    aes(x = UMAP_1, y=UMAP_2, color = Sample))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    scale_color_manual(values = colors_samples)+
    guides(colour = guide_legend(override.aes = list(size=6),
			       ncol=ncol, byrow=FALSE)) +
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size),
      #plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_blank(),  
      legend.text=element_text(size=font_size.legend.text),
      legend.position=legend.position,
      legend.spacing.x=unit(0.5,'cm'),
      legend.spacing.y=unit(0.5,'cm')
  )

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  if (is.null(title)) {
        title <- sprintf("n = %s", format(nrow(df), nsmall=0, big.mark=",", scientific=F))
  }
  if (!is.null(title)) gg <- gg + ggtitle(title)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_samples%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_umap_samples












# print_umap_clusters
# usage:
# gg <- print_umap_clusters(list_sort, ncol=2, width=9, height=7, str_condition=str_condition)
print_umap_clusters <- function(list_sort, pattern_cluster.type_removal=NULL, title=NULL, ncol=2, font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-1, xlim=NULL, ylim=NULL, width=9, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  df <- subset(list_sort$df, !is.na(cluster.new))
  nv_colors <- list_sort$nv_color_cluster_type_ordered

  if (!is.null(pattern_cluster.type_removal)) {
        f <- grepl(pattern_cluster.type_removal, df$cluster.type)
        df <- df[!f,]
        f <- grepl(pattern_cluster.type_removal, names(nv_colors))
        nv_colors <- nv_colors[!f]
  } # if

  gg <- ggplot(df,
    aes(x = UMAP_1, y=UMAP_2, color = cluster.new))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    scale_color_manual(values = unname(nv_colors))+
    guides(colour = guide_legend(title="cluster",
			override.aes = list(size=6),
			ncol=ncol, bycol=TRUE)) +
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size),
      #plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_text(size=font_size),  
      legend.text=element_text(size=font_size.legend.text),
      #legend.key.size = unit(0.1, "cm")
  )

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  if (is.null(title)) {
        title <- sprintf("n = %s", format(nrow(df), nsmall=0, big.mark=",", scientific=F))
  }
  if (!is.null(title)) gg <- gg + ggtitle(title)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_clusters%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_umap_clusters








# print_umap_cluster_labels
#
# input:
#   list_sort:
#   pattern_cluster.type_removal:
#
# usage:
# gg <- print_umap_cluster_labels(list_sort, ncol=2, width=9, height=7, str_condition=str_condition)
print_umap_cluster_labels <- function(list_sort, pattern_cluster.type_removal=NULL, title=NULL, ncol=2, font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-1, xlim=NULL, ylim=NULL, width=9, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {


  df <- subset(list_sort$df, !is.na(cluster.new))
  nv_colors <- list_sort$nv_color_cluster_type_ordered

  if (!is.null(pattern_cluster.type_removal)) {
	f <- grepl(pattern_cluster.type_removal, df$cluster.type)
	df <- df[!f,]
	f <- grepl(pattern_cluster.type_removal, names(nv_colors))
	nv_colors <- nv_colors[!f]
  } # if

  gg <- ggplot(df,
    aes(x = UMAP_1, y=UMAP_2, color = cluster.new))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    scale_color_manual(values = unname(nv_colors))+
    guides(colour = guide_legend(title="cluster",
			override.aes = list(size=6),
			ncol=2, bycol=TRUE)) +
    theme_classic()+
    theme(,
      plot.title=element_text(size=font_size),
      #plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_text(size=font_size),  
      legend.text=element_text(size=font_size.legend.text),
      #legend.key.size = unit(0.1, "cm")
  )

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  #+NoLegend()

  gg <- LabelClusters(gg,
	    id="cluster.new",
	    color="black",
	    repel = T,
	    size=6)

  if (is.null(title)) {
        title <- sprintf("n = %s", format(nrow(df), nsmall=0, big.mark=",", scientific=F))
  }
  if (!is.null(title)) gg <- gg + ggtitle(title)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_cluster_labels%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg


} # print_umap_cluster_labels









# print_umap_cluster_types
#
# input:
#   list_sort:
#   pattern_cluster.type_removal:
#   pattern_cell.type_removal:
#
# usage:
# gg <- print_umap_cluster_types(list_sort, ncol=1, width=10, height=7, str_condition=str_condition)
print_umap_cluster_types <- function(list_sort, col_cluster.type="cluster.type", col_cluster_type="cluster_type", pattern_cluster.type_removal=NULL, pattern_cell.type_removal=NULL, title=NULL, ncol=1, font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-1, xlim=NULL, ylim=NULL, width=10, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  # unique(list_sort$df$cluster_type) # check if there is any NA
  # unique(list_sort$nv_color_cell_type)

  # convert_cell_type_names
  list_sort$df$cluster_type <- convert_cell_type_names(list_sort$df$cluster_type)
  names(list_sort$nv_color_cell_type) <- convert_cell_type_names(names(list_sort$nv_color_cell_type))
  list_sort$df[, col_cluster_type] <- convert_cell_type_names(list_sort$df[, col_cluster_type])

  df <- subset(list_sort$df, !is.na(cluster.new))
  if (!is.null(pattern_cluster.type_removal)) {
	f <- grepl(pattern_cluster.type_removal, df[, col_cluster.type])
	df <- df[!f,]
  } # if

  nv_colors <- list_sort$nv_color_cell_type
  if (!is.null(pattern_cell.type_removal)) {
	f <- grepl(pattern_cell.type_removal, names(nv_colors))
	nv_colors <- nv_colors[!f]
  } # if

  str_legend_title <- gsub("_", " ", col_cluster_type)
  gg <- ggplot(df,
    aes_string(x = "UMAP_1", y = "UMAP_2", color = col_cluster_type))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    scale_color_manual(values = nv_colors)+
    guides(colour = guide_legend(title = str_legend_title,
			override.aes = list(size=6),
			ncol=1, bycol=TRUE)) +
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size),
      #plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_text(size=font_size),  
      legend.text=element_text(size=font_size.legend.text),
      #legend.key.size = unit(0.1, "cm")
  )

  #+NoLegend()

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  gg <- LabelClusters(gg,
	    id="cluster_type",      
	    color="black",
	    repel = T,
	    size=6)

  if (is.null(title)) {
        title <- sprintf("n = %s", format(nrow(df), nsmall=0, big.mark=",", scientific=F))
  }
  if (!is.null(title)) gg <- gg + ggtitle(title)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_cluster_types%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_umap_cluster_types






# print_umap_meta.data
# usage:
# gg <- print_umap_meta.datav(rna, list_sort, "CNV.value", width=8, height=7, str_condition=str_condition)
print_umap_meta.data <- function(rna, list_sort, str_column, pattern_cluster.type_removal=NULL, f_scale=FALSE, color_name="", color_limits=NULL, na.value = "grey90", colors = c("#0000ff", "#ffffcc", "#ffffcc", "#ffffcc", "#ff0000"), color_values=c(0, 0.45, 0.5, 0.55, 1), title=NULL, font_size=font_size_sc_clusters, font_size.legend.title=font_size_sc_clusters-3, font_size.legend.text=font_size_sc_clusters-3, legend.key.width = unit(0.8, "cm"), legend.key.height=unit(1.5, "cm"), xlim=NULL, ylim=NULL, width=8, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {


  # update list_sort
  idx <- match(rownames(list_sort$df), rownames(rna@meta.data))
  f <- !is.na(idx); idx <- idx[f]
  list_sort$df[f, str_column] <- rna@meta.data[idx, str_column]

  df <- subset(list_sort$df, !is.na(cluster.new))

  switch(str_column,
	"treatment"={
		color_name <- ""
		colors <- colorRampPalette(c("#0000ff", "#ffffcc", "#ff0000"))(100)
		if (is.null(color_limits)) {
			color_limits <- c(0.25, 0.75)
		}
	},
	"CNV.value"={
		color_name <- "sCNV"
		colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
		color_values <- NULL
		if (is.null(color_limits)) {
			color_limits <- c(0, 0.4)
		}
	},
	{}
  ) # switch

  gg <- ggplot(df,
    aes_string(x = "UMAP_1", y = "UMAP_2", color = str_column))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    scale_colour_gradientn(name = color_name,
                 limits = color_limits,
                 oob = scales::squish,
                 na.value = na.value,
                 colours = colors,
                 values = color_values)+
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size),
      #plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_text(size=font_size.legend.title),
      legend.text=element_text(size=font_size.legend.text),
      legend.key.width = legend.key.width,
      legend.key.height = legend.key.height
    )

  #+NoLegend()

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  gg <- LabelClusters(gg,
	    id="cluster.new",
	    color="black",
	    repel = T,
	    size=6)

  if (is.null(title)) {
        title <- sprintf("n = %s", format(nrow(df), nsmall=0, big.mark=",", scientific=F))
  }
  if (!is.null(title)) gg <- gg + ggtitle(title)

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/umap_%s_%s%s.%s", figure_format, str_condition, tolower(str_column), fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_umap_meta.data







# print_umap_cell_cycle
# usage:
# gg <- print_umap_cell_cycle(list_sort, ncol=2, width=9, height=7, str_condition=str_condition)
print_umap_cell_cycle <- function(list_sort, sample=NULL, pattern_sample=NULL, group_name=NULL, ncol=2, font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-1, xlim=NULL, ylim=NULL, width=9, height=7, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  df <- list_sort$df

  group_name_ <- ""
  if (!is.null(sample)) {
	f <- (df$Sample == sample)
	df <- df[f,]
	group_name_ <- sample
  } else if (!is.null(pattern_sample)) {
	f <- grepl(pattern_sample, df$Sample)
	df <- df[f,]
	group_name_ <- gsub("[\\^\\*\\$]", "", pattern_sample)
  }

  if (is.null(group_name)) group_name <- group_name_
  if (nchar(group_name) == 0) group_name <- "others"

  tb <- table(df$phase)
  nv <- sprintf("%s (%s, %s%%)", names(tb), tb, round(tb / sum(tb),2))
  names(nv) <- names(tb)

  # update df$phase
  # cell cycle: G1 - S - G2 - M
  df$phase <- factor(df$phase, level=c("G1", "S", "G2M"))

  colors <- brewer.pal(length(unique(df$phase)), "Set1")
  gg <- ggplot(df,
    aes(x = UMAP_1, y=UMAP_2, color = phase))+
    { if (!f_rasterize) geom_point(size = .1) } +
    { if (f_rasterize) rasterise(geom_point(size = .1)) } +
    xlab(TeX(r'($UMAP_1$)'))+
    ylab(TeX(r'($UMAP_2$)'))+
    ggtitle(group_name)+
    scale_color_manual(labels = nv, values = colors)+
    guides(colour = guide_legend(title="phase",
			override.aes = list(size=6),
			ncol=1, bycol=TRUE)) +
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      axis.title.y=element_text(size=font_size),
      legend.title=element_blank(),  
      legend.text=element_text(size=font_size.legend.text),
      legend.position="bottom"
      #legend.key.size = unit(0.1, "cm")
    )
    #+NoLegend()

  if (is.null(xlim) && all(df$UMAP_1 > -min_umap_range[1] & df$UMAP_1 < min_umap_range[1])) xlim <- c(-min_umap_range[1], min_umap_range[1])
  if (is.null(ylim) && all(df$UMAP_2 > -min_umap_range[2] & df$UMAP_2 < min_umap_range[2])) ylim <- c(-min_umap_range[2], min_umap_range[2])

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  if (!is.null(str_condition)) {
	group_name <- gsub(" ", "_", group_name)
	filename <- sprintf("%s/umap_%s_%s_cell_cycle%s.%s", figure_format, str_condition, group_name, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg


} # print_umap_cell_cycle









### featureplot


# print_featureplot_markers
# usage:
# gg <- print_featureplot_markers(rna, ncol=4, width=12, height=9, str_condition=str_condition)
# options(repr.plot.width=12, repr.plot.height=9)
# gg
print_featureplot_markers <- function(rna, genes=NULL, min.cutoff=0, max.cutoff=10, ncol=4, xlim=NULL, ylim=NULL, width=12, height=9, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  genes_ <- c("ESR1", "PGR", "ERBB2", "FOXA1")
  # genes_basal <- c("KRT5", "KRT6A", "KRT14", "KRT17")
  genes_basal <- c("KRT5", "KRT14", "KRT17")
  genes_luminal <- c("KRT8", "KRT18", "KRT19")
  # https://www.nature.com/articles/s41588-021-00911-1
  genes_epi <- c("EPCAM")
  if ("PScore" %in% colnames(rna@meta.data)) {
	genes_proliferation <- c()  
	meta_proliferation <- c("PScore")  
  } else {
	genes_proliferation <- c("MKI67")  
	meta_proliferation <- c()
  }
  genes_tcells <- c("CD3D")  
  genes_myeloid_cells <- c("CD68")  
  genes_bcells <- c("MS4A1")  
  genes_plasmablasts <- c("JCHAIN")  
  genes_endothelial <- c("PECAM1")  
  genes_mesenchymal <- c("PDGFRB") # fibroblasts/perivascular-like cells

  if (is.null(genes)) {

	if (grepl("-bc|+bc|tnbc|-breast", args$cancer_type)) {
		# brca1-mut-bc, er+bc, her2+bc, male-bc, normal-breast, tr-bc (tamoxifen resistance breast cancer)
		genes <- c(genes_, genes_basal, genes_luminal, genes_proliferation)
		#genes <- c(genes_epi, genes_proliferation, genes_tcells, genes_myeloid_cells, genes_bcells, genes_plasmablasts, genes_endothelial, genes_mesenchymal)
		meta <- c(meta_proliferation)
	} else if (grepl("-oc", args$cancer_type)) {
		genes <- c(igenes_basal, genes_luminal, genes_proliferation)
		meta <- c(meta_proliferation)
	} else {
		stop(sprintf("no reference data defined for %s", args$cancer_type))
	}

	f <- genes %in% rownames(rna)
	genes <- genes[f]

	f <- meta %in% colnames(rna@meta.data)
	meta <- meta[f]

	genes <- c(genes, meta)

  } # if


  # min.cutoff, max.cutoff
  n_genes <- length(genes)
  min.cutoff <- rep(min.cutoff, n_genes) 
  max.cutoff <- rep(max.cutoff, n_genes)
  names(max.cutoff) <- genes

  if ("PScore" %in% names(max.cutoff)) {
	max_pscore <- max(rna$PScore, na.rm=T)
	if (max_pscore > 1.5) max_pscore <- 1.0
	else max_pscore <- 0.75
	max.cutoff["PScore"] <- max_pscore
  }

  # https://www.rdocumentation.org/packages/Seurat/versions/4.1.0/topics/FeaturePlot
  plots <- FeaturePlot(rna,
		features = genes,
		dims = c(1, 2),
		cells = NULL,
		#cols = if (blend) {     c("lightgrey", "#ff0000", "#00ff00") } else {    c("lightgrey", "blue") },
		pt.size = 2, # default=NULL
		order = FALSE,
		min.cutoff = min.cutoff,
		max.cutoff = max.cutoff,
		reduction = str_umap_reduction, 
		split.by = NULL,
		keep.scale = "feature",
 		shape.by = NULL,
		slot = "data",
		blend = FALSE, # default=FALSE, Scale and blend expression values to visualize coexpression of two features, Blending feature plots only works with two features
		blend.threshold = 0.5, # default=0.5, The color cutoff from weak signal to strong signal; ranges from 0 to 1.
		label = FALSE,
		label.size = 4,  # default=4
		repel = FALSE,
		ncol = ncol,
		coord.fixed = FALSE,
		by.col = TRUE,
		sort.cell = NULL,
		interactive = FALSE,
		combine = FALSE, # default=TRUE
		raster = f_rasterize # default=NULL, Rasterizing points since number of points exceeds 100,000.  To disable this behavior set `raster=FALSE`
	)


  gg <- patchwork::wrap_plots(plots, ncol = ncol) * xlab(TeX(r'($UMAP_1$)')) * ylab(TeX(r'($UMAP_2$)'))

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/featureplot_%s_markers%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_featureplot_markers






# featureplot_enrichment_analysis
# usage:
# gg <- featureplot_enrichment_analysis(rna, list_ea, str_umap_reduction, str_condition_tmp)
#
featureplot_enrichment_analysis <- function(rna, list_ea, str_umap_reduction="umap", str_condition=NULL, pattern_gene_removal="^MT-|^L[0-9]+$", ncol = 2, width_screen=9, width_file=9) {

  if (is.null(list_ea)) {
	return(NULL)
  }

  df_up <- list_ea[["df_up"]]
  df_dn <- list_ea[["df_dn"]]

  genes_up <- rownames(df_up)
  genes_dn <- rownames(df_dn)

  if (!is.null(pattern_gene_removal)) {
    genes_up <- genes_up[!grepl(pattern_gene_removal, genes_up)]
    genes_dn <- genes_dn[!grepl(pattern_gene_removal, genes_dn)]
  }

  genes <- c(head(genes_up,2), rev(head(genes_dn,2)))
  if (length(genes) == 0) {
	return(NULL)
  }

  
  plots <- FeaturePlot(rna, features = genes,
		 min.cutoff = 0, max.cutoff = 10,
		 reduction = str_umap_reduction,
		 ncol = ncol,
		 raster = f_rasterize, # Rasterizing points since number of points exceeds 100,000.  To disable this behavior set `raster=FALSE`
		 combine = FALSE)

  gg <- patchwork::wrap_plots(plots, ncol = ncol) * xlab(TeX(r'($UMAP_1$)')) * ylab(TeX(r'($UMAP_2$)'))

  if (!is.null(gg)) {
    height <- 2*2*round(length(genes)/2)
    options(repr.plot.width=width_screen, repr.plot.height=height)
    print(gg)

    if (!is.null(str_condition)) {
	# convert_cell_type_names
	str_condition <- convert_cell_type_names(str_condition)
	str_condition <- gsub(" ", "_", str_condition)
	ggsave(sprintf("%s/featureplot_%s.%s", figure_format, str_condition, figure_format), width = width_file, height = height, plot=gg)
    }
  }

  gg

} # featureplot_enrichment_analysis

















### boxplots





# print_boxplot_total_cnv
# input:
#   list_sort: list of output sort_cluster_members()
#   col_cluster: {cluster, cluster.type}
print_boxplot_total_cnv <- function(list_sort, col_cluster="cluster", clusters_include=NULL, str_column_of_meta_data_cluster_=str_column_of_meta_data_cluster, font_size=font_size_sc_clusters, width=NULL, height=NULL, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {


  df <- keep_most_common_cell_type_for_each_cluster(list_sort$df, args, str_column_of_meta_data_cluster_, "cluster.type", "cell_type", th_ratio=-1, min_ncell=-1, n_log=0)

  switch(col_cluster,
	"cluster"={
  		df <- df %>% dplyr::mutate(cluster = factor(cluster, levels = rev(unique(list_sort$vec_cluster_num))))
		if (is.null(width)) width <- 4
	},
	"cluster.type"={
  		df <- df %>% dplyr::mutate(cluster = factor(cluster.type, levels = rev(unique(list_sort$vec_cluster_name))))
		if (is.null(width)) width <- 6
	},
	{}
   ) # switch


  if (!is.null(clusters_include)) {
	# select clusters
	f <- df[,col_cluster] %in% clusters_include
	df <- df[f,]
  }

  if (is.null(height)) {
	height <- length(unique(df$cluster))*0.4
	#cat(sprintf("\twidth=%g, height=%g\n", width, height))
  }

  # convert_cluster.type
  df$cluster <- convert_cluster.type(df$cluster)
  names(list_sort$nv_color_cluster_type_ordered) <- convert_cluster.type(names(list_sort$nv_color_cluster_type_ordered))

  gg <- ggplot(df,
     aes(x=cluster, y=CNV.value, fill=cluster)) +
     geom_boxplot()+
     coord_flip()+
     scale_fill_manual(values = rev(list_sort$nv_color_cluster_type_ordered))+
     theme_classic()+
     theme(plot.title=element_text(size=font_size, face = "bold"),
       axis.text.x=element_text(size=font_size),
       axis.text.y=element_text(size=font_size, lineheight=0.9),
       axis.title.x=element_text(size=font_size),
       #axis.title.y=element_text(size=font_size),
       axis.title.y=element_blank(),
       plot.margin = margin(t=0,r=1,b=0,l=0, "cm") ) +
     NoLegend()

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/boxplot_%s_total_cnv%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_boxplot_total_cnv
 













### barplots


# print_barplot_cell_type_percent
# usage:
# gg <- print_barplot_cell_type_percent(df_cell_type, width=9, height=8, str_condition=str_condition)
print_barplot_cell_type_percent <- function(df_cell_type, cell.types_others=NULL, cell.types_merge=NULL, cell.types_include=NULL, cell.types_exclude=NULL, font_size=14, legend.position="right", legend.title="cell type", legend.font_size=14, legend.ncol=1, width=9, height=8, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  colnames(df_cell_type) <- convert_cell_type_names(colnames(df_cell_type))
  # cell.type_others
  if (!is.null(cell.types_others)) {
  	idx <- which(colnames(df_cell_type) %in% cell.types_others)
	if (length(idx) > 0) {
		if (!"Others" %in% colnames(df_cell_type)) {
			df_cell_type$Others <- 0
 		}
		if (nrow(df_cell_type) < 2) {
			df_cell_type$Others <- df_cell_type$Others + sum(df_cell_type[,idx], na.rm=T)
		} else {
			df_cell_type$Others <- df_cell_type$Others + rowSums(df_cell_type[,idx,drop=F], na.rm=T)
		}
		df_cell_type <- df_cell_type[,-idx]
	}
  }

  # cell.type_merge
  if (!is.null(cell.types_merge)) {
	for (cell.type1 in cell.types_merge) {
		if (!cell.type1 %in% names(list_pattern_cells)) next
		cols <- colnames(df_cell_type)
		idx <- grep(list_pattern_cells[[cell.type1]], cols)
		if (length(idx) < 1) next
		if (!cell.type1 %in% cols) {
			df_cell_type[[cell.type1]] <- 0
		}
		if (nrow(df_cell_type) < 2) {
			df_cell_type[[cell.type1]] <- df_cell_type[[cell.type1]] + sum(df_cell_type[,idx], na.rm=T)
		} else {
			df_cell_type[[cell.type1]] <- df_cell_type[[cell.type1]] + rowSums(df_cell_type[,idx,drop=F], na.rm=T)
		}
		df_cell_type <- df_cell_type[,-idx]
	}
  } # if

  # df_cell_type_percent
  cols <- 2:ncol(df_cell_type)
  vec_sum <- apply(df_cell_type[,cols,drop=F], 1, sum, na.rm=T)
  #vec_sum

  df_cell_type_percent <- df_cell_type
  df_cell_type_percent[,cols] <- 100*df_cell_type[,cols] / vec_sum

  df_gg <- df_cell_type_percent %>% tidyr::pivot_longer(!Sample, names_to = "cell_types", values_to = "percent") %>%
	    as.data.frame %>% 
	    dplyr::rename(sample = Sample)

  # convert_cell_type_names
  #df_gg$cell_types <- convert_cell_type_names(df_gg$cell_types)
  names(nv_color) <- convert_cell_type_names(names(nv_color))

  f <- names(nv_color) %in% unique(df_gg$cell_types)
  df_gg$cell_types <- factor(df_gg$cell_types, levels=names(nv_color)[f])

  # cell.types_include
  if (!is.null(cell.types_include)) {
	lvl <- levels(df_gg$cell_types)
	f <- df_gg$cell_types %in% cell.types_include
	df_gg <- subset(df_gg, f)
	df_gg$cell_types <- droplevels(df_gg$cell_types, lvl[!lvl %in% cell.types_include])
  }

  # cell.type_exclude
  if (!is.null(cell.types_exclude)) {
	f <- df_gg$cell_types %in% cell.types_exclude
	df_gg <- subset(df_gg, !f)
	df_gg$cell_types <- droplevels(df_gg$cell_types, cell.types_exclude)
  }

  #head(df_gg)

  f <- names(nv_color) %in% levels(df_gg$cell_types)
  gg <- ggplot(data = df_gg, mapping = aes(x = sample, y = percent, fill = cell_types )) + 
	geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(labels = names(nv_color)[f],
			  values= nv_color[f]) +
	xlab("") + ylab("percent of cells")
  gg <- gg + theme_bw() +
	   theme(plot.title=element_text(size=font_size),
		axis.text.x=element_text(size=font_size-4, angle=45, hjust=1),
		axis.text.y=element_text(size=font_size),
		axis.title.x=element_text(size=font_size),
		axis.title.y=element_text(size=font_size),
		legend.title=element_text(size=font_size),
		legend.text=element_text(size=legend.font_size),
		legend.position=legend.position)

  gg <- gg + guides(fill=guide_legend(title=legend.title, ncol=legend.ncol))

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/barplot_%s_cell_type_percent%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg


} # print_barplot_cell_type_percent





# print_barplot_cell_numbers_for_cell_types
# usage:
# gg <- print_barplot_cell_numbers_for_cell_types(df_cell_type, width=9, height=8, str_condition=str_condition)
print_barplot_cell_numbers_for_cell_types <- function(df_cell_type, cell.types_others=NULL, cell.types_merge=NULL, cell.types_include=NULL, cell.types_exclude=NULL, font_size=14, legend.position="right", legend.title="cell type", legend.font_size=14, legend.ncol=1, width=9, height=8, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  colnames(df_cell_type) <- convert_cell_type_names(colnames(df_cell_type))
  # cell.type_others
  if (!is.null(cell.types_others)) {
  	idx <- which(colnames(df_cell_type) %in% cell.types_others)
	if (length(idx) > 0) {
		if (!"Others" %in% colnames(df_cell_type)) {
			df_cell_type$Others <- 0
 		}
		if (nrow(df_cell_type) < 2) {
			df_cell_type$Others <- df_cell_type$Others + sum(df_cell_type[,idx], na.rm=T)
		} else {
			df_cell_type$Others <- df_cell_type$Others + rowSums(df_cell_type[,idx,drop=F], na.rm=T)
		}
		df_cell_type <- df_cell_type[,-idx]
	}
  }

  # cell.type_merge
  if (!is.null(cell.types_merge)) {
	for (cell.type1 in cell.types_merge) {
		if (!cell.type1 %in% names(list_pattern_cells)) next
		cols <- colnames(df_cell_type)
		idx <- grep(list_pattern_cells[[cell.type1]], cols)
		if (length(idx) < 1) next
		if (!cell.type1 %in% cols) {
			df_cell_type[[cell.type1]] <- 0
		}
		if (nrow(df_cell_type) < 2) {
			df_cell_type[[cell.type1]] <- df_cell_type[[cell.type1]] + sum(df_cell_type[,idx], na.rm=T)
		} else {
			df_cell_type[[cell.type1]] <- df_cell_type[[cell.type1]] + rowSums(df_cell_type[,idx,drop=F], na.rm=T)
		}
		df_cell_type <- df_cell_type[,-idx]
	}
  } # if

  df_gg <- df_cell_type %>% tidyr::pivot_longer(!Sample,
		names_to = "cell_types",
		values_to = "cell_numbers") %>%
	    as.data.frame %>% 
	    dplyr::rename(sample = Sample)

  # convert_cell_type_names
  #df_gg$cell_types <- convert_cell_type_names(df_gg$cell_types)
  names(nv_color) <- convert_cell_type_names(names(nv_color))

  f <- names(nv_color) %in% unique(df_gg$cell_types)
  df_gg$cell_types <- factor(df_gg$cell_types, levels=names(nv_color)[f])

  # cell.types_include
  if (!is.null(cell.types_include)) {
	lvl <- levels(df_gg$cell_types)
	f <- df_gg$cell_types %in% cell.types_include
	df_gg <- subset(df_gg, f)
	df_gg$cell_types <- droplevels(df_gg$cell_types, lvl[!lvl %in% cell.types_include])
  }

  # cell.type_exclude
  if (!is.null(cell.types_exclude)) {
	f <- df_gg$cell_types %in% cell.types_exclude
	df_gg <- subset(df_gg, !f)
	df_gg$cell_types <- droplevels(df_gg$cell_types, cell.types_exclude)
  }

  #head(df_gg)

  f <- names(nv_color) %in% levels(df_gg$cell_types)
  gg <- ggplot(data = df_gg, mapping = aes(x = sample, y = cell_numbers, fill = cell_types )) + 
	geom_bar(position = "stack", stat = "identity") +
	scale_fill_manual(labels = names(nv_color)[f],
			  values= nv_color[f]) +
	xlab("") + ylab("number of cells")
  gg <- gg + theme_bw() +
	   theme(plot.title=element_text(size=font_size),
		axis.text.x=element_text(size=font_size-4, angle=45, hjust=1),
		axis.text.y=element_text(size=font_size),
		axis.title.x=element_text(size=font_size),
		axis.title.y=element_text(size=font_size),
		legend.title=element_text(size=font_size),
		legend.text=element_text(size=legend.font_size),
		legend.position=legend.position)

  gg <- gg + guides(fill=guide_legend(title=legend.title, ncol=legend.ncol))

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/barplot_%s_cell_numbers_for_cell_types%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_barplot_cell_numbers_for_cell_types






# print_barplot_ratio_of_cnv_pos
#
# input:
#
# usage:
# gg <- print_barplot_ratio_of_cnv_pos(list_sort, width=6, height=12, str_condition=str_condition)
print_barplot_ratio_of_cnv_pos <- function(list_sort, col_cluster="cluster", str_column_of_meta_data_cluster_=str_column_of_meta_data_cluster, font_size=font_size_sc_clusters, width=NULL, height=NULL, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  df <- keep_most_common_cell_type_for_each_cluster(list_sort$df, args, str_column_of_meta_data_cluster_, "cluster.type", "cell_type", th_ratio=-1, min_ncell=-1, n_log=0)

  md <- df %>% as.data.table
  md$CNV.cancer <- ifelse(md$CNV.Pos == "11", 1, 0)

  dt.rna2 <- md[, .N, by = c("CNV.cancer", col_cluster)] %>% data.table::dcast(., sprintf("CNV.cancer ~ %s", col_cluster), value.var = "N")

  df <- as.data.frame(t(dt.rna2))
  df$cluster1 <- rownames(df)

  switch(col_cluster,
	"cluster"={
  		df$cluster1 <- factor(df$cluster, levels = rev(unique(list_sort$vec_cluster_num)))
		if (is.null(width)) width <- 4
	},
	"cluster.type"={
  		df$cluster1 <- factor(df$cluster, levels = rev(unique(list_sort$vec_cluster_name)))
		if (is.null(width)) width <- 6
	},
	{}
  ) # switch


  cols <- df[1,]
  df <- df[-1,]
  if (ncol(df) == 2) {
    if (cols[1] == "1") {
	# all cells are CNA+
	colnames(df) <- c("cnv.pos.1", "cluster")
	df$sum <- df[,1]
    } else if (cols[1] == "0") {
	# all cells are CNA-
	colnames(df) <- c("cnv.pos.0", "cluster")
	df$sum <- 0
    }
  } else if (ncol(df) < 4) {
    # cols <- c(0, 1, NA)
    colnames(df) <- c("cnv.pos.0", "cnv.pos.1", "cluster")
    df$sum <- rowSums(as.matrix(df[,1:2]), na.rm=T)
  } else {
    # cols <- c(NA, 0, 1, NA)
    colnames(df) <- c("cnv.pos.na", "cnv.pos.0", "cnv.pos.1", "cluster")
    df$sum <- rowSums(as.matrix(df[,1:3]), na.rm=T)
  } # if

  df$ratio_cnv.pos.1 <- df[,"cnv.pos.1"] / df$sum
  df$ratio_cnv.pos.1[is.na(df$ratio_cnv.pos.1)] <- 0

  #print(head(df))
  #print(class(df$cluster))
  #print(levels(df$cluster))
  #print(length(unique(df$cluster)))
  #print(length(list_sort$nv_color_cluster_type_ordered))

  if (is.null(height)) {
	height <- length(unique(df$cluster))*0.4
	#cat(sprintf("\twidth=%g, height=%g\n", width, height))
  }

  # convert_cluster.type
  df$cluster <- convert_cluster.type(df$cluster)
  names(list_sort$nv_color_cluster_type_ordered) <- convert_cluster.type(names(list_sort$nv_color_cluster_type_ordered))

  gg <- ggplot(df,
    aes(x=cluster, y=ratio_cnv.pos.1, fill=cluster)) +
    geom_bar(stat="identity", width=0.5) +
    coord_flip() + theme_classic() +
    scale_fill_manual(values = rev(list_sort$nv_color_cluster_type_ordered))+
    theme(plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      #axis.title.y=element_text(size=font_size)) +
      axis.title.y=element_blank() ) +
    ylab("ratio of CNV.Pos=11") +
    NoLegend()

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/barplot_%s_cnvpos%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_barplot_ratio_of_cnv_pos



# print_barplot_sample_proportion_per_subcluster
print_barplot_sample_proportion_per_subcluster <- function(list_sort, sample_colors, col_cluster="cluster", ncol=1, font_size=font_size_sc_clusters, font_size.legend.text=font_size_sc_clusters-3, width=NULL, height=NULL, str_condition=NULL, filename=NULL, fname_appendix="_rna", f_display=TRUE) {

  df <- list_sort$df %>% group_by_at(col_cluster) %>% 
	dplyr::count(Sample)
  colnames(df) <- c("Cluster", "Sample", "Cells")

  switch(col_cluster,
	"cluster"={
  		rna_levels <- unique(list_sort$vec_cluster_num)
		if (is.null(width)) width <- 6
	},
	"cluster.type"={
  		rna_levels <- unique(list_sort$vec_cluster_name)
		if (is.null(width)) width <- 8
	},
	{}
  ) # switch

  if (is.null(height)) {
	height <- length(unique(df$Cluster))*0.4
	height <- height + length(unique(df$Sample))*0.2 / ncol
	#cat(sprintf("\twidth=%g, height=%g\n", width, height))
  }

  df <- df %>% dplyr::mutate(cluster.type = factor(Cluster, levels = rna_levels))
 
  # convert_cluster.type
  df$cluster.type <- convert_cluster.type(df$cluster.type)

  # reorder cluster factor levels to group by cell type 
  gg <- ggplot(df, aes(fill=Sample, y=Cells, x=forcats::fct_rev(cluster.type))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip() + xlab("Clusters") + ylab("# of cells")+
    scale_fill_manual(values = sample_colors)+
    theme_classic()+
    theme(
      plot.title=element_text(size=font_size, face = "bold"),
      axis.text.x=element_text(size=font_size),
      axis.text.y=element_text(size=font_size, lineheight=0.9),
      axis.title.x=element_text(size=font_size),
      #axis.title.y=element_text(size=font_size),
      axis.title.y=element_blank(),
      legend.title=element_text(size=font_size),  
      legend.text=element_text(size=font_size.legend.text),
      #legend.key.size = unit(0.1, "cm")
      legend.position="bottom"
    ) +
    guides(fill=guide_legend(title="", ncol=ncol))

  if (!is.null(str_condition)) {
	filename <- sprintf("%s/barplot_%s_cluster_type_prop%s.%s", figure_format, str_condition, fname_appendix, figure_format)
  }

  if (!is.null(filename)) {
    ggsave(filename, width = width, height = height, plot = gg)
  }

  if (f_display) {
    options(repr.plot.width=width, repr.plot.height=height)
    print(gg)
  }

  gg

} # print_barplot_sample_proportion_per_subcluster








### browser track





# plot_browser_track
#
# requirement:
# write permission of ./archr_output/{ArrowFiles, Plots}
# chmod g+w -R /datastore/nextgenout5/share/labs/francolab/hyunsoo.kim/sc-atac-seq/male-bc/run-20211030/output_male-bc/archr_output/{ArrowFiles,Plots}
#
# input:
#   atac: 
#   atac.for_comparison: with cancer-specific peaks
#   df_markpeaks:
#     list.markerpeaks <- getMarkers(markerpeaks.for_comparison, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1.0"); df_markpeaks <- as.data.frame(list.markerpeaks[[1]])
#   marker_genes:
#   list_sort_atac: output of sort_cluster_members()
#     list_sort_atac <- sort_cluster_members(atac, args, col_cluster_types = col_cluster_types, str_umap_reduction = str_umap_reduction, f_merge_immune_cell = FALSE)
#   grangeslist_features: GRangesList
#     df_encode.all <- read.delim("./reference/genome_annotation/GRCh38-cCREs.bed.gz", header=F); colnames(df_encode.all)[1:3] <- c("seqnames","start","end"); gr_encode <<- makeGRangesFromDataFrame(df_encode.all); gr_hmec <- readRDS("reference/normal_cell_lines_breast/HMEC_H3K27ac_peaks.rds"); grangeslist_features <- GRangesList(Peaks = getPeakSet(atac), HMEC = gr_hmec, Encode = gr_encode)
#   
plot_browser_track <- function(atac, atac.for_comparison, df_markpeaks, marker_genes, list_sort_atac, grangeslist_features, str_upstream_downstream=NULL, usegroups=NULL, nv_cluster.type_conversion=NULL, f_debug=FALSE) {


  # consider only cluster types used in atac.for_comparison.
  peakset.for_comparison <- getPeakSet(atac.for_comparison)
  if (is.null(usegroups)) {
	usegroups <- unique(names(peakset.for_comparison))
  }

  geneAnnotation <- getGeneAnnotation(atac)
  nv_color_cluster_type <- list_sort_atac$nv_color_cluster_type
  #names(nv_color_cluster_type) <- gsub("Epithelial cells", "Epithelial", names(nv_color_cluster_type))
  f <- grepl("Epi", names(nv_color_cluster_type))
  nv_color_cluster_type[f] <- "forestgreen"
  if (!is.null(usegroups)) {
	nv_color_cluster_type <- nv_color_cluster_type[names(nv_color_cluster_type) %in% usegroups]
  }

  # sort usegroups by nv_color_cluster_type
  usegroups <- names(nv_color_cluster_type)
    
  if (!is.null(nv_cluster.type_conversion)) {
	# update nv_color_cluster_type
	cluster_types <- names(nv_color_cluster_type)
	idx <- match(cluster_types, names(nv_cluster.type_conversion))
	f <- !is.na(idx); idx <- idx[f]
	cluster_types[f] <- nv_cluster.type_conversion[idx]
	names(nv_color_cluster_type) <- cluster_types
  }

  # nv_color_cancer_normal
  nv_color_cancer_normal <- nv_color_cluster_type
  f.tumor <- grepl("Tumor", names(nv_color_cancer_normal))
  nv_color_cancer_normal[f.tumor] <- "gray80"
  nv_color_cancer_normal[!f.tumor] <- "gray90"

  # when you want to modify archr_browser.R
  source("r/archr/archr_functions_modified.R")

  loops_p2g <- getPeak2GeneLinks(atac.for_comparison,
		    corCutOff = 0.45,
		    FDRCutOff = 1e-12,
		    #PValCutOff = 1e-12,
		    #RawPValCutOff = 1e-12,
		    peakName = df_markpeaks$peakName,
		    varCutOffATAC = NULL,
		    varCutOffRNA = NULL,
		    resolution = 1,
		    returnLoops = TRUE)

  p2g <- loops_p2g[[1]]

  loops_coa <- getCoAccessibility(
	 ArchRProj = atac.for_comparison,
	 corCutOff = 0.70,  # ArchR v0.9.5 default=0.5
	 resolution = 1,    # ArchR v0.9.5 default=1000
	 returnLoops = TRUE)
  coa <- loops_coa[[1]]

  loops <- list()
  loops$Peak2GeneLinks <- p2g
  loops$CoAccessibility <- coa


  p2g.df <- NULL
  #p2g.df <- p2g.df.sub
  #p2g.df <- p2g.df.sub.for_comparison
  #p2g.df <- p2g.df.sub.plot.cancer_specific

  for (marker_gene in marker_genes) {
    
	if (is.null(p2g.df)) {
		idx <- which(p2g$geneName == marker_gene)
		p2g_gene <- p2g[idx]  
	} else { 
		# peak2gene links  
		idx <- which(p2g.df$geneName == marker_gene)
		#display(p2g.df[idx,])
		p2g_gene <- peakset.for_comparison[p2g.df[idx, "idxATAC"]]
		#display(p2g_gene)
	}
    
	if (length(p2g_gene) == 0) {
      
		cat(sprintf("%s: p2g_gene is empty\n", marker_gene))

		upstream <- 40000
		downstream <- 40000
      
	} else {
		 
		idx <- which(mcols(geneAnnotation$genes)$symbol == marker_gene)
		gr_gene <- geneAnnotation$genes[idx]
		if (f_debug) {
			display(gr_gene)
		}
		s <- start(gr_gene)
		e <- end(gr_gene)
		 
		isMinus <- BiocGenerics::which(strand(gr_gene) == "-")
 		#isOther <- BiocGenerics::which(strand(gr_gene) != "-")
		f_neg_strand <- (length(isMinus) > 0)   
      
 		if (f_neg_strand) {
			# s------p2g------e
			#	s--gene--e
			#		 | 
			upstream <- (e - min(start(p2g_gene)))
			#		 s------p2g------e
			#	s--gene--e
			#		  |
			downstream <- (max(end(p2g_gene)) - e)
		} else {
			# positive strand 
			# s------p2g------e
			#		 s--gene--e
			#		  |
			upstream <- (s - min(start(p2g_gene)))
			#		  s------p2g------e
			#		  s--gene--e
			#		  |
			downstream <- max(end(p2g_gene)) - s
		} # if



		if (upstream <= 0 || (f_neg_strand && upstream < (e-s+1))) {
			#  [add space] gene       p2g
			total_bp <- (e-s+1) + downstream
			if (f_neg_strand) {
 				#	     s--e
				# [add space] gene       p2g
				#		s <- resize(gr_gene, 1, "start")
				upstream <- round(total_bp/10)+(e-s+1)  
			} else {
				#	     s--e
				# [add space] gene       p2g
				#	     s <- resize(gr_gene, 1, "start")
				upstream <- round(total_bp/10)  
			}
		} # if
      
		if (downstream <= 0 || (!f_neg_strand && downstream < (e-s+1))) {
			#   p2g       gene [add space]
 			total_bp <- upstream + (e-s+1)
			total_bp <- (e-s+1) + downstream
			if (f_neg_strand) {
				#	   s--e
				# p2g       gene [add space]
				#	      s <- resize(gr_gene, 1, "start")
				downstream <- round(total_bp/10)  
			} else {
				#	   s--e
				# p2g       gene [add space]
				#	   s <- resize(gr_gene, 1, "start")
				downstream <- round(total_bp/10)+(e-s+1)  
			}
		} # if
    
		upstream <- max(upstream, 1e3)
		downstream <- max(downstream, 1e3)  

		upstream <- round(upstream*1.2)
		downstream <- round(downstream*1.2)
      
	} # if


	if (!is.null(str_upstream_downstream)) {
		vec_upstream_downstream <- strsplit(str_upstream_downstream, "[=, ]")[[1]]
		idx_loc <- match(marker_gene, vec_upstream_downstream)
		if (!is.na(idx_loc)) {
			upstream <- round(as.numeric(vec_upstream_downstream[idx_loc+1]))
			downstream <- round(as.numeric(vec_upstream_downstream[idx_loc+2]))
		}
	} # if
    
	if (f_debug) {
		display_html(sprintf("upstream: %d, downstream: %d", upstream, downstream))
	}
    
	idx <- which(p2g %over% p2g_gene)
	#idx <- which(p2g %over% p2g_gene | p2g %over% gr_gene)
	#display(length(idx) > 0)  
	p2g_ov <- p2g[idx]
	#display(p2g_ov)
  
	idx <- which(coa %over% p2g_gene)
	#display(length(idx) > 0)  
	coa_ov <- coa[idx]
	#display(coa_ov)
    
	# loops_gene
	loops_gene <- loops  
	idx <- which(p2g$geneName == marker_gene)
	loops_gene[["Peak2GeneLinks"]] <- p2g[idx]
    
	f_out <- tryCatch({
		p <- plotBrowserTrack(
			ArchRProj = atac,
			groupBy = "predictedGroup_ArchR",
			useGroups = usegroups,
 			geneSymbol = marker_gene,
 			upstream = upstream,
			downstream = downstream,
			loops = loops_gene,
			features = grangeslist_features,
			pal = nv_color_cluster_type,
			pal_strip = nv_color_cancer_normal,
			bulktracks_scale = 1.21,
			bulktracks_hjust = 0.082,
			bulktracks_vjust = 0, 
			logFile = createLogFile("plotBrowserTrack", logDir="log")
		) # plotBrowserTrack
		TRUE
	}, error = function(e) {
		cat(sprintf('%s',e))
		FALSE
	}, finally = {
	}) # try
      
	if (!f_out) next
		  
	filename_pdf <- sprintf("archr_peak2gene_%s.pdf", marker_gene)
	display_html(filename_pdf)  
    
	out <- capture.output( suppressMessages(
		plotPDF(plotList = p,
			name = filename_pdf,
			ArchRProj = atac,
			addDOC = T, width = 8, height = 5.0)
	) ) # out
    
  } # for

  marker_genes


} # plot_browser_track





### heatmap






# heatmap_enrichment_analysis
#
# input:
#   mtx: data matrix
#   list_markers: output of find_markers()
#   list_ea: execute_enrichment_analysis()
#   th_log2fc: abs(log2fc) > th_log2fc
#   th_padj: pajd < th_padj
#   min_pct: pct >= min_pct
#   pattern_gene_removal: (default="^MT-|^L[0-9]+$")
#   max_up: (default=20)
#   max_dn: (default=20)
#   n_sampling: (default=-1) the number columns for sampling. heatmap cannot be drawn when the number of columns is too large.
#   type_heatmap: {["plot_heatmap"], "get_heatmap"}
#   cluster_rows: 
#   list_top_annotation_legend_param:
#   heatmap_legend_side: (defualt="right")
#
# usage:
# list_markers <- find_markers()
# list_ea <- execute_enrichment_analysis(list_markers$markers)
# list_out <- heatmap_enrichment_analysis(list_markers, list_ea)
heatmap_enrichment_analysis <- function(list_markers, list_ea, str_condition=NULL, col_log2fc="avg_log2FC", col_pvalue="p_val", col_padj="p_val_adj", th_log2fc=0, th_padj=1.0, min_pct=0, pattern_gene_removal="^MT-|^L[0-9]+$", max_up=20, max_dn=20, n_sampling=-1, type_heatmap="plot_heatmap", cluster_rows=FALSE, list_top_annotation_legend_param=NULL, heatmap_legend_side="right", width=4, height=4) {

  if (is.null(list_ea)) {
	return(NULL)
  }

  df_up <- list_ea[["df_up"]]
  df_dn <- list_ea[["df_dn"]]

  idx_up <- which((df_up[,col_log2fc] > th_log2fc) & (df_up[,col_padj] < th_padj))
  df_up <- df_up[idx_up,,drop=F]
  df_up <- df_up[order(df_up[,col_padj]),]

  idx_dn <- which((df_dn[,col_log2fc] < -th_log2fc) & (df_dn[,col_padj] < th_padj))
  df_dn <- df_dn[idx_dn,,drop=F]
  df_dn <- df_dn[order(df_dn[,col_padj]),]

  f <- (df_up$pct.1 >= min_pct) & (df_up$pct.2 >= min_pct)
  df_up <- df_up[f,]
  f <- (df_dn$pct.1 >= min_pct) & (df_dn$pct.2 >= min_pct)
  df_dn <- df_dn[f,]

  genes_up <- rownames(df_up)
  genes_dn <- rownames(df_dn)

  if (!is.null(pattern_gene_removal)) {
    genes_up <- genes_up[!grepl(pattern_gene_removal, genes_up)]
    genes_dn <- genes_dn[!grepl(pattern_gene_removal, genes_dn)]
  }

  genes <- c(head(genes_up, max_up), rev(head(genes_dn, max_dn)))
  if (length(genes) == 0) {
	return(NULL)
  }

  mtx <- list_markers$mtx
  idx1 <- list_markers$idx1
  idx_ref <- list_markers$idx_ref

  if (n_sampling > 0) {
    # sampling
    n_idx1 <- length(idx1)
    n_sampling <- min(n_sampling, n_idx1)
    if (n_idx1 > n_sampling) {
	idx <- sample.int(n_idx1, n_sampling)
	idx1 <- idx1[idx]
	#cat(sprintf("\tn_idx1=%d\n", length(idx1)))
    } 

    n_idx_ref <- length(idx_ref)
    n_sampling <- round(n_sampling * (n_idx_ref / n_idx1))
    n_sampling <- min(n_sampling, n_idx_ref)
    if (n_idx_ref > n_sampling) {
	idx <- sample.int(n_idx_ref, n_sampling)
	idx_ref <- idx_ref[idx]
	#cat(sprintf("\tn_idx_ref=%d\n", length(idx_ref)))
    } 
  } # if

  idx <- c(idx_ref, idx1)

  group_name1 <- list_markers$group_name1
  group_name_ref <- list_markers$group_name_ref
  df_top_annot <- data.frame(type=c(rep(group_name_ref, length(idx_ref)), rep(group_name1, length(idx1))))
  rownames(df_top_annot) <- colnames(mtx)[idx]
  list_top_annotation_col <- list(type=map_row_annot_color(df_top_annot[,1]))
  if (is.null(list_top_annotation_legend_param)) {
    list_top_annotation_legend_param=list(type=list(title="",
	title_gp=gpar(fontsize=9), labels_gp=gpar(fontsize=9),
	grid_height=unit(3,"mm"), ncol=1))
  }

  if (!any(genes %in% rownames(mtx))) {
	cat(sprintf("no genes\n"))
	return(NULL)
  }

  if (length(idx) == 0) {
	cat(sprintf("no matched column\n"))
	return(NULL)
  }


  switch(type_heatmap,
	"plot_heatmap"={

		mtx <- as.matrix(mtx[genes, idx, drop=F])
		mtx_score <- t(scale(t(mtx), center=T, scale=T))
		list_out <- plot_heatmap(mtx_score,
				markers = genes,
				sort_var = c("type"),
				n = 8,
				anno_var = c("type"),
				anno_colors = NULL,
				anno_name = c("type"),
				df_col_info = df_top_annot,
				list_annotation_legend_param = list_top_annotation_legend_param,
				heatmap_legend_param = list(direction = "horizontal", legend_width = unit(4.25, "cm"), title = "Expression"),
				simple_anno_size = unit(0.75, "cm"),
				hm_limit = c(-2, 0, 2),
				#hm_colors = c("purple","black","yellow"),
				hm_colors = c("#ff00ff","black","#ffff00"),
				#hm_colors = c("#4575b4","white","#d73027"),
				#hm_colors = c("#0000ff", "#ffffcc","#ff0000"),
				row_font_size = 12)

		heatmap_legend_side <- "bottom"
		nchar_max <- 0
		for (j in 1:ncol(df_top_annot)) {
			nchar_max <- max(nchar_max, max(nchar(unique(df_top_annot[,j]))))
		}
		if (nchar_max > 30) {
			width <- width + (nchar_max-30)*0.2
		}
		if (length(genes) < 30) {
			height <- max(4, length(genes)*0.25 + 1.5)
		} else {
			height <- length(genes)*0.20 + 1.5
		}
	},
	"get_heatmap"={
		list_out <- get_heatmap(mtx[genes, idx, drop=F],
				type_heatmap="overview",
				show_column_names=FALSE,
				df_top_annot=df_top_annot,
				list_top_annotation_col=list_top_annotation_col,
				list_top_annotation_legend_param=list_top_annotation_legend_param,
				cluster_rows=cluster_rows,
				clustering_distance_rows="euclidean",
				clustering_method_rows="complete")
	},
	{}
  ) # switch

  # convert_cell_type_names
  str_condition <- convert_cell_type_names(str_condition)
  str_condition <- gsub(" ", "_", str_condition)

  filename_figure <- sprintf("heatmap_%s_%s",str_condition, list_out$type_col_short)
  log_txt(sprintf("filename_figure: %s\n", filename_figure))
  print_figure(list_out$list_ht, width=width, height=height,
    file=filename_figure,
    heatmap_legend_side=heatmap_legend_side,
    figure_format_=figure_format,
    resolution=300)

  list_out

} # heatmap_enrichment_analysis




