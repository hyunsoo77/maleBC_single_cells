#
# utilities_for_sc_analyses.R
# author: H. Kim
# date craeted: 2021, Oct.
# date last modified: 2021, Dec.
#
# dependencies:
# jupyter_message.R
#
# contents:
# convert_cell_type_names()
# convert_cluster.type()
# get_column_name_for_cell_types()
# get_vec_cell_types()
# set_vec_cell_types()
# select_cells_with_high_counts()
# select_cluster.types()
# keep_most_common_cell_type_for_each_cluster()
# calculate_sig_scores()
# inspect_gene_expr_distributions()
# get_diet_seurat_obj()
#
# get_med_mad()
# file.remove_when_file.exists()
# gzip_when_file.exists()
#
# comment:
#
# reference:
#
#




# patterns to select cell types

# pattern_epi <- "Epithelial|Epi\\.|LEp|BEp|Keratinocytes"
# pattern_epi_fibro <- "Ciliated|Epithelial|Epithelia|Epi\\.|LEp|BEp|Basal|Her2E|Lum[AB]|Fibroblast"


pattern_epi <- "(^|-| )Epithelial|(^|-| )Epi\\.|(^|-| )LEp|(^|-| )BEp|(^|-| )Keratinocytes|(^|-| )Basal|(^|-| )CLow|(^|-| )Her2[E]*|(^|-| )Lum[AB]|(^|-| )Normal-like|(^|-| )NBL|Myoepithelial|Luminal Progenitors|Mature Luminal|Cancer [Cc]ycling"


# normal
pattern_normal_epi <- "(^|Normal |Non-tumor )Epithelial|(^|-| )Epi\\. Non-tumor|(^|-| )LEp|(^|-| )BEp|(^|-| )Keratinocytes|Myoepithelial|Luminal Progenitors|Mature Luminal"
pattern_normal_epi_basal <- "(^|-| )BEp|(^|-| )LEp_prog|(^|-| )LEp_secretory|Myoepithelial|Luminal Progenitors"
pattern_normal_epi_luminal <- "(^|-| )LEp$|(^|-| )LEp_hormone|Mature Luminal"


# tumor
pattern_tumor_epi <- "(Cancer |Tumor )Epithelial|Epithelial [Cc]ycling|(^|-| )Epi. Tumor|(^|-| )Basal|(^|-| )CLow|(^|-| )Her2[E]*|(^|-| )Lum[AB]|(^|-| )Normal-like|(^|-| )NBL|Cancer [Cc]ycling"
pattern_tumor_basal <- "(^|-| )Basal"
pattern_tumor_luminal <- "(^|-| )Lum[AB]|(^|-| )Her2[E]*|(^|-| )Normal-like|(^|-| )NBL"


# stromal
pattern_stromal_cells <- "(^|-)Stromal|(^|-)Fibroblasts|(^|-)Smooth muscle|(^|-)Endothelial cells|(^|-)Pericytes|(^|-)Adipocytes"


# immune
pattern_immune_cells <- "(^|-)Immune|(^|-)T-cells|(^|-)CD4\\+ T-cells|(^|-)CD8\\+ T-cells|(^|-)Treg-cells|(^|-)NK cells|(^|-)NKT cells|(^|-)Monocytes|(^|-)DC|(^|-)Macrophages|(^|-)Mast cells|(^|-)B-cells"
pattern_immune_plus_endothelial_cells <- paste(pattern_immune_cells, "(^|-)Endothelial cells", sep="|")
pattern_immune_plus_stromal_cells <- paste(pattern_immune_cells, pattern_stromal_cells, sep="|")


list_pattern_cells <- list()
list_pattern_cells[["Immune cells"]] <- pattern_immune_cells
list_pattern_cells[["T-cells"]] <- c("(^|-)T-cells|(^|-)CD4\\+ T-cells|(^|-)CD8\\+ T-cells|(^|-)Treg-cells|(^|-)NK cells|(^|-)NKT cells")
list_pattern_cells[["B-cells"]] <- c("(^|-)B-cells")
list_pattern_cells[["Myeloid cells"]] <- c("(^|-)Myeloid|(^|-)Monocytes|(^|-)DCs|(^|-)Macrophages")


# usage:
# celltypes <- c('Myoepithelial', 'Luminal Progenitors', 'Mature Luminal', 'Plasmablasts', 'Cancer Cycling', 'Cancer Her2 SC', 'Cancer LumB SC', 'Cancer Basal SC', 'Cycling PVL', 'Cancer LumA SC')
# grepl(pattern_normal_epi, celltypes)
# grepl(pattern_tumor_epi, celltypes)


pattern_stromal_cells_for_replacement <- gsub("\\(\\^\\|\\-\\)", "", pattern_stromal_cells)
pattern_immune_cells_for_replacement <- gsub("\\(\\^\\|\\-\\)", "", pattern_immune_cells)


# for epithelial tumors
# these global variables can be changed in get_cancer_type_specific_info()
pattern_normal_cells <- pattern_immune_cells
pattern_tumor_cells <- pattern_epi







# nv_cell_type_conversion_table
if (!exists("nv_cell_type_conversion_table")) {
	nv_cell_type_conversion_table <- c("LEp_prog"="LEp_prog", "LEp_secretory"="LEp_secretory", "LEp"="LEp", "LEp_hormone"="LEp_hormone", "BEp"="BEp", "BEp_myo"="BEp_myo", "Normal"="Normal", "Normal-like"="Normal-like", "NBL"="NBL", "Basal"="Basal", "CLow"="CLow", "Her2E"="Her2E", "LumA"="LumA" ,"LumB"="LumB")
} # if




# update_args_cancer_type_standard
update_args_cancer_type_standard <- function(args) {

  if (!is.null(args$cancer_type_standard) && (nchar(args$cancer_type_standard) > 0)) return(args)

  if (grepl("brca|-bc|+bc|tnbc|-breast", args$cancer_type)) {
	# brca1-mut-bc, er+bc, her2+bc, male-bc, normal-breast, tr-bc (tamoxifen resistance breast cancer)
	args$cancer_type_standard <- "bc"

  } else if (grepl("-oc", args$cancer_type)) {

	args$cancer_type_standard <- "oc"

  } else if (grepl("-atll", args$cancer_type)) {

	args$cancer_type_standard <- "atll"

  } # if

  args

} # update_args_cancer_type_standard











# convert_cell_type_names
convert_cell_type_names <- function(vec_cell_types, nv_converion_table=nv_cell_type_conversion_table, nv_converion_pattern=NULL) {

	if (is.null(vec_cell_types)) return(NULL)

	if (!is.null(nv_converion_table)) {
		if (is.factor(vec_cell_types)) {
			vec_levels <- levels(vec_cell_types)
			f <- vec_levels %in% names(nv_converion_table)
			vec_levels[f] <- nv_converion_table[vec_levels[f]]
			levels(vec_cell_types) <- vec_levels

		} else {
			f <- vec_cell_types %in% names(nv_converion_table)
			vec_cell_types[f] <- nv_converion_table[vec_cell_types[f]]
		}
		return(vec_cell_types)
	} # if


	names_ <- names(nv_converion_pattern)
	idx <- which(names_ != as.vector(nv_converion_pattern))

        if (length(idx) == 0) return(vec_cell_types)

	if (is.factor(vec_cell_types)) {
		vec_levels <- levels(vec_cell_types)
		for (x in idx) {
	                vec_levels <- gsub(names_[x], nv_converion_pattern[[x]], vec_levels)
		}
		levels(vec_cell_types) <- vec_levels
	} else {
		for (x in idx) {
			vec_cell_types <- gsub(names_[x], nv_converion_pattern[[x]], vec_cell_types)
		}		
	}

	vec_cell_types

} # convert_cell_type_names





# convert_cluster.type
#
# usage:
# nv_cell_type_conversion_table <- c("Epi. Non-tumor"="Epi. CNA-"); str <- convert_cluster.type("Epi. Non-tumor")
convert_cluster.type <- function(obj, nv_converion_table=nv_cell_type_conversion_table) {


	if ("factor" %in% class(obj)) {
		vec <- levels(obj)
		mtx_ <- str_match(vec, "([0-9]+)-(.*)$")
		f <- mtx_[,3] %in% names(nv_cell_type_conversion_table)

		mtx_[f,3] <- nv_cell_type_conversion_table[mtx_[f,3]]
		mtx_[,3] <- paste0(mtx_[,2], "-", mtx_[,3])

		levels_out <- mtx_[,3]
		idx <- match(levels_out, names(nv_cell_type_conversion_table))
		f <- !is.na(idx); idx <- idx[f]
		levels_out[f] <- nv_cell_type_conversion_table[idx]

		levels(obj) <- levels_out
	} else {
		# 9-Epi. Non-tumor
		mtx_ <- str_match(obj, "([0-9]+)-(.*)$")
		f <- mtx_[,3] %in% names(nv_cell_type_conversion_table)

		# 9-Epi. CNA-
		mtx_[f,3] <- nv_cell_type_conversion_table[mtx_[f,3]]
		mtx_[,3] <- paste0(mtx_[,2], "-", mtx_[,3])
		obj_out <- mtx_[,3]

		idx <- match(obj_out, names(nv_cell_type_conversion_table))
		f <- !is.na(idx); idx <- idx[f]
		obj_out[f] <- nv_cell_type_conversion_table[idx]

		obj <- obj_out
	} # if

	
	obj

} # convert_cluster.type




# get_column_name_for_cell_types
#
# input:
#    obj: seurat obj or archrproject
#    args$method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", "singler_blueprint_encode"}
#
# output:
#    column name: {"SingleR.HTAPP_toolbox", "SingleR.HPCA", "SingleR.BED", ...}
#
get_column_name_for_cell_types <- function(obj, args) {

  switch(class(obj),
	"Seurat"={
		df_meta.data <- obj@meta.data
	},
	"ArchRProject"={
		df_meta.data <- obj@cellColData
	},
	{
		stop(sprintf("unknown class: %s", class(obj)))
	}
  ) # swtich

  if ("cell.type" %in% colnames(df_meta.data)) {

	col <- "cell.type"

  } else {
	switch(args$method_to_identify_cell_types,

		"singler_htapp_toolbox"={
               		 col <- "SingleR.HTAPP_toolbox"
		},

		"singler_hpca_cellidx"={
			col <- "SingleR.HPCA"
		},

		"singler_blueprint_encode"={
			col <- "SingleR.BED"
		},
		{
			col <- args$method_to_identify_cell_types
		}
	) # switch
  } # if

  col

} # get_column_name_for_cell_types









# get_vec_cell_types
#
# input:
#    args$method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", "singler_blueprint_encode"}
#
# usage:
# get_vec_cell_types(rna, args)
get_vec_cell_types <- function(rna, args) {

  col_cell_types <- get_column_name_for_cell_types(rna, args)

  nv_cell_types <- NULL

  if (col_cell_types %in% colnames(rna@meta.data)) {
        # mimic rna$SingleR.BED
        nv_cell_types <- rna@meta.data[,col_cell_types]
        names(nv_cell_types) <- rownames(rna@meta.data)
  } # if

  nv_cell_types

} # get_vec_cell_types









# set_vec_cell_types
#
# input:
#    args$method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", "singler_blueprint_encode"}
#
set_vec_cell_types <- function(rna, vec_cell_types, args) {

  col_cell_types <- get_column_name_for_cell_types(rna, args)

  # update rna@meta.data column of col_cell_types
  rna@meta.data[, col_cell_types] <- vec_cell_types

  rna

} # set_vec_cell_types












# select_cells_with_high_counts
#
# input:
#   col_cell_types: column name of cell types
#
# usage:
# col_cell_types <- get_column_name_for_cell_types(rna, args)
# idx_cells <- select_cells_with_high_counts(rna, col_cell_types, max_n_cells = args$max_n_cells, min_n_cells_per_cell_type = 20, assay = NULL, slot = "counts")
# idx_cells <- select_cells_with_high_counts(obj[[i]], "SingleR", max_n_cells = args$max_n_cells, min_n_cells_per_cell_type = 20)
select_cells_with_high_counts <- function(rna, col_cell_types, max_n_cells = 2000, min_n_cells_per_cell_type = 20, assay = NULL, slot = "data") {

  # mtx, df_meta
  mtx <- GetAssayData(object = rna, assay = assay, slot = slot)
  df_meta <- rna@meta.data[, col_cell_types, drop=F]
  colnames(df_meta) <- c("cell_type")

  n_cells <- ncol(mtx)
  if (n_cells <= max_n_cells) {
	idx_cells <- 1:n_cells
	return(idx_cells)
  } # if

  # select cell types
  tb <- df_meta %>% dplyr::group_by(cell_type) %>% dplyr::count()


  # check normal vs. tumor ratio
  f.normal <- grepl(pattern_normal_epi, tb$cell_type)
  f.tumor <- grepl(pattern_tumor_epi, tb$cell_type)
  n.normal <- sum(tb[f.normal,2])
  n.tumor <- sum(tb[f.tumor,2])
  if (n.normal > 0 && n.tumor > 0) {
	if (n.tumor/n.normal < 0.20) {
		# the number of tumor cells is less than 20% of # of normal cells.
		# e.g. 52BC3L n.tumor/n.normal=0.13271
		# exclude normal cells in order to focus on tumor cells.
		#tb <- tb[!f.normal,]
		# reduce the number of normal cells
		# so that n.normal ~ n.tumor*0.10
  		tb[f.normal,2] <- floor(tb[f.normal,2]*(n.tumor/n.normal)*0.10)
		n_cells <- sum(tb[,2])
	}
  } # if

  tb_max_n_cells <- tb
  tb_max_n_cells[,2] <- floor(tb[,2]*(max_n_cells/n_cells))
  tb_max_n_cells <- tb_max_n_cells[tb_max_n_cells[,2] >= min_n_cells_per_cell_type, ,drop=T]

  nnz_cells <- diff(mtx@p)

  idx_cells <- c()
  for (i in 1:nrow(tb_max_n_cells)) {
    if (is.na(tb_max_n_cells[[i,1]])) next
    idx_cell_type1 <- which(df_meta[,1] == tb_max_n_cells[[i,1]])
    idx_order <- order(nnz_cells[idx_cell_type1], decreasing=T)
    idx_cell_type1 <- idx_cell_type1[idx_order]

    # select cells with large number of expressed genes
    # tb[,2] can be larger than actual number of cells for a cell type when the total number of cells is smaller than max_n_cells=2000. (see slurm_er+bc_GSM4909300_ER-MH0032.out)
    n_cell_type1 <- min(length(idx_cell_type1), tb_max_n_cells[[i,2]])
    idx_cell_type1 <- idx_cell_type1[1:n_cell_type1]

    idx_cells <- c(idx_cells, idx_cell_type1)
  } # for

  idx_cells

} # select_cells_with_high_counts










# select_cluster.types
#
# input:
#   cluster.type: vector (e.g. '0-Epithelial cells', ... )
#   pattern_cell_type: {"pattern_epi", "pattern_normal_epi", "pattern_tumor_epi", ...}
#
# output:
#   list_out$cluster.type
#   list_out$cluster_num
#   list_out$cluster_name
#
# comment:
# this is not used yet.
#
# usage:
# list_cluster.type_tumor_epi <- select_cluster.types(list_sort$vec_cluster_name, pattern_tumor_epi)
select_cluster.types <- function(cluster.type, pattern_cell_type = NULL) {

  list_out <- list()
  list_out$cluster.type <- cluster.type

  if (!is.null(pattern_cell_type)) {
    idx <- grep(pattern_cell_type, cluster.type)
    list_out$cluster.type <- cluster.type[idx]
  }

  list_out$cluster_num <- gsub("-.*", "", list_out$cluster.type)
  list_out$cluster_name <- gsub("^[0-9]+-", "", list_out$cluster.type)

  list_out

} # select_cluster.types











# keep_most_common_cell_type_for_each_cluster
#
# input:
#   obj:
#   args:
#   col_seurat_cluster: (e.g. sprintf("RNA_snn_res.%g", args$seurat_resolution) or sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse="_"))
#   col_cluster_types: column name of cluster types (e.g. cluster.type)
#   col_cell_types: column name of cell types (e.g. SingleR.BED)
#   th_ratio: include clusters when (n_most_common_cell_type/ncells_cluster.type1 > th_ratio. use -1 to include all clusters.
#   min_ncells: include clusters when n_most_common_cell_type >= min_ncell. use -1 to include all clusters.
#   pattern_cluster.type_include: e.g. "B-cells|Mast cells"
#
# usage:
# rna <- keep_most_common_cell_type_for_each_cluster(rna, args, col_cluster.type, th_ratio=0.2, min_ncell=50, n_log=1)
# df <- keep_most_common_cell_type_for_each_cluster(list_sort$df, args, "cluster.type", "label_cell_type", th_ratio=0, min_ncell=0 , n_log=0)
keep_most_common_cell_type_for_each_cluster <- function(obj, args, col_seurat_cluster, col_cluster_types, col_cell_types=NULL, th_ratio=0.2, min_ncell=50, pattern_cluster.type_include=NULL, n_log=0) {

  if (n_log > 0) {
	cat(sprintf("\tkeep_most_common_cell_type_for_each_cluster(%s)\n", col_cluster_types))
  }

  idx_cell <- c()

  switch(class(obj),
	"Seurat"={
		rna <- obj
  		if (is.null(col_cell_types)) {
  			col_cell_types <- get_column_name_for_cell_types(rna, args)
  		}

		cluster.type <- rna@meta.data[, col_cluster_types]
		vec_cell_type <- rna@meta.data[, col_cell_types]
	},
	"data.frame"={
		df <- obj
		cluster.type <- df[, col_cluster_types]
		vec_cell_type <- df[, col_cell_types]
	},
	{}
  ) # switch

  cluster.types <- naturalsort::naturalsort(unique(cluster.type))
  for (cluster.type1 in cluster.types) {
	idx <- which(cluster.type == cluster.type1)
	ncells_cluster.type1 <- length(idx)

	str_cell_type <- gsub("^[0-9]+-", "", cluster.type1)
	idx_s <- which(vec_cell_type[idx] == str_cell_type)
	idx_most_common_cell_type <- idx[idx_s]
	n_most_common_cell_type <- length(idx_most_common_cell_type)

	str_log <- sprintf("\t\t%s %d/%d=%g", cluster.type1, n_most_common_cell_type, ncells_cluster.type1, n_most_common_cell_type/ncells_cluster.type1)

	if (n_most_common_cell_type == 0) {
		# strange, but let's pass this step
		str_log <- sprintf("%s -- abnormally included", str_log)
		idx_cell <- c(idx_cell, idx)
	} else if ( !is.null(pattern_cluster.type_include) && grepl(pattern_cluster.type_include, cluster.type1) ) {
		str_log <- sprintf("%s -- included", str_log)
		idx_cell <- c(idx_cell, idx)
	} else if ((n_most_common_cell_type/ncells_cluster.type1 > th_ratio) && 
	    	   (n_most_common_cell_type >= min_ncell))  {
		idx_cell <- c(idx_cell, idx_most_common_cell_type)
	} else {
		str_log <- sprintf("%s -- removed", str_log)
	}

	if (n_log > 0) {
		cat(sprintf("%s\n", str_log))
	}
  } # for

 
  switch(class(obj),
	"Seurat"={
		rna <- rna[,idx_cell]

		# get_most_common_cell_type_for_each_cluster
		rna@meta.data[,col_cluster_types] <- get_most_common_cell_type_for_each_cluster(rna, args, col_seurat_cluster, col_cell_types)

		obj <- rna
	},
	"data.frame"={

		obj <- df[idx_cell,]
	},
	{}
  ) # switch

  obj
 
} # keep_most_common_cell_type_for_each_cluster











# calculate_sig_scores
# 
# input:
#   rna: seurat obj
#   args:
#   df_sig:
#   type_cell_selection: {"all", "epi", "nonepi"}
#   assay: (default="RNA")
#   slot: (default="data")
#
# output:
#   finalmt: data.frame (sigs x cells)
#     1st row: scores for the first sig
#     2nd row: scores for the second sig
#     ...
#
# reference:
# https://www.nature.com/articles/s41588-021-00911-1
# 
# usage:
# df_sig <- data.frame(proliferation=c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C")) https://www.nature.com/articles/s41588-021-00911-1#ref-CR33
# proliferation_score <- calculate_sig_scores(rna, args, df_sig, "epi")
calculate_sig_scores <- function(rna, args, df_sig, type_cell_selection, assay = "RNA", slot = "data", col_cell_types="cell.type") {

	#if ("col_cell_types" %in% names(args)) {
	#	col_cell_types <- args$col_cell_types
	#}

	switch(type_cell_selection,
		"all"={
			idx_cells <- 1:ncol(rna)
		},
		"epi"={
			idx_cells <- grep(pattern_epi, rna@meta.data[, col_cell_types])
		},
		"nonepi"={
			f <- !grepl(pattern_epi, rna@meta.data[, col_cell_types])
			idx_cells <- which(f)
		},
		{}
	) # switch

	if (length(idx_cells) == 0) {
		return(NULL)
	}
	
	# signature score calculations
	switch(slot,
		"data"={
			# Averaged normalized expression of 11 genes34 (BIRC5, CCNB1, CDC20, NUF2, CEP55, NDC80, MKI67, PTTG1, RRM2, TYMS and UBE2C), independent of the SCSubtype gene lists, was used to compute the proliferation score.  https://www.nature.com/articles/s41588-021-00911-1
			tocalc <- GetAssayData(object = rna, assay = assay, slot = "data")
			switch(type_cell_selection,
				"all"={	 tocalc <- tocalc },
				{ tocalc <- tocalc[,idx_cells,drop=F] }
			) # switch
		},
		"scale.data"={
			#tocalc <- as.data.frame(tumors.combined@assays$RNA@scale.data)
			switch(type_cell_selection,
				"all"={ rna.tmp <- rna },
				{ rna.tmp <- rna[,idx_cells] }
			) # switch
			rna.tmp <- FindVariableFeatures(object = rna.tmp, assay = assay, selection.method = "vst", nfeatures = 3000)
			features <- union(VariableFeatures(rna.tmp, assay = assay), setdiff(unlist(df_sig), c(NA, "")))
			rna.tmp <- ScaleData(object = rna.tmp, features = features, assay = assay)
			tocalc <- GetAssayData(object = rna.tmp, assay = assay, slot = "scale.data")
		},
		{}
	) # switch

	outdat <- matrix(0,
			nrow = ncol(df_sig),
			ncol = ncol(tocalc),
			dimnames = list(colnames(df_sig),
			colnames(tocalc)))

	#for(i in 1:nrow(sigdat)) {
	for(i in 1:ncol(df_sig)) {
		#sigdat[i,!is.na(sigdat[i,])]->module
		#row <- as.character(unlist(module))
	  	genes <- as.character(df_sig[,i])
		genes <- unique(genes[genes != ""])
		genes <- which(rownames(tocalc) %in% genes)
		if (length(genes) < 1) next
		temp <- apply(tocalc[genes,,drop=F], 2,
			function(x){
				mean(as.numeric(x), na.rm=TRUE)
			} )
		outdat[i,] <- as.numeric(temp)
	} # for

	final <- outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
	final <- as.data.frame(final)
	is.num <- sapply(final, is.numeric)
	final[is.num] <- lapply(final[is.num], round, 4)
	finalm <- as.matrix(final)
	finalmt <- as.data.frame(t(finalm))

	#center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
	#		get_average <- function(v) sum(v * row.w)/sum(row.w)
	#		average <- apply(x, 2, get_average)
	#		sweep(x, 2, average)
	#}
	#finalm.sweep.t <- center_sweep(finalmt)
	#finalnames <- colnames(finalm.sweep.t)[max.col(finalm.sweep.t, ties.method="first")]

	finalmt

} # calculate_sig_scores




# inspect_gene_expr_distributions
#
# input:
#    nv_samples: e.g. c("Cancer LumA SC"="CID4290A", "Cancer LumB SC"="CID4535", "Cancer Her2 SC"="CID3921", "Cancer Basal SC"="CID4515")
#
# usage:
# df <- inspect_gene_expr_distributions(rna, "KRT5", n_log=1)
# df <- inspect_gene_expr_distributions(rna, "ERBB2", col_cell_types="celltype_minor", n_log=1)
# df <- inspect_gene_expr_distributions(rna, "ERBB2", col_cell_types="celltype_minor", nv_samples=c("LumA"="CID4290A", "LumB"="CID4535", "Her2"="CID3921", "Basal"="CID4515"), n_log=1)
inspect_gene_expr_distributions <- function(rna, gene, assay = "RNA", slot = "data", col_cell_types="cell.type", subtypes=c("LumA", "LumB", "Her2", "Basal", "CLow", "NBL", "BEp", "LEp_prog", "LEp"), col_sample="Sample", nv_samples=NULL, probs=c(seq(0.1,0.9,by=0.1), 1.0), n_log=0) {


  mtx <- GetAssayData(object = rna, assay = "RNA", slot = "data")

  df <- data.frame()
  subtypes_ <- c()
  for (subtype in subtypes) {
	f <- grepl(subtype, rna@meta.data[,col_cell_types])
	if (length(which(f)) == 0) next
	if (!is.null(col_sample) && !is.null(nv_samples)) {
		f <- f & (rna@meta.data[,col_sample] == nv_samples[subtype])
	}
	if (gene %in% rownames(mtx)) {
		df <- rbind(df, quantile(mtx[gene, which(f)], probs=probs, na.rm = TRUE))
	} else {
		df <- rbind(df, quantile(rna@meta.data[which(f), gene], probs=probs, na.rm = TRUE))
	}
	subtypes_ <- c(subtypes_, subtype)
  } # for

  if (nrow(df) == 0) return(df)

  colnames(df) <- sprintf("p%.1f", probs)
  rownames(df) <- subtypes_

  if (n_log > 0) {
	log_txt(sprintf("\t\t%s\n", gene))
	log_obj(round(df,2), tab=3)
  }

  df

} # inspect_gene_expr_distributions








# get_diet_seurat_obj
# 
# input:
#    rna: seurat obj
#    args:
#
# output:
#    rna: seurat obj without data/scale.data slots 
#
# comment:
#    reproduce data/scale.data slots
#    rna <- SetAssayData( object = rna, slot = 'data', new.data = log1p(x = GetAssayData(object = rna, slot = 'counts')) )
#    rna <- ScaleData(rna, features = rownames(rna) )
#
# usage:
# rna_for_saverds <- get_diet_seurat_obj(rna, args)
# path_rds <- sprintf("%s/%s_rna_passedpc1checks.rds", dir_seurat_obj, cancer_type)
# cat(sprintf("\tsaveRDS(rna, '%s')\n", path_rds))
# saveRDS(rna_for_saverds, path_rds)
get_diet_seurat_obj <- function(rna, args, f_force_diet=FALSE, dimreducs=c('pca', 'umap'), diet_level=NULL) {

  if (is.null(diet_level)) {
	diet_level <- args$diet_seurat_level
  }

  if (args$f_diet_seurat || f_force_diet) {

    str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)
    str_column_of_meta_data_harmony <- sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse="_"))
    str_column_of_meta_data_cluster_multik <- "RNA_multik"
    str_column_of_meta_data_harmony_multik <- sprintf("RNA_harmony_th.%s_multik", paste(args$harmony_theta, collapse="_"))
  
    # select columns for reducing memory/disk usage
    cols_epi_type <- c("epi_krt_epcam", "epi_normal_tumor", "normal_epi_type", "tumor_epi_type", "epi_type")

    # cell type
    #cols_cell_type <- c("celltype_garnett.major", "celltype_garnett.major.ext", "celltype_garnett.minor", "celltype_garnett.minor.ext", "celltype_garnett.subset", "celltype_garnett.subset.ext", "celltype_td.major", "celltype_td.major.score.max", "celltype_td.minor", "celltype_td.minor.score.max", "celltype_td.subset", "celltype_td.subset.score.max", "celltype_markercount", "celltype.cycling")
    cols_cell_type <- c("celltype_garnett.major", "celltype_td.major", "celltype_markercount", "celltype.cycling")

    # scores
    cols_score <- c("PScore")

    switch(as.character(diet_level),
	"0"={
		cols <- colnames(rna@meta.data)
	},
	"1"={
		cols <- c("Sample", "dataset", "species", "tumor.type", "donor", "sample.type", "technology", "condition", "treatment", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase", str_column_of_meta_data_cluster, "seurat_clusters", cols_epi_type, "SingleR.HTAPP_toolbox", "SingleR.HPCA", "SingleR.BED", cols_cell_type, cols_score, "cell.type", "CNV.value", "CNV.corr", "CNV.Pos", "CNV.type", "CNV.genes", "cluster.type", str_column_of_meta_data_harmony, "cluster.type.harmony", str_column_of_meta_data_cluster_multik, str_column_of_meta_data_harmony_multik)
		# https://www.rdocumentation.org/packages/Seurat/versions/4.0.0/topics/DietSeurat
		# include dimention reduction
		#if ("harmony" %in% names(rna)) {
		#	dimreducs <- union(dimreducs, "harmony")
		#}
		if ("umapharmony" %in% names(rna)) {
			dimreducs <- union(dimreducs, "umapharmony")
		}
	},
	"2"={
		# make_sc-rna-seq_seurat_obj.R
		cols <- c("Sample", "nCount_RNA", "nFeature_RNA", str_column_of_meta_data_cluster, "CNV.value", "CNV.corr", "CNV.Pos", "CNV.type", "cell.type")
	},
	{}
    ) # switch	

    vec_cols <- intersect(cols, colnames(rna@meta.data))
    rna@meta.data <- rna@meta.data[, vec_cols]

    # usage:
    # rna.df <- Embeddings(object = rna[["umap"]]) %>% as.data.frame
    # rna.df: cell x harmony (1:50) 
   
    rna <- DietSeurat(rna,
                counts = TRUE,
                data = TRUE, # data = FALSE currently not supported
                scale.data = FALSE,
                features = NULL,
                assays = NULL,
                #dimreducs = c('pca','tsne','FItSNE','umap')),
                dimreducs = dimreducs,
                graphs = NULL
    )

  } # if

  rna

} # get_diet_seurat_obj









##### other utilities



### stats



# get_med_mad
# https://github.com/davismcc/archive-scater/blob/master/R/qc.R
get_med_mad <- function(metric, nmads = 5, type = c("both", "lower", "higher"), log = FALSE, min.diff = NA, n_log=0, nround=1) {

    if (log) {
        metric <- log10(metric)
    }

    cur.med <- median(metric, na.rm = TRUE)
    cur.mad <- mad(metric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min.diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val 
    lower.limit <- cur.med - diff.val 

    if (n_log > 0) {
	switch(type,
		"lower"={
			p <- sprintf("\t\t\tmed - nmad*mad: %%.%df - %%.%df = %%.%df\n", nround, nround, nround)
			cat(sprintf(p, cur.med, diff.val, lower.limit))
		},
		"higher"={
			p <- sprintf("\t\t\tmed + nmad*mad: %%.%df + %%.%df = %%.%df\n", nround, nround, nround)
			cat(sprintf(p, cur.med, diff.val, upper.limit))
		},
		{}
	) # switch
    }

    list_out <- list()
    list_out$med <- cur.med
    list_out$mad <- cur.mad
    list_out$upper.limit <- upper.limit
    list_out$lower.limit <- lower.limit

} # get_med_mad






# summarize_meta_data
#
# input:
#   obj:
#   str_row:
#   str_colum:
#   type_summarization: {["dt"], "for_loop"}
#   pattern_remove_row: e.g. "[0-9]+-"
#   pattern_remove_col: e.g. "-.*$"
#
# usage:
# dt <- summarize_meta_data(rna, "Sample", col_seurat_cluster)
# dt <- summarize_meta_data(rna, "SingleR", col_seurat_cluster)
# dt <- summarize_meta_data(rna, col_cluster_types, "Sample")
# dt <- summarize_meta_data(rna, col_cluster_types, "Sample", pattern_remove_row="[0-9]+-")
# dt <- summarize_meta_data(atac, "predictedGroup_ArchR", "Sample")
summarize_meta_data <- function(obj, str_row, str_column, type_summarization="dt", pattern_remove_row=NULL, pattern_remove_col=NULL, tab=0, n_log=1) {

  switch(class(obj)[1],
	"Seurat"={
		md <- obj@meta.data %>% as.data.table
	},
	"ArchRProject"={
		md <- obj@cellColData %>% as.data.table
	},
	{}
  ) # switch

  if (!is.null(pattern_remove_row)) {
	md[, str_row] <- gsub(pattern_remove_row, "", md[, str_row])
  }

  if (!is.null(pattern_remove_col)) {
	md[, str_col] <- gsub(pattern_remove_col, "", md[, str_col])
  }

  switch(type_summarization,
	"dt"={
		# count the number of cells per unique combinations of "Sample" and "seurat_clusters"
		#dt.n <- md[, .N, by = c(str_row, str_column)]
		# columns: str_row, str_column, N
		#total_num <- sum(dt.n$N)

		# with additional casting after the counting
		dt <- md[, .N, by = c(str_row, str_column)] %>% data.table::dcast(., sprintf("%s ~ %s", str_row, str_column), value.var = "N") %>% as.data.table

		dt[, `:=`(SUM = rowSums(.SD, na.rm=T)), .SDcols=colnames(dt)[-1]]
		#dt[, `:=`(MIN = matrixStats::rowMins(as.matrix(.SD), na.rm=T), MAX = matrixStats::rowMaxs(as.matrix(.SD), na.rm=T), AVG = rowMeans(.SD, na.rm=T), SUM = rowSums(.SD, na.rm=T)), .SDcols=colnames(dt)[-1]]

		dt <- rbind(dt,
			cbind(setnames(data.table(c("SUM")), str_row), dt[, lapply(.SD, function(x) list(sum=sum(x, na.rm=TRUE))), .SDcols=2:ncol(dt)]))
			#cbind(setnames(data.table(c("MIN", "MAX", "AVG", "SUM")), str_row), dt[, lapply(.SD, function(x) list(min=min(x, na.rm=TRUE), max=max(x, na.rm=TRUE), mean=mean(x, na.rm=TRUE), sum=sum(x, na.rm=TRUE))), .SDcols=2:ncol(dt)]))

		out <- as.data.frame(dt)
		if (n_log > 0) {
			log_obj(out, tab=tab)
		}

	},
	"for_loop"={
		for (str_row1 in naturalsort::naturalsort(unique(md[,str_row]))) {

			cat(sprintf("\t%s\n", str_row1))
			idx_row1 <- which(md[,str_row] == str_row1)
			for (str_col1 in levels(factor(md[,str_column]))) {
				idx_col1 <- which(md[idx_row1, str_column] == str_col1)
				cat(sprintf("\t\t%s: %d\n", str_col1, length(idx_col1)))
			} # for

		} # for
		out <- TRUE
	},
	{}
  ) # switch

  out


} # summarize_meta_data












### files

# file.remove_when_file.exists
file.remove_when_file.exists <- function(filename) {

  out <- FALSE
  if (file.exists(filename)) {
        out <- file.remove(filename)
  }

  out

} # file.remove_when_file.exists





# gzip_when_file.exists
gzip_when_file.exists <- function(filename) {

  require(R.utils)

  out <- FALSE
  if (file.exists(filename)) {
        out <- gzip(filename, overwrite=TRUE)
  }

  out

} # gzip_when_file.exists











