# seurat_heatmap.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, June
#
# content:
# set_colors()
# plot_heatmap()
#



suppressPackageStartupMessages(library(Scillus))



# are_colors
# reference:
# https://github.com/xmc811/Scillus
are_colors <- function(x) {
        sapply(x, function(X) {
                tryCatch(is.matrix(col2rgb(X)),
                         error = function(e) FALSE)
        })
} # are_colors


# set_colors
# modified by H. Kim
# date last modified: 2022, June
#
# reference:
# https://github.com/xmc811/Scillus
# /disk2/github/scillus/Scillus/R/helpers.R
set_colors <- function(pal, n) {

	if (all(pal %in% rownames(brewer.pal.info))) {
		num <- c()
		for (i in seq(length(pal))) {
			num[i] <- brewer.pal.info[pal[i],][[1]]
		}
		full_pal <- do.call(c, purrr::map2(.x = num, .y = pal, .f = brewer.pal))
	} else if (all(are_colors(pal))) {
		full_pal <- pal
	} else {
		stop('Incorrect palette setup. Please input valid RColorBrewer palette names or color names.')
	}

	if (n <= length(full_pal)) {
		return(full_pal[1:n])
	} else {
		warning("Number of colors required exceeds palette capacity. RdYlBu spectrum will be used instead.", immediate. = TRUE)
		return(colorRampPalette(brewer.pal(11, "RdYlBu"))(n))
	}

} # set_colors



# get_top_genes_with_df_markers
# modified by H. Kim
# date last modified: 2022, June
#
# reference:
# https://github.com/xmc811/Scillus
# /disk2/github/scillus/Scillus/R/helpers.R
get_top_genes_with_df_markers <- function(dataset, df_markers, n, assay=NULL) {

	# begin of modification by H. Kim
        #int_features <- rownames(dataset@assays$integrated@scale.data)
	mtx <- GetAssayData(object = dataset, assay = assay, slot = "scale.data")
	int_features <- rownames(mtx)
	# end of modification

        df <- df_markers %>%
                filter(.data$gene %in% int_features) %>%
                arrange(desc(.data$avg_log2FC)) %>%
                group_by(.data$cluster) %>%
                filter(row_number() <= n) %>%
                arrange(.data$cluster)

        return(df$gene)

} # get_top_genes_with_df_markers






# plot_heatmap
# modified by H. Kim
# date last modified: 2022, June
#
# input:
#   dataset: matrix or seurat obj
#   markers: data.frame or characters
#   sort_var: column names of dataset@meta.data
#   anno_var: column names of dataset@meta.data
#   anno_colors: list of (named) vectors for anno_var
#   anno_name: annotation legend names
#   df_col_info: cells x var data.frame for column information when class(dataset)=matrix
#   list_annotation_legend_param: e.g. list(type=list(title="type", grid_height=unit(1, "cm"), ncol=2)),
#   heatmap_legend_param:
#   hm_limit: heatmap color limit
#   hm_colors: heatmap colors
#
# output:
#   list_out: list
#      list_ht:
#
# reference:
# https://github.com/xmc811/Scillus
# https://scillus.netlify.app/vignettes/plotting.html
# /disk2/github/scillus/Scillus/R/visualization.R
plot_heatmap <- function(dataset,
	markers,
	sort_var = c('seurat_clusters'),
	n = 8,
	anno_var,
	anno_colors,
	# begin of addition by H. Kim
	anno_name,
	df_col_info = NULL,
	list_annotation_legend_param = NULL,
	heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm"), title = "Expression"),
	simple_anno_size = unit(1, "cm"),
	# end of addition
	hm_limit = c(-2, 0, 2),
	hm_colors = c("#4575b4","white","#d73027"),
	row_font_size = 12) {




  # begin of modification by H. Kim

  # mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
  #mat <- mat[match(genes, rownames(mat)),]
  #anno <- dataset@meta.data %>%
  #	rownames_to_column(var = "barcode") %>%
  #	arrange(!!!syms(sort_var))

  class_dataset <- class(dataset)[1]
  if (class_dataset == "dgCMatrix") {
	dataset <- as.matrix(dataset)
	class_dataset <- "matrix"
  }
  switch(class_dataset,
        "matrix"={
                mat <- dataset
		genes <- markers
		mat <- mat[match(genes, rownames(mat)),,drop=F]
		anno <- df_col_info %>%
			tibble::rownames_to_column(var = "barcode") %>%
			arrange(!!!syms(sort_var))
		if (is.null(anno_colors)) {
			anno_colors <- list()
			for (anno_var1 in anno_var) {
				anno_colors[[anno_var1]] <- map_row_annot_color(df_col_info[,anno_var1])
			}
		}
        },
        "Seurat"={
		assay <- DefaultAssay(dataset)
		mat <- GetAssayData(object = dataset, assay = assay, slot = "scale.data")
		if (nrow(mat) == 0) {
			genes <- markers
			if (is.data.frame(markers)) genes <- markers$gene
			dataset <- ScaleData(object = dataset, features = genes, assay = assay, verbose = FALSE)
			mat <- GetAssayData(object = dataset, assay = assay, slot = "scale.data")
		} # if

		if (is.data.frame(markers)) {
			#genes <- get_top_genes(dataset, markers, n)
			genes <- get_top_genes_with_df_markers(dataset, markers, n, assay = assay)
		} else if (is.character(markers)) {
			genes <- markers
		} else {
			stop('Incorrect input of markers')
 		} # if

		mat <- mat[match(genes, rownames(mat)),]
		anno <- dataset@meta.data %>%
			tibble::rownames_to_column(var = "barcode") %>%
			arrange(!!!syms(sort_var))
	},
	{}
  ) # switch
  # end of modification




  mat <- t(mat)
  mat <- mat[match(anno$barcode, rownames(mat)),,drop=F]
  mat <- t(mat)

  annos <- list()
  for (i in seq_along(1:length(anno_var))) {

	err_msg <- paste('Incorrect specification for annotation colors for', anno_var[i])
	value <- anno[[anno_var[i]]]
	if (is.numeric(value)) {

		if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 'qual'])) {
			n <- brewer.pal.info[anno_colors[[i]],]['maxcolors'][[1]]
			pal <- brewer.pal(n = n, name = anno_colors[[i]])
			col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), c(pal[2], pal[(n+1)/2], pal[n-1]))
		} else if (length(anno_colors[[i]]) == 3 & all(are_colors(anno_colors[[i]]))) {
			col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), anno_colors[[i]])
		} else {
			stop(err_msg)
		} # if

		# begin of modificaiton by H. Kim
		#ha <- HeatmapAnnotation(a = anno[[anno_var[i]]], col = list(a = col_fun), border = TRUE, annotation_label = anno_var[i])
		ha <- HeatmapAnnotation(a = anno[[anno_var[i]]], col = list(a = col_fun), border = TRUE, annotation_label = anno_name[i])
		# end of modification

	} else {

		l <- levels(factor(anno[[anno_var[i]]]))
		if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
			col <- set_colors(anno_colors[[i]], length(l))
		} else if (length(anno_colors[[i]]) >= length(l) & all(are_colors(anno_colors[[i]]))) {
			col <- anno_colors[[i]]
		} else {
			stop(err_msg)
		} # if

		names(col) <- l
		col <- col[!is.na(names(col))]
		col <- list(a = col)

		# begin of modificaiton by H. Kim
		#ha <- HeatmapAnnotation(a = anno[[anno_var[i]]], col = col, border = TRUE, annotation_label = anno_var[i])
		names(col) <- "df"
		ha <- HeatmapAnnotation(df = anno[[anno_var[i]]], col = col, border = TRUE, annotation_label = anno_name[i], annotation_legend_param=list_annotation_legend_param[[anno_name[i]]], simple_anno_size = simple_anno_size)
		# end of modification
	} # if

	# begin of modification by H. Kim
	#names(ha) <- anno_var[i]
	names(ha) <- anno_name[i]
	# end of modification
	annos[[i]] <- ha

  } # for


  annos <- do.call(c, annos)
  annos@gap <- rep(unit(1,"mm"), length(annos))

  #m <- matrix(nrow = 0, ncol = ncol(mat))
  #colnames(m) <- colnames(mat)
  #m <- mat[1:2,]
  m <- mat
  ht <- Heatmap(m, 
	  cluster_rows = FALSE,
	  cluster_columns = FALSE,
	  # begin of modification by H. Kim
	  #heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm"), title = "Expression"),
	  heatmap_legend_param = heatmap_legend_param,
	  # end of modification
	  col = colorRamp2(hm_limit, hm_colors),
	  show_column_names = FALSE,
	  row_names_side = "left",
	  row_names_gp = gpar(fontsize = row_font_size),
	  top_annotation = annos)

  # begin of modification by H. Kim
  #ht
  list_out <- list()
  list_out$type_col_short <- "zscore"
  list_out$list_ht <- ht

  list_out
  # end of modification

} # plot_heatmap












