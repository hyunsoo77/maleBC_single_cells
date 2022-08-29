#
# jupyter_common_heatmap.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Mar.
#
# content:
# get_heatmap()
# get_heatmap_row_annotation()
# heatmap_pca()
# 




### heatmap

ch_ver <- packageVersion("ComplexHeatmap")
if (ch_ver >= 2.0) {
  # set global options for heatmaps
  # https://rdrr.io/github/jokergoo/ComplexHeatmap/man/ht_opt.html
  # `use_raster` is automatically set to TRUE for a matrix with more than 2000 columns You can control `use_raster` argument by explicitly setting TRUE/FALSE to it.
  ht_opt$message = FALSE
}







# get_heatmap
# input:
#   clustering_distance_rows: it can be a pre-defined character which is in ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"). It can also be a function. If the function has one argument, the input argument should be a matrix and the returned value should be a dist object. If the function has two arguments, the input arguments are two vectors and the function calculates distance between these two vectors.
#   clustering_method_rows: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
get_heatmap <- function(df, genes=NULL, type_col="zscore", col=NULL, platform="rnaseq", type_heatmap="details", row_pattern=NULL, column_pattern=NULL, row_title=NULL, row_title_side="left",row_title_gp=gpar(fontsize=10), column_title=NULL, column_title_side="top", column_title_gp=gpar(fontsize=10), show_column_names=T, show_row_names=T, show_heatmap_legend=T, heatmap_legend_param=list(title=type_col, title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=9)), df_top_annot=NULL, list_top_annotation_col=NULL, list_top_annotation_legend_param=NULL, left_annotation=NULL, right_annotation=NULL, fontsize_row_names=6, fontsize_column_names=5, cluster_rows=T, clustering_distance_rows="euclidean", clustering_method_rows="complete", cluster_columns=F, clustering_distance_columns="euclidean", clustering_method_columns="complete", sample_names=NULL, column_name_angle=NULL, ...) {


  ch_ver <- packageVersion("ComplexHeatmap")

  # special treatment
  colnames(df) <- mgsub::mgsub(colnames(df),
         c("tss", "_plusminus_", "_minus_", "_plus_", "five_prime_utr|5UTR", "CDS", "three_prime_utr|3UTR"),
         c("TSS", " ±", " -", " +", "5'UTR", "CDS", "3'UTR"))
  if (!is.null(row_pattern)) df <- df[grepl(row_pattern, rownames(df), perl=T), , drop=F]
  if (!is.null(column_pattern)) df <- df[, grepl(column_pattern, colnames(df),perl=T), drop=F]

  top_annotation <- NULL
  top_annotation_height <- NULL
  if (!is.null(df_top_annot)) {
    if (ch_ver < 2.0) {
      top_annotation=HeatmapAnnotation(df=df_top_annot, col=list_top_annotation_col, annotation_legend_param=list_top_annotation_legend_param, annotation_name_gp=gpar(fontsize=fontsize_column_names))
      top_annotation_height <- top_annotation@size
    } else {
      top_annotation=HeatmapAnnotation(df=df_top_annot, col=list_top_annotation_col, annotation_legend_param=list_top_annotation_legend_param, annotation_name_gp=gpar(fontsize=fontsize_column_names),
         annotation_height = NULL, # $ v4.x
         annotation_width = NULL,  # $ v4.x
         height = NULL,
         width = NULL
      )
      top_annotation_height <- top_annotation@anno_size # R v4.x
    }
  } # if

  bottom_annotation <- NULL
  bottom_annotation_height <- NULL
  if ((!is.null(column_name_angle)) && (show_column_names)) {
    # default: angle=90, bottom to top
    just="right"
    bottom_annotation_height <- max_text_width(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    if (column_name_angle == -90) {
       # angle=-90, top to bottom
       just='left'
    } else if (column_name_angle == 0) { just='top';
       bottom_annotation_height <- max_text_height(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    }
    if (ch_ver < 2.0) {
      bottom_annotation=HeatmapAnnotation(cn=anno_text(colnames(df),'column', gp=gpar(fontsize=fontsize_column_names), rot=column_name_angle, just=just, offset=unit(1,"npc")) )
    } else {
      bottom_annotation=HeatmapAnnotation(cn=anno_text(colnames(df),'column', gp=gpar(fontsize=fontsize_column_names), rot=column_name_angle, just=just, location=unit(1,"npc")),
         annotation_height = bottom_annotation_height,
         annotation_width = NULL,
         height = NULL,
         width = NULL
      )
    } # if
    show_column_names <- F
  } # if

  type_col_short <- tolower(gsub('[()]','',type_col))
  if (!is.null(genes)) {
    f <- !is.na(genes) & genes %in% rownames(df)
    genes <- genes[f]
    df <- df[genes,,drop=F]
    if (nrow(df) == 0) {
	return(list(type_col_short=type_col_short, list_ht=NULL))
    }
  } else if (nrow(df) > 5000) {
    return(list(type_col_short=type_col_short, list_ht=NULL))
  }

  # df_score
  df_score <- df

  switch(type_col_short,
    'zscore' = {
      if (class(df)[1] == "dgCMatrix") {
	df_score <- singleCellTK::computeZScore(df)
	#class(df_score) = "DelayedMatrix"
	#df_score <- as(df_score, "dgCMatrix")
	df_score <- as.matrix(df_score)
      } else {
        df_score <- t(scale(t(df_score), center=T, scale=T))
      }
      col_score <- colorRamp2(c(-2,0,2),c("#0000ff", "#ffffcc","#ff0000")) },   
    'log2cpm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },
    'log2fpkm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },         
    'log2fc'={
      col_score <- colorRamp2(c(-5,0,5),c("#0000ff","#ffffcc","#ff0000")) },
    '-log10pvalue'={
      type_col <- expression('-log'[10]*'(pvalue)')
      col_score <- colorRamp2(c(0,4),c("#ffffcc","#ff0000")) },
    'dlog10p'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },   
    'gsva'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },   
    'gsea_nes'={
      col_score <- colorRamp2(c(-5,0,5),c("#0000ff","#ffffcc","#ff0000")) },
    'ssgsea'={
      col_score <- colorRamp2(c(-0.5,0,0.5),c("#0000ff","#ffffcc","#ff0000")) },
    {}
  )
 
  if (!is.null(col)) col_score <- col
    
  switch(type_heatmap,
    "overview"={
       rect_gp=gpar(col=NA)
     },
     "details"={
       rect_gp=gpar(col="white", lty=1, lwd=1)
     },
     {}
  )
    
  if (!is.null(sample_names)) colnames(df_score) <- sample_names
  if (ch_ver < 2.0) {
    ht11=Heatmap(df_score, col=col_score,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=heatmap_legend_param,
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation, top_annotation_height=top_annotation_height,
        bottom_annotation=bottom_annotation,
        bottom_annotation_height=bottom_annotation_height,
        show_heatmap_legend=show_heatmap_legend,...)

  } else {

    ht11=Heatmap(df_score, col=col_score,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=heatmap_legend_param,
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation,
        # top_annotation_height` is removed. Set the height directly in `HeatmapAnnotation()  https://github.com/jokergoo/ComplexHeatmap/blob/master/R/Heatmap-class.R
	left_annotation=left_annotation,
	right_annotation=right_annotation,
        bottom_annotation=bottom_annotation,
        # bottom_annotation_height` is removed. Set the height directly in `HeatmapAnnotation()  https://github.com/jokergoo/ComplexHeatmap/blob/master/R/Heatmap-class.R
        show_heatmap_legend=show_heatmap_legend,...)
  } # if

  result <- list()
  result$type_col_short <- type_col_short
  result$list_ht <- ht11
  return(result)
} # get_heatmap



# get_heatmap_row_annotation
# input:
#   clustering_distance_rows: it can be a pre-defined character which is in ("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"). It can also be a function. If the function has one argument, the input argument should be a matrix and the returned value should be a dist object. If the function has two arguments, the input arguments are two vectors and the function calculates distance between these two vectors.
#   clustering_method_rows: the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
get_heatmap_row_annotation <- function(df, vec_color_group, nv_genes=NULL, vec_gene_type=NULL, type_col="zscore", nv_row_annot=NULL, platform="rnaseq", type_heatmap="details", row_pattern=NULL, column_pattern=NULL, row_title=NULL, row_title_side="left", row_title_gp=gpar(fontsize=10), column_title=NULL, column_title_side="top", column_title_gp=gpar(fontsize=10), show_column_names=T, show_row_names=T, show_heatmap_legend=T, heatmap_legend_param=list(title=type_col, title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=9)), row_annot_bar_width=unit(0.1,"in"), row_annot_max_len=15, row_annot_f_conv=F, row_annot_f_remove_first_word=F,row_annot_offset=-0.1, nv_sym_up=NULL, nv_sym_dn=NULL, df_top_annot=NULL, list_top_annotation_col=NULL, list_top_annotation_legend_param=NULL, left_annotation=NULL, right_annotation=NULL, fontsize_row_names=6, fontsize_column_names=5, fontsize_row_annot=8, cluster_rows=T, clustering_distance_rows="euclidean", clustering_method_rows="complete", cluster_columns=F, clustering_distance_columns="euclidean", clustering_method_columns="complete", sample_names=NULL,...) {

  ch_ver <- packageVersion("ComplexHeatmap")

  # special treatment
  colnames(df) <- mgsub::mgsub(colnames(df),
	 c("tss", "_plusminus_", "_minus_", "_plus_", "five_prime_utr|5UTR", "CDS", "three_prime_utr|3UTR"),
	 c("TSS", " ±", " -", " +", "5'UTR", "CDS", "3'UTR"))
  if (!is.null(row_pattern)) df <- df[grepl(row_pattern, rownames(df), perl=T), , drop=F]
  if (!is.null(column_pattern)) df <- df[, grepl(column_pattern, colnames(df), perl=T), drop=F]
  
  
  top_annotation <- NULL
  top_annotation_height <- NULL
  if (!is.null(df_top_annot)) {
    # column annotation
    if (ch_ver < 2.0) {
      top_annotation=HeatmapAnnotation(df=df_top_annot, col=list_top_annotation_col, annotation_legend_param=list_top_annotation_legend_param, annotation_name_gp=gpar(fontsize=fontsize_column_names))
      top_annotation_height <- top_annotation@size
    } else {
      top_annotation=HeatmapAnnotation(df=df_top_annot, col=list_top_annotation_col, annotation_legend_param=list_top_annotation_legend_param, annotation_name_gp=gpar(fontsize=fontsize_column_names),
         annotation_height = NULL,
         annotation_width = NULL,
         height = NULL,
         width = NULL
      )
      top_annotation_height <- top_annotation@anno_size
    }
  } # if


  bottom_annotation <- NULL
  bottom_annotation_height <- NULL
  if ((!is.null(column_name_angle)) && (show_column_names)) {
    # default: angle=90, bottom to top
    just="right"
    bottom_annotation_height <- max_text_width(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    if (column_name_angle == -90) {
       # angle=-90, top to bottom
       just='left'
    } else if (column_name_angle == 0) { just='top';
       bottom_annotation_height <- max_text_height(colnames(df), gp=gpar(fontsize=fontsize_column_names))
    }
    if (ch_ver < 2.0) {
      bottom_annotation=HeatmapAnnotation(cn=anno_text(colnames(df),'column', gp=gpar(fontsize=fontsize_column_names), rot=column_name_angle, just=just, offset=unit(1,"npc"))
      )
    } else {
      bottom_annotation=HeatmapAnnotation(cn=anno_text(colnames(df),'column', gp=gpar(fontsize=fontsize_column_names), rot=column_name_angle, just=just, location=unit(1,"npc")),
         annotation_height = bottom_annotation_height,
         annotation_width = NULL,
         height = NULL,
         width = NULL
      )
    } # if
    show_column_names <- F
  } # if

  type_col_short <- tolower(gsub('[()]','',type_col))
  if (!is.null(nv_genes)) {
    f <- !is.na(nv_genes) & nv_genes %in% rownames(df)
    nv_genes <- nv_genes[f]
    df <- df[nv_genes,,drop=F]
    if (nrow(df) == 0) {
	return(list(type_col_short=type_col_short, list_ht=NULL))
    }
    if (!is.null(vec_gene_type)) vec_gene_type <- vec_gene_type[f]
  } else if (nrow(df) > 5000) {
    return(list(type_col_short=type_col_short, list_ht=NULL))
  }

  # colors for row annotation
  if (!is.null(vec_gene_type)) {
    # gene groups
    df_annot <- data.frame(group=vec_gene_type, row.names=rownames(df))
    group_u <- unique(vec_gene_type); n_group <- length(group_u)
    col_group <- map_row_annot_color(group_u)
    vec_split <- vec_gene_type
    row_title_gp <- gpar(fontsize=0) # suppress the original row title
  } else {
    vec_gene_type <- names(nv_genes)
    df_annot <- data.frame(group=vec_gene_type, row.names=rownames(df))
    col_group_up <- map_row_annot_color(unique(names(nv_sym_up)),'up')
    col_group_dn <- map_row_annot_color(unique(names(rev(nv_sym_dn))),'dn')
    col_group <- c(col_group_up, col_group_dn)
    vec_split <- NULL
  }
    
  # row annotation
  ha12 <- rowAnnotation(df=df_annot, col=list(group=col_group),
    annotation_legend_param=list( group=list(title="gene set",
      title_gp=gpar(fontsize=8), labels_gp=gpar(fontsize=5))),
        width=row_annot_bar_width, show_legend=F)
   
  # df_score
  df_score <- df    
    
  switch(type_col_short,
    'zscore' = {
      if (class(df)[1] == "dgCMatrix") {
	df_score <- singleCellTK::computeZScore(df)
	#class(df_score) = "DelayedMatrix"
	#df_score <- as(df_score, "dgCMatrix")
	df_score <- as.matrix(df_score)
      } else {
	df_score <- t(scale(t(df_score),center=T,scale=T))
      }
      col_score <- colorRamp2(c(-2,0,2),c("#0000ff", "#ffffcc","#ff0000")) },
    'log2cpm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },
    'log2fpkm'={
      col_score <- colorRamp2(c(0,4,10),c("#0000ff","#ffffcc","#ff0000")) },            
    'log2fc'={
      col_score <- colorRamp2(c(-5,0,5),c("#0000ff","#ffffcc","#ff0000")) },
    '-log10pvalue'={
      type_col <- expression('-log'[10]*'(pvalue)')
      col_score <- colorRamp2(c(0,4),c("#ffffcc","#ff0000")) },
    'dlog10p'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },
    'gsva'={
      col_score <- colorRamp2(c(-4,0,4),c("#0000ff", "#ffffcc","#ff0000")) },
    'gsea_nes'={
      col_score <- colorRamp2(c(-5,0,5),c("#0000ff","#ffffcc","#ff0000")) },
    {}
  )
    
  switch(type_heatmap,
    "overview"={
       rect_gp=gpar(col=NA)
     },
     "details"={
       rect_gp=gpar(col="white", lty=1, lwd=1)
     },
     {}
  )

  if (!is.null(sample_names)) colnames(df_score) <- sample_names
  if (ch_ver < 2.0) {
    ht11 <- Heatmap(df_score, col=col_score,
        split=vec_split,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=heatmap_legend_param,
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation,
        bottom_annotation=bottom_annotation,
        bottom_annotation_height=bottom_annotation_height,
        show_heatmap_legend=show_heatmap_legend,...)
  } else {
    ht11 <- Heatmap(df_score, col=col_score,
        split=vec_split,
        cluster_rows=cluster_rows,
        clustering_distance_rows=clustering_distance_rows,
        clustering_method_rows=clustering_method_rows,
        cluster_columns=cluster_columns,
        clustering_distance_columns=clustering_distance_columns,
        clustering_method_columns=clustering_method_columns,
        show_column_names=show_column_names, show_row_names=show_row_names,
        row_names_gp=gpar(fontsize=fontsize_row_names),
        column_names_gp=gpar(fontsize=fontsize_column_names),
        rect_gp=rect_gp,
        heatmap_legend_param=heatmap_legend_param,
        row_title=row_title,
        row_title_side=row_title_side,
        row_title_gp=row_title_gp,
        row_title_rot = 90,
        column_title=column_title,
        column_title_side="top",
        column_title_gp=column_title_gp,
        top_annotation=top_annotation,
        # top_annotation_height` is removed. Set the height directly in `HeatmapAnnotation()  https://github.com/jokergoo/ComplexHeatmap/blob/master/R/Heatmap-class.R
	left_annotation=left_annotation,
	right_annotation=right_annotation,
        bottom_annotation=bottom_annotation,
        show_heatmap_legend=show_heatmap_legend,...)
  } # if

    
  # row annotation text
  ht_split <- ht11
  if (is.null(vec_split)) ht_split <- NULL

  # row annotation
  row_annot_text <- rep('',1,length(vec_gene_type))
  idx <- match(nv_row_annot, rownames(df_score))
  row_annot_text[idx] <- names(nv_row_annot)

  ha13 <- rowAnnotation(
         id=anno_text(row_annot_text,'row', gp=gpar(fontsize=fontsize_row_annot),
                rot=0, just="left", offset=row_annot_offset),
         #width=unit(0.5,"cm")+max_text_width(row_annot_text, just="left", gp=gpar(fontsize=fontsize_row_annot)),
         width=max_text_width(row_annot_text, just="left", gp=gpar(fontsize=fontsize_row_annot)))

  # result
  result <- list()
  result$type_col_short <- type_col_short
  result$list_ht <- ht11 + ha12 + ha13
  return(result)
    
} # get_heatmap_row_annotation




# heatmap_pca
# input:
#   list_pca:
#     df_loading: df_loading <- as.data.frame(prcomp_obj$rotation)
#   str_pc: {'PC1','PC2','PC3',...}
#   df_log2cpm:
#   vec_color_group:
#   n_pos:
#   n_neg:
#   fname_table:
#   f_format:
# output:
#   list_out$df_loading_sorted
# usage:
# list_out <- heatmap_pca(df_loading, 'PC2', df_log2cpm, vec_color_group=project$condition, n_pos=10, n_neg=10, cluster_rows = F, column_title='', column_annot_ncol=3, fontsize_row_names=6, fontsize_column_names=6, column_name_angle=45)
# print_figure(list_out$list_ht, width = 3.25, height = 3.25, file = sprintf("heatmap.pca2_%s", list_out$type_col_short))
# df <- cbind(list_out$df_dn[,'PC2',drop=F], df_log2cpm[rownames(list_out$df_dn),])
# cols <- colnames(df)
# df[,cols] <- sapply(cols, function(x) { formatC(df[,x], digits=2) })
# write.table(df, file = "table/pca_pc2_170224+190122.txt", row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
#
#
heatmap_pca <- function(list_pca, str_pc, df_log2cpm, vec_color_group, n_pos=10, n_neg=10, fname_table=NULL, f_formatc=T, ...) {

  if (is.list(list_pca)) {
    df_loading <- list_pca$df_loading
  } else {
    df_loading <- list_pca
  }

  df_loading_sorted <- df_loading[order(df_loading[,str_pc], decreasing=T),]
  pc_pos <- head(rownames(df_loading_sorted), n_pos)
  pc_neg <- tail(rownames(df_loading_sorted), n_neg)

  list_out <- get_heatmap(df_log2cpm, vec_color_group = vec_color_group,
      genes=c(pc_pos, pc_neg), type_col = "zscore", ...)
  list_out$df_loading_sorted <- df_loading_sorted

  if (!is.null(fname_table)) {
    df <- cbind(df_loading_sorted[,str_pc,drop=F], df_log2cpm[rownames(df_loading_sorted),])
    if (f_formatc) {
      cols <- colnames(df)
      df[,cols] <- sapply(cols, function(x) { formatC(df[,x], digits=2) })
    }
    write.table(df, file = fname_table, row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  }

  return(list_out)

} # heatmap_pca





