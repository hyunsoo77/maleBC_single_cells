#
# find_markers.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2021, Dec.
#
# content:
#
# comment:
# called by make_sc-rna-seq_seurat_obj.R
#
# debug:
#
#


th_log2fc_global <- 0.25 # for seurat_findmarkers_enricher
th_padj_global <- 0.01

min.pct_global <- 0.5 # 0.5 for speed, default=0.1, matt=0.25, only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
min.diff.pct_global <- 0.25 # 0.25 for speed, default=-Inf, only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
max.cells.per.ident_global <- 200 # 200 for speed, default=Inf, # Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)


# find_all_markers
# 
# input:
#    args$dir_seurat_obj
#    assay: {NULL, RNA, integrated}
#    col_cluster_types: {cluster.type, cluster.type.harmony}
#    sample_id: 
#    f_save_rds:
#
# output:
#    ${dir_seurat_obj}/wilcox_degs/${cancer_type}_${col_cluster_types}_wilcox_degs.rds
#
# usage:
# markers <- find_all_markers(rna, args, col_cluster_types="cluster.type")
find_all_markers <- function(rna, args, assay=NULL, col_cluster_types="cluster.type", sample_id=NULL, th_log2fc=th_log2fc_global, th_padj=th_padj_global, min.pct=min.pct_global, min.diff.pct=min.diff.pct_global, max.cells.per.ident=max.cells.per.ident_global, n_log=0, f_save_rds=FALSE, dir_rds="./rds") {

  # perform DEGs analysis with cell type annotated clusters 
  if (n_log > 0) {
    cat(sprintf("\tFindAllMarkers for %s\n", col_cluster_types))
    cat(sprintf("\t\tmin.pct=%g\n", min.pct))
    cat(sprintf("\t\tmin.diff.pct=%g\n", min.diff.pct))
    cat(sprintf("\t\tmax.cells.per.ident=%g\n", max.cells.per.ident))

    cat(sprintf("\tFindAllMarkers for %s\n", col_cluster_types), file=stderr())
    cat(sprintf("\t\tmin.pct=%g\n", min.pct), file=stderr())
    cat(sprintf("\t\tmin.diff.pct=%g\n", min.diff.pct), file=stderr())
    cat(sprintf("\t\tmax.cells.per.ident=%g\n", max.cells.per.ident), file=stderr())
  } # if

  Idents(rna) <- col_cluster_types

  # https://satijalab.org/seurat/reference/findallmarkers
  # Finds markers (differentially expressed genes) for each of the identity classes in a dataset
  markers <- FindAllMarkers(object = rna,
    assay = assay,
    features = NULL,
    logfc.threshold = th_log2fc, # default
    test.use = "wilcox",
    slot = "data", # Slot to pull data from; note that if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts"
    min.pct = min.pct, # 0.5 for spped, default=0.1, matt=0.25, only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
    min.diff.pct = min.diff.pct, # 0.25 for speed, default=-Inf, only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
    node = NULL,
    verbose = TRUE,
    only.pos = FALSE,
    max.cells.per.ident = max.cells.per.ident, # 200 for speed, default=Inf, # Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
    random.seed = 1,
    latent.vars = NULL,
    min.cells.feature = 3, # Minimum number of cells expressing the feature in at least one of the two groups, currently only used for poisson and negative binomial tests
    min.cells.group = 3, # Minimum number of cells in one of the groups
    pseudocount.use = 1,
    mean.fxn = NULL,
    fc.name = NULL,
    base = 2,
    return.thresh = 0.01, # Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
    densify = FALSE
  ) # FindAllMarkers

  markers <- markers[markers$p_val_adj < th_padj, ]

  if (f_save_rds) {
    #saveRDS(markers, "./wilcox_DEGs.rds")
    dir.create(sprintf("%s/wilcox_degs", dir_rds), showWarnings = FALSE, recursive = TRUE)
    if (is.null(sample_id)) {
      path_rds <- sprintf("%s/wilcox_degs/%s_%s_wilcox_degs.rds", dir_rds, cancer_type, col_cluster_types)
    } else {
      path_rds <- sprintf("%s/wilcox_degs/%s_wilcox_degs.rds", dir_rds, sample_id)
    }
    cat(sprintf("\t\tsaveRDS(markers, '%s')\n", path_rds))
    saveRDS(markers, path_rds)
  }

  markers
  

} # find_all_markers








# find_markers
#
# input:
#   sample1:
#   sample_ref:
#   pattern_sample1:
#   pattern_sample_ref:
#   pattern_cnv1: 11, 01, 10, or 00, e.g. "11"
#   pattern_cnv_ref: 11, 01, 10 or 00, e.g. "01|10|00"
#   pattern_cell_type1:
#   pattern_cell_type_ref:
#   str_column_of_meta_data_cluster: e.g. RNA_snn_res.0.7, RNA_harmony_th.0
#   clusters1: e.g. c(5, 6, 12)
#   clusters_ref: e.g. c(0, 2, 3)
#   group_name1:
#   group_name_ref:
#   col_cluster_types: (e.g. cluster.type)
#   col_cell_types: (e.g. cell.type)
#
# output:
#   list_markers: list
#	markers: data.frame of columns {p_val, avg_log2FC, pct.1, pct.2, p_val_adj}
#   	sample1, sample_ref
#	group_name1, group_name_ref
#	idx1, idx_ref: idx of cells
#
# usage:
# list_out <- find_markers(rna, sample1=sample1, sample_ref=sample_ref, pattern_cell_type1=pattern_epi, pattern_cell_type_ref=pattern_epi, col_cell_types="cell.type")
# markers <- find_markers(rna, sample1=sample1, sample_ref=sample_ref, pattern_cell_type1="Fibroblasts", pattern_cell_type_ref="Fibroblasts", col_cell_types="cell.type")
find_markers <- function(rna, assay=NULL, slot="data", idx1=NULL, idx_ref=NULL, sample1=NULL, sample_ref=NULL, pattern_sample1=NULL, pattern_sample_ref=NULL, pattern_cnv1=NULL, pattern_cnv_ref=NULL, cell_type1=NULL, cell_type_ref=NULL, pattern_cell_type1=NULL, pattern_cell_type_ref=NULL, str_column_of_meta_data_cluster=NULL, clusters1=NULL, clusters_ref=NULL, genes=NULL, str_cond1=NULL, str_cond_ref=NULL, group_name1=NULL, group_name_ref=NULL, col_cluster_types=NULL, col_cell_types="cell.type", th_log2fc=th_log2fc_global, th_padj=th_padj_global, min.pct=min.pct_global, min.diff.pct=min.diff.pct_global, max.cells.per.ident=max.cells.per.ident_global, method_dge="seurat_findmarkers_enricher", str_condition=NULL, n_log=0) {

  if (n_log > 10) {
    log_txt(sprintf("\tFindMarkers\n"))
    log_txt(sprintf("\t\tth_log2fc=%g\n", th_log2fc))
    log_txt(sprintf("\t\tth_padj=%g\n", th_padj))
    log_txt(sprintf("\t\tmin.pct=%g\n", min.pct))
    log_txt(sprintf("\t\tmin.diff.pct=%g\n", min.diff.pct))
    log_txt(sprintf("\t\tmax.cells.per.ident=%g\n", max.cells.per.ident))
  } # if


  if (is.null(assay)) {
	assay <- DefaultAssay(rna)
  }
  if (n_log > 0) {
	log_txt(sprintf("\t\tassay=%s, slot=%s\n", assay, slot))
  }
  mtx <- GetAssayData(object = rna, assay = assay, slot = slot)

  # convert_cell_type_names
  group_name1 <- convert_cell_type_names(group_name1)
  group_name_ref <- convert_cell_type_names(group_name_ref)

  # idx1
  group_name1_ <- ""
  if (is.null(idx1)) {
	f <- rep(TRUE, nrow(rna@meta.data))
	if (!is.null(sample1)) {
		f <- f & (sample1 == rna@meta.data$Sample)
		group_name1_ <- sample1
	} else if (!is.null(pattern_sample1)) {
		f <- f & grepl(pattern_sample1, rna@meta.data$Sample, perl=T)
		group_name1_ <- gsub("[\\^\\*\\$]", "", pattern_sample1)
	} else if (!is.null(pattern_sample_ref)) {
		# pattern_sample1=NULL
		f <- f & !grepl(pattern_sample_ref, rna@meta.data$Sample, perl=T)
		group_name1_ <- paste("no", gsub("[\\^\\*\\$]", "", pattern_sample_ref))
	}

	if (!is.null(pattern_cnv1)) {
		f <- f & grepl(pattern_cnv1, rna@meta.data$CNV.Pos)
		group_name1_ <- paste(group_name1_, "CNV")
	}

	if (!is.null(cell_type1)) {
		if (!is.null(col_cluster_types)) {
			f <- f & (cell_type1 == rna@meta.data[,col_cluster_types])
		} 
		if (!is.null(col_cell_types)) {
			f <- f & (cell_type1 == rna@meta.data[,col_cell_types])
		}
		group_name1_ <- paste(group_name1_, cell_type1)
	} else if (!is.null(pattern_cell_type1)) {
		if (!is.null(col_cluster_types)) {
			f <- f & grepl(pattern_cell_type1, rna@meta.data[,col_cluster_types])
		} 
		if (!is.null(col_cell_types)) {
			f <- f & grepl(pattern_cell_type1, rna@meta.data[,col_cell_types])
		}
		group_name1_ <- paste(group_name1_, gsub("[\\^\\*\\$]", "", pattern_cell_type1))
	}

	if (!is.null(clusters1)) {
		if (is.null(str_column_of_meta_data_cluster)) {
			log_txt(sprintf("is.null(str_column_of_meta_data_cluster)=TRUE"))
		}
		f <- f & (as.vector(rna@meta.data[,str_column_of_meta_data_cluster]) %in% clusters1)
		group_name1_ <- paste(group_name1_, "clusters")
	}
	
	if (!is.null(genes) && !is.null(str_cond1)) {
		genes <- genes[genes %in% rownames(mtx)]
		if (length(genes) > 0) {
			vec <- Matrix::colMeans(mtx[genes,,drop=F])
			# str_cond1 <- "vec == 0", "vec < quantile(vec, 0.2)"
			f <- f & eval(parse(text=str_cond1))
			group_name1_ <- paste(group_name1_, genes[1])
		}
	}

	idx1 <- which(f)
  } # if
  if (is.null(group_name1)) group_name1 <- group_name1_
  

  # idx_ref
  group_name_ref_ <- ""
  if (is.null(idx_ref)) {
	f <- rep(TRUE, nrow(rna@meta.data))
	if (!is.null(sample_ref)) {
		f <- f & (sample_ref ==  rna@meta.data$Sample)
  		group_name_ref_ <- sample_ref
	} else if (!is.null(pattern_sample_ref)) {
		f <- f & grepl(pattern_sample_ref, rna@meta.data$Sample, perl=T)
		group_name_ref_ <- gsub("[\\^\\*\\$]", "", pattern_sample_ref)
	} else if (!is.null(pattern_sample1)) {
		# pattern_sample_ref=NULL
		f <- f & !grepl(pattern_sample1, rna@meta.data$Sample, perl=T)
		group_name_ref_ <- paste("no", gsub("[\\^\\*\\$]", "", pattern_sample1))
	}

	if (!is.null(pattern_cnv_ref)) {
		f <- f & grepl(pattern_cnv_ref, rna@meta.data$CNV.Pos)
		group_name_ref_ <- paste(group_name_ref_, "CNV")
	}

	if (!is.null(cell_type_ref)) {
		if (!is.null(col_cluster_types)) {
			f <- f & (cell_type_ref == rna@meta.data[,col_cluster_types])
		}
		if (!is.null(col_cell_types)) {
			f <- f & (cell_type_ref == rna@meta.data[,col_cell_types])
		}
		group_name_ref_ <- paste(group_name_ref_, cell_type_ref)
	} else if (!is.null(cell_type1)) {
		if (!is.null(col_cluster_types)) {
			f <- f & (cell_type1 != rna@meta.data[,col_cluster_types])
		} 
		if (!is.null(col_cell_types)) {
			f <- f & (cell_type1 != rna@meta.data[,col_cell_types])
		}
		group_name_ref_ <- paste("no", gsub("[\\^\\*\\$]", "", cell_type1))
	}

	if (!is.null(pattern_cell_type_ref)) {
		if (!is.null(col_cluster_types)) {
			f <- f & grepl(pattern_cell_type_ref, rna@meta.data[,col_cluster_types])
		}
		if (!is.null(col_cell_types)) {
			f <- f & grepl(pattern_cell_type_ref, rna@meta.data[,col_cell_types])
		}
		group_name_ref_ <- paste(group_name_ref_, gsub("[\\^\\*\\$]", "", pattern_cell_type_ref))
	}

	if (!is.null(clusters_ref)) {
		if (is.null(str_column_of_meta_data_cluster)) {
			log_txt(sprintf("is.null(str_column_of_meta_data_cluster)=TRUE"))
		}
		f <- f & (as.vector(rna@meta.data[,str_column_of_meta_data_cluster]) %in% clusters_ref)
		group_name_ref_ <- paste(group_name_ref_, "clusters")
	}

	if (!is.null(genes) && !is.null(str_cond_ref)) {
		genes <- genes[genes %in% rownames(mtx)]
		if (length(genes) > 0) {
			vec <- Matrix::colMeans(mtx[genes,,drop=F])
			# str_cond1 <- "vec > 0", "vec > quantile(vec, 0.8)"
			f <- f & eval(parse(text=str_cond_ref))
			group_name_ref_ <- paste(group_name_ref_, genes[1])
		}
	}

	idx_ref <- which(f)
  } # if
  if (is.null(group_name_ref)) group_name_ref <- group_name_ref_


  n_idx1 <- length(idx1)
  n_idx_ref <- length(idx_ref)
  if (n_log > 0) {
	log_txt(sprintf("\t\t%s n_idx1: %d, samples: %s\n", group_name1, n_idx1, paste(unique(rna$Sample[idx1]), collapse=", ")))
	log_txt(sprintf("\t\t%s n_idx_ref: %d, samples: %s\n", group_name_ref, n_idx_ref, paste(unique(rna$Sample[idx_ref]), collapse=", ")))
  }

  if (n_idx1 < 3) {
	# why < 3?
	# Error in ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2, : Cell group 1 has fewer than 3 cells
	return(NULL)
  }

  if (length(idx_ref) < 3) {
	# why < 3?
	# Error in ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2, : Cell group 1 has fewer than 3 cells
	return(NULL)
  }


  ident.1 <- rownames(rna@meta.data[idx1,])
  ident.ref <- rownames(rna@meta.data[idx_ref,])


  switch(method_dge,

        "seurat_findmarkers_enricher"={

		# https://satijalab.org/seurat/reference/findmarkers
		markers <- FindMarkers(
			rna,
			ident.1 = ident.1,
			ident.2 = ident.ref,
 			group.by = NULL,
			subset.ident = NULL,
			assay = NULL,
			slot = "data",
 			eduction = NULL, # Reduction to use in differential expression testing - will test for DE on cell embeddings

			features = NULL,
			logfc.threshold = th_log2fc,
			test.use = "wilcox",
			min.pct = min.pct, # 0.5 for speed, default=0.1, matt=0.25, only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
			min.diff.pct = min.diff.pct, # 0.25 for speed, default=-Inf, only test genes that show a minimum difference in the fraction of detection between the two groups. Set to -Inf by default
			max.cells.per.ident = max.cells.per.ident # 200 for speed, default=Inf, # Down sample each identity class to a max number. Default is no downsampling. Not activated by default (set to Inf)
		) # FindMarkers

		if (nrow(markers) == 0) {
			log_txt(sprintf("%s no markers found\n", group_name1))
			return(NULL)
		}

		idxsort <- order(markers$p_val)
		if (n_log > 100) {
			log_obj(head(markers[idxsort,]), tab=2)
		}

		markers <- markers[markers$p_val_adj < th_padj, ]
		# {p_val, avg_log2FC, pct.1, pct.2, p_val_adj}
	},

	"presto_wilcoxauc_gsea"={

                rna[["group_name"]] <- NA
                rna@meta.data[idx1, "group_name"] <- "grp1"
                rna@meta.data[idx_ref, "group_name"] <- "ref"

                markers <- presto::wilcoxauc(rna,
                        group_by = "group_name",
                        assay = "data", # name of feature matrix slot (e.g. 'data' or 'logcounts').
                        groups_use = c("ref", "grp1"), # must have at least 2 groups defined.
                        seurat_assay = DefaultAssay(rna),
                        verbose = TRUE
                ) # wilcoxauc
                # c("feature", "group", "avgExpr", "logFC", "statistic", "auc", "pval", "padj", "pct_in", "pct_out")
                # *auc* - area under the receiver operator curve.
                # *padj* - Benjamini-Hochberg adjusted p value.
                # *pct_in* - Percent of observations in the group with non-zero feature value.
                # *pct_out* - Percent of observations out of the group with non-zero feature value.

		if (nrow(markers) == 0) {
			log_txt(sprintf("%s no markers found\n", group_name1))
			return(NULL)
		}

                idx <- which(markers$group == "grp1")
                markers <- markers[idx, ]

                # sort by padj
                idx_order <- order(markers$padj)
                markers <- markers[idx_order, ]

                markers <- markers %>% remove_rownames %>% column_to_rownames(var="feature") %>% dplyr::select(-c("group"))

		if (n_log > 100) {
			log_obj(head(markers), tab=2)
		}

        },

	{}
  ) # switch


  # update information
  syms <- rownames(markers)
  if (length(syms) == 0) {
	markers$expr_min <- numeric(0)
	#markers$expr_med <- numeric(0)
	markers$expr_mean <- numeric(0)
	markers$expr_max <- numeric(0)
  } else if (length(syms) < 2) {
	markers$expr_min <- min(mtx[syms, ident.1], na.rm = T)
	#markers$expr_med <- median(mtx[syms, ident.1], na.rm = T)
	markers$expr_mean <- mean(mtx[syms, ident.1], na.rm = T)
	markers$expr_max <- max(mtx[syms, ident.1], na.rm = T)
  } else {
	markers$expr_min <- sparseMatrixStats::rowMins(mtx[syms, ident.1], na.rm = T)
	#markers$expr_med <- sparseMatrixStats::rowMedians(mtx[syms, ident.1], na.rm = T)
	markers$expr_mean <- rowMeans(mtx[syms, ident.1], na.rm = T)
	markers$expr_max <- sparseMatrixStats::rowMaxs(mtx[syms, ident.1], na.rm = T)
  }

  if (nrow(markers) == 0) {
	log_txt(sprintf("\t\t%s nrow(markers)=0\n", group_name1))
  }

  if (!is.null(str_condition)) {
	filename_dge <- sprintf("tsv/%s.tsv", str_condition)
	log_txt(sprintf("\t\twrite.table(markers, file='%s')\n", filename_dge))
	write.table(markers, file=filename_dge, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
  }


  list_out <- list()
  list_out$markers <- markers
  list_out$sample1 <- sample1
  list_out$sample_ref <- sample_ref
  list_out$group_name1 <- group_name1
  list_out$group_name_ref <- group_name_ref
  list_out$mtx <- mtx
  list_out$idx1 <- idx1
  list_out$idx_ref <- idx_ref

  list_out


} # find_markers






