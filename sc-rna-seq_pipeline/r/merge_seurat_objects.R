#
# merge_seurat_objects.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Mar.
#
# content:
# merge_seurat_objects()
#
#
# 



# merge_seurat_objects
#
# input:
#   filenames_rds
#   args:
#	min_n_cells: minimum number of cells per object
#
# output:
#   list_out$rna:
#   list_out$args:
#	batch_vars: {dataset, tumor.type, donor, species, condition, treatment, Sample}
#	harmony_vars: same as batch_vars
#   list_out$df_sample.meta:
#	columns: nGenes, nCells, nCounts, median_nCount, median_nFeature
#

merge_seurat_objects <- function(filenames_rds, args) {


  cat(sprintf("------------------------------------\n"))
  cat(sprintf("\tmerge_seurat_objects\n\n"))

  cancer_type <- args$cancer_type

  if (is.null(args$filename_to_select_barcodes)) {
	args$filename_to_select_barcodes <- ""
  }
  if (is.null(args$pattern_to_select_cell_types)) {
	args$pattern_to_select_cell_types <- ""
  }
  if (is.null(args$samples_to_exclude)) {
	args$samples_to_exclude <- ""
  }



  obj <- list()
  vec_obj_n_cells <- rep(0, length(filenames_rds))
  vec_obj_sample_ids <- rep("", length(filenames_rds))


  ### loop for each filename_rds
  for (i in 1:length(filenames_rds)) {
  
    cat(sprintf("\t\tread %s\n", filenames_rds[i]))
    #obj[[i]] <- readRDS(filenames_rds[i])
    obj[[i]] <- tryCatch({
		  readRDS(filenames_rds[i])
	  }, error = function(e) {
		  cat(sprintf("\t\t\t#error-handler-code"))
		  cat(sprintf("\t\t\t%s",e))
		  return(NULL)
	  }, finally = {
    })
    
    if (is.null(obj[[i]])) {
  	next
    }
  
    # get_diet_seurat_obj
    obj[[i]] <- get_diet_seurat_obj(obj[[i]], args, f_force_diet=TRUE, dimreducs=NULL)
  

  
  
    # obj[[i]]$cell.type
    vec_cell_type <- get_vec_cell_types(obj[[i]], args)

    switch(args$method_to_identify_cell_types,
	  "singler_htapp_toolbox"={
  		# A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors, Slyper et al., Nat. Medicine, 2020. Single-cell analysis of tumors is rapidly expanding, including the launch of a Human Tumor Atlas Network (HTAPP) as part of the Cancer Moonshot7.  https://www.nature.com/articles/s41591-020-0844-1, Finally, we provide a website that displays a comprehensive analysis summary for each sample tested (https://tumor-toolbox.broadinstitute.org).
    		if ("SingleR.HTAPP_toolbox" %in% names(obj[[i]]@meta.data)) {

  			# Macrophage, Epithelial cell, T cell, B cell, NK cell
    			#vec_cell_type <- obj[[i]]$SingleR.HTAPP_toolbox

  			# Epithelial cell --> Epithelial cells, NK cell --> NK cells, T cell --> T-cells, B cell --> B-cells
  			vec_cell_type <- gsub("cell$", "cells", vec_cell_type)
  			vec_cell_type <- gsub("phage$", "phages", vec_cell_type)
  			vec_cell_type <- gsub("cyte$", "cytes", vec_cell_type)
  			vec_cell_type <- gsub("T cells", "T-cells", vec_cell_type)
  			vec_cell_type <- gsub("B cells", "B-cells", vec_cell_type)
  			vec_cell_type <- gsub("^DC$", "DCs", vec_cell_type)
    			obj[[i]]$cell.type <- vec_cell_type

  		} else {

    			obj[[i]]$cell.type <- obj[[i]]$SingleR.cancer_type
  		}
	  },
  
	  "singler_hpca_celldex"={
  		# Human Primary Cell Atlas Data (microarray)
  		# Epithelial_cells, NK_cell, T_cells, B_cells, Macrophage, Monocyte, DC
  		#vec_cell_type <- obj[[i]]$SingleR.HPCA

  		# Epithelial_cells --> Epithelial cells, NK_cell --> NK cells, T_cells --> T-cells, B_cells --> B-cells
  		vec_cell_type <- gsub("_", " ", vec_cell_type)
  		vec_cell_type <- gsub("cell$", "cells", vec_cell_type)
  		vec_cell_type <- gsub("phage$", "phages", vec_cell_type)
  		vec_cell_type <- gsub("cyte$", "cytes", vec_cell_type)
  		vec_cell_type <- gsub("T cells", "T-cells", vec_cell_type)
  		vec_cell_type <- gsub("B cells", "B-cells", vec_cell_type)
  		vec_cell_type <- gsub("^DC$", "DCs", vec_cell_type)
    		obj[[i]]$cell.type <- vec_cell_type
	  },
  
	  "singler_blueprint_encode"={
  		# BluePrint Encode (bulk RNA-seq)
  		# Epithelial cells, NK cells, CD8+ T-cells, B-cells, Macrophages, Monocytes, DC
    		#vec_cell_type <- obj[[i]]$SingleR.BED

  		vec_cell_type <- gsub("^DC$", "DCs", vec_cell_type)
    		obj[[i]]$cell.type <- vec_cell_type
	  },
	  {}
    ) # switch
  
  
  
    filename_rds <- basename(filenames_rds[i])
  
    # parsing filename_rds for sample_id, tumor.type, donor, sample.type
    tumor.type <- "unknown"; donor <- "unknown"; sample.type <- "unknown"

    types_parsing_rds_filename <- strsplit(args$type_parsing_rds_filename, ",")[[1]]
    for (type_parsing_rds_filename in types_parsing_rds_filename) {
  
      switch(type_parsing_rds_filename,
  	"2nd_item_after_parsing_with_underbar"={
    		# male-bc_4CC61L_sc-rna-seq_sample_seurat_obj.rds --> 4CC61L
  		items <- strsplit(filename_rds, split="_")[[1]]
  		sample_id <- items[2]
  	},
  	"gsm"={
    		# er+bc_GSM4909296_ER-MH0001_sc-rna-seq_sample_seurat_obj.rds --> GSM4909296
    		# normal-breast_GSM4909255_N-N280-Epi_sc-rna-seq_sample_seurat_obj.rds --> GSM4909255
  		sample_id <- gsub(".*_(GSM[0-9]+)_.*", "\\1", filename_rds)
  		if (sample_id == filename_rds) {
  			sample_id <- gsub(sprintf("%s_(.*)_sc-rna-seq_sample_seurat_obj.rds", cancer_type), "\\1", filename_rds)
  		}
  	},
  	"gsub"={
    		# normal-breast_49758L_12hrSUS_E2_sc-rna-seq_sample_seurat_obj.rds --> 49758L_12hrSUS_E2
  		cancer_type_tmp <- cancer_type
  		if (nchar(args$cancer_type_for_parsing_rds_filename) > 0) {
  			cancer_type_tmp <- args$cancer_type_for_parsing_rds_filename
  		}
  		items <- strsplit(cancer_type_tmp, ",")[[1]]
  		for (item in items) {
  			item <- gsub("\\+","\\\\\\+", item)
  			pattern <- sprintf("%s_(.*)_sc-rna-seq_sample_seurat_obj.rds", item)
  			sample_id <- gsub(pattern, "\\1", filename_rds)
  			if (sample_id != filename_rds) break
  		} # for
  	},
  	"unc-male-bc"={
    		# male-bc_4CC61L_sc-rna-seq_sample_seurat_obj.rds --> 4CC61L
  		items <- strsplit(filename_rds, split="_")[[1]]
  		sample_id <- items[2]
    		tumor.type <- "ER"
  		donor <- sample_id
  		sample.type <- ""
  	},
  	"GSE161529"={
  		# er+bc_GSM4909296_ER-MH0001_sc-rna-seq_sample_seurat_obj.rds --> tumor.type: ER, donor: MH0001, sample.type: ""
  		#filename_rds <- list.files("/datastore/nextgenout5/share/labs/francolab/hyunsoo.kim/sc-rna-seq/dataset/GSE161529/rds")
  		mtx <- str_match(filename_rds, "(.*)_(GSM[0-9]+)_(B1|ER|mER|HER2|N|TN|TN-B1)-((?:[A-Z0-9]+|Tum[0-9]+)(?:-7C|-9C)*)(-T|-T[0-9]+|-LN|-Total|-Epi)*_sc-rna-seq_sample")
		  if (is.na(mtx[1,1])) {
			  next
  		} else {
  			# B1, ER, mER
		  	tumor.type <- mtx[1,4]
  			# KCF0894, MH0033, PM0360, N280, NF, NE, N1105, N1B
  			donor <- mtx[1,5]
  			# NA, T3, T, LN, Total, Epi
  			sample.type <- ""
  			if (is.na(mtx[1,6])) {
  				sample_id <- sprintf("%s_%s", tumor.type, donor)
  			} else {
  				sample.type <- gsub("-", "", mtx[1,6])
  				if (sample.type == "LN") {
  					# ER-LN
  					tumor.type <- sprintf("%s-LN", mtx[1,4])
  					sample_id <- sprintf("%s_%s", tumor.type, donor)
  				} else {
  					sample_id <- sprintf("%s_%s_%s", tumor.type, donor, sample.type)
  				}
  			} 
  			break
  		} # if
  	},
  	{}
      ) # switch
    } # for 


    # apply sample_id
    cat(sprintf("\t\tsample_id: %s\n", sample_id))
    vec_obj_sample_ids[i] <- sample_id
    obj[[i]]$Sample <- sample_id

    # change barcodes
    obj[[i]] <- RenameCells(obj[[i]], new.names=paste0(sample_id, "#", colnames(obj[[i]])))
  
    # batch_vars
    batch_vars <- c()
  
    # dataset
    if (is.null(args$type_parsing_rds_filename_for_dataset)) {
	args$type_parsing_rds_filename_for_dataset <- ""
    }
    switch(args$type_parsing_rds_filename_for_dataset,
	"bc-scrnaseq"={
		obj[[i]]$dataset <- "bc-scrnaseq"
		batch_vars <- c(batch_vars, "dataset")
	},
  	"GSE161529"={
		obj[[i]]$dataset <- "GSE161529"
  		obj[[i]]$tumor.type <- tumor.type
  		obj[[i]]$donor <- donor
  		obj[[i]]$sample.type <- sample.type
		  batch_vars <- c(batch_vars, "tumor.type", "donor")
  	},
	{
  		if (nchar(args$type_parsing_rds_filename_for_dataset) > 0) {
		  	obj[[i]]$dataset <- args$type_parsing_rds_filename_for_dataset
		  	batch_vars <- c(batch_vars, "dataset")
  		}
  	}
    ) # switch
  
    # species
    if (is.null(args$type_parsing_rds_filename_for_species)) {
	args$type_parsing_rds_filename_for_species <- ""
    }
    switch(args$type_parsing_rds_filename_for_species,
	  "human"={
		  obj[[i]]$species <- "human"
		  batch_vars <- c(batch_vars, "species")
	  },
	  {
  		if (nchar(args$type_parsing_rds_filename_for_species) > 0) {
		  	obj[[i]]$species <- args$type_parsing_rds_filename_for_species
		  	batch_vars <- c(batch_vars, "species")
  		}
  	}
    ) # switch
  
    # tumor.type
    if (is.null(args$type_parsing_rds_filename_for_tumor.type)) {
	args$type_parsing_rds_filename_for_tumor.type <- ""
    }
    switch(args$type_parsing_rds_filename_for_tumor.type,
	  "1st_item_after_parsing_with_underbar"={
		  items <- strsplit(filename_rds, split="_")[[1]]
		  obj[[i]]$tumor.type <- items[1]
		  batch_vars <- c(batch_vars, "tumor.type")
	  },
	  {
  		if (nchar(args$type_parsing_rds_filename_for_tumor.type) > 0) {
  			obj[[i]]$tumor.type <- tumor.type
  			batch_vars <- c(batch_vars, "tumor.type")
  		}
  	}
    ) # switch
  
    # donor
    if (is.null(args$type_parsing_rds_filename_for_donor)) {
	args$type_parsing_rds_filename_for_donor <- ""
    }
    switch(args$type_parsing_rds_filename_for_donor,
	"2nd_item_after_parsing_with_underbar"={
		items <- strsplit(filename_rds, split="_")[[1]]
  		obj[[i]]$donor <- items[2]
  		batch_vars <- c(batch_vars, "donor")
	},
	"3rd_item_after_parsing_with_underbar"={
		items <- strsplit(filename_rds, split="_")[[1]]
  		obj[[i]]$donor <- items[3]
  		batch_vars <- c(batch_vars, "donor")
	},
  	{
  		if (nchar(args$type_parsing_rds_filename_for_donor) > 0) {
  			obj[[i]]$donor <- donor
  			batch_vars <- c(batch_vars, "donor")
  		}
  	}
    ) # switch	
  
    # technology
    if (is.null(args$type_parsing_rds_filename_for_technology)) {
	args$type_parsing_rds_filename_for_technology <- ""
    }
    switch(args$type_parsing_rds_filename_for_technology,
	  "not_implemented"={
	  },
	  {}
    ) # switch
  
    # condition
    if (is.null(args$type_parsing_rds_filename_for_condition)) {
 	args$type_parsing_rds_filename_for_condition <- ""
    }
    switch(args$type_parsing_rds_filename_for_condition,
	  "hrsus"={
  		mtx <- str_match(filename_rds, ".*([0-9]+hr[A-Z]+).*")
  		if (is.na(mtx[1,1])) {
		  	obj[[i]]$condition <- "unspecified"
  		} else {
		  	obj[[i]]$condition <- mtx[1,2]
  		}
  		batch_vars <- c(batch_vars, "condition")
  	},
	  {}
    ) # switch
  
    # treatment
    if (is.null(args$type_parsing_rds_filename_for_treatment)) {
	args$type_parsing_rds_filename_for_treatment <- ""
    }
    switch(args$type_parsing_rds_filename_for_treatment,
	"2nd_item_after_parsing_with_underbar"={
		items <- strsplit(filename_rds, split="_")[[1]]
		obj[[i]]$treatment <- items[2]
		batch_vars <- c(batch_vars, "treatment")
	},
	"3rd_item_after_parsing_with_underbar"={
		items <- strsplit(filename_rds, split="_")[[1]]
		obj[[i]]$treatment <- items[3]
		batch_vars <- c(batch_vars, "treatment")
	},
	"tam"={
  		mtx <- str_match(filename_rds, ".*(TAM).*")
  		if (is.na(mtx[1,1])) {
		  	obj[[i]]$treatment <- "no"
  		} else {
		  	obj[[i]]$treatment <- mtx[1,2]
  		}
  		batch_vars <- c(batch_vars, "treatment")
  	},
	  {}
    ) # switch
   
    batch_vars <- unique(batch_vars)
  
    if (length(batch_vars) == 0) {
  	batch_vars <- "Sample"
    } # if




    # determine n_cells
    n_cells <- length(colnames(obj[[i]]))
    cat(sprintf("\t\t# of cells: %d\n", n_cells))


    if (nchar(args$filename_to_select_barcodes) > 0) {

	# select cells with args$filename_to_select_barcodes
	df <- read.table(args$filename_to_select_barcodes, header=F, comment="")
	idx_cells <- match(df[,1], colnames(obj[[i]]))
	f <- !is.na(idx_cells); idx_cells <- idx_cells[f]
	obj[[i]] <- obj[[i]][, idx_cells]

	n_cells <- length(colnames(obj[[i]]))
	cat(sprintf("\t\t# of cells after selecting cell types with %s: %d\n", args$filename_to_select_barcodes, n_cells))

    } # if


    if (nchar(args$pattern_to_select_cell_types) > 0) {
  
	# select cells with args$pattern_to_select_cell_types
	switch(args$pattern_to_select_cell_types,
		"pattern_epi"={
    			idx_cell_type <- grep(pattern_epi, obj[[i]]$cell.type)
		},
		"pattern_normal_epi"={
    			idx_cell_type <- grep(pattern_normal_epi, obj[[i]]$cell.type)
		},
		"pattern_tumor_epi"={
    			idx_cell_type <- grep(pattern_tumor_epi, obj[[i]]$cell.type)
		},
		"pattern_stromal_cells"={
    			idx_cell_type <- grep(pattern_stromal_cells, obj[[i]]$cell.type)
		},
		"pattern_immune_cells"={
    			idx_cell_type <- grep(pattern_immune_cells, obj[[i]]$cell.type)
		},
		{
    			idx_cell_type <- grep(args$pattern_to_select_cell_types, obj[[i]]$cell.type)
		}
   	) # switch
	obj[[i]] <- obj[[i]][, idx_cell_type]

	n_cells <- length(colnames(obj[[i]]))
	cat(sprintf("\t\t# of cells after selecting cell types with %s: %d\n", args$pattern_to_select_cell_types, n_cells))

    } # if
  



    if (n_cells > args$max_n_cells) {
  
	# select cells with high counts when n_cells is larger than args$max_n_cells
	idx_cells <- select_cells_with_high_counts(obj[[i]], "cell.type", max_n_cells = args$max_n_cells, min_n_cells_per_cell_type = 20)
	obj[[i]] <- obj[[i]][, idx_cells]
  
	n_cells <- length(colnames(obj[[i]]))
	cat(sprintf("\t\t# of cells after applying max_n_cells: %d\n", n_cells))
  
    } # if
  
  
    vec_obj_n_cells[i] <- n_cells


  
  
  } # for i

  ### end of loop








  # store batch_vars
  args$batch_vars <- batch_vars
  args$harmony_vars <- batch_vars
  
  
  
  
  # exclude nonexistant items
  obj <- obj[!sapply(obj, is.null)]
  



  # exclude objs with too small number of cells
  idx_obj_excluded_snc <- which(vec_obj_n_cells < args$min_n_cells)
  if (length(idx_obj_excluded_snc) > 0) {

	cat(sprintf("\t\t\texclude %s due to n_cells=%d < %d\n", vec_obj_sample_ids[idx_obj_excluded_snc], vec_obj_n_cells[idx_obj_excluded_snc], args$min_n_cells)) 

  } # if


  # exclude objs with args$samples_to_exclude
  idx_obj_excluded_ste <- c()
  if (nchar(args$samples_to_exclude) > 0) {
	vec_samples_to_exclude <- strsplit(args$samples_to_exclude, ",")[[1]]

	idx_obj_excluded_ste <- match(vec_samples_to_exclude, vec_obj_sample_ids)
	f <- !is.na(idx_obj_excluded_ste)
	idx_obj_excluded_ste <- idx_obj_excluded_ste[f]
	if (length(idx_obj_excluded_ste) > 0) {
		cat(sprintf("\t\t\texclude %s since it is a sample to exclude\n", vec_obj_sample_ids[idx_obj_excluded_ste]))
	}

  } # if
  
  idx_obj_excluded <- Reduce(union, list(idx_obj_excluded_snc, idx_obj_excluded_ste))
  if (length(idx_obj_excluded) > 0) {
	obj[idx_obj_excluded] <- NULL
	vec_obj_sample_ids <- vec_obj_sample_ids[-idx_obj_excluded]
  } # if
  


  
  # df_sample.meta
  df_sample.meta <- data.frame()
  for (i in 1:length(obj)) {
  
    sample <- unique(obj[[i]]$Sample)
    mtx <- GetAssayData(object = obj[[i]], assay="RNA", slot = "counts")
    df_sample.meta[sample, "nGenes"] <- nrow(mtx)
    df_sample.meta[sample, "nCells"] <- ncol(mtx) # PostDoubletNumCells
    df_sample.meta[sample, "nCounts"] <- sum(sum(mtx, na.rm = TRUE)) # sum(rna$nCount_RNA)
    df_sample.meta[sample, "median_nCount"] <- median(obj[[i]]$nCount_RNA, na.rm = TRUE)
    df_sample.meta[sample, "median_nFeature"] <- median(obj[[i]]$nFeature_RNA, na.rm = TRUE)

    if ("tumor.type" %in% colnames(obj[[1]]@meta.data)) {
	tumor.type <- unique(obj[[i]]$tumor.type)
	df_sample.meta[sample, "tumor.type"] <- tumor.type
    }

    if ("donor" %in% colnames(obj[[1]]@meta.data)) {
	donor <- unique(obj[[i]]$donor)
	df_sample.meta[sample, "donor"] <- donor
    }
  
  } # for
  
  #file_name_sample.meta <- sprintf("tsv/%s_sc-rna-seq_pipeline_summary.tsv", cancer_type)
  #write.table(df_sample.meta, file_name_sample.meta, sep = "\t", row.names = TRUE, col.names = NA, quote=FALSE)
  
  print(df_sample.meta)
  
  
  
  
  
  if (length(obj) > 1) {
  
    # merge Seurat objects
    # https://satijalab.org/seurat/articles/merge_vignette.html
    cat(sprintf("\t\tmerge Seurat objects\n"))
    #rna <- merge(x = endo_3533EL, y = list(ovar_3BAE2L, ovar_3E5CFL))
    rna <- merge(x = obj[[1]], y = obj[2:length(obj)])
  
  } else {
  
    rna <- obj[[1]]
  
  } # if
  
  
  
  # print information of merged seurat obj
  #
  # e.g.)
  # An object of class Seurat
  # 17704 features across 6468 samples within 1 assay
  # Active assay: RNA (17704 features, 0 variable features)
  
  cat(sprintf("\t\tmerged Seurat object:\n"))
  print(rna)
  
  n_cells <- length(colnames(rna))
  cat(sprintf("\t\t# of cells: %d\n", n_cells))
  
  
  list_out <- list()
  list_out$rna <- rna
  list_out$args <- args
  list_out$df_sample.meta <- df_sample.meta

  list_out




} # merge_seurat_objects







