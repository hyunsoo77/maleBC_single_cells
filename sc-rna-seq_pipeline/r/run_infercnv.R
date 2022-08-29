#
# run_infercnv.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2021, Dec.
#
# contents:
# select_barcodes_reference_and_observation()
# get_infercnv_obj()
# run_infercnv()
# update_seurat_obj_with_infercnv()
#
# comment:
# called by make_sc-rna-seq_seurat_obj.R
#
# reference:
# https://bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html
# https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
# https://github-wiki-see.page/m/broadinstitute/inferCNV/wiki/Output-Files
#
# functions:
# https://rdrr.io/bioc/infercnv/man/CreateInfercnvObject.html
# https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
#
# 


# for pattern_epi
source("./r/utilities_for_sc_analyses.R")



file_name_gene_order_default <- "reference/genome_annotation/Homo_sapiens.GRCh38.86.symbol.txt"




# get_infercnv_obj
#
# input:
#   colname_group: {predoublet.idents, postdoublet.idents}
#
# usage:
# infercnv_obj <- get_infercnv_obj(rna, reference.clusters, "predoublet.idents", args, file_name_gene_order)
# 
get_infercnv_obj <- function(rna, reference.clusters, colname_group, args, file_name_gene_order=file_name_gene_order_default, file_name_annotation=NULL) {

  # create gene order file
  if (!file.exists(file_name_gene_order)) {
	GRCH38.annotations <- "reference/genome_annotation/Homo_sapiens.GRCh38.86.txt"
	df_gtf <- read.delim(GRCH38.annotations, header = F)
	df_gtf <- convert.symbol(df_gtf)
	#head(df_gtf)
	write.table(df_gtf, file_name_gene_order, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
  
  # select refernce barcodes and observation barcodes
  list_infercnv_input_barcodes <- select_barcodes_reference_and_observation(rna, reference.clusters, colname_group, args)

  if (length(list_infercnv_input_barcodes$barcodes_reference) == 0) {
	# e.g. no immune cells
	# no infercnv error even if there is no reference cells
	#message_stop <- sprintf("\n\n%s\nstop due to empty barcodes_reference\n%s\n\n",
	#     "------------------------------------",
	#     "------------------------------------")
	#cat(message_stop)
        #stop(message_stop)

	args$method_to_identify_reference_clusters <- "pattern_immune_cells"
	cat(sprintf("\tmethod_to_identify_reference_clusters (changed due to n_reference_cells=0): %s\n", args$method_to_identify_reference_clusters))

  } # if

  if (length(list_infercnv_input_barcodes$barcodes_observation) == 0) {
	# slurm_er+bc_4D8C8L_12hrSUS_E2.out   # of cells after post QC: 1848
	# slurm_er+bc_4D8C8L_12hrSUS_E2_TAM.out # of cells after post QC: 778
	message_stop <- sprintf("\n\n%s\nstop due to empty barcodes_observation\n%s\n\n",
	      "--------------------------------------",
	      "--------------------------------------")
	cat(message_stop)
	stop(message_stop)
  } # if
  

  # make_infercnv_input_mtx
  list_out <- make_infercnv_input_mtx(rna, list_infercnv_input_barcodes, args)

  counts_matrix <- list_out$counts_matrix
  df_infercnv_annot <- list_out$df_infercnv_annot
  ref_group_names <- list_out$ref_group_names

  if (length(list_out$barcodes_observation) == 0) {
	return(NULL)
  }

  if (is.null(file_name_annotation)) {
	sample_id <- sprintf("%s_%s", args$cancer_type, args$sample_id)
	file_name_annotation <- sprintf("%s/infercnv_input_barcode_group_%s.tsv", dir_tsv, sample_id)
  }

  #colnames(df_infercnv_annot) <- c("V1","V2")
  rownames(df_infercnv_annot) <- NULL
  colnames(df_infercnv_annot) <- NULL
  #write.table(df_infercnv_annot,"./sample_annotation_file_inferCNV.txt",sep = "\t",row.names = FALSE)
  #df_infercnv_annot <- read.delim("./sample_annotation_file_inferCNV.txt",header = F)
  write.table(df_infercnv_annot, file_name_annotation, sep = "\t", row.names = FALSE)
  #df_infercnv_annot <- read.delim(file_name_annotation, header = F)


  # https://github.com/broadinstitute/infercnv/blob/master/R/inferCNV.R
  # https://rdrr.io/bioc/infercnv/man/CreateInfercnvObject.html
  infercnv_obj <- CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
		annotations_file=file_name_annotation,
		gene_order_file=file_name_gene_order,
		ref_group_names=ref_group_names, # a vector containing the classifications of the reference (normal) cells to use for infering cnv
		delim = "\t",
		max_cells_per_group = NULL,
		min_max_counts_per_cell = c(100, +Inf),
		chr_exclude = c("chrX", "chrY", "chrM", "chrMT") # Because in most cases chrX, chrY and MT have very few genes, we remove them from the analysis altogether with the chr_exclude=c('chrX', 'chrY', 'chrM') option when creating the infercnv object.  https://github.com/broadinstitute/infercnv/issues/290
  ) # CreateInfercnvObject

  infercnv_obj

} # get_infercnv_obj













# select_barcodes_reference_and_observation
# comment:
# called by get_infercnv_obj()
select_barcodes_reference_and_observation <- function(rna, reference.clusters, colname_group, args) {

  cat(sprintf("\t[select_barcodes_reference_and_observation()]\n"))

  list_out <- list()

  vec_barcode <- rownames(rna@meta.data)
  vec_cluster <- rna@meta.data[, colname_group]
  vec_cell_types <- get_vec_cell_types(rna, args)

  clusters <- unique(vec_cluster)
  nv_clusters_cell_type <- rep("unknown", length(clusters))
  names(nv_clusters_cell_type) <- clusters

  #ref_group_names <- sprintf("immune.%s", reference.clusters),
  ref_group_names <- sprintf("%s%s", args$prefix_reference_clusters, reference.clusters)

  # barcodes_observation
  cat(sprintf("\tmethod_barcodes_observation: %s\n", args$method_barcodes_observation))
  barcodes_observation <- c()
  observation.clusters <- c()
  for (cluster in clusters) {

    idx_cells_cluster <- which(Idents(rna) == cluster)
    n_cells_cluster <- length(idx_cells_cluster)
    tb_cluster <- sort(table(vec_cell_types[idx_cells_cluster]), decreasing=T)

    # consider only obsrvation clusters
    if (cluster %in% ref_group_names) {
	idx_epi <- grep(pattern_epi, vec_cell_types)
	idx_epi_cell_cluster <- intersect(idx_cells_cluster, idx_epi)  
	n_epi_cell_cluster <- length(idx_epi_cell_cluster)
	if (n_epi_cell_cluster > 0) {
    		cat(sprintf("\treference cluster: %s\n", cluster))
		cat(sprintf("\t\tratio_epi_cells: %d/%d=%.2f\n", n_epi_cell_cluster, n_cells_cluster, n_epi_cell_cluster/n_cells_cluster))
		cat(sprintf("\t\tmost frequent cell.type: %s: %d\n", names(tb_cluster[1]), tb_cluster[1]))
	}
	next
    } # if

    cat(sprintf("\tobservation cluster: %s\n", cluster))
    cat(sprintf("\t\tn_cells: %d\n", n_cells_cluster))
    
    # apply additional filter
    idx_cells_cluster.filtered <- idx_cells_cluster

    f_apply_additional_filter <- TRUE
    if (f_apply_additional_filter) {
    	switch(args$method_barcodes_observation,

		"epi_cells_in_epi_clusters"={
			# warning: normal epithelial cells can be included.
			# cancer cells will be selected by correlation with the average vector of top copy number altered cells.
			idx_cells_cluster.filtered <- idx_cells_cluster
		},

		"epi_krt_epcam"={
			# warning: normal epithelial cells can be included.
			# cancer cells will be selected by correlation with the average vector of top copy number altered cells.
			idx_epi <- grep("^epi", rna$epi_krt_epcam)
			idx_cells_cluster.filtered <- intersect(idx_cells_cluster, idx_epi)  
		},

		"tumor_epithelial_cells"={
			# use rna$epi_normal_tumor set by identify_cell_types()
			# rna$epi_normal_tumor: {nonepi, epi, normal_epi}
			# calling tumor epithelial cells from non-normal epithelial cells tends to call many false positive tumor_epi, especially for normal samples. these tumor epithelial cells needs to be filtered by copy number alteration or other somatic mutations.
			idx_epi <- grep("^epi", rna$epi_normal_tumor)
			idx_cells_cluster.filtered <- intersect(idx_cells_cluster, idx_epi)  
		},

		"epi_singler"={
			idx_epi <- grep(pattern_epi, vec_cell_types)
			idx_cells_cluster.filtered <- intersect(idx_cells_cluster, idx_epi)  
		},
		{ }
	) # switch

	cat(sprintf("\t\tn_cells after applying additional filter: %d\n", length(idx_cells_cluster.filtered)))
    } # if

    if (length(idx_cells_cluster.filtered) <= 0) {
	log_obj(tb_cluster, tab=2)
	next
    }

    tb_cluster.filtered <- sort(table(vec_cell_types[idx_cells_cluster.filtered]), decreasing=T)
    if (length(names(tb_cluster.filtered)) == 0) {
	log_obj(tb_cluster, tab=2)
	next
    }
    nv_clusters_cell_type[cluster] <- names(tb_cluster.filtered)[1]

    cat(sprintf("\t\tmost frequent cell.type: %s: %d\n", names(tb_cluster.filtered)[1], tb_cluster.filtered[1]))

    if (length(idx_cells_cluster.filtered) <  args$min_num_cells_per_cluster) {
    	cat(sprintf("\t\tn_cells: %d < min_num_cells_per_cluster=%d\n", length(idx_cells_cluster.filtered), args$min_num_cells_per_cluster))
	#log_obj(tb_cluster, tab=2)
	next
    }
    
    if (args$f_infercnv_observations_only_epithelial) {

      if (!grepl(pattern_epi, names(tb_cluster.filtered)[1])) {
    	cat(sprintf("\t\tmost frequent cell.type is not epi\n"))
	#log_obj(tb_cluster, tab=2)
	next
      }

      observation.clusters <- c(observation.clusters, cluster)
      idx_epi <- grep(pattern_epi, vec_cell_types[idx_cells_cluster.filtered])
      if (length(idx_epi) > 0) {
          idx_epi_cell_cluster <- idx_cells_cluster.filtered[idx_epi]
    	  cat(sprintf("\t\tn_epi_cells: %d\n", length(idx_epi_cell_cluster)))
          barcodes_observation <- c(barcodes_observation, vec_barcode[idx_epi_cell_cluster])
      }

    } else {

      # observations: epithelial, endothelial, fibroblast, keratinocytes
      observation.clusters <- c(observation.clusters, cluster)
      barcodes_observation <- c(barcodes_observation, vec_barcode[idx_cells_cluster.filtered])

    } # if

  } # for

  epi.clusters <- names(nv_clusters_cell_type)[grepl(pattern_epi, nv_clusters_cell_type)]
  cat(sprintf("\tepi.clusters: %d/%d, %s\n", length(epi.clusters), length(clusters), paste(epi.clusters, collapse=", ")))
  cat(sprintf("\tobservation.clusters: %d/%d, %s\n", length(observation.clusters), length(clusters), paste(observation.clusters, collapse=", ")))
  cat(sprintf("\t\tn_observation_cells: %d\n", length(barcodes_observation)))


  # barcodes_reference, ref_group_names
  cat(sprintf("\tmethod_barcodes_reference: %s\n", args$method_barcodes_reference))

  barcodes_reference <- c()

  switch(args$method_barcodes_reference,

        "ref_group_names"={
          idx_ref <- which(vec_cluster %in% ref_group_names)
          barcodes_reference <- vec_barcode[idx_ref]
        },

        "nonepi_krt_epcam"={
          # other cells except epithelial cells
	  # e.g. immune.4 & non-epithelial cells defined by weaker expression of Keratins and Epcam
          idx_ref <- which( (vec_cluster %in% ref_group_names) & grepl("^nonepi", rna$epi_krt_epcam) )
          barcodes_reference <- vec_barcode[idx_ref]
        },

	"normal_epithelial_cells"={
          idx_ref <- which( (vec_cluster %in% epi.clusters) & grepl("normal_epi", rna$epi_normal_tumor) )
    	  cat(sprintf("\t\tn_normal_epi_cells: %d\n", length(idx_ref)))
          if (length(idx_ref) > 0) {
                levels(rna@meta.data[, colname_group]) <- unique(c(levels(rna@meta.data[, colname_group]), "normal_epi.99"))
                barcodes_normal_epi <- vec_barcode[idx_ref]
                rna@meta.data[barcodes_normal_epi, colname_group] <- "normal_epi.99"
                ref_group_names <- c(ref_group_names, "normal_epi.99")
                barcodes_reference <- barcodes_normal_epi
		#print(barcodes_normal_epi)
	  } else {
          	idx_ref <- (vec_cluster %in% ref_group_names)
         	barcodes_reference <- vec_barcode[idx_ref]
          }
	},

	"nonepi_singler"={
          idx_ref <- which( (vec_cluster %in% ref_group_names) & !grepl(pattern_epi, vec_cell_types) )
          barcodes_reference <- vec_barcode[idx_ref]
	},

        "add_normal_epithelial_cells"={
          idx_ref <- which(vec_cluster %in% ref_group_names)
          barcodes_reference <- vec_barcode[idx_ref]

	  # add normal_epi
          idx_ref <- which( (vec_cluster %in% epi.clusters) & grepl("normal_epi", rna$epi_normal_tumor) )
    	  cat(sprintf("\t\tn_normal_epi_cells: %d\n", length(idx_ref)))
	  if (length(idx_ref) >= args$min_num_cells_per_cluster) {
		levels(rna@meta.data[, colname_group]) <- unique(c(levels(rna@meta.data[, colname_group]), "normal_epi.99"))
          	barcodes_normal_epi <- vec_barcode[idx_ref]
          	rna@meta.data[barcodes_normal_epi, colname_group] <- "normal_epi.99"
	  	ref_group_names <- c(ref_group_names, "normal_epi.99")
	  	barcodes_reference <- c(barcodes_reference, barcodes_normal_epi)
	  }
        },
        {
        }
  ) # switch

  # update ref_group_names 
  ref_group_names <- ref_group_names[ ref_group_names %in% unique(rna@meta.data[barcodes_reference, colname_group]) ]
  cat(sprintf("\tref_group_names: %d/%d, %s\n", length(ref_group_names), length(clusters), paste(ref_group_names, collapse=", ")))
  cat(sprintf("\t\tn_reference_cells: %d\n", length(barcodes_reference)))


  list_out$colname_group <- colname_group
  list_out$ref_group_names <- ref_group_names
  list_out$observation.clusters <- observation.clusters
  list_out$barcodes_reference <- barcodes_reference
  list_out$barcodes_observation <- barcodes_observation
  
  list_out


} # select_barcodes_reference_and_observation















# make_infercnv_input_mtx
#
# input:
#   rna: seurat obj
#   list_infercnv_input_barcodes: output of select_barcodes_reference_and_observation()
#
# comment:
# called by get_infercnv_obj()
#
make_infercnv_input_mtx <- function(rna, list_infercnv_input_barcodes, args) {

  colname_group <- list_infercnv_input_barcodes$colname_group
  ref_group_names <- list_infercnv_input_barcodes$ref_group_names
  observation.clusters <- list_infercnv_input_barcodes$observation.clusters
  barcodes_reference <- list_infercnv_input_barcodes$barcodes_reference
  barcodes_observation <- list_infercnv_input_barcodes$barcodes_observation

  n_reference_cells <- length(barcodes_reference)
  counts_matrix <- GetAssayData(rna, assay="RNA", slot="counts")

  # to generate cnv information for all cells in observation.clusters
  #f_select_cells_with_clusters <- TRUE

  # to generate cnv information for tumor cells with cnv
  f_select_cells_with_clusters <- FALSE

  # special treatment for pattern_cell_type
  # update pattern_cell_type, ref_group_name, f_search_ref_inside_obs
  f_search_ref_inside_obs <- FALSE
  switch(args$method_to_identify_reference_clusters,

    "pattern_immune_cells"={
	args$method_to_identify_reference_clusters <- "pattern_cell_type"
	pattern_cell_type <- pattern_immune_cells
	ref_group_name <- "Immune cells"

    },

    "normal-like"={
	args$method_to_identify_reference_clusters <- "pattern_cell_type"
	pattern_cell_type <- "Normal-like"
	ref_group_name <- "Normal-like"
	f_search_ref_inside_obs <- TRUE
    },

    {}

  ) # switch
  


  # update counts_matrix, df_infercnv_annot, ref_group_names
  switch(args$method_to_identify_reference_clusters,

    "normal_sample_epi"={
    
	# read_rds_normal_samples
	rna_normal <- read_rds_normal_samples(args)

	# mtx_counts_normal, barcodes_reference
	idx_gene <- match(rownames(counts_matrix), rownames(rna_normal))
	f.gene <- !is.na(idx_gene); idx_gene <- idx_gene[f.gene]
	if (is.null(args$epi_type_as_infercnv_ref)) {
		idx_cell <- grep(pattern_normal_epi, rna_normal$cell.type)
	} else {
		switch(args$epi_type_as_infercnv_ref,
			"basal"={
				idx_cell <- grep(pattern_normal_epi_basal, rna_normal$cell.type)
			},
			"luminal"={
				idx_cell <- grep(pattern_normal_epi_luminal, rna_normal$cell.type)
			},
			{}
		) # swtich
	} # if

	rna_normal <- rna_normal[idx_gene, idx_cell]
  	mtx_counts_normal <- GetAssayData(rna_normal, assay="RNA", slot="counts")
	colnames(mtx_counts_normal) <- paste0(colnames(mtx_counts_normal), "n")
	barcodes_reference <- colnames(mtx_counts_normal)

        # f.cells.observation, barcodes_observation
	if (f_select_cells_with_clusters) {
		f.cells.observation <- (rna@meta.data[, colname_group] %in% observation.clusters)
	} else {
		f.cells.observation <- (rownames(rna@meta.data) %in% barcodes_observation)
	}

	if (!is.null(args$epi_type_as_infercnv_ref)) {
		vec_cell_types <- get_vec_cell_types(rna, args)
		f.cells.basal <- grepl(paste(pattern_tumor_basal, pattern_normal_epi_basal, sep="|"), vec_cell_types)
		f.cells.luminal <- grepl(paste(pattern_tumor_luminal, pattern_normal_epi_luminal, sep="|"), vec_cell_types)
		switch(args$epi_type_as_infercnv_ref,
			"basal"={
				f.cells.observation <- f.cells.observation & f.cells.basal
			},
			"luminal"={
				if (any(f.cells.luminal)) {
					f.cells.observation <- f.cells.observation & f.cells.luminal
				} else if (!any(f.cells.basal) & !any(f.cells.luminal)) {
					# try to use --method_to_update_cell_types cancer_epithelial_cell_types --method_to_identify_subtypes brca_pam50_genefu
					# check all epi. cells. just for debugging
					f.cells.observation <- f.cells.observation & grepl(pattern_epi, vec_cell_types)
				} else {
					# skip
				} # if
			},
			{}
		) # swtich
	} # if

	barcodes_observation <- colnames(counts_matrix)[f.cells.observation]
	counts_matrix <- counts_matrix[f.gene, f.cells.observation]
	counts_matrix <- cbind(counts_matrix, mtx_counts_normal)
        df_infercnv_annot <- data.frame(colnames(counts_matrix),
		c(rna@meta.data[f.cells.observation, colname_group], rna_normal$cell.type))

	# update ref_group_names
	ref_group_names <- unique(rna_normal$cell.type)

    },

    "pattern_cell_type"={

	# f.normal.cells
	vec_cell_types <- get_vec_cell_types(rna, args)
        f.normal.cells <- grepl(pattern_cell_type, vec_cell_types)
	f.tumor.cells <- grepl(pattern_tumor_cells, vec_cell_types)

	# f.cells.reference
	if (f_select_cells_with_clusters) {
		f.cells.reference <- (rna@meta.data[, colname_group] %in% ref_group_names)
		f.cells.observation <- (rna@meta.data[, colname_group] %in% observation.clusters)
	} else {
		f.cells.reference <- (rownames(rna@meta.data) %in% barcodes_reference)
		f.cells.observation <- (rownames(rna@meta.data) %in% barcodes_observation)
	}

	if (f_search_ref_inside_obs) {
		# select reference inside observation.clusters
		# e.g. pattern_cell_type = "Normal-like"
		f.cells.reference <- f.normal.cells & f.cells.observation
	} else if (any(f.cells.reference)) {
		# select normal cells inside reference.clusters
		# e.g. pattern_cell_type = pattern_immune cells
		f.cells.reference <- f.normal.cells & f.cells.reference
	} else {
		# no reference.cluster
		cat(sprintf("\t\tselect any normal cells inside observation.clusters\n"))
		f.cells.reference <- f.normal.cells
	}

	if (any(f.cells.observation)) {
		# select tumor cells inside observation.clusters
		#f.cells.observation <- f.tumor.cells & f.cells.observation
		# to generate cnv information for all cells in observation.clusters
		f.cells.observation <- f.cells.observation
	} else {
		# no observation.cluster
		cat(sprintf("\t\tselect any tumor cells inside reference.clusters\n"))
		f.cells.observation <- f.tumor.cells
	}

	if (!any(f.cells.reference)) {
		# no infercnv error even if there is no reference cells
		# --method_to_identify_reference_clusters normal-like --method_to_update_cell_types cancer_epithelial_cell_types --method_to_identify_subtypes brca_pam50_genefu
		#message_stop <- sprintf("\n\n%s\nno normal cells with pattern_cell_type=%s.\n%s\n\n",
	     	#	 "--------------------------------------", pattern_cell_type
	     	#	 "--------------------------------------")
		#cat(message_stop)
		#stop(message_stop)
		log_obj(vec_cell_types, tab=2)
		ref_group_name <- NULL
	}

	# vec_group
	vec_group <- as.character(rna@meta.data[, colname_group])
	vec_group[f.cells.reference] <- ref_group_name

	barcodes_reference <- colnames(counts_matrix)[f.cells.reference]
	barcodes_observation <- colnames(counts_matrix)[f.cells.observation]
	f.cells <- (f.cells.reference | f.cells.observation)
	counts_matrix <- counts_matrix[, f.cells]
        df_infercnv_annot <- data.frame(barcode=colnames(counts_matrix),
				 	group=vec_group[f.cells])
	
	# update ref_group_names
	ref_group_names <- ref_group_name

    },

    "low_total_cnv"={

	if (f_select_cells_with_clusters) {
		f.cells.observation <- (rna@meta.data[, colname_group] %in% observation.clusters)
	} else {
		f.cells.observation <- (rownames(rna@meta.data) %in% barcodes_observation)
	}

	barcodes_observation <- colnames(counts_matrix)[f.cells.observation]
        counts_matrix <- counts_matrix[, f.cells.observation]
	vec_group <- as.character(rna@meta.data[f.cells.observation, colname_group])

	vec_total_cnv <- rna$CNV.value
        vec_total_cnv <- vec_total_cnv[f.cells.observation]
	th_total_cnv <- quantile(vec_total_cnv, probs=c(0.20, 0.50, 0.80))
        th_total_cnv <- min(th_total_cnv, 0.1)
        idx.low_cna.cells <- which(vec_total_cnv < th_total_cnv)
	vec_group[idx.low_cna.cells] <- "Low CNV"

	barcodes_reference <- barcodes_observation[idx.low_cna.cells]
	barcodes_observation <- barcodes_observation[-idx.low_cna.cells]
        df_infercnv_annot <- data.frame(barcode=colnames(counts_matrix),
				 	group=vec_group)
	
	# update ref_group_names
	ref_group_names <- "Low CNV"

    },

    {
  	# make inferCNV input files
	if (f_select_cells_with_clusters) {
		# to generate cnv information for all cells in observation.clusters
		clusters <- c(ref_group_names, observation.clusters)
		f.cells <- (rna@meta.data[, colname_group] %in% clusters)
	} else {
		barcodes <- c(barcodes_reference, barcodes_observation)
		f.cells <- (rownames(rna@meta.data) %in% barcodes)
 	}

  	df_infercnv_annot <- data.frame(
		barcode=rownames(rna@meta.data)[f.cells],
		group=rna@meta.data[f.cells, colname_group])
    }

  ) # switch


  n_reference_cells_modified <- length(barcodes_reference)
  n_observation_cells_modified <- length(barcodes_observation)
  if (n_reference_cells_modified > 0) {
	cat(sprintf("\tmethod_to_identify_reference_clusters: %s\n", args$method_to_identify_reference_clusters))
	cat(sprintf("\tobservation.clusters: %d, %s\n", length(observation.clusters), paste(observation.clusters, collapse=", ")))
	cat(sprintf("\t\tn_observation_cells: %d\n", n_observation_cells_modified))
	cat(sprintf("\tref_group_names: %d, %s\n", length(ref_group_names), paste(ref_group_names, collapse=", ")))
	cat(sprintf("\t\tn_reference_cells: %d\n", n_reference_cells_modified))
  } # if


  list_out <- list()
  list_out$counts_matrix <- counts_matrix
  list_out$df_infercnv_annot <- df_infercnv_annot
  list_out$ref_group_names <- ref_group_names
  list_out$observation.clusters <- observation.clusters
  list_out$barcodes_reference <- barcodes_reference
  list_out$barcodes_observation <- barcodes_observation

  list_out


} # make_infercnv_input_mtx














# run_infercnv
#
# input:
#   args$type_infercnv_argset: {["vignettes"], "wiki", "matt", "default", "none"}
#    
# output files:
# https://github-wiki-see.page/m/broadinstitute/inferCNV/wiki/Output-Files
# Note, the current inferCNV software generates a lot of output files, and only a small subset are of greatest interest. In the future, outputs will be better organized to facilitate identification of the most useful outputs.
# 1) infercnv.preliminary.png : the preliminary inferCNV view (prior to denoising or HMM prediction)
# 2) infercnv.png : the final heatmap generated by inferCNV with denoising methods applied.
# 3) infercnv.references.txt : the 'normal' cell matrix data values.
# 4) infercnv.observations.txt : the tumor cell matrix data values
# 5) infercnv.observation_groupings.txt : group memberships for the tumor cells as clustered.
# 6) infercnv.observations_dendrogram.txt : the newick formatted dendrogram for the tumor cells that matches the heatmap.
#
run_infercnv <- function(infercnv_obj, dir_output, args) {

  dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

  switch(args$type_infercnv_argset,

    "vignettes"={
      # https://bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html
      # https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
      out <- capture.output(infercnv_obj <- infercnv::run(infercnv_obj,
                #cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                cutoff=args$infercnv_cutoff,
                min_cells_per_gene=3,
                #out_dir="./output_dir_CNV_predoublet",  # dir is auto-created for storing outputs
                out_dir=dir_output,  # dir is auto-created for storing outputs
                cluster_by_groups=TRUE, # cluster
                scale_data=FALSE,
                HMM=FALSE,
                analysis_mode="samples",
                denoise=TRUE,
		no_prelim_plot=TRUE,
		png_res=60,
                num_threads=8,
		plot_steps = FALSE,
		resume_mode = FALSE,
		save_rds = FALSE
      ) )
    },

    "wiki"={
      # https://github.com/broadinstitute/inferCNV/wiki/Running-InferCNV
      # https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
      out <- capture.output(infercnv_obj <- infercnv::run(infercnv_obj,
		#cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                cutoff=args$infercnv_cutoff,
		min_cells_per_gene=3,
		#out_dir="./output_dir_CNV_predoublet",  # dir is auto-created for storing outputs
		out_dir=dir_output,  # dir is auto-created for storing outputs
		cluster_by_groups=TRUE, # cluster
		scale_data=FALSE,
		HMM=TRUE,
		HMM_type="i6",
		analysis_mode="samples",
		BayesMaxPNormal=args$bayes_max_pnormal,
		denoise=TRUE,
		no_prelim_plot=TRUE,
		num_threads=8,
		plot_steps = FALSE,
		resume_mode = FALSE,
		save_rds = FALSE
      ) )
    },

    "matt"={
      # https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
      out <- capture.output(infercnv_obj <- infercnv::run(infercnv_obj,
		#cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                cutoff=args$infercnv_cutoff,
		min_cells_per_gene=10,
		 #out_dir="./output_dir_CNV_predoublet",  # dir is auto-created for storing outputs
		out_dir=dir_output,  # dir is auto-created for storing outputs
		cluster_by_groups=TRUE, # cluster
		scale_data=TRUE,
		HMM=TRUE,
		HMM_type="i6",
		analysis_mode="samples",
		BayesMaxPNormal=args$bayes_max_pnormal,
		denoise=TRUE,
		num_threads=8,
		plot_steps = FALSE,
		resume_mode = FALSE,
		save_rds = FALSE
      ) )
    },

    "default"={
      # https://rdrr.io/github/broadinstitute/inferCNV/man/run.html
      out <- capture.output(infercnv_obj <- infercnv::run(infercnv_obj,
		cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
		min_cells_per_gene = 3,
		out_dir=dir_output,  # dir is auto-created for storing outputs
		window_length = 101, # Length of the window for the moving average (smoothing). Should be an odd integer. (default: 101)
		#smooth_method = c("pyramidinal", "runmeans", "coordinates"),
		smooth_method = "pyramidinal",
		num_ref_groups = NULL,
		ref_subtract_use_mean_bounds = TRUE,
		cluster_by_groups = FALSE, # If observations are defined according to groups (ie. patients), each group of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
		cluster_references = TRUE, # Whether to cluster references within their annotations or not. (dendrogram not displayed) (default: TRUE)
		k_obs_groups = 1,
		hclust_method = "ward.D2",
		max_centered_threshold = 3,
		scale_data = FALSE, # perform Z-scaling of logtransformed data (default: FALSE). This may be turned on if you have very different kinds of data for your normal and tumor samples. For example, you need to use GTEx representative normal expression profiles rather than being able to leverage normal single cell data that goes with your experiment.
		HMM = FALSE,
		HMM_transition_prob = 1e-06,
		#HMM_report_by = c("subcluster", "consensus", "cell"),
		HMM_report_by = "subcluster",
		#HMM_type = c("i6", "i3"),
		HMM_type = "i6",
		HMM_i3_pval = 0.05,
		HMM_i3_use_KS = TRUE,
		#BayesMaxPNormal = 0.5,
		BayesMaxPNormal = args$bayes_max_pnormal,
		sim_method = "meanvar",
		sim_foreground = FALSE,
		reassignCNVs = TRUE,
		#analysis_mode = c("samples", "subclusters", "cells"),
		analysis_mode = "samples",
		#tumor_subcluster_partition_method = c("leiden", "random_trees", "qnorm", "pheight", "qgamma", "shc"),
		tumor_subcluster_partition_method = "leiden",
		tumor_subcluster_pval = 0.1,
		k_nn = 30,
		leiden_resolution = 1,
		denoise = FALSE,
		noise_filter = NA,
		sd_amplifier = 1.5,
		noise_logistic = FALSE,
		outlier_method_bound = "average_bound",
		outlier_lower_bound = NA,
		outlier_upper_bound = NA,
		final_scale_limits = NULL,
		final_center_val = NULL,
		debug = FALSE,
		num_threads = 8,
		plot_steps = FALSE, # if true, saves infercnv objects and plots data at the intermediate steps.
		resume_mode = TRUE, # leverage pre-computed and stored infercnv objects where possible. (default=TRUE)
		png_res = 300,
		plot_probabilities = TRUE,
		save_rds = TRUE, # Whether to save the current step object results as an .rds file (default: TRUE)
		save_final_rds = TRUE,
		diagnostics = FALSE,
		remove_genes_at_chr_ends = FALSE,
		prune_outliers = FALSE,
		mask_nonDE_genes = FALSE,
		mask_nonDE_pval = 0.05,
		test.use = "wilcoxon",
		require_DE_all_normals = "any",
		hspike_aggregate_normals = FALSE,
		no_plot = FALSE,
		no_prelim_plot = FALSE,
		output_format = "png",
		useRaster = TRUE,
		up_to_step = 100
      ) )
    },
    {
    }
  ) # switch

  # post-processing
  out <- file.remove_when_file.exists(sprintf("%s/expr.infercnv.dat", dir_output))
  out <- gzip_when_file.exists(sprintf("%s/infercnv.observations.txt", dir_output))
  out <- gzip_when_file.exists(sprintf("%s/infercnv.references.txt", dir_output))

  infercnv_obj

} # run_infercnv











# get_infercnv_obj_old
get_infercnv_obj_old <- function(counts_matrix, reference.clusters, file_name_annotation, file_name_gene_order=file_name_gene_order_default) {

    num.reference.clusters = length(reference.clusters)

    # create the infercnv object
    if ( num.reference.clusters == 1) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
             annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=paste0("immune.", reference.clusters[1]))
      
    } else if (num.reference.clusters == 2) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
             annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.", reference.clusters[1]),
			paste0("immune.", reference.clusters[2])))
      
    } else if ( num.reference.clusters == 3 ) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
             annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3])))

    } else if (num.reference.clusters == 4) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
             annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4])))

    } else if (num.reference.clusters == 5) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
               annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5])))

    } else if (num.reference.clusters == 6) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
               annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5]),
			paste0("immune.",reference.clusters[6])))

    }else if (num.reference.clusters == 7) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
               annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5]),
			paste0("immune.",reference.clusters[6]),
			paste0("immune.",reference.clusters[7])))
    }else if (num.reference.clusters == 8) {
      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
             #annotations_file="./sample_annotation_file_inferCNV.txt",
               annotations_file=file_name_annotation,
             gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5]),
			paste0("immune.",reference.clusters[6]),
			paste0("immune.",reference.clusters[7]),
			paste0("immune.",reference.clusters[8])))

    }else if (num.reference.clusters == 9) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
                #annotations_file="./sample_annotation_file_inferCNV.txt",
                annotations_file=file_name_annotation,
                gene_order_file=file_name_gene_order,
                ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5]),
			paste0("immune.",reference.clusters[6]),
			paste0("immune.",reference.clusters[7]),
			paste0("immune.",reference.clusters[8]),
			paste0("immune.",reference.clusters[9])))

    } else if (num.reference.clusters == 10) {

      infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(counts_matrix),
           #annotations_file="./sample_annotation_file_inferCNV.txt",
           annotations_file=file_name_annotation,
           gene_order_file=file_name_gene_order,
             ref_group_names=c(paste0("immune.",reference.clusters[1]),
			paste0("immune.",reference.clusters[2]),
			paste0("immune.",reference.clusters[3]),
			paste0("immune.",reference.clusters[4]),
			paste0("immune.",reference.clusters[5]),
			paste0("immune.",reference.clusters[6]),
			paste0("immune.",reference.clusters[7]),
			paste0("immune.",reference.clusters[8]),
			paste0("immune.",reference.clusters[9]),
			paste0("immune.",reference.clusters[10])))
    }


    infercnv_obj

} # get_infercnv_obj_old


























# update_seurat_obj_with_infercnv
#
# input:
#   infercnv_objs: not used currently, but for future usage.
#   colname_group: {predoublet.idents, postdoublet.idents}
#   args$method_to_update_seurat_obj_with_infercnv: {"peroulab", "cna_cor_per_cluster", "default"}
#   dirs_output: vector of dir_output
#
# output:
#   rna: Seurat object
#     CNV.value: sum of CNV values
#     CNV.Pos: 11=CNV, (01,10,00)=no CNV
#
update_seurat_obj_with_infercnv <- function(rna, infercnv_objs, colname_group, dirs_output, args) {

  if (is.null(dirs_output)) {
	return(rna)
  }

  switch(args$method_to_update_seurat_obj_with_infercnv,
    "peroulab"={
	rna <- update_seurat_obj_with_infercnv_wo_hmm(rna, infercnv_objs, colname_group, dirs_output, args)
     },
    "cna_cor_per_cluster"={
	rna <- update_seurat_obj_with_infercnv_wo_hmm_cna_cor_per_cluster(rna, infercnv_objs, colname_group, dirs_output, args)
     },
     {
	rna <- update_seurat_obj_with_infercnv_hmm(rna, infercnv_objs, colname_group, dirs_output, args)
     }

  ) # switch

  rna

} # update_seurat_obj_with_infercnv


















# update_seurat_obj_with_infercnv_wo_hmm
#
# correlation to itself
#
# input:
#   infercnv_objs: not used currently, but for future usage.
#   colname_group: {predoublet.idents, postdoublet.idents}
#   dirs_output: vector of dir_output
#
# output:
#   rna: seurat object
#     rna@meta.data[, "CNV.value"]
#     rna@meta.data[, "CNV.corr"]
#     rna@meta.data[, "CNV.Pos"]
#     rna@meta.data[, "CNV.type"]
#     rna@meta.data[, "CNV.genes"]
#
# reference:
# correlationforinfercnv_peroulab.R
# https://www.sciencedirect.com/science/article/pii/S0092867419306877
update_seurat_obj_with_infercnv_wo_hmm <- function(rna, infercnv_objs, colname_group, dirs_output, args) {


  suppressPackageStartupMessages(library(ICIKendallTau))

  rna@meta.data[, "CNV.value"] <- NA
  rna@meta.data[, "CNV.corr"] <- NA
  rna@meta.data[, "CNV.Pos"] <- NA

  vec_cluster <- as.character(rna@meta.data[, colname_group])
  vec_barcode <- rownames(rna@meta.data)
  vec_cell_types <- get_vec_cell_types(rna, args)

  counts <- data.frame()
  observations_metadata <- data.frame()

  
  cat(sprintf("\t[update_seurat_obj_with_infercnv_wo_hmm()]\n"))
  df_infercnv_obs <- NULL
  for (dir_output in dirs_output) {
	# 1) Read in your CNV values file for the TUMOR
	filename_obs <- sprintf("%s/infercnv.observations.txt.gz", dir_output)
	cat(sprintf("\t\tread %s\n", filename_obs))
	df_infercnv_obs_ <- read.table(filename_obs, sep="", header=TRUE) # data.frame of genes x cells, this takes time since it is a big data.frame (e.g. 7802 x 7092)
	if (is.null(df_infercnv_obs)) {
		df_infercnv_obs <- df_infercnv_obs_
	} else {
		#df_infercnv_obs <- cbind(df_infercnv_obs, df_infercnv_obs_)
		df_infercnv_obs <- merge(df_infercnv_obs, df_infercnv_obs_, by="row.names", all=TRUE)
		rownames(df_infercnv_obs) <- df_infercnv_obs$Row.names
		df_infercnv_obs <- df_infercnv_obs[2:length(df_infercnv_obs)]
	}
  } # for


  # CNV.genes
  rna <- update_seurat_obj_with_df_infercnv_obs(rna, df_infercnv_obs, args)


  # CNV.value, CNV.Pos
  infercnv_output <- as.data.frame(t(df_infercnv_obs)) # transpose to make data.frame of cells x genes
  # AGGGTGATCCATACTT.1 --> AGGGTGATCCATACTT-1
  rownames(infercnv_output) <- gsub(".([0-9]+)$", "-\\1", rownames(infercnv_output))
  barcodes_infercnv_output <- rownames(infercnv_output)


	if (args$f_infercnv_observations_only_epithelial) {
        	idx_epi <- grep(pattern_epi, vec_cell_types)
        	n_epi_cells <- length(idx_epi)
		cat(sprintf("\t\tn_epi_cells=%d\n", n_epi_cells))
		if (n_epi_cells < args$min_num_cells_per_cluster) {
			next
		}
		barcodes_epi <- vec_barcode[idx_epi]
	} else {
		barcodes_epi <- vec_barcode
 	} # if

	# infercnv_output_epi
	barcodes_epi <- intersect(barcodes_epi, barcodes_infercnv_output)
	cat(sprintf("\t\tn_barcodes_epi_clusters=%d\n", length(barcodes_epi)))
	infercnv_output_epi <- infercnv_output[barcodes_epi, ]

	# rescale all the CNV values for the TUMOR
	scaled_df <- as.data.frame(scales::rescale(as.matrix(infercnv_output_epi), c(-1,1))) #normalizing

	# correlate values to the top 5% of cells with high CNV values in the TUMOR
	CNV.values <- apply(scaled_df, 1, function(y) {
		#y[is.na(y)] <- 0
		#scaled_y <- rescale(y, c(-1, 1))
		return(mean(y^2))
	}) #remove the ngatives

	CNV.value_df <- data.frame(
		row.names = names(CNV.values),
		CNV.value = CNV.values
	)

	CNV_order <- order(CNV.value_df$CNV.value, decreasing=T)

	#ordering dataframe based on previous ordering
	ordered_CNV.values  <- data.frame(
		row.names = rownames(CNV.value_df)[CNV_order],
		CNV.value = CNV.value_df[CNV_order,]
	)

	# just extracting the top 0.05
	top_cancer <- head(ordered_CNV.values, max(nrow(ordered_CNV.values)*0.05, 2))
	top_cancer_CNV_average <- apply(infercnv_output_epi[rownames(top_cancer),], 2, mean)
	#top_cancer_CNV_average <- apply(as.data.frame(t(observ)), 2, mean)

	# correlation_df
	type_cor <- "ici_kt"
	switch(type_cor,
		"cor.test.kendall"={
			cancer_correlations <- apply(infercnv_output_epi, 1, function(x) {
				if (length(unique(as.numeric(x))) == 1) {
					cor_result <- data.frame(CNV.corr="no_CNVs_recorded", CNV.corr.p="no_CNVs_recorded")
				} else {
					cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
					cor_result <- data.frame(CNV.corr, CNV.corr.p)
				}
				return(cor_result)
			})
			correlation_df <- do.call("rbind", cancer_correlations)
		},

		"ici_kt"={
			cancer_correlations <- apply(infercnv_output_epi, 1, function(x) {
				if (length(unique(as.numeric(x))) == 1) {
					cor_result <- data.frame(CNV.corr="no_CNVs_recorded", CNV.corr.p="no_CNVs_recorded")
				} else {
					nv_cor <- ici_kt(as.numeric(x), top_cancer_CNV_average)
					cor_result <- data.frame(CNV.corr=nv_cor[["tau"]], CNV.corr.p=nv_cor[["pvalue"]])
				}
				return(cor_result)
			})
			correlation_df <- do.call("rbind", cancer_correlations)
		},

		"ici_kendalltau"={

			#suppressPackageStartupMessages(library(purrr))
			#suppressPackageStartupMessages(library(furrr))

			f_unique_rows <- apply(infercnv_output_epi, 1, function(x) { (length(unique(as.numeric(x))) == 1) })

			# ici_kendalltau
			list_out <- ici_kendalltau(rbind(top_cancer_CNV_average, as.matrix(infercnv_output_epi[!f_unique_rows,])), global_na = NA, perspective = "global", scale_max = FALSE, diag_good = TRUE, check_timing = FALSE)

			# correlation_df
			correlation_df <- infercnv_output_epi[,1,drop=F]
			colnames(correlation_df) <- c("CNV.corr")
			correlation_df[!f_unique_rows,1] <- list_out$cor[2:nrow(list_out$cor),1]
			correlation_df$CNV.corr.p <- NA

			correlation_df[f_unique_rows,1] <- "no_CNVs_recorded"
			correlation_df[f_unique_rows,2] <- "no_CNVs_recorded"
			
		},

		{}
	) # switch

	observations_metadata <- cbind(CNV.value_df, correlation_df)

	list_th <- determine_th_cna_value_corr(observations_metadata, args)
	th_cna_value <- list_th$th_cna_value
	th_cna_corr <- list_th$th_cna_corr


	#f_cancer <- (  observations_metadata$CNV.value > 0.02 & observations_metadata$CNV.corr > 0.4  ) # Aatish's criteria
 	#observations_metadata$Cancer.Cell <- ifelse(f_cancer, "Cancer","Normal")

	idx_11 <- which(  observations_metadata$CNV.value > th_cna_value & observations_metadata$CNV.corr > th_cna_corr  ) # Aatish's criteria
	cat(sprintf("\t\tn_11=%d\n", length(idx_11)))
	if (length(idx_11) > 0) {
		observations_metadata[idx_11, "CNV.Pos"] <- "11"
		observations_metadata[idx_11, "CNV.type"] <- "cancer"
	}

	idx_01 <- which(  observations_metadata$CNV.value <= th_cna_value & observations_metadata$CNV.corr > th_cna_corr  )
	cat(sprintf("\t\tn_01=%d\n", length(idx_01)))
	if (length(idx_01) > 0) {
		observations_metadata[idx_01, "CNV.Pos"] <- "01"
		observations_metadata[idx_01, "CNV.type"] <- "unassigned"
	}

	idx_10 <- which(  observations_metadata$CNV.value > th_cna_value & observations_metadata$CNV.corr <= th_cna_corr  )
	cat(sprintf("\t\tn_10=%d\n", length(idx_10)))
	if (length(idx_10) > 0) {
		observations_metadata[idx_10, "CNV.Pos"] <- "10"
		observations_metadata[idx_10, "CNV.type"] <- "unassigned"
	}
	
	idx_00 <- which(  observations_metadata$CNV.value <= th_cna_value & observations_metadata$CNV.corr <= th_cna_corr  )
	cat(sprintf("\t\tn_00=%d\n", length(idx_00)))
	if (length(idx_00) > 0) {
		observations_metadata[idx_00, "CNV.Pos"] <- "00"
		observations_metadata[idx_00, "CNV.type"] <- "normal"
	}

	observations_metadata$CNV.type <- factor(observations_metadata$CNV.type, level=c("cancer", "unassigned", "normal"))

	# update rna@meta.data
	barcodes_observations <- rownames(observations_metadata)
        cols <- c("CNV.value", "CNV.corr", "CNV.corr.p", "CNV.Pos", "CNV.type")
	rna@meta.data[barcodes_observations, cols] <- observations_metadata[, cols]
  

	# plot scatter plot
	
	observations_metadata_3cols <- observations_metadata[, c("CNV.value", "CNV.corr", "CNV.type")]
	idx <- which(observations_metadata_3cols$CNV.corr == "no_CNVs_recorded")
	if (length(idx) > 0) {
		observations_metadata_3cols <- observations_metadata_3cols[-idx,,drop=F]
		observations_metadata_3cols[,1:2] <- apply(observations_metadata_3cols[,1:2], 2, as.numeric)
	} # if

	gg <- ggplot(observations_metadata_3cols, aes(x=CNV.value, y=CNV.corr, color=CNV.type))+
		geom_point(size=2.0, alpha=0.5) +
		geom_vline(xintercept = th_cna_value, linetype="dashed") +
		geom_hline(yintercept = th_cna_corr, linetype="dashed") +
		#scale_color_manual(values=c("11"="red", "01"="grey", "10"="grey", "00"="blue"), labels=c("cancer", "unassigned", "unassigned", "normal")) +
		scale_color_manual(values=c("cancer"="red", "unassigned"="grey", "normal"="blue")) +
		ggtitle(sprintf("inferCNV CNV values vs. correlation")) +
		xlab("CNV values") +
		ylab("kendall tau with top 5% CN altered cells") +
		theme_bw()

	
	filename_plot <- sprintf("%s/scatterplot_infercnv_cna_vs_cor.pdf", dir_output)
	cat(sprintf("\t\tggsave %s\n", filename_plot))
	ggsave(filename_plot, width = 6, height = 5, plot=gg)



  write.table(observations_metadata, sprintf("%s/infercnv.observations_tumor_correlationvalues_toitself.txt", dir_output), sep="\t", quote=FALSE)

  rna

} # update_seurat_obj_with_infercnv_wo_hmm

























# update_seurat_obj_with_infercnv_wo_hmm_cna_cor_per_cluster
#
# correlation to itself
#
# input:
#   infercnv_objs: not used currently, but for future usage.
#   colname_group: {predoublet.idents, postdoublet.idents}
#   dirs_output: vector of dir_output
#
# output:
#   rna: seurat object
#     rna@meta.data[, "CNV.value"]
#     rna@meta.data[, "CNV.corr"]
#     rna@meta.data[, "CNV.Pos"]
#     rna@meta.data[, "CNV.type"]
#     rna@meta.data[, "CNV.genes"]
#
# reference:
# correlationforinfercnv_peroulab.R
# https://www.sciencedirect.com/science/article/pii/S0092867419306877
update_seurat_obj_with_infercnv_wo_hmm_cna_cor_per_cluster <- function(rna, infercnv_objs, colname_group, dirs_output, args) {

  suppressPackageStartupMessages(library(ICIKendallTau))

  rna@meta.data[, "CNV.value"] <- NA
  rna@meta.data[, "CNV.corr"] <- NA
  rna@meta.data[, "CNV.Pos"] <- NA

  vec_cluster <- as.character(rna@meta.data[, colname_group])
  vec_barcode <- rownames(rna@meta.data)
  vec_cell_types <- get_vec_cell_types(rna, args)

  df_counts <- data.frame()
  observations_metadata <- data.frame()

  
  cat(sprintf("\t[update_seurat_obj_with_infercnv_wo_hmm()]\n"))
  df_infercnv_obs <- NULL
  for (dir_output in dirs_output) {
	# 1) Read in your CNV values file for the TUMOR
	filename_obs <- sprintf("%s/infercnv.observations.txt.gz", dir_output)
	cat(sprintf("\t\tread %s\n", filename_obs))
	df_infercnv_obs_ <- read.table(filename_obs, sep="", header=TRUE) # data.frame of genes x cells, this takes time since it is a big data.frame (e.g. 7802 x 7092)
	if (is.null(df_infercnv_obs)) {
		df_infercnv_obs <- df_infercnv_obs_
	} else {
		#df_infercnv_obs <- cbind(df_infercnv_obs, df_infercnv_obs_)
		df_infercnv_obs <- merge(df_infercnv_obs, df_infercnv_obs_, by="row.names", all=TRUE)
		rownames(df_infercnv_obs) <- df_infercnv_obs$Row.names
		df_infercnv_obs <- df_infercnv_obs[2:length(df_infercnv_obs)]
	}
  } # for


  # CNV.genes
  rna <- update_seurat_obj_with_df_infercnv_obs(rna, df_infercnv_obs, args)


  # CNV.value, CNV.Pos
  infercnv_output <- as.data.frame(t(df_infercnv_obs)) # transpose to make data.frame of cells x genes
  # AGGGTGATCCATACTT.1 --> AGGGTGATCCATACTT-1
  rownames(infercnv_output) <- gsub(".([0-9]+)$", "-\\1", rownames(infercnv_output))
  barcodes_infercnv_output <- rownames(infercnv_output)

  observation.clusters <- unique(rna@meta.data[rownames(infercnv_output), colname_group])
  for (cluster in observation.clusters) {

	cat(sprintf("\tcluster=%s\n", cluster))
	if (args$f_infercnv_observations_only_epithelial) {
        	idx_cluster <- which((vec_cluster == cluster) & (grepl(pattern_epi, vec_cell_types)))
        	n_epi_cells <- length(idx_cluster)
		cat(sprintf("\t\tn_epi_cells=%d\n", n_epi_cells))
		if (n_epi_cells < args$min_num_cells_per_cluster) {
			next
		}
	} else {
        	idx_cluster <- which(vec_cluster == cluster)
		cat(sprintf("\t\tn_cells=%d\n", length(idx_cluster)))
        	if (length(idx_cluster) < 1) {
			next
		}
 	} # if

	barcodes_cluster <- vec_barcode[idx_cluster]

	# infercnv_output_cluster
	barcodes_cluster <- intersect(barcodes_cluster, barcodes_infercnv_output)
	infercnv_output_cluster <- infercnv_output[barcodes_cluster, ]

	# rescale all the CNV values for the TUMOR
	scaled_df <- as.data.frame(scales::rescale(as.matrix(infercnv_output_cluster), c(-1,1))) #normalizing

	# correlate values to the top 5% of cells with high CNV values in the TUMOR
	CNV.values <- apply(scaled_df, 1, function(y) {
		#y[is.na(y)] <- 0
		#scaled_y <- rescale(y, c(-1, 1))
		return(mean(y^2))
	}) #remove the ngatives

	CNV.value_df <- data.frame(
		row.names = names(CNV.values),
		CNV.value = CNV.values
	)

	CNV_order <- order(CNV.value_df$CNV.value, decreasing=T)

	#ordering dataframe based on previous ordering
	ordered_CNV.values  <- data.frame(
		row.names = rownames(CNV.value_df)[CNV_order],
		CNV.value = CNV.value_df[CNV_order,]
	)

	# just extracting the top 0.05
	top_cancer <- head(ordered_CNV.values, max(nrow(ordered_CNV.values)*0.05, 2))
	top_cancer_CNV_average <- apply(infercnv_output_cluster[rownames(top_cancer),], 2, mean)
	#top_cancer_CNV_average <- apply(as.data.frame(t(observ)), 2, mean)

	# correlation_df
	type_cor <- "ici_kt"
	switch(type_cor,
		"cor.test.kendall"={
			cancer_correlations <- apply(infercnv_output_cluster, 1, function(x) {
				if (length(unique(as.numeric(x))) == 1) {
					cor_result <- data.frame(CNV.corr="no_CNVs_recorded",CNV.corr.p="no_CNVs_recorded")
				} else {
					cor <- cor.test(as.numeric(x), top_cancer_CNV_average, method = "kendall")
					cor_result <- data.frame(CNV.corr, CNV.corr.p)
				}
				return(cor_result)
			})
			correlation_df <- do.call("rbind", cancer_correlations)
		},

		"ici_kt"={
			cancer_correlations <- apply(infercnv_output_cluster, 1, function(x) {
				if (length(unique(as.numeric(x))) == 1) {
					cor_result <- data.frame(CNV.corr="no_CNVs_recorded", CNV.corr.p="no_CNVs_recorded")
				} else {
					nv_cor <- ici_kt(as.numeric(x), top_cancer_CNV_average)
					cor_result <- data.frame(CNV.corr=nv_cor[["tau"]], CNV.corr.p=nv_cor[["pvalue"]])
				}
				return(cor_result)
			})
			correlation_df <- do.call("rbind", cancer_correlations)
		},

		"ici_kendalltau"={

			#suppressPackageStartupMessages(library(purrr))
			#suppressPackageStartupMessages(library(furrr))

			f_unique_rows <- apply(infercnv_output_cluster, 1, function(x) { (length(unique(as.numeric(x))) == 1) })

			# ici_kendalltau
			list_out <- ici_kendalltau(rbind(top_cancer_CNV_average, as.matrix(infercnv_output_cluster[!f_unique_rows,])), global_na = NA, perspective = "global", scale_max = FALSE, diag_good = TRUE, check_timing = FALSE)

			# correlation_df
			correlation_df <- infercnv_output_cluster[,1,drop=F]
			colnames(correlation_df) <- c("CNV.corr")
			correlation_df[!f_unique_rows,1] <- list_out$cor[2:nrow(list_out$cor),1]
			correlation_df$CNV.corr.p <- NA

			correlation_df[f_unique_rows,1] <- "no_CNVs_recorded"
			correlation_df[f_unique_rows,2] <- "no_CNVs_recorded"
			
		},

		{}
	) # switch

	observations_metadata_cluster <- cbind(CNV.value_df, correlation_df)
	observations_metadata_cluster$cluster <- cluster

	list_th <- determine_th_cna_value_corr(observations_metadata_cluster, args)
	th_cna_value <- list_th$th_cna_value
	th_cna_corr <- list_th$th_cna_corr


	#f_cancer <- (  observations_metadata_cluster$CNV.value > 0.02 & observations_metadata_cluster$CNV.corr > 0.4  ) # Aatish's criteria
 	#observations_metadata_cluster$Cancer.Cell <- ifelse(f_cancer, "Cancer","Normal")

	idx_11 <- which(  observations_metadata_cluster$CNV.value > th_cna_value & observations_metadata_cluster$CNV.corr > th_cna_corr  ) # Aatish's criteria
	cat(sprintf("\t\tn_11=%d\n", length(idx_11)))
	if (length(idx_11) > 0) {
		observations_metadata_cluster[idx_11, "CNV.Pos"] <- "11"
		observations_metadata_cluster[idx_11, "CNV.type"] <- "cancer"
	}

	idx_01 <- which(  observations_metadata_cluster$CNV.value <= th_cna_value & observations_metadata_cluster$CNV.corr > th_cna_corr  )
	cat(sprintf("\t\tn_01=%d\n", length(idx_01)))
	if (length(idx_01) > 0) {
		observations_metadata_cluster[idx_01, "CNV.Pos"] <- "01"
		observations_metadata_cluster[idx_01, "CNV.type"] <- "unassigned"
	}

	idx_10 <- which(  observations_metadata_cluster$CNV.value > th_cna_value & observations_metadata_cluster$CNV.corr <= th_cna_corr  )
	cat(sprintf("\t\tn_10=%d\n", length(idx_10)))
	if (length(idx_10) > 0) {
		observations_metadata_cluster[idx_10, "CNV.Pos"] <- "10"
		observations_metadata_cluster[idx_10, "CNV.type"] <- "unassigned"
	}
	
	idx_00 <- which(  observations_metadata_cluster$CNV.value <= th_cna_value & observations_metadata_cluster$CNV.corr <= th_cna_corr  )
	cat(sprintf("\t\tn_00=%d\n", length(idx_00)))
	if (length(idx_00) > 0) {
		observations_metadata_cluster[idx_00, "CNV.Pos"] <- "00"
		observations_metadata_cluster[idx_00, "CNV.type"] <- "normal"
	}

	observations_metadata_cluster$CNV.type <- factor(observations_metadata_cluster$CNV.type, level=c("cancer", "unassigned", "normal"))

	# update rna@meta.data
	barcodes_observations <- rownames(observations_metadata_cluster)
        cols <- c("CNV.value", "CNV.corr", "CNV.corr.p", "CNV.Pos", "CNV.type")
	rna@meta.data[barcodes_observations, cols] <- observations_metadata_cluster[, cols]
  
	df_counts <- rbind( df_counts, c(cluster, dim(observations_metadata_cluster)[1], length(idx_11), length(idx_01), length(idx_10), length(idx_00)) )
	observations_metadata <- rbind(observations_metadata, observations_metadata_cluster)


	# plot scatter plot

	observations_metadata_cluster_3cols <- observations_metadata_cluster[, c("CNV.value", "CNV.corr", "CNV.type")]
	idx <- which(observations_metadata_cluster_3cols$CNV.corr == "no_CNVs_recorded")
	if (length(idx) > 0) {
		observations_metadata_cluster_3cols <- observations_metadata_cluster_3cols[-idx,,drop=F]
		observations_metadata_cluster_3cols[,1:2] <- apply(observations_metadata_cluster_3cols[,1:2], 2, as.numeric)
	} # if

	gg <- ggplot(observations_metadata_cluster_3cols, aes(x=CNV.value, y=CNV.corr, color=CNV.type)) +
		geom_point(size=2.0, alpha=0.5) +
		geom_vline(xintercept = th_cna_value, linetype="dashed") +
		geom_hline(yintercept = th_cna_corr, linetype="dashed") +
		#scale_color_manual(values=c("11"="red", "01"="grey", "10"="grey", "00"="blue"), labels=c("cancer", "unassigned", "unassigned", "normal")) +
		scale_color_manual(values=c("cancer"="red", "unassigned"="grey", "normal"="blue")) +
		ggtitle(sprintf("inferCNV for cluster %s", cluster)) +
		xlab("CNV values") +
		ylab("kendall tau with top 5% CN altered cells") +
		theme_bw()

	filename_plot <- sprintf("%s/scatterplot_infercnv_cna_vs_cor_for_cluster%s.pdf", dir_output, cluster)
	cat(sprintf("\t\tggsave %s\n", filename_plot))
	ggsave(filename_plot, width = 6, height = 5, plot=gg)

  } # for






  if (nrow(df_counts) > 0) {
	df_counts[,3:6] <- apply(df_counts[,3:6], 2, as.numeric)
	df_counts$Total <- rowSums(df_counts[,3:6], na.rm = TRUE)
	colnames(df_counts) <- c("Cluster", "Epithelial", "11", "01", "10", "00", "Total")
	write.table(df_counts, sprintf("%s/infercnv.observations_tumor_correlationvalues_toitself_summary.txt", dir_output), sep="\t", quote=FALSE)
  }

  write.table(observations_metadata, sprintf("%s/infercnv.observations_tumor_correlationvalues_toitself.txt", dir_output), sep="\t", quote=FALSE)

  rna




} # update_seurat_obj_with_infercnv_wo_hmm_cna_cor_per_cluster




















# determine_th_cna_value_corr
#
# input:
#   observations_metadata
#   args$method_to_determine_th_cna_value_corr
#
# output:
#   list_th$th_cna_value
#   list_th$th_cna_corr
#
# comment:
# called by update_seurat_obj_with_infercnv_wo_hmm()
# called by update_seurat_obj_with_infercnv_wo_hmm_cna_cor_per_cluster()
#
# reference:
# https://www.nature.com/articles/s41588-021-00911-1
determine_th_cna_value_corr <- function(observations_metadata, args) {

	if (args$method_to_determine_th_cna_value_corr == "fixed") {
		th_cna_value <- args$th_cna_value
		th_cna_corr <- args$th_cna_corr
		cat(sprintf("\t\tth_cna_value: %.3f\n", th_cna_value))
		cat(sprintf("\t\tth_cna_corr: %.3f\n", th_cna_corr))

		list_th <- list()
		list_th$th_cna_value <- th_cna_value
		list_th$th_cna_corr <- th_cna_corr

		return(list_th)
	} # if




	# automatic determination of th_cna_value and th_cna_corr
	# https://www.nature.com/articles/s41588-021-00911-1

	suppressPackageStartupMessages(library(fpc))
	suppressPackageStartupMessages(library(cluster))


	observations_metadata_2cols <- observations_metadata[, c("CNV.value", "CNV.corr")]

	idx <- which(observations_metadata_2cols$CNV.corr == "no_CNVs_recorded")
	if (length(idx) > 0) {
		observations_metadata_2cols <- observations_metadata_2cols[-idx,,drop=F]
		mtx_observations_metadata <- apply(observations_metadata_2cols, 2, as.numeric)
		rownames(mtx_observations_metadata) <- rownames(observations_metadata_2cols)
		observations_metadata_2cols <- mtx_observations_metadata
	} # if

	observations_metadata_scaled <- scale(observations_metadata_2cols)

	criterion <- "asw"
	usepam <- TRUE
	f_use_multiasw_for_large_data <- FALSE
	if (f_use_multiasw_for_large_data && length(observations_metadata_scaled) > 2000) {
		# https://search.r-project.org/CRAN/refmans/fpc/html/pamk.html
		# "multiasw" is better for larger data sets, use larger ns then.
		criterion <- "multiasw"
		usepam <- FALSE
	} else {
		set.seed(1)
	} # if

	# https://www.rdocumentation.org/packages/fpc/versions/2.2-9/topics/pamk
	list_pamk_out <- pamk(observations_metadata_scaled, krange=1:4,
		criterion = criterion, # {asw, multiasw, ch}, average silhouette width
		usepam = usepam, # logical. If TRUE, pam is used, otherwise clara (recommended for large datasets with 2,000 or more observations; dissimilarity matrices can not be used with clara).
          	scaling = FALSE,
		alpha = 0.001,
		#diss = inherits(data, "dist"),
          	critout = FALSE,
		ns = 10, # passed on to distcritmulti if criterion="multiasw".
		seed = 1 # default=NULL, passed on to distcritmulti if criterion="multiasw".
	  ) # pamk

	cat(sprintf("\t\tpamk # of clusters with criterion=%s usepam=%s: %d\n", criterion, usepam, list_pamk_out$nc))
	# https://www.rdocumentation.org/packages/cluster/versions/2.1.3/topics/pam
	# Partitioning (clustering) of the data into k clusters “around medoids”, a more robust version of K-means.
	set.seed(1)
	pam_obj <- pam(observations_metadata_scaled, list_pamk_out$nc)
	cn_upper <- which.max(pam_obj$medoids[,2])
	cn_lower <- which.min(pam_obj$medoids[,2])

	barcodes_upper <- names(pam_obj$clustering)[pam_obj$clustering == cn_upper]
	idx_upper <- match(barcodes_upper, rownames(observations_metadata_2cols))
	if (any(is.na(idx_upper))) {
		stop("NAs in idx_upper")
 	}

	vec_cna_value <- na.omit(observations_metadata_2cols[idx_upper, "CNV.value"])
	vec_cna_corr <- na.omit(observations_metadata_2cols[idx_upper, "CNV.corr"])

	switch(args$method_to_determine_th_cna_value_corr,
		"ng_2021"={
			if (list_pamk_out$nc > 1) {
				# Thresholds defining normal and neoplastic cells were set at 2 cluster s.d. to the left and 1.5 s.d. below the first cancer cluster means.
				th_cna_value <- mean(vec_cna_value) - 2*sd(vec_cna_value)
				th_cna_corr <- mean(vec_cna_corr) - 1.5*sd(vec_cna_corr)
			} else {
				# For tumors where partitioning around medoids could not define more than 1 cluster, the thresholds were set at 1 s.d. to the left and 1.25 s.d. below the cluster means.
				th_cna_value <- mean(vec_cna_value) - 1*sd(vec_cna_value)
				th_cna_corr <- mean(vec_cna_corr) - 1.25*sd(vec_cna_corr)
			} # if
		},
		"ng_2021_and_thresholding"={
			if (list_pamk_out$nc > 1) {
				# Thresholds defining normal and neoplastic cells were set at 2 cluster s.d. to the left and 1.5 s.d. below the first cancer cluster means.
				th_cna_value <- mean(vec_cna_value) - 2*sd(vec_cna_value)
				th_cna_corr <- mean(vec_cna_corr) - 1.5*sd(vec_cna_corr)
			} else {
				# For tumors where partitioning around medoids could not define more than 1 cluster, the thresholds were set at 1 s.d. to the left and 1.25 s.d. below the cluster means.
				th_cna_value <- mean(vec_cna_value) - 1*sd(vec_cna_value)
				th_cna_corr <- mean(vec_cna_corr) - 1.25*sd(vec_cna_corr)
			} # if
			cat(sprintf("\t\tth_cna_value: %.3f\n", th_cna_value))
			cat(sprintf("\t\tth_cna_corr: %.3f\n", th_cna_corr))

			# apply maximum to exclude normal cells.
			cat(sprintf("\t\tapply maximum to exclude normal cells.\n"))
			# when cna_value < 0.01, no cna
			th_cna_value <- max(th_cna_value, 0.01)
			# when cna_corr < 0.2, no corr
			th_cna_corr <- max(th_cna_corr, 0.2)

			cat(sprintf("\t\tth_cna_value: %.3f\n", th_cna_value))
			cat(sprintf("\t\tth_cna_corr: %.3f\n", th_cna_corr))

			# apply minimum to include more cells.
			cat(sprintf("\t\tapply minimum to include more cells.\n"))
			# when th_cna_value > 0.05, use th_cna_value=0.05 since it is large enough.
			th_cna_value <- min(th_cna_value, 0.05)
			# when th_cna_corr > 0.4, use th_cna_corr=0.4 since it is large enough.
			th_cna_corr <- min(th_cna_corr, 0.4)
		},
		{}
	) # switch

	cat(sprintf("\t\tth_cna_value: %.3f\n", th_cna_value))
	cat(sprintf("\t\tth_cna_corr: %.3f\n", th_cna_corr))

	list_th <- list()
	list_th$th_cna_value <- th_cna_value
	list_th$th_cna_corr <- th_cna_corr

	list_th


} # determine_th_cna_value_corr





















# update_seurat_obj_with_df_infercnv_obs
update_seurat_obj_with_df_infercnv_obs <- function(rna, df_infercnv_obs, args) {


  # CNV gene names
  # args <- list()
  if (is.null(args$min_amp) || (args$min_amp < 0)) {
    col_numeric <- which( sapply(df_infercnv_obs, is.numeric ) )
    mtx_qt <- sapply( col_numeric, function( y ) {
      quantile( x = unlist( df_infercnv_obs[,  y ] ),
            c(.01, .99),
            na.rm = TRUE )
    })
    vec <- quantile(mtx_qt[1,], c(0.01, 0.05)); args$max_del <- vec[1]
    vec <- quantile(mtx_qt[2,], c(0.95, 0.99)); args$min_amp <- vec[2]
  }

  cat(sprintf("\tmin_amp: %.2f\n", args$min_amp))
  cat(sprintf("\tmax_del: %.2f\n", args$max_del))

  # remove genes amp/del for a small number of cells
  vec_n_cna <- apply(df_infercnv_obs, 1, function(x) { 
		n.del <- length(which(x <= args$max_del))
		n.amp <- length(which(x >= args$min_amp))
		if (n.del > 0 && n.amp > 0) return(0)
		n.del + n.amp
	})


  f.cna_genes <- (vec_n_cna > 100) # CNV in more than 100 cells
  if (exists("features_for_pca")) {
	# exclude genes of cell type biomarkers
	f.cna_genes <- f.cna_genes & (!rownames(df_infercnv_obs) %in% features_for_pca)
  }
  df_infercnv_cna <- df_infercnv_obs[f.cna_genes,,drop=F]

  
  rna@meta.data[, "CNV.genes"] <- ""
  x <- apply(df_infercnv_cna, 2, function(x) {
		vec_cna <- paste0(sort(rownames(df_infercnv_cna)[x >= args$min_amp]), " amp")
  		if (!is.null(args$list_cancer_type_specific_info$df_cna)) {
			vec_cna <- intersect(vec_cna, args$list_cancer_type_specific_info$df_cna$cna) # select amp
		}
		paste(vec_cna, collapse=", ")
	}); x[x  == " amp"] <- ""
  
  idx <- match(gsub(".([0-9]+)$", "-\\1", names(x)), colnames(rna))
  rna@meta.data[idx, "CNV.genes"] <- x

  x <- apply(df_infercnv_cna, 2, function(x) {
		vec_cna <- paste0(sort(rownames(df_infercnv_cna)[x <= args$max_del]), " del")
  		if (!is.null(args$list_cancer_type_specific_info$df_cna)) {
			vec_cna <- intersect(vec_cna, args$list_cancer_type_specific_info$df_cna$cna) # select del
		}
		paste(vec_cna, collapse=", ")
	}); x[x  == " del"] <- ""

  idx <- match(gsub(".([0-9]+)$", "-\\1", names(x)), colnames(rna))
  f.amp <- (nchar(rna@meta.data[idx, "CNV.genes"]) > 0)
  rna@meta.data[idx[f.amp], "CNV.genes"] <- sprintf("%s, %s", rna@meta.data[idx[f.amp], "CNV.genes"], x[f.amp]) # gene1 amp, gene2 del
  rna@meta.data[idx[!f.amp], "CNV.genes"] <- x[!f.amp] # gene2 del
  # tb <- table(rna$CNV.genes); table(tb)
 
  # frequent amp/del
  vec <- unlist(sapply(rna$CNV.genes, function(x) {
	strsplit(x, ", ")
  }) )
  tb <- table(vec)
  tb <- tb[order(tb, decreasing=T)]
  tb_amp <- tb[grepl("amp", names(tb))]
  tb_del <- tb[grepl("del", names(tb))]
  
  freq_amp <- names(head(tb_amp,5))
  freq_del <- names(head(tb_del,5))
  cat(sprintf("\tfrequent amp: %s\n", paste(freq_amp, collapse=", ")))
  cat(sprintf("\tfrequent del: %s\n", paste(freq_del, collapse=", ")))


  rna

} # update_seurat_obj_with_df_infercnv_obs














# update_seurat_obj_with_infercnv_hmm
#
# input:
#   infercnv_objs: not used currently, but for future usage.
#   colname_group: {predoublet.idents, postdoublet.idents}
#   dirs_output: vector of dir_output
#
# usage:
# rna <- update_seurat_obj_with_infercnv_hmm(rna, infercnv_objs, dirs_output, args)
update_seurat_obj_with_infercnv_hmm <- function(rna, infercnv_objs, colname_group, dirs_output, args) {

    cat(sprintf("\t[update_seurat_obj_with_infercnv_hmm()]\n"))

    # not implemented yet for dirs_output
    dir_output <- dirs_output

    #regions <- read.delim(sprintf("./output_dir_CNV_predoublet/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_%.1f.pred_cnv_regions.dat", args$bayes_max_pnormal))
    #probs <- read.delim("./output_dir_CNV_predoublet/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat")

    regions <- read.delim(sprintf("%s/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_%.1f.pred_cnv_regions.dat", dir_output, args$bayes_max_pnormal))
    probs <- read.delim(sprintf("%s/BayesNetOutput.HMMi6.hmm_mode-samples/CNV_State_Probabilities.dat", dir_output))

    # posterior probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV’s identified by inferCNV’s HMM  https://bioconductor.org/packages/devel/bioc/manuals/infercnv/man/infercnv.pdf
    # https://github-wiki-see.page/m/broadinstitute/infercnv/wiki/Infercnv-i6-HMM-type
    # State 1: 0x: complete loss
    # State 2: 0.5x: loss of one copy
    # State 3: 1x: neutral
    # State 4: 1.5x: addition of one copy
    # State 5: 2x: addition of two copies
    # State 6: 3x: essentially a placeholder for >2x copies but modeled as 3x.
    # CNV regions identified by the HMM are filtered out if the CNV region's posterior probability of being normal exceeds a specified threshold. This combats possibility of miss identified CNV's by removing CNV's that are most likely to be normal and not a true CNV events. By default this threshold is set to 0.5, given this any CNV region that has a posterior probability of being of a normal state greater than 0.5 is relabeled as "normal" and no longer considered an identified CNV region.  https://github-wiki-see.page/m/broadinstitute/infercnv/wiki/Infercnv-i6-HMM-type
    # If you want to use the HMM predictions, the regions are defined in the 17_HMM_predHMMi6.*.pred_cnv_regions.dat file with coordinates, HMM state (1 = 2 copies loss, 2 = 1 copy loss, 3 = neutral, 4 = 1 copy gain, 5 = 2 copies gain, 6 = 3+ copies gain), and subcluster id (17_HMM_predHMMi6.*.cell_groupings gives the subclusters to cells information).  https://github.com/broadinstitute/infercnv/issues/210
    probs <- as.data.frame(t(probs[3,]))

    colnames(probs) <- "Prob.Normal"
    probs <- dplyr::filter(probs, Prob.Normal < args$min_prob_normal)
    cnvs <- rownames(probs)
    cnvs <- gsub("\\.","-",cnvs)

    regions <- regions[regions$cnv_name %in% cnvs, ]
    cnv.groups <- sub("\\..*", "", regions$cell_group_name)

    # begin of modification by H. Kim
    #rna$CNV.Pos <- ifelse(as.character(rna$predoublet.idents) %in% cnv.groups, 1, 0)
    vec_group <- as.character(rna@meta.data[,colname_group])
    #rna$CNV.Pos <- ifelse(vec_group %in% cnv.groups, 1, 0)
    rna$CNV.Pos <- ifelse(vec_group %in% cnv.groups, "11", "00")

    cnv.freq <- data.frame(table(regions$cell_group_name))
    cnv.freq$Var1 <- sub("\\..*", "", cnv.freq$Var1)

    rownames(cnv.freq) <- cnv.freq$Var1
    rna$CNV.value <- ifelse(vec_group %in% cnv.freq$Var1, cnv.freq[vec_group, "Freq"], 0)
    # end of modification

    rna

} # update_seurat_obj_with_infercnv_hmm



















# read_rds_normal_samples
read_rds_normal_samples <- function(args) {


  # filenames_rds, args_normal
  if (grepl("-bc|+bc|tnbc|-breast", args$cancer_type)) {

     dir_rds_normal <- "/datastore/nextgenout5/share/labs/spanheimer_lab/hyunsoo.kim/sc-rna-seq/tr-bc/run-20220121/rds_normal-breast"
     samples <- c("49758L")
     filenames_rds <- sprintf("%s/normal-breast_%s_sc-rna-seq_sample_seurat_obj.rds", dir_rds_normal, samples)
     args_normal <- args
     args_normal$type_parsing_rds_filename <- "gsub"
     args_normal$cancer_type_for_parsing_rds_filename <- "normal-breast"

  } else if (grepl("-oc", args$cancer_type)) {

    stop(sprintf("no reference data defined for %s", args$cancer_type))

  } else {
 
    stop(sprintf("no reference data defined for %s", args$cancer_type))

  }

  list_out <- merge_seurat_objects(filenames_rds, args_normal)
  rna_normal <- list_out$rna

  rna_normal

} # read_rds_normal_samples



# batch_correction_with_normal_samples
batch_correction_with_normal_samples <- function(rna, rna_normal, args) {


  list_out <- list()
  list_out$rna <- rna
  list_out$rna_normal <- rna_normal

} # batch_correction_with_normal_samples



