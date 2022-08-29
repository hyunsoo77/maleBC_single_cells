#
# find_clusters_after_doublet_removal.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
#
# output:
#    rna: Seurat obj
#    rna@meta.data 
#      RNA_snn_res.{args$seurat_resolution}
#      CNV.Pos: 11=CNV, (01,10,00)=no CNV
#      Total_CNVs: sum of CNV values
#
# called by make_sc-rna-seq_seurat_obj.R
#
# reference:
#
#

# find_clusters_after_doublet_removal
find_clusters_after_doublet_removal <- function(rna, args) {





  cancer_type <- args$cancer_type
  dir_seurat_obj <- args$dir_seurat_obj

  #str_column_of_meta_data_cluster <- "RNA_snn_res.0.7"
  str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)

  #SAMPLE.ID = "ovar_3E5CFL"
  sample_id <- sprintf("%s_%s", cancer_type, args$sample_id)

  
  cat(sprintf("\tRunUMAP\n"))
  # https://satijalab.org/seurat/reference/runumap
  rna <- RunUMAP(rna,
                reduction = "pca",
                dims = args$dimsToUse,
                n.neighbors = args$umap_n_neighbors, # default=30L,
                metric = args$umap_metric, # default="cosine",
                min.dist = args$umap_min_dist, # default=0.3
		seed.use = 42L,
		verbose = FALSE
        )

  cat(sprintf("\tFindNeighbors\n"))
  # https://satijalab.org/seurat/reference/findneighbors
  rna <- FindNeighbors(rna,
		 reduction = "pca",
		 dims = args$dimsToUse,
		 verbose=F)

  cat(sprintf("\tFindClusters\n"))
  # https://satijalab.org/seurat/reference/findclusters
  # Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
  rna <- FindClusters(rna,
		graph.name = NULL,
		modularity.fxn = 1,
		initial.membership = NULL,
		node.sizes = NULL,
		resolution = args$seurat_resolution,
		method = "matrix", # Method for running leiden (defaults to matrix which is fast for small datasets). Enable method = "igraph" to avoid casting large data to a dense matrix.
		algorithm = 1, # Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
		n.start = 10,
		n.iter = 10,
		random.seed = args$seed.seurat_clustering, # default=0
		group.singletons = TRUE,
		temp.file.location = NULL,
		edge.file.name = NULL,
		verbose=F)
  # FindClusters() Returns a Seurat object where the idents have been updated with new cluster info; latest clustering results will be stored in object metadata under 'seurat_clusters'. Note that 'seurat_clusters' will be overwritten everytime FindClusters is run.  https://satijalab.org/seurat/reference/findclusters
  # The clustering results are stored in meta.data under ASSAY_snn_res.X  https://satijalab.org/seurat/reference/findclusters

  #Idents(rna) <- "RNA_snn_res.0.7"
  Idents(rna) <- str_column_of_meta_data_cluster
  
  # identify cell types
  rna <- identify_cell_types_with_multiple_methods(rna, args)


  if ((args$method_to_identify_reference_clusters == "none") || (args$type_infercnv_argset == "none")) {

	rna$CNV.value <- NA
	rna$CNV.corr <- NA

	# consider all epithelial cells as CNV+
	if (grepl("normal", cancer_type)) {
		rna$CNV.Pos <- "00"
		rna$CNV.type <- "normal"
	} else {
		rna$CNV.Pos <- "11"
		rna$CNV.type <- "cancer"
	}

  } else {

	# identify reference clusters
	reference.clusters <- identify_reference_clusters(rna, args)
	rna <- update_rna_idents_with_reference_clusters(rna, reference.clusters, "postdoublet.idents", args)
      
	if (args$f_save_infercnv_rds) {
		rna_for_saverds <- get_diet_seurat_obj(rna, args)
		path_rds <- sprintf("%s/%s_rna_postdoublet_preinfercnv.rds", dir_seurat_obj, sample_id)
		cat(sprintf("\tsaveRDS(rna, '%s')\n", path_rds))
		saveRDS(rna_for_saverds, path_rds)
	}
  
	# perform infercnv operations to reveal cnv signal
	cat(sprintf("------------------------------------\n"))
	cat(sprintf("\tPerform infercnv operations to reveal CNV signal\n\n"))


	infercnv_objs <- list()
	dirs_output <- c()
	switch(args$method_to_identify_reference_clusters,
		"normal_sample_epi_types"={
			args$method_to_identify_reference_clusters <- "normal_sample_epi"
			epi_types_as_infercnv_ref <- c("basal", "luminal")
			for (epi_type_as_infercnv_ref in epi_types_as_infercnv_ref) {
  				cat(sprintf("------------------------------------\n"))
  				cat(sprintf("\tPerform infercnv with %s\n\n", epi_type_as_infercnv_ref))
				args$epi_type_as_infercnv_ref <- epi_type_as_infercnv_ref
				infercnv_obj <- get_infercnv_obj(rna, reference.clusters, "postdoublet.idents", args)
				if (!is.null(infercnv_obj)) {
					dir_output <- sprintf("%s/infercnv/%s_cnv_postdoublet_%s", args$dir_output, sample_id, epi_type_as_infercnv_ref)
					cat(sprintf("\tdir_output: %s\n\n", dir_output))
					infercnv_obj <- run_infercnv(infercnv_obj, dir_output, args)
				
					infercnv_objs[[epi_type_as_infercnv_ref]] <- infercnv_obj
					dirs_output <- c(dirs_output, dir_output)
				}
			} # for

		},
		{
			infercnv_obj <- get_infercnv_obj(rna, reference.clusters, "postdoublet.idents", args)
			if (!is.null(infercnv_obj)) {
				dir_output <- sprintf("%s/infercnv/%s_cnv_postdoublet", args$dir_output, sample_id)
				cat(sprintf("\tdir_output: %s\n\n", dir_output))
				infercnv_obj <- run_infercnv(infercnv_obj, dir_output, args)
				infercnv_objs <- infercnv_obj
				dirs_output <- dir_output
			}
		}
	) # switch

	rna <- update_seurat_obj_with_infercnv(rna, infercnv_objs, "postdoublet.idents", dirs_output, args)

	if (grepl("low_total_cnv", args$method_to_identify_reference_clusters)) {

		cat(sprintf("------------------------------------\n"))
		cat(sprintf("\tPerform infercnv with low Total_CNVs\n\n"))

		args$method_to_identify_reference_clusters <- "low_total_cnv"

		#out_unlink <- unlink(dir_output, recursive = TRUE, force = FALSE, expand = TRUE)	
		infercnv_obj <- get_infercnv_obj(rna, reference.clusters, "postdoublet.idents", args)
		if (!is.null(infercnv_obj)) {
			dir_output <- sprintf("%s/infercnv/%s_cnv", args$dir_output, sample_id)
			cat(sprintf("\tdir_output: %s\n\n", dir_output))

			infercnv_obj <- run_infercnv(infercnv_obj, dir_output, args)
			rna <- update_seurat_obj_with_infercnv(rna, infercnv_obj, "postdoublet.idents", dir_output, args)
		}

	}
  
	if (args$f_save_infercnv_rds) {

		rna_for_saverds <- get_diet_seurat_obj(rna, args)
		path_rds <- sprintf("%s/%s_rna_postdoublet_postinferecnv.rds", dir_seurat_obj, sample_id)
		cat(sprintf("\tsaveRDS(rna, '%s')\n", path_rds))
		saveRDS(rna_for_saverds, path_rds)

	}

  } # if


  rna

} # function


