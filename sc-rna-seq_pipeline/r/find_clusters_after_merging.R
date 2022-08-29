#
# find_clusters_after_merging.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
#
# description:
#
# (1) without batch effect correction
# assay="RNA", slot="data"
# reduction="pca" --> umap_reduction="umap"
#
# (2) batch effect correction with harmony
# assay="RNA", slot="data"
# reduction="harmony" --> umap_reduction="umapharmonoy"
#
# (3) batch effect correction with seurat integration
# assay="integrated", slot="data"
# 	reduction="pca" --> umap_reduction="umap"
# assay="RNA", slot="data"
# 	reduction="harmony" --> umap_reduction="umapharmonoy"
#
# output:
#    rna: Seurat obj
#    rna@meta.data 
#      RNA_snn_res.{args$seurat_resolution}
#      RNA_harmony_th.{paste(args$harmony_theta, collapse="_")}
#
# reference:
#













# find_clusters_after_merging
#
# input:
#   args$dimsToUse
#
# comment:
# called by make_sc-rna-seq_merged_seurat_obj.R
find_clusters_after_merging <- function(rna, args) {




  cancer_type <- args$cancer_type
  dir_seurat_obj <- args$dir_seurat_obj

  pattern_seurat_batch_effect_correction <- "cca|rpca|rlsi"
  #str_column_of_meta_data_cluster <- "RNA_snn_res.0.7"
  str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)
  str_column_of_meta_data_harmony <- sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse="_"))
  if (args$f_multik) {
	str_column_of_meta_data_cluster_multik <- "RNA_multik"
	str_column_of_meta_data_harmony_multik <- sprintf("RNA_harmony_th.%s_multik", paste(args$harmony_theta, collapse="_"))
  }

  sample_id <- sprintf("%s_merged", cancer_type)
  

  str_assay_harmony <- DefaultAssay(rna)
  if (grepl(pattern_seurat_batch_effect_correction, args$method_batch_effect_correction)) {
	str_assay_harmony <- "RNA"
  }


  # https://satijalab.org/seurat/reference/runumap
  cat(sprintf("\tRunUMAP reduction=pca to make embedding=umap\n"))
  rna <- RunUMAP(rna,
                dims = args$dimsToUse,
                reduction = "pca",
		features = NULL,
		graph = NULL,
		assay = DefaultAssay(rna),
		nn.name = NULL,
		slot = "data",
		umap.method = "uwot",
		reduction.model = NULL,
		return.model = FALSE,
                n.neighbors = args$umap_n_neighbors, # default=30L
		n.components = 2L,
                metric = args$umap_metric, # default="cosine"
		n.epochs = NULL,
		learning.rate = 1,
                min.dist = args$umap_min_dist, # default=0.3
		spread = 1,
		set.op.mix.ratio = 1,
		local.connectivity = 1L,
		repulsion.strength = 1,
		negative.sample.rate = 5L,
		a = NULL,
		b = NULL,
		uwot.sgd = FALSE,
		seed.use = 42L,
		metric.kwds = NULL,
		angular.rp.forest = FALSE,
		densmap = FALSE,
		dens.lambda = 2,
		dens.frac = 0.3,
		dens.var.shift = 0.1,
		verbose = FALSE,
                reduction.name = "umap",
                reduction.key = "UMAP_"
	)


  cat(sprintf("\tFindNeighbors\n"))
  rna <- FindNeighbors(rna,
		reduction = "pca",
		dims = args$dimsToUse,
		verbose=F)


  cat(sprintf("\tFindClusters\n"))
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

  rna@meta.data[,str_column_of_meta_data_cluster] <- rna$seurat_clusters 




  # multik
  if (args$f_multik) {

		list_multik <- run_multik(rna, args, filename_rds=sprintf("%s_multik_40.rds", cancer_type))
		if (is.null(list_multik$mtx_clusters)) {
			Idents(rna) <- str_column_of_meta_data_cluster
		} else {
			rna@meta.data[,str_column_of_meta_data_cluster_multik] <- list_multik$mtx_clusters[,1]
			Idents(rna) <- str_column_of_meta_data_cluster_multik
		}
  } else {
		#Idents(rna) <- "RNA_snn_res.0.7"
		Idents(rna) <- str_column_of_meta_data_cluster
  } # if











  if (grepl("harmony", args$method_batch_effect_correction) && ("harmony" %in% names(rna))) {

		# reduction="harmony" was already created by make_sc-rna-seq_merged_seurat_obj.R
        	# build umapharmony with reduction="harmony"
    		cat(sprintf("\tRunUMAP reduction=harmony to make embedding=umapharmony\n"))
        	rna <- RunUMAP(rna,
			reduction = "harmony",
			dims = args$dimsToUse,
			n.neighbors = args$umap_n_neighbors, # default=30L
			metric = args$umap_metric, # default="cosine"
			min.dist = args$umap_min_dist, # default=0.3
			verbose = FALSE,
			reduction.name = "umapharmony",
			reduction.key = "umapharmony_" # Keys should be one or more alphanumeric characters followed by an underscore, setting key from umap_harmony_ to umapharmony_
		)

		# batch corrected with seurat integratedata()
		# clusters were defined by assay="integrated", slot="data"
		# umapharmony was genreated by run_harmony()

		# add clusters with assay="RNA", slot="data", reduction="harmony"
		rna.tmp <- rna
		DefaultAssay(rna.tmp) <- "RNA"

      		cat(sprintf("\t\tFindNeighbors assay=RNA\n"))
		rna.tmp <- FindNeighbors(rna.tmp, reduction = "harmony", dims = args$dimsToUse, verbose=F)
      		cat(sprintf("\t\tFindClusters assay=RNA\n"))
		rna.tmp <- FindClusters(rna.tmp, resolution = args$seurat_resolution, verbose=F)
		rna@meta.data[,str_column_of_meta_data_harmony] <- rna.tmp$seurat_clusters

		# multik
		if (args$f_multik) {
			list_multik <- run_multik(rna.tmp, args, f_harmony=T, filename_rds=sprintf("%s_multik_1h.rds", cancer_type))
			if (!is.null(list_multik$mtx_clusters)) {
				rna@meta.data[,str_column_of_meta_data_harmony_multik] <- list_multik$mtx_clusters[,1]
			}
		}
		remove(list = c("rna.tmp"))

  } # if

  rna

} # function





