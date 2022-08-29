#
# find_clusters_before_doublet_removal.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
#
# output:
#    rna: Seurat obj
#    rna@meta.data 
#      RNA_snn_res.{args$seurat_resolution}
#
# called by make_sc-rna-seq_seurat_obj.R
#
# reference:
#
  
# find_clusters_before_doublet_removal
find_clusters_before_doublet_removal <- function(rna, args) {
  


  cancer_type <- args$cancer_type
  dir_seurat_obj <- args$dir_seurat_obj

  #str_column_of_meta_data_cluster <- "RNA_snn_res.0.7"
  str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)

  #SAMPLE.ID = "ovar_3E5CFL"
  sample_id <- sprintf("%s_%s", args$cancer_type, args$sample_id)

  
  cat(sprintf("\tRunUMAP\n"))
  rna <- RunUMAP(rna,
                reduction = "pca",
                dims = args$dimsToUse,
                n.neighbors = args$umap_n_neighbors, # default=30L,
                metric = args$umap_metric, # default="cosine",
                min.dist = args$umap_min_dist, # default=0.3
		verbose = FALSE
        )

  cat(sprintf("\tFindNeighbors\n"))
  rna <- FindNeighbors(rna, reduction = "pca", dims = args$dimsToUse, verbose=F)

  cat(sprintf("\tFindClusters\n"))
  rna <- FindClusters(rna, resolution = args$seurat_resolution, verbose=F)
  # FindClusters() Returns a Seurat object where the idents have been updated with new cluster info; latest clustering results will be stored in object metadata under 'seurat_clusters'. Note that 'seurat_clusters' will be overwritten everytime FindClusters is run.  https://satijalab.org/seurat/reference/findclusters
  # The clustering results are stored in meta.data under ASSAY_snn_res.X  https://satijalab.org/seurat/reference/findclusters

  #Idents(rna) <- "RNA_snn_res.0.7"
  Idents(rna) <- str_column_of_meta_data_cluster
  
  if (args$f_save_infercnv_rds) {
      rna_for_saverds <- get_diet_seurat_obj(rna, args)
      path_rds <- sprintf("%s/%s_rna_predoublet.rds", dir_seurat_obj, sample_id)
      cat(sprintf("\tsaveRDS(rna, '%s')\n", path_rds))
      saveRDS(rna_for_saverds, path_rds)
  }

  rna


} # function




