#
# batch_correction_harmony.R
# author: H. Kim
# date created: 2021, Dec.
# date last modified: 2021, Dec.
#
# contents:
#   run_harmony()
#   save_symphony_reference_rds()
#
#
#
#





# run_harmony
#
# input:
#   args$cancer_type: (e.g. "normal-bc")
#   args$dimsToUse: (e.g. 1:30)
#   args$vars.to.regress:
#   args$seed.harmony:
#   args$harmony_vars: (e.g. "donor")
#   args$harmony_theta:
#   args$harmony_lambda:
#
# output:
#   rna: seurat obj
#     harmony: reduction name
#
# test:
# library(Seurat); library(harmony); source("../r/batch_correction_harmony.R"); args <- list(); args$batch_vars <- "Sample"; args$seed.harmony <- 51; args$harmony_theta <- 0; args$umap_n_neighbors <- 30; args$umap_metric <- "cosine"; min.dist = args$umap_min_dist <- 0.3; args$umap_min_dist <- 0.3; rna1 <- run_harmony(rna, args)
run_harmony <- function(rna, args, f_preprocess_for_assay_integrated=TRUE, f_force_preprocessing=FALSE, vars.to.regress=NULL) {


  if (is.null(vars.to.regress)) {
	vars.to.regress <- args$vars.to.regress
  }

  cat(sprintf("\tharmony_vars=%s\n", paste(args$harmony_vars, collapse=",")))

  vec_batch <- do.call(paste, c(rna@meta.data[,args$harmony_vars,drop=F], sep=","))
  batch_names  <- names(table(vec_batch))
  n_batches <- length(batch_names)
  rna$batch_key <- vec_batch

  if (n_batches < 2) {
	return(rna)
  } # if

  str_assay <- DefaultAssay(rna)
  str_reduction_pca <- "pca"

  if (((str_assay == "integrated") && f_preprocess_for_assay_integrated) || f_force_preprocessing) {

	# assay="integrated", slot="data", "scale.data"
	# make assay="RNA", slot="data", "scale.data" here.

	# https://satijalab.org/seurat/reference/findvariablefeatures
	cat(sprintf("\tFindVariableFeatures for assay=RNA\n"))
        rna <- FindVariableFeatures(object = rna,
		assay = "RNA",
		selection.method = "vst",
		loess.span = 0.3,
		clip.max = "auto",
		num.bin = 20,
		binning.method = "equal_width",
		nfeatures = 2000,
		mean.cutoff = c(0.1, 8),
		dispersion.cutoff = c(1, Inf),
		verbose = FALSE)

	# https://satijalab.org/seurat/reference/scaledata
	if (is.null(vars.to.regress)) {
		cat(sprintf("\tScaleData assay=RNA\n"))
	} else {
		cat(sprintf("\tScaleData assay=RNA with vars.to.regress=%s\n", paste(vars.to.regress, collapse=", ")))
	}
	rna <- ScaleData(rna, features = NULL, assay = "RNA", vars.to.regress = vars.to.regress, verbose = FALSE)

	# https://satijalab.org/seurat/reference/runpca
	cat(sprintf("\tRunPCA for assay=RNA to make pca_rna\n"))

	npcs <- min(50, ncol(rna)-1)

	rna <- RunPCA(rna,
	  assay = "RNA",
	  features = NULL, # default=NULL, Features to compute PCA on. If features=NULL, PCA will be run using the variable features for the Assay. Note that the features must be present in the scaled data. Any requested features that are not scaled or have 0 variance will be dropped, and the PCA will be run using the remaining features.
	  npcs = npcs,
	  rev.pca = FALSE,
	  weight.by.var = TRUE,
	  verbose = FALSE,
	  ndims.print = 1:5, # default=1:5
	  nfeatures.print = 30, # default=30
	  reduction.name = "pca_rna",
	  reduction.key = "pc_", # this should be different with "PC_"
	  seed.use = 42 # default
	)

	str_reduction_pca <- "pca_rna"
	str_assay <- "RNA"

  } # if

  # set.seed for RunHarmony()  https://github.com/immunogenomics/harmony/issues/13
  set.seed(args$seed.harmony)

  # compute Harmony with theta
  # https://github.com/satijalab/seurat-wrappers/blob/master/docs/harmony.md
  # https://github.com/immunogenomics/harmony/blob/master/R/RunHarmony.R
  cat(sprintf("\tRunHarmony assay=%s\n", str_assay))
  cat(sprintf("\t\tharmony_vars=%s\n", paste(args$harmony_vars, collapse=",")))
  cat(sprintf("\t\tharmony_theta=%s\n", paste(args$harmony_theta, collapse=",")))
  cat(sprintf("\t\tharmony_lambda=%s\n", paste(args$harmony_lambda, collapse=",")))
  rna <- RunHarmony(object = rna,
		group.by.vars = args$harmony_vars,
		reduction = str_reduction_pca, # default="pca"
		dims.use = args$dimsToUse,
		theta = args$harmony_theta, # Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
		lambda = args$harmony_lambda, # Ridge regression penalty parameter. Specify for each variable in group.by.vars. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.
		sigma = 0.1, # Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.
		nclust = NULL, # Number of clusters in model. nclust=1 equivalent to simple linear regression.
		tau = 0, # Protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster.
		block.size = 0.05, # What proportion of cells to update during clustering. Between 0 to 1, default 0.05. Larger values may be faster but less accurate
		max.iter.harmony = 10,
		max.iter.cluster = 20,
		epsilon.cluster = 1e-05,
		epsilon.harmony = 1e-04,
		plot_convergence = FALSE,
		verbose = FALSE,
		reference_values = NULL,
		reduction.save = "harmony",
		assay.use = str_assay,
		project.dim = TRUE
	)

  rna

} # run_harmony








# save_symphony_reference_rds
#
# input:
#   args$seed.symphony
#   args$method_sctransform_vst
#   args$harmony_vars: (e.g. "donor")
#   args$harmony_theta:
#   args:dir_seurat_obj
#
save_symphony_reference_rds <- function(rna, args) {


  f_save <- FALSE

  if (args$seed.symphony > 0) {

	# symphony

	# https://rdrr.io/github/immunogenomics/symphony/src/vignettes/utils_seurat.R
	default_assay <- DefaultAssay(rna)
	switch(default_assay,
		"RNA"={
			ref_exp <- GetAssayData(object = rna, assay = "RNA", slot = "data")
		},
		"SCT"={
			ref_exp <- GetAssayData(object = rna, assay = "SCT", slot = "scale.data")
		},
		"integrated"={
			if (nchar(args$method_sctransform_vst) == 0) {
				ref_exp <- GetAssayData(object = rna, assay = "integrated", slot = "data")
			} else {
				# integration on datasets normalized with SCTransform
				ref_exp <- GetAssayData(object = rna, assay = "integrated", slot = "data")
			}
		},
		{}
	) # switch

	# https://github.com/immunogenomics/symphony/blob/main/vignettes/pbmcs_tutorial.ipynb/
	cat(sprintf("\tsymphony::buildReference\n"))
	cat(sprintf("\t\tharmony_vars=%s\n", paste(args$harmony_vars, collapse=",")))
	cat(sprintf("\t\tharmony_theta=%s\n", paste(args$harmony_theta, collapse=",")))
	symphony_reference <- symphony::buildReference(
            ref_exp,
            rna@meta.data,
            vars = args$harmony_vars,  # variables to integrate over
            K = 100,    # number of Harmony clusters
            verbose = TRUE, 
            do_umap = TRUE, 
            do_normalize = FALSE,
            vargenes_method = 'vst', # method for variable gene selection ('vst' or 'mvp')
            vargenes_groups = NULL, # metadata column specifying groups for variable gene selection (e.g. "donor")
            topn = 1000,  # number of variable genes to choose per group
	    tau = 0,	# tau parameter for harmony
	    theta = args$harmony_theta,	# theta parameter(s) for harmony
            save_uwot_path = NULL,   # file path to save uwot model
            d = 20,      # number of PCs
	    additional_genes = NULL,
            umap_min_dist = 0.1,
	    seed = args$seed.symphony
	) # symphony::buildReference

	# save reference (modify with your desired output path)
	dir.create(sprintf("%s/symphony_reference", args$dir_seurat_obj), showWarnings = FALSE, recursive = TRUE)

	path_rds <- sprintf("%s/symphony_reference/%s_symphony_reference.rds", args$dir_seurat_obj, cancer_type)
	cat(sprintf("\t\tsaveRDS(symphony_reference, '%s')\n", path_rds))
	saveRDS(symphony_reference, path_rds)

	f_save <- TRUE

  } # if


  f_save

} # save_symphony_reference_rds








