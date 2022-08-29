#
# batch_correction_seurat.R
# author: H. Kim
# date created: 2021, Dec.
# date last modified: 2021, Dec.
#
# contents:
# batch_correction_with_seurat_lognormalize()
# batch_correction_with_seurat_sctransform()
# update_seurat_obj_after_batch_correction_seurat()
#







# batch_correction_with_seurat_lognormalize
#
# input:
#   args$batch_vars: (e.g. "tumor.type,donor") note: the separator is ',' when there are multiple batch vars
#   args$batch_keys_for_reference: (e.g. B1/KCF0894,ER/MH0025) note: the separator is ',' when there are multiple references. each reference needs to contain the information of batch_vars with the separator of '/'. when there are multiple sampels matched with a reference kay (e.g. tumor.type/donor), the first match will be selected.
#   args$cancer_type: (e.g. "normal-bc")
#   args$dimsToUse: (e.g. 1:30)
#   args$method_integration: {["cca"], "rpca", "rlsi"}
#     use "rpca" for large sample sizes  https://github.com/satijalab/seurat/issues/3650
#   args$type_integration_anchor_features: {["all"], "2000"}
#   args$k.anchor: default=5
#
batch_correction_with_seurat_lognormalize <- function(rna, args, vars.to.regress=NULL) {


    if (is.null(vars.to.regress)) {
	vars.to.regress <- args$vars.to.regress
    }

    cat(sprintf("\tbatch_vars=%s\n", paste(args$batch_vars, collapse=",")))

    vec_batch <- do.call(paste, c(rna@meta.data[,args$batch_vars,drop=F], sep="/"))
    batch_names  <- names(table(vec_batch))
    n_batches <- length(batch_names)
    rna$batch_key <- vec_batch

    if (n_batches < 2) {
	return(rna)
    } # if

    # https://satijalab.org/seurat/reference/splitobject
    cat(sprintf("\tSplitObject\n"))
    list_batch <- SplitObject(rna, split.by="batch_key")

    #cat(sprintf("\tNormalizeData and FindVariableFeatures\n"))
    cat(sprintf("\tNormalizeData\n"))

    #vec_samples <- rep("unknown", n_batches)
    vec_batch_key <- rep("unknown", n_batches)
    for (i in 1:n_batches) {

	#vec_samples[i] <- unique(list_batch[[i]]$Sample)
	vec_batch_key[i] <- unique(list_batch[[i]]$batch_key)

	# https://satijalab.org/seurat/reference/normalizedata
    	cat(sprintf("\t\tNormalizeData batch=%d\n", i))
	list_batch[[i]] <- NormalizeData(object = list_batch[[i]],
	 	assay = NULL,
		normalization.method = "LogNormalize",
		scale.factor = 10000,
		margin = 1,
		verbose = FALSE)

	# https://satijalab.org/seurat/reference/findvariablefeatures
    	cat(sprintf("\t\tFindVariableFeatures batch=%d\n", i))
	list_batch[[i]] <- FindVariableFeatures(object = list_batch[[i]],
		assay = NULL,
		selection.method = "vst", # default = "vst"
		loess.span = 0.3,
		clip.max = "auto",
		num.bin = 20,
		binning.method = "equal_width",
		nfeatures = 2000, # default = 2000
		mean.cutoff = c(0.1, 8),
		dispersion.cutoff = c(1, Inf),
		verbose = FALSE)

	# cell cycle after normalization
	list_batch[[i]] <- identify_cell_types(list_batch[[i]], "cell_cycle", args, n_log=0)

    } # for

    idx_reference <- NULL
    if (args$method_integration == "rpca") {

	# RPCA
	cat(sprintf("\tIntegration with reciprocal PCA\n"))
	# https://satijalab.org/seurat/articles/integration_rpca.html
	# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features

	# https://satijalab.org/seurat/reference/selectintegrationfeatures
	#features <- SelectIntegrationFeatures(object.list = list_batch, nfeatures = 2000, assay = NULL, verbose = TRUE, fvf.nfeatures = 2000)
	#features <- rownames(rna)
	features <- NULL
	if (is.null(vars.to.regress)) {
		cat(sprintf("\t\tScaleData\n"))
	} else {
		cat(sprintf("\t\tScaleData with vars.to.regress=%s\n", paste(vars.to.regress, collapse=", ")))
	}
	list_batch <- lapply(X = list_batch, FUN = function(x) {
		#x <- ScaleData(x, features = features, verbose = FALSE)
		x <- ScaleData(x, features = features, vars.to.regress = vars.to.regress, verbose = FALSE)
		npcs <- min(50, ncol(x)-1)
		x <- RunPCA(x, features = features, npcs = npcs, verbose = FALSE)
	})
	vec_batch_key_u <- unique(vec_batch_key)
	idx_reference <- match(vec_batch_key_u, vec_batch_key)
	
	if (nchar(args$batch_keys_for_reference) > 0) {
		vec_batch_key_u <- strsplit(args$batch_keys_for_reference, ",")[[1]]
		#idx_reference <- match(vec_batch_key_u, vec_samples)
		idx_reference <- match(vec_batch_key_u, vec_batch_key)
	}
	#cat(sprintf("\t\treference: %s\n", paste(vec_samples[idx_reference], collapse=", ")))
	cat(sprintf("\t\treference: %s\n", paste(vec_batch_key[idx_reference], collapse=", ")))
    } # if

    # https://github.com/cellgeni/batchbench/blob/master/bin/seurat3_method.R
    # prevent small datasets from not having enough neighbors (k) to use when filtering anchors 
    if (any(sapply(list_batch, ncol) < 200)) {
	k_filter <- (min(sapply(list_batch, ncol)))
    } else {
	k_filter <- 200
    }





    # anchor.features, features.to.integrate
    switch(args$type_integration_anchor_features,
	"all"={
		nfeatures <- nrow(rna)
		fvf.nfeatures <- nrow(rna)
		anchor.features <- nrow(rna)
		features.to.integrate = rownames(rna)
	},
	{
		nfeatures <- as.numeric(args$type_integration_anchor_features) # default=2000
		if (is.na(nfeatures)) {
			nfeatures <- nrow(rna)
			fvf.nfeatures <- nrow(rna)
			anchor.features <- nrow(rna)
			features.to.integrate <- rownames(rna)
		} else {
			fvf.nfeatures <- nfeatures
			anchor.features <- nfeatures
			features.to.integrate <- NULL
		}
	}
    ) # switch


    
    # https://satijalab.org/seurat/reference/selectintegrationfeatures
    #cat(sprintf("\tSelectIntegrationFeatures\n"))
    #features <- SelectIntegrationFeatures(object.list = list_batch,
#	nfeatures = nfeatures,
#	assay = NULL,
#	verbose = FALSE,
#	fvf.nfeatures = fvf.nfeatures )

    # https://satijalab.org/seurat/reference/findintegrationanchors
    cat(sprintf("\tFindIntegrationAnchors with %s\n", args$method_integration))



    anchors <- FindIntegrationAnchors(object.list = list_batch,
	assay = NULL,
	reference = idx_reference, # default=NULL, A vector specifying the object/s to be used as a reference during integration. If NULL (default), all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in object.list will be used as references. When using a set of specified references, anchors are first found between each query and each reference. The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference.
	anchor.features = anchor.features, # default=2000 or feature from SelectIntegrationFeatures, Can be either: (1) A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding (2) A vector of features to be used as input to the anchor finding process
	scale = TRUE,
	normalization.method = "LogNormalize", # c("LogNormalize", "SCT")
	sct.clip.range = NULL,
	reduction = args$method_integration, # Dimensional reduction to perform when finding anchors. Can be one of: cca: Canonical correlation analysis, rpca: Reciprocal PCA, rlsi: Reciprocal LSI
	l2.norm = TRUE,
	dims = args$dimsToUse,
	k.anchor = args$k.anchor, # default=5, increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default.  https://satijalab.org/seurat/articles/integration_rpca.html
	k.filter = k_filter, # default=200
	k.score = 30,
	max.features = 200,
	nn.method = "annoy",
	n.trees = 50,
	eps = 0,
	verbose = FALSE )




    # https://satijalab.org/seurat/reference/integratedata
    cat(sprintf("\tIntegrateData\n"))
    rna <- IntegrateData(anchorset = anchors, 
	new.assay.name = "integrated",
	normalization.method = "LogNormalize",
	features = NULL, # Vector of features to use when computing the PCA to determine the weights. Only set if you want a different set from those used in the anchor finding process
	features.to.integrate = features.to.integrate, # default=NULL, Vector of features to integrate. By default, will use the features used in anchor finding.
	dims = args$dimsToUse,
	k.weight = args$k.weight, # default=100, Number of neighbors to consider when weighting anchors
	weight.reduction = NULL,
	sd.weight = 1,
	sample.tree = NULL,
	preserve.order = FALSE,
	eps = 0,
	verbose = FALSE )



    # specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
    DefaultAssay(rna) <- "integrated"

    #remove(list=c("list_batch"))

    rna



} # batch_correction_with_seurat_lognormalize

























# batch_correction_with_seurat_sctransform
#
# input:
#   args$batch_vars: (e.g. "donor")  
#   args$cancer_type: (e.g. "normal-bc")
#   args$dimsToUse: (e.g. 1:30)
#   args$method_integration: {["cca"], "rpca", "rlsi"}
#     use "rpca" for large sample sizes  https://github.com/satijalab/seurat/issues/3650
#   args$type_integration_anchor_features: {["all"], "2000"}
#   args$k.anchor: default=5
#
batch_correction_with_seurat_sctransform <- function(rna, args, vars.to.regress=NULL) {

    if (is.null(vars.to.regress)) {
	vars.to.regress <- args$vars.to.regress
    }

    cat(sprintf("\tbatch_vars=%s\n", paste(args$batch_vars, collapse=",")))

    vec_batch <- do.call(paste, c(rna@meta.data[,args$batch_vars,drop=F], sep=","))
    batch_names  <- names(table(vec_batch))
    n_batches <- length(batch_names)
    rna$batch_key <- vec_batch

    if (n_batches < 2) {
	return(rna)
    } # if

    # https://satijalab.org/seurat/reference/splitobject
    cat(sprintf("\tSplitObject\n"))
    list_batch <- SplitObject(rna, split.by="batch_key")

    cat(sprintf("\tSCTransform\n"))

    vec_batch_key <- rep("unknown", n_batches)
    for (i in 1:n_batches) {

	vec_batch_key[i] <- unique(list_batch[[i]]$batch_key)

	# https://satijalab.org/seurat/reference/sctransform
	if (is.null(vars.to.regress)) {
    		cat(sprintf("\t\tSCTransform batch=%d\n", i))
	} else {
		cat(sprintf("\t\tSCTransform batch=%d with vars.to.regress=%s\n", paste(vars.to.regress, collapse=", ")))
	}
	# SCTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
 	list_batch[[i]] <- SCTransform(list_batch[[i]],
		 method = args$method_sctransform_vst, # default method for sctransform::vst() = "poisson".  https://rdrr.io/cran/sctransform/man/vst.html
		 vars.to.regress = vars.to.regress,
		 verbose = FALSE)

    } # for

    idx_reference <- NULL
    if (args$method_integration == "rpca") {

	# RPCA
	cat(sprintf("\tIntegration with reciprocal PCA\n"))
	# https://satijalab.org/seurat/articles/integration_rpca.html
	# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features

	# https://satijalab.org/seurat/reference/selectintegrationfeatures
	#features <- SelectIntegrationFeatures(object.list = list_batch, nfeatures = 2000, assay = NULL, verbose = TRUE, fvf.nfeatures = 2000)
	#features <- rownames(rna)
	features <- NULL

	list_batch <- lapply(X = list_batch, FUN = function(x) {
		# do not scale data since SCTransform already scaled data.
		#x <- ScaleData(x, features = features, verbose = FALSE)
		#x <- ScaleData(x, features = features, vars.to.regress = vars.to.regress, verbose = FALSE)
		npcs <- min(50, ncol(x)-1)
		x <- RunPCA(x, features = features, npcs = npcs, verbose = FALSE)
	})
	vec_batch_key_u <- unique(vec_batch_key)
	idx_reference <- match(vec_batch_key_u, vec_batch_key)

	if (nchar(args$batch_keys_for_reference) > 0) {
		vec_batch_key_u <- strsplit(args$batch_keys_for_reference, ",")[[1]]
		idx_reference <- match(vec_batch_key_u, vec_batch_key)
	}
	cat(sprintf("\t\treference: %s\n", paste(vec_batch_key[idx_reference], collapse=", ")))

    } # if

    # https://github.com/cellgeni/batchbench/blob/master/bin/seurat3_method.R
    # prevent small datasets from not having enough neighbors (k) to use when filtering anchors 
    if (any(sapply(list_batch, ncol) < 200)) {
	k_filter <- (min(sapply(list_batch, ncol)))
    } else {
	k_filter <- 200
    }




    # anchor.features, features.to.integrate
    switch(args$type_integration_anchor_features,
	"all"={
		nfeatures <- nrow(rna)
		fvf.nfeatures <- nrow(rna)
		anchor.features <- nrow(rna)
		features.to.integrate = rownames(rna)
	},
	{
		nfeatures <- as.numeric(args$type_integration_anchor_features) # default=2000
		if (is.na(nfeatures)) {
			nfeatures <- nrow(rna)
			fvf.nfeatures <- nrow(rna)
			anchor.features <- nrow(rna)
			features.to.integrate <- rownames(rna)
		} else {
			fvf.nfeatures <- nfeatures
			anchor.features <- nfeatures
			features.to.integrate <- NULL
		}
	}
    ) # switch







    # https://satijalab.org/seurat/reference/selectintegrationfeatures
    #cat(sprintf("\tSelectIntegrationFeatures\n"))
    #features <- SelectIntegrationFeatures(object.list = list_batch,
#	nfeatures = nfeatures,
#	assay = NULL,
#	verbose = FALSE,
#	fvf.nfeatures = fvf.nfeatures )






    # https://satijalab.org/seurat/reference/prepsctintegration
    cat(sprintf("\tPrepSCTIntegration\n"))
    list_batch <- PrepSCTIntegration(object.list = list_batch,
	assay = NULL,
	anchor.features = anchor.features, # default=2000, features from SelectIntegrationFeatures() Can be either: (1) A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding (2) A vector of features to be used as input to the anchor finding process
	sct.clip.range = NULL,
	verbose = FALSE )




    # https://satijalab.org/seurat/reference/findintegrationanchors
    cat(sprintf("\tFindIntegrationAnchors with %s\n", args$method_integration))
    anchors <- FindIntegrationAnchors(object.list = list_batch,
	assay = NULL,
	reference = idx_reference, # default=NULL, A vector specifying the object/s to be used as a reference during integration. If NULL (default), all pairwise anchors are found (no reference/s). If not NULL, the corresponding objects in object.list will be used as references. When using a set of specified references, anchors are first found between each query and each reference. The references are then integrated through pairwise integration. Each query is then mapped to the integrated reference.
	anchor.features = anchor.features, # default=2000, features from SelectIntegrationFeatures() Can be either: (1) A numeric value. This will call SelectIntegrationFeatures to select the provided number of features to be used in anchor finding (2) A vector of features to be used as input to the anchor finding process
	scale = TRUE,
	normalization.method = "SCT", # c("LogNormalize", "SCT")
	sct.clip.range = NULL,
	reduction = args$method_integration, # c("cca", "rpca", "rlsi")
	l2.norm = TRUE,
	dims = args$dimsToUse,
	k.anchor = args$k.anchor, # default=5, increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default.  https://satijalab.org/seurat/articles/integration_rpca.html
	k.filter = k_filter, # default=200
	k.score = 30,
	max.features = 200,
	nn.method = "annoy",
	n.trees = 50,
	eps = 0,
	verbose = FALSE )




    # https://satijalab.org/seurat/reference/integratedata
    cat(sprintf("\tIntegrateData\n"))
    rna <- IntegrateData(anchorset = anchors, 
	new.assay.name = "integrated",
	normalization.method = "SCT",
	features = NULL, # Vector of features to use when computing the PCA to determine the weights. Only set if you want a different set from those used in the anchor finding process
	features.to.integrate = features.to.integrate, # default=NULL, Vector of features to integrate. By default, will use the features used in anchor finding.
	dims = args$dimsToUse,
	k.weight = 100,
	weight.reduction = NULL,
	sd.weight = 1,
	sample.tree = NULL,
	preserve.order = FALSE,
	eps = 0,
	verbose = FALSE )





    # specify that we will perform downstream analysis on the corrected data note that the original unmodified data still resides in the 'RNA' assay
    DefaultAssay(rna) <- "integrated"

    #remove(list=c("list_batch"))

    rna

} # batch_correction_with_seurat_sctransform




















# update_seurat_obj_after_batch_correction_seurat
# comment:
# set data element=0 when count=0
update_seurat_obj_after_batch_correction_seurat <- function(rna) {


  if (!"integrated" %in% names(rna)) {
	return(rna)
  } # if

  mtx <- GetAssayData(object = rna, assay = "integrated", slot = "data")
  mtx_counts <- GetAssayData(object = rna, assay = "RNA", slot = "counts")

  mtx <- as(mtx, "dgTMatrix")
  mtx_counts <- as(mtx_counts, "dgTMatrix")

  dt_mtx <- data.table(i=mtx@i, j=mtx@j, x=mtx@x, key=c("i", "j"))
  dt_mtx_counts <- data.table(i=mtx_counts@i, j=mtx_counts@j)

  idx_nz <- dt_mtx[.(dt_mtx_counts[,1], dt_mtx_counts[,2]), which=TRUE]
  f <- !is.na(idx_nz); idx_nz <- idx_nz[f]
  dt_mtx <- dt_mtx[idx_nz,]

  mtx@i <- dt_mtx[["i"]]
  mtx@j <- dt_mtx[["j"]]
  mtx@x <- dt_mtx[["x"]]

  rna <- SetAssayData(object = rna, assay = "integrated", slot = "data", new.data=as(mtx, "dgCMatrix"))
  
  rna

} # update_seurat_obj_after_batch_correction_seurat



