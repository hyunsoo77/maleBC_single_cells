#
# seurat_reproducible.R
# modified by H. Kim
# last modified on 2022, Jul.
#
# content:
# RegressOutMatrix()
# RunPCA.default()
# RunCCA.default()
# 
#
# reference:
# https://github.com/lyc-1995/seurat_reproducible/blob/main/R/modify_seurat.R





# RegressOutMatrix
# comment: no modification so far.
# reference:
# copied from https://github.com/lyc-1995/seurat_reproducible/blob/main/R/modify_seurat.R
RegressOutMatrix <- function(
  data.expr,
  latent.data = NULL,
  features.regress = NULL,
  model.use = NULL,
  use.umi = FALSE,
  verbose = TRUE
) {
  # Do we bypass regression and simply return data.expr?
  bypass <- vapply(
    X = list(latent.data, model.use),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  )
  if (any(bypass)) {
    return(data.expr)
  }
  # Check model.use
  possible.models <- c("linear", "poisson", "negbinom")
  if (!model.use %in% possible.models) {
    stop(paste(
      model.use,
      "is not a valid model. Please use one the following:",
      paste0(possible.models, collapse = ", ")
    ))
  }
  # Check features.regress
  if (is.null(x = features.regress)) {
    features.regress <- 1:nrow(x = data.expr)
  }
  if (is.character(x = features.regress)) {
    features.regress <- intersect(x = features.regress, y = rownames(x = data.expr))
    if (length(x = features.regress) == 0) {
      stop("Cannot use features that are beyond the scope of data.expr")
    }
  } else if (max(features.regress) > nrow(x = data.expr)) {
    stop("Cannot use features that are beyond the scope of data.expr")
  }
  # Check dataset dimensions
  if (nrow(x = latent.data) != ncol(x = data.expr)) {
    stop("Uneven number of cells between latent data and expression data")
  }
  use.umi <- ifelse(test = model.use != 'linear', yes = TRUE, no = use.umi)
  # Create formula for regression
  vars.to.regress <- colnames(x = latent.data)
  fmla <- paste('GENE ~', paste(vars.to.regress, collapse = '+'))
  fmla <- as.formula(object = fmla)
  if (model.use == "linear") {
    # In this code, we'll repeatedly regress different Y against the same X
    # (latent.data) in order to calculate residuals.  Rather that repeatedly
    # call lm to do this, we'll avoid recalculating the QR decomposition for the
    # latent.data matrix each time by reusing it after calculating it once
    # regression.mat <- cbind(latent.data, data.expr[1,])
    # colnames(regression.mat) <- c(colnames(x = latent.data), "GENE")
    # qr <- lm(fmla, data = regression.mat, qr = TRUE)$qr
    # rm(regression.mat)
    fmla <- paste(' ~', paste(vars.to.regress, collapse = '+'))
    fmla <- as.formula(object = fmla)
    mod <- model.matrix(fmla, data = latent.data)
    qr <- qr(mod, tol = 1e-09, LAPACK = TRUE)
    attr(qr, "useLAPACK") <- NULL
  }
  # Make results matrix
  data.resid <- matrix(
    nrow = nrow(x = data.expr),
    ncol = ncol(x = data.expr)
  )
  if (verbose) {
    pb <- txtProgressBar(char = '=', style = 3, file = stderr())
  }
  for (i in 1:length(x = features.regress)) {
    x <- features.regress[i]
    regression.mat <- cbind(latent.data, data.expr[x, ])
    colnames(x = regression.mat) <- c(vars.to.regress, 'GENE')
    regression.mat <- switch(
      EXPR = model.use,
      'linear' = qr.resid(qr = qr, y = data.expr[x,]),
      'poisson' = residuals(object = glm(
        formula = fmla,
        family = 'poisson',
        data = regression.mat),
        type = 'pearson'
      ),
      'negbinom' = NBResiduals(
        fmla = fmla,
        regression.mat = regression.mat,
        gene = x
      )
    )
    data.resid[i, ] <- regression.mat
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / length(x = features.regress))
    }
  }
  if (verbose) {
    close(con = pb)
  }
  if (use.umi) {
    data.resid <- log1p(x = Sweep(
      x = data.resid,
      MARGIN = 1,
      STATS = apply(X = data.resid, MARGIN = 1, FUN = min),
      FUN = '-'
    ))
  }
  dimnames(x = data.resid) <- dimnames(x = data.expr)
  return(data.resid)
} # RegressOutMatrix

environment(RegressOutMatrix) <- asNamespace("Seurat")
assignInNamespace("RegressOutMatrix", RegressOutMatrix, ns = "Seurat")









suppressPackageStartupMessages(library(RSpectra)) # for RSpectra::svds()

# RunPCA.default
# comment: no modification so far.
# reference:
# copied from https://github.com/lyc-1995/seurat_reproducible/blob/main/R/modify_seurat.R
RunPCA.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "PC_",
  seed.use = 42,
  approx = TRUE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  object <- object[order(rownames(x = object)), ]
  gc(reset = TRUE)
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    # pca.results <- irlba(A = object, nv = npcs, ...)
    pca.results <- RSpectra::svds(A = object, k = npcs, ...)
    total.variance <- sum(RowVar(x = t(x = object)))
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    } else{
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  }
  else {
    total.variance <- sum(RowVar(x = object))
    if (approx) {
      npcs <- min(npcs, nrow(x = object) - 1)
      # pca.results <- irlba(A = t(x = object), nv = npcs, ...)
      pca.results <- RSpectra::svds(A = t(x = object), k = npcs, ...)
      feature.loadings <- pca.results$v
      sdev <- pca.results$d/sqrt(max(1, ncol(object) - 1))
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
      } else {
        cell.embeddings <- pca.results$u
      }
    } else {
      npcs <- min(npcs, nrow(x = object))
      pca.results <- prcomp(x = t(object), rank. = npcs, ...)
      feature.loadings <- pca.results$rotation
      sdev <- pca.results$sdev
      if (weight.by.var) {
        cell.embeddings <- pca.results$x
      } else {
        cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = object) - 1))
      }
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key,
    misc = list(total.variance = total.variance)
  )
  if (verbose) {
    msg <- capture.output(print(
      x = reduction.data,
      dims = ndims.print,
      nfeatures = nfeatures.print
    ))
    message(paste(msg, collapse = '\n'))
  }
  return(reduction.data)
} # RunPCA.default

environment(RunPCA.default) <- asNamespace("Seurat")
assignInNamespace("RunPCA.default", RunPCA.default, ns = "Seurat")







# RunCCA.default
RunCCA.default <- function(
  object1,
  object2,
  standardize = TRUE,
  num.cc = 20,
  seed.use = 42,
  verbose = FALSE,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  cells1 <- colnames(x = object1)
  cells2 <- colnames(x = object2)
  if (standardize) {
    object1 <- Standardize(mat = object1, display_progress = FALSE)
    object2 <- Standardize(mat = object2, display_progress = FALSE)
  }
  mat3 <- crossprod(x = object1, y = object2)
  # begin of modification by H. Kim
  #cca.svd <- irlba(A = mat3, nv = num.cc)
  cca.svd <- RSpectra::svds(A = mat3, k = num.cc)
  # end of modification
  cca.data <- rbind(cca.svd$u, cca.svd$v)
  colnames(x = cca.data) <- paste0("CC", 1:num.cc)
  rownames(cca.data) <- c(cells1, cells2)
  cca.data <- apply(
    X = cca.data,
    MARGIN = 2,
    FUN = function(x) {
      if (sign(x[1]) == -1) {
        x <- x * -1
      }
      return(x)
    }
  )

  return(list(ccv = cca.data, d = cca.svd$d))

} # RunCCA.default

environment(RunCCA.default) <- asNamespace("Seurat")
assignInNamespace("RunCCA.default", RunCCA.default, ns = "Seurat")





# RunSPCA.default
RunSPCA.default <- function(
  object,
  assay = NULL,
  npcs = 50,
  reduction.key = "SPC_",
  graph = NULL,
  verbose = FALSE,
  seed.use = 42,
  ...
) {
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  npcs <- min(npcs, nrow(x = object) - 1)
  if (verbose) {
    message("Computing sPCA transformation")
  }
  HSIC <- object %*% graph %*% t(x = object)
  # begin of modification by H. Kim
  #pca.results <- irlba(A = HSIC, nv = npcs)
  pca.results <- RSpectra::svds(A = HSIC, k = npcs)
  # end of modification
  feature.loadings <- pca.results$u
  rownames(x = feature.loadings) <- rownames(x = object)
  cell.embeddings <- t(object) %*% feature.loadings
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings) <-
    paste0(reduction.key, 1:ncol(x = cell.embeddings))
  sdev <- pca.results$d / sqrt(max(1, nrow(x = HSIC) - 1))
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    stdev = sdev,
    key = reduction.key
  )

  return(reduction.data)

} # RunSPCA.default

environment(RunSPCA.default) <- asNamespace("Seurat")
assignInNamespace("RunSPCA.default", RunSPCA.default, ns = "Seurat")







