#
# seurat_utilities.R
# modified by H. Kim
# last modified on 2021, Dec.
# 
#
# reference:
# seurat_v4.0.5/R/utilities.R
#




# AddModuleScore_modified
# origianl: seurat_v4.0.5/R/utilities.R/AddModuleScore()

#' Calculate module scores for feature expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single
#' cell level, subtracted by the aggregated expression of control feature sets.
#' All analyzed features are binned based on averaged expression, and the
#' control features are randomly selected from each bin.
#'
#' @param object Seurat object
#' @param features A list of vectors of features for expression programs; each
#' entry should be a vector of feature names
#' @param pool List of features to check expression levels against, defaults to
#' \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all
#' analyzed features
#' @param ctrl Number of control features selected from the same bin per
#' analyzed feature
#' @param k Use feature clusters returned from DoKMeans
#' @param assay Name of assay to use
#' @param name Name for the expression programs; will append a number to the
#' end for each entry in \code{features} (eg. if \code{features} has three
#' programs, the results will be stored as \code{name1}, \code{name2},
#' \code{name3}, respectively)
#' @param seed Set a random seed. If NULL, seed is not set.
#' @param search Search for symbol synonyms for features in \code{features} that
#' don't match features in \code{object}? Searches the HGNC's gene names
#' database; see \code{\link{UpdateSymbolList}} for more details
#' @param ... Extra parameters passed to \code{\link{UpdateSymbolList}}
#'
#' @return Returns a Seurat object with module scores added to object meta data;
#' each module is stored as \code{name#} for each module program present in
#' \code{features}
#'
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#' @concept utilities
#'
#' @examples
#' \dontrun{
#' data("pbmc_small")
#' cd_features <- list(c(
#'   'CD79B',
#'   'CD79A',
#'   'CD19',
#'   'CD180',
#'   'CD200',
#'   'CD3D',
#'   'CD2',
#'   'CD3E',
#'   'CD7',
#'   'CD8A',
#'   'CD14',
#'   'CD1C',
#'   'CD68',
#'   'CD9',
#'   'CD247'
#' ))
#' pbmc_small <- AddModuleScore(
#'   object = pbmc_small,
#'   features = cd_features,
#'   ctrl = 5,
#'   name = 'CD_Features'
#' )
#' head(x = pbmc_small[])
#' }
#'
AddModuleScore_modified <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  slot = "data", # H. Kim added
  name = 'Cluster',
  seed = 1,
  search = FALSE,
  ...
) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay

  # H. Kim modified
  #assay.data <- GetAssayData(object = object)
  assay.data <- GetAssayData(object = object, slot = slot)
  # end of modification

  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(
      X = features,
      FUN = function(x) {
        missing.features <- setdiff(x = x, y = rownames(x = object))
        if (length(x = missing.features) > 0) {
          warning(
            "The following features are not present in the object: ",
            paste(missing.features, collapse = ", "),
            ifelse(
              test = search,
              yes = ", attempting to find updated synonyms",
              no = ", not searching for symbol synonyms"
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          if (search) {
            tryCatch(
              expr = {
                updated.features <- UpdateSymbolList(symbols = missing.features, ...)
                names(x = updated.features) <- missing.features
                for (miss in names(x = updated.features)) {
                  index <- which(x == miss)
                  x[index] <- updated.features[miss]
                }
              },
              error = function(...) {
                warning(
                  "Could not reach HGNC's gene names database",
                  call. = FALSE,
                  immediate. = TRUE
                )
              }
            )
            missing.features <- setdiff(x = x, y = rownames(x = object))
            if (length(x = missing.features) > 0) {
              warning(
                "The following features are still not present in the object: ",
                paste(missing.features, collapse = ", "),
                call. = FALSE,
                immediate. = TRUE
              )
            }
          }
        }
        return(intersect(x = x, y = rownames(x = object)))
      }
    )
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }

  # begin of modification by H. Kim
  #pool <- pool %||% rownames(x = object)
  pool <- pool %||% rownames(assay.data)
  # end of modification

  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e30, n = nbin, labels = FALSE, right = FALSE)
  #data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)
  ctrl.scores <- matrix(
    data = numeric(length = 1L),
    nrow = length(x = ctrl.use),
    ncol = ncol(x = object)
  )
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
  }
  features.scores <- matrix(
    data = numeric(length = 1L),
    nrow = cluster.length,
    ncol = ncol(x = object)
  )
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use)
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}





# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}




