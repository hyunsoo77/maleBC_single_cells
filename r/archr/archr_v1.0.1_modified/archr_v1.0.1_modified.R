#
# archr_v1.0.1_modified.R
# author: H. Kim
# date created: 2022, Apr.
# date last modified: 2022, Apr.
#
# requirement:
#   /home/hkim77/github/archr/v1.0.1/R
#   /home/hkim77/github/archr/ArchR.dpeerlab/R
#   
#   global variables:
#     args$archr_package: {["ArchR"], "ArchR.dpeerlab"}
#     dir_tmp: the location of tmp directory
#
# content:
#   IntegrativeAnalysis.R
#	addPeak2GeneLinks()
#	.getNullCorrelations()
#
# usage:
# source("./archr_v1.0.1_modified/archr_v1.0.1_modified.R")
# source("./archr_v0.9.5_modified/filterDoublets_modified.R")
# source("./archr_v1.0.1_modified/GgplotUtils.R")
# 
#




if ("package:ArchR.dpeerlab" %in% (search())) {

  dir_archr.src <- "./r/archr/ArchR.dpeerlab/R"

} else {

  dir_archr.src <- "./r/archr/ArchR-1.0.1/R"

}


# load utiltities for modified functions
source(sprintf("%s/ArrowRead.R", dir_archr.src))
#	.getGroupMatrix
source(sprintf("%s/ArrowUtils.R", dir_archr.src))
#	.getFeatureDF
source(sprintf("%s/ArchRHeatmap.R", dir_archr.src))
source(sprintf("%s/HiddenUtils.R", dir_archr.src))
#	.computeKNN
#	.getQuantiles
#	.safelapply
#       .tempfile
source(sprintf("%s/LoggerUtils.R", dir_archr.src))
source(sprintf("%s/ValidationUtils.R", dir_archr.src))









##### IntegrativeAnalysis.R




#' Add Peak2GeneLinks to an ArchRProject
#' 
#' This function will add peak-to-gene links to a given ArchRProject
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param reducedDims The name of the `reducedDims` object (i.e. "IterativeLSI") to retrieve from the designated `ArchRProject`.
#' @param dimsToUse A vector containing the dimensions from the `reducedDims` object to use in clustering.
#' @param scaleDims A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing
#' the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific
#' biases since it is over-weighting latent PCs. If set to `NULL` this will scale the dimensions based on the value of `scaleDims` when the
#' `reducedDims` were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
#' @param corCutOff A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a
#' correlation to sequencing depth that is greater than the `corCutOff`, it will be excluded from analysis.
#' @param cellsToUse A character vector of cellNames to compute coAccessibility on if desired to run on a subset of the total cells.
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param addEmpiricalPval Add empirical p-values based on randomly correlating peaks and genes not on the same seqname.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addPeak2GeneLinks <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  cellsToUse = NULL,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  addEmpiricalPval = FALSE,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = knnIteration, name = "knnIteration", valid = c("integer"))
  .validInput(input = overlapCutoff, name = "overlapCutoff", valid = c("numeric"))
  .validInput(input = maxDist, name = "maxDist", valid = c("integer"))
  .validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  .validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeak2GeneLinks Input-Parameters", logFile = logFile)

  .logDiffTime(main="Getting Available Matrices", t1=tstart, verbose=verbose, logFile=logFile)
  AvailableMatrices <- getAvailableMatrices(ArchRProj)

  if("PeakMatrix" %ni% AvailableMatrices){
    stop("PeakMatrix not in AvailableMatrices")
  }

  if(useMatrix %ni% AvailableMatrices){
    stop(paste0(useMatrix, " not in AvailableMatrices"))
  }

  ArrowFiles <- getArrowFiles(ArchRProj)

  tstart <- Sys.time()

  dfAll <- .safelapply(seq_along(ArrowFiles), function(x){
    cNx <- paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames")))
    pSx <- tryCatch({
      h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
    }, error = function(e){
      if(getArchRVerbose()) message("No predictionScore found. Continuing without predictionScore!")
      rep(9999999, length(cNx))
    })
    DataFrame(
      cellNames = cNx,
      predictionScore = pSx
    )
  }, threads = threads) %>% Reduce("rbind", .)

  .logDiffTime(
    sprintf("Filtered Low Prediction Score Cells (%s of %s, %s)", 
    sum(dfAll[,2] < predictionCutoff), 
    nrow(dfAll), 
    round(sum(dfAll[,2] < predictionCutoff) / nrow(dfAll), 3)
    ), t1=tstart, verbose=verbose, logFile=logFile)

  keep <- sum(dfAll[,2] >= predictionCutoff) / nrow(dfAll)
  dfAll <- dfAll[which(dfAll[,2] > predictionCutoff),]

  set.seed(seed)

  #Get Peak Set
  peakSet <- getPeakSet(ArchRProj)
  .logThis(peakSet, "peakSet", logFile = logFile)

  #Gene Info
  geneSet <- .getFeatureDF(ArrowFiles, useMatrix, threads = threads)
  geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, width = 1), name = geneSet$name, idx = geneSet$idx)
  .logThis(geneStart, "geneStart", logFile = logFile)

  #Get Reduced Dims
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, corCutOff = corCutOff, dimsToUse = dimsToUse)
  if(!is.null(cellsToUse)){
    rD <- rD[cellsToUse, ,drop=FALSE]
  }

  #Subsample
  idx <- sample(seq_len(nrow(rD)), knnIteration, replace = !nrow(rD) >= knnIteration)

  #KNN Matrix
  .logDiffTime(main="Computing KNN", t1=tstart, verbose=verbose, logFile=logFile)
  knnObj <- .computeKNN(data = rD, query = rD[idx,], k = k)

  #Determin Overlap
  .logDiffTime(main="Identifying Non-Overlapping KNN pairs", t1=tstart, verbose=verbose, logFile=logFile)

  # begin of modification by H. Kim
  #keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  if ("package:ArchR.dpeerlab" %in% (search())) {
    keepKnn <- ArchR.dpeerlab:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  } else {
    keepKnn <- ArchR:::determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  } # if
  # end of modification

  #Keep Above Cutoff
  knnObj <- knnObj[keepKnn==0,]
  .logDiffTime(paste0("Identified ", nrow(knnObj), " Groupings!"), t1=tstart, verbose=verbose, logFile=logFile)

  #Convert To Names List
  knnObj <- lapply(seq_len(nrow(knnObj)), function(x){
    rownames(rD)[knnObj[x, ]]
  }) %>% SimpleList

  #Check Chromosomes
  chri <- gtools::mixedsort(unique(paste0(seqnames(peakSet))))
  chrj <- gtools::mixedsort(unique(paste0(seqnames(geneStart))))
  chrij <- intersect(chri, chrj)

  #Features
  geneDF <- mcols(geneStart)
  peakDF <- mcols(peakSet)
  geneDF$seqnames <- seqnames(geneStart)
  peakDF$seqnames <- seqnames(peakSet)

  #Group Matrix RNA
  .logDiffTime(main="Getting Group RNA Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatRNA <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = geneDF, 
    groupList = knnObj, 
    useMatrix = useMatrix,
    threads = threads,
    verbose = FALSE
  )
  rawMatRNA <- groupMatRNA
  .logThis(groupMatRNA, "groupMatRNA", logFile = logFile)

  #Group Matrix ATAC
  .logDiffTime(main="Getting Group ATAC Matrix", t1=tstart, verbose=verbose, logFile=logFile)
  groupMatATAC <- .getGroupMatrix(
    ArrowFiles = getArrowFiles(ArchRProj), 
    featureDF = peakDF, 
    groupList = knnObj, 
    useMatrix = "PeakMatrix",
    threads = threads,
    verbose = FALSE
  )
  rawMatATAC <- groupMatATAC
  .logThis(groupMatATAC, "groupMatATAC", logFile = logFile)

  .logDiffTime(main="Normalizing Group Matrices", t1=tstart, verbose=verbose, logFile=logFile)

  groupMatRNA <- t(t(groupMatRNA) / colSums(groupMatRNA)) * scaleTo
  groupMatATAC <- t(t(groupMatATAC) / colSums(groupMatATAC)) * scaleTo

  if(log2Norm){
    groupMatRNA  <- log2(groupMatRNA + 1)
    groupMatATAC <- log2(groupMatATAC + 1)    
  }

  names(geneStart) <- NULL

  seRNA <- SummarizedExperiment(
    # begin of modification by H. Kim
    # reduce disk storage
    #assays = SimpleList(RNA = groupMatRNA, RawRNA = rawMatRNA), 
    assays = SimpleList(RNA = groupMatRNA), 
    # end of modification
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)

  names(peakSet) <- NULL

  seATAC <- SummarizedExperiment(
    # begin of modification by H. Kim
    # reduce disk storage
    #assays = SimpleList(ATAC = groupMatATAC, RawATAC = rawMatATAC), 
    assays = SimpleList(ATAC = groupMatATAC), 
    # end of modification
    rowRanges = peakSet
  )
  metadata(seATAC)$KNNList <- knnObj
  .logThis(seATAC, "seATAC", logFile = logFile)

  rm(groupMatRNA, groupMatATAC)
  gc()

  #Overlaps
  .logDiffTime(main="Finding Peak Gene Pairings", t1=tstart, verbose=verbose, logFile=logFile)
  o <- DataFrame(
    findOverlaps(
      .suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
      resize(rowRanges(seATAC), 1, "center"), 
      ignore.strand = TRUE
    )
  )

  #Get Distance from Fixed point A B 
  o$distance <- distance(rowRanges(seRNA)[o[,1]] , rowRanges(seATAC)[o[,2]] )
  colnames(o) <- c("B", "A", "distance")

  #Null Correlations
  if(addEmpiricalPval){
    .logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
    nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
  }

  .logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)

  # begin of modification by H. Kim
  #o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  if ("package:ArchR.dpeerlab" %in% (search())) {
    o$Correlation <- ArchR.dpeerlab:::rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  } else {
  o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  } # if
  # end of modification

  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((pmax(1-o$Correlation^2, 0.00000000000000001, na.rm = TRUE))/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")  
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart

  if(addEmpiricalPval){
    out$EmpPval <- 2*pnorm(-abs(((out$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
    out$EmpFDR <- p.adjust(out$EmpPval, method = "fdr")
  }

  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  .safeSaveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  .safeSaveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA

  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out

  .logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  .endLogging(logFile = logFile)

  ArchRProj

} # addPeak2GeneLinks









# .getNullCorrelations
# comment:
# no modification
.getNullCorrelations <- function(seA, seB, o, n){

  o$seq <- seqnames(seA)[o$A]

  nullCor <- lapply(seq_along(unique(o$seq)), function(i){

    #Get chr from olist
    chri <- unique(o$seq)[i]
    #message(chri, " ", appendLF = FALSE)

    #Randomly get n seA
    id <- which(as.character(seqnames(seA)) != chri)
    if(length(id) > n){
      transAidx <- sample(id, n)
    }else{
      transAidx <- id
    }

    #Calculate Correlations
    grid <- expand.grid(transAidx, unique(o[o$seq==chri,]$B))

    idxA <- unique(grid[,1])
    idxB <- unique(grid[,2])

    seSubA <- seA[idxA]
    seSubB <- seB[idxB]

    grid[,3] <- match(grid[,1], idxA)
    grid[,4] <- match(grid[,2], idxB)

    colnames(grid) <- c("A", "B")
   
    # begin of modification by H. Kim
    #out <- rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    if ("package:ArchR.dpeerlab" %in% (search())) {
	out <- ArchR.dpeerlab:::rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    } else {
	out <- ArchR:::rowCorCpp(grid[,3], grid[,4], assay(seSubA), assay(seSubB))
    } # if
    # end of modification

    out <- na.omit(out)

    return(out)

  }) %>% SimpleList
  #message("")

  summaryDF <- lapply(nullCor, function(x){
    data.frame(mean = mean(x), sd = sd(x), median = median(x), n = length(x))
  }) %>% Reduce("rbind",.)

  return(list(summaryDF, unlist(nullCor)))

} # .getNullCorrelations
















##### HiddenUtils.R



# begin of modification by H. Kim
#.tempfile <- function(pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE){
.tempfile_modified <- function(pattern = "tmp", tmpdir = dir_tmp, fileext = "", addDOC = TRUE){
# end of modification

  dir.create(tmpdir, showWarnings = FALSE)

  if(addDOC){
    doc <- paste0("-Date-", Sys.Date(), "_Time-", gsub(":","-", stringr::str_split(Sys.time(), pattern=" ",simplify=TRUE)[1,2]))
  }else{
    doc <- ""
  }

  tempfile(pattern = paste0(pattern, "-"), tmpdir = tmpdir, fileext = paste0(doc, fileext))

} # .tempfile


# override/replace a function in a package
R.utils::reassignInPackage(".tempfile", pkgName=args$archr_package, .tempfile_modified)



