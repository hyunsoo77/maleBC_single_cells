#
# Archr_Peak_RawPval.R
# author: Matt Regner
# created on 2021 Nov.
#
# content:
# addPeak2GeneLinks_RawPval()
#
#
# reference:
#
#







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
#' @param k The number of k-nearest neighbors to use for creating single-cell groups for correlation analyses.
#' @param knnIteration The number of k-nearest neighbor groupings to test for passing the supplied `overlapCutoff`.
#' @param overlapCutoff The maximum allowable overlap between the current group and all previous groups to permit the current
#' group be added to the group list during k-nearest neighbor calculations.
#' @param maxDist The maximum allowable distance in basepairs between two peaks to consider for co-accessibility.
#' @param scaleTo The total insertion counts from the designated group of single cells is summed across all relevant peak regions
#' from the `peakSet` of the `ArchRProject` and normalized to the total depth provided by `scaleTo`.
#' @param log2Norm A boolean value indicating whether to log2 transform the single-cell groups prior to computing co-accessibility correlations.
#' @param predictionCutoff A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
#' @param seed A number to be used as the seed for random number generation required in knn determination. It is recommended
#' to keep track of the seed used so that you can reproduce results downstream.
#' @param threads The number of threads to be used for parallel computing.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
addPeak2GeneLinks_RawPval <- function(
  ArchRProj = NULL,
  reducedDims = "IterativeLSI",
  useMatrix = "GeneIntegrationMatrix",
  dimsToUse = 1:30,
  scaleDims = NULL,
  cellsToUse = NULL,
  corCutOff = 0.75,
  k = 100, 
  knnIteration = 500, 
  overlapCutoff = 0.8, 
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.4,
  seed = 1, 
  threads = max(floor(getArchRThreads() / 2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks")
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  .validInput(input = dimsToUse, name = "dimsToUse", valid = c("numeric", "null"))
  .validInput(input = cellsToUse, name = "cellsToUse", valid = c("character", "null"))
  .validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric", "null"))
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
    DataFrame(
      cellNames = paste0(names(ArrowFiles)[x], "#", h5read(ArrowFiles[x], paste0(useMatrix, "/Info/CellNames"))),
      predictionScore = h5read(ArrowFiles[x], paste0(useMatrix, "/Info/predictionScore"))
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
  keepKnn <- determineOverlapCpp(knnObj, floor(overlapCutoff * k))
  
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
    assays = SimpleList(RNA = groupMatRNA), 
    rowRanges = geneStart
  )
  metadata(seRNA)$KNNList <- knnObj
  .logThis(seRNA, "seRNA", logFile = logFile)
  
  names(peakSet) <- NULL
  
  seATAC <- SummarizedExperiment(
    assays = SimpleList(ATAC = groupMatATAC), 
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
  .logDiffTime(main="Computing Background Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  nullCor <- .getNullCorrelations(seATAC, seRNA, o, 1000)
  
  .logDiffTime(main="Computing Correlations", t1=tstart, verbose=verbose, logFile=logFile)
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), assay(seATAC), assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  o$TStat <- (o$Correlation / sqrt((1-o$Correlation^2)/(ncol(seATAC)-2))) #T-statistic P-value
  o$Pval <- 2*pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  o$EmpPval <- 2*pnorm(-abs(((o$Correlation - mean(nullCor[[2]])) / sd(nullCor[[2]]))))
  o$EmpFDR <- p.adjust(o$EmpPval, method = "fdr")
  
  out <- o[, c("A", "B", "Correlation","Pval", "FDR", "VarAssayA", "VarAssayB","EmpPval","EmpFDR")]
  # begin of modification by H. Kim
  #colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "RawPVal","FDR", "VarQATAC", "VarQRNA","EmpPval","EmpFDR")
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "RawPval","FDR", "VarQATAC", "VarQRNA","EmpPval","EmpFDR")
  # end of modification
  mcols(peakSet) <- NULL
  names(peakSet) <- NULL
  metadata(out)$peakSet <- peakSet
  metadata(out)$geneSet <- geneStart
  
  #Save Group Matrices
  dir.create(file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks"), showWarnings = FALSE)
  outATAC <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seATAC-Group-KNN.rds")
  saveRDS(seATAC, outATAC, compress = FALSE)
  outRNA <- file.path(getOutputDirectory(ArchRProj), "Peak2GeneLinks", "seRNA-Group-KNN.rds")
  saveRDS(seRNA, outRNA, compress = FALSE)
  metadata(out)$seATAC <- outATAC
  metadata(out)$seRNA <- outRNA
  
  metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
  
  .logDiffTime(main="Completed Peak2Gene Correlations!", t1=tstart, verbose=verbose, logFile=logFile)
  .endLogging(logFile = logFile)
  
  ArchRProj
  
} # addPeak2GeneLinks



