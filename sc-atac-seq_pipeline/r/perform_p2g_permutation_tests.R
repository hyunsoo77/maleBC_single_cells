#
# perform_p2g_permutation_tests.R
# original author: Matt Regner
# modified by H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
# 
#
# usage:
# called by analyze_cancer_specific_p2g.R
#
# reference:
#
# 0. turorials
# https://www.archrproject.com/bookdown
# 1. papers
# https://www.nature.com/articles/s41588-021-00790-6
# 2. codes
# https://github.com/RegnerM2015/scENDO_scOVAR_2020
#




source("./r/archr/archr_v0.9.5_modified/Archr_Peak_Null_Permute.R")



# Permuated p2g links 
####################################################################
#proj.archr <- readRDS("./final_archr_proj_archrGS.rds")
proj.archr <- readRDS(fname_final_archrgs_rds)

cat(sprintf("\n------------------------------------\n"))
cat(sprintf("addPermPeak2GeneLinks loop 100\n"))

store <- as.numeric(0)
for (i in 1:100) {

  # see archr_v0.9.5_modified/Archr_Peak_Null_Permute.R
  proj.null <- addPermPeak2GeneLinks(
    ArchRProj = proj.archr ,
    reducedDims = reducedDims, # default="LSI_ATAC"
    useMatrix = "GeneIntegrationMatrix_ArchR",
    dimsToUse = dimsToUse,
    scaleDims = NULL,
    corCutOff = 0.75,
    k = 100,
    knnIteration = 500,
    overlapCutoff = 0.8,
    maxDist = 250000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    predictionCutoff = 0.5,
    seed = 1,
    threads = max(floor(getArchRThreads()/2), 1),
    verbose = TRUE,
    logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log, useLogs=TRUE)
  ) 

  p2geneDF <- metadata(proj.null@peakSet)$Peak2GeneLinks
  p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
  p2g.df.null <- as.data.frame(p2geneDF)
  p2g.df.null <- p2g.df.null[complete.cases(p2g.df.null),]
  p2g.null.sub <- dplyr::filter(p2g.df.null, RawPVal <= 1e-12)
  store[i] <- nrow(p2g.null.sub)

} # for



#Example histograms
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("addPermPeak2GeneLinks\n"))

proj.null <- addPermPeak2GeneLinks(
  ArchRProj = proj.archr ,
  reducedDims = reducedDims, # default="LSI_ATAC"
  useMatrix = "GeneIntegrationMatrix_ArchR",
  dimsToUse = dimsToUse,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 250000,
  scaleTo = 10^4,
  log2Norm = TRUE,
  predictionCutoff = 0.5,
  seed = 1,
  threads = max(floor(getArchRThreads()/2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log, useLogs=TRUE)
) 

store.prop <- numeric(0)
# All/Peak2GeneLinks/{seATAC-Group-KNN.rds,  seRNA-Group-KNN.rds}
test <- readRDS(sprintf("%s/Peak2GeneLinks/seATAC-Group-KNN.rds", dir_archr))
test <- test@metadata$KNNList@listData
for ( i in 1:length(test)){
  
  test[[i]] <- gsub("\\#.*","",test[[i]])
  num <- max(table(test[[i]]))
  store.prop[i] <- num/100
} # for

saveRDS(store.prop, sprintf("%s/store_knn_proportions-null.rds", dir_rds))


pdf(sprintf("%s/hist_patient_purity_per_cell_aggregate-null.pdf", dir_pdf), width = 5,height = 3.5)
hist(store.prop, main="Distribtion of patient purity per cell aggregate")
dev.off()


p2geneDF <- metadata(proj.null@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2g.df.null <- as.data.frame(p2geneDF)
p2g.df.null <- p2g.df.null[complete.cases(p2g.df.null),]
#saveRDS(p2g.df.null, sprintf("%s/p2g.df.null.rds", dir_rds))




pdf(sprintf("%s/hist_p2g_correlation-null.pdf", dir_pdf), width = 5,height = 3.5)
hist(p2g.df.null$Correlation,col = "lightblue",main = "Histogram of null peak-to-gene correlations",xlab = "Correlation")
dev.off()

pdf(sprintf("%s/hist_p2g_pval-null.pdf", dir_pdf), width = 5,height = 3.5)
hist(p2g.df.null$RawPVal,col="lightblue",main = "Histogram null peak-to-gene p-values",xlab = "p-value")
abline(v=0.01,col = "red")
dev.off()




# Compute eFDR for alpha 1e-12
print(median(store) / nrow(p2g.df.obs[p2g.df.obs$RawPVal <= 1e-12,]))
#saveRDS(store, sprintf("%s/store_null_tests.rds", dir_rds))








