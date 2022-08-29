library(ArchR)
library(SummarizedExperiment)
library(HDF5Array)
library(Matrix)
library(stringr)

setwd("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/scATAC-seq_Processing/All_ATAC_v2")


files <- list.files(pattern = ".arrow")


for (i in 1:length(files)){
  sumExp <- getMatrixFromArrow(ArrowFile = files[i],useMatrix = "TileMatrix",binarize = T)
  name <- str_remove(files[i],".arrow")

  saveRDS(sumExp,file = paste0("/datastore/nextgenout5/share/labs/francolab/scENDO_scOVAR_Proj/GEO_Submission/Processed_Data/",name,"_scATAC_Tile_Matrix.rds"))
}

