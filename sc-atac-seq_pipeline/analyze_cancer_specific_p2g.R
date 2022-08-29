#!/usr/bin/env Rscript
#
# analyze_cancer_specific_p2g.R
# analyze cancer specific p2g links
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Apr.
#
# requirement:
#   ArchR v1.0.1
#	or modified packages of ArchR v1.0.1
#		1) ArchR.dpeerlab
#
# options:
#   --n_cores: (default=20) the number of cores for parallel computing
#   --dir_output: (default="./output_p2g")
#   --dir_log: (default="./output_p2g/log")
#
#   --dir_seurat_obj: (default="")
#   --max_dimstouse: (default=30)
#
#   --dir_archr_output: (default="./output/archr_output")
#   --archr_package: {["ArchR"], "ArchR.dpeerlab"}
#
#   --harmony_theta: theta for computing harmony
#
#   --permutation_tests: perform permutation tests (obsoleted)
#   --permutation_tests_firstquart: perfrom permutation tests for firstquart (obsoleted)
#   --seed_kmeans: (default=1) seed for kmeans clustering in plotPeak2GeneHeatmap.distal()
#   --make_tfsee_input: make tfsee input (obsoleted)
#
#   --exclude_encode_peaks: (default=FALSE)
#   --exclude_cluster_epi_unassigned: exclude cluster of epi. unassigned
#   --exclude_cluster_highly_overlapped_with_normal_peaks: (default=FALSE)
#   --max_ratio_overlapped_with_normal_peaks: (default=0.6)
#   --exclude_cluster_low_nfrags: (default=FALSE)
#   --subset_archrproject_force_to_update: force to update subset archrproject dir
#
# input:
#   cancer_type: (e.g. ovar)
#
# output:
#   ./out_cancer_specific_p2g
#       all_p2g_observed.rds
#
# usage:
# ./analyze_cancer_specific_p2g.R male-bc
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

suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser(description="hs script",python_cmd="python")
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--n_log", default=1, type="integer", help="n_log")
parser$add_argument("--n_debug", default=0, type="integer", help="n_debug")
parser$add_argument("-l","--n_log_level", default=0, type="integer", help="n_log_level")
parser$add_argument("-nc", "--n_cores", default=20, type="double", help="the number of cores for parallel computing")

# arguments for input/output
parser$add_argument("-do", "--dir_output", default="./output_p2g", type="character", help="direcotry for output")
parser$add_argument("-dl", "--dir_log", default="./output_p2g/log", type="character", help="direcotry for log")
#parser$add_argument("-ff", "--figure_format", default="pdf", type="character", help="figure format") # coded to print figures with pdf format (pdf and png whenever it is possible since png is helpful for lighter PPT slides).


# arguments for seurat
parser$add_argument("-dr", "--dir_seurat_obj", default="", type="character", help="direcotry for seurat objects")
parser$add_argument("-md", "--max_dimstouse", default=30, type="integer", help="max dimsToUse")



# arguments for archr
parser$add_argument("-da", "--dir_archr_output", default="./output/archr_output", type="character", help="ArchR output directory")
parser$add_argument("-ar", "--archr_package", default="ArchR", type="character", help="ArchR package name")




# arguments for harmony
parser$add_argument("-ht", "--harmony_theta", default="-1", type="character", help="theta for computing harmony")

# arguments for permutation tests
parser$add_argument("-p2g", "--addpeak2genelinks_rawpval", dest="f_addpeak2genelinks_rawpval", action="store_true", default=FALSE, help="use addPeakGgeneLinks_RawPval")
parser$add_argument("-pt", "--permutation_tests", dest="f_permutation_tests", action="store_true", default=FALSE, help="perfrom p2g permutation tests")
parser$add_argument("-ptfq", "--permutation_tests_firstquart", dest="f_permutation_tests_firstquart", action="store_true", default=FALSE, help="perfrom p2g permutation tests for the first quartile")
parser$add_argument("-p2gks", "--seed_kmeans", default=1, type="double", help="seed number for plotPeak2GeneHeatmap.distal(seed_kmeans)")



# arguments for cancer specific
# figure1. cnv plot for sc-rna-seq clusters, select epithelial clusters of which median total cnvs is larger than 0.
parser$add_argument("-cc", "--atac_clusters_cancer", default="", type="character", help="cancer cluster numbers (e.g. \"0,2,6,10,12\")")
parser$add_argument("-eep", "--exclude_encode_peaks", dest="f_exclude_encode_peaks", action="store_true", default=FALSE, help="exclude encode peaks")

parser$add_argument("-ect", "--exclude_cluster_epi_unassigned", dest="f_exclude_cluster_epi_unassigned", action="store_true", default=FALSE, help="exclude cluster of epi. unassigned")
parser$add_argument("-echonp", "--exclude_cluster_highly_overlapped_with_normal_peaks", dest="f_exclude_cluster_highly_overlapped_with_normal_peaks", action="store_true", default=FALSE, help="exclude cluster highly overlapped with normal_peaks")
parser$add_argument("-mronp", "--max_ratio_overlapped_with_normal_peaks", default=0.6, type="double", help="max_ratio_overlapped_with_normal_peaks")

parser$add_argument("-eclnf", "--exclude_cluster_low_nfrags", dest="f_exclude_cluster_low_nfrags", action="store_true", default=FALSE, help="exclude cluster low nfrags")
parser$add_argument("-mnf", "--min_nfrags", default=-1, type="double", help="min number of fragments for a cluster")
parser$add_argument("-saftu", "--subset_archrproject_force_to_update", dest="f_subset_archrproject_force_to_update", action="store_true", default=FALSE, help="force to update subset_archrproject dir")

parser$add_argument("-mti", "--make_tfsee_input", dest="f_make_tfsee_input", action="store_true", default=FALSE, help="make TFSEE input")

# arguments for others
parser$add_argument("-cts","--cancer_type_standard", default="", type="character", help="standard cancer type")


parser$add_argument("cancer_type", nargs=1, help="cancer type")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

dir_output <- args$dir_output
#figure_format <- args$figure_format
dir_seurat_obj <- args$dir_seurat_obj
cancer_type <- args$cancer_type




cat(sprintf("------------------------------------\n"))
#cat(sprintf("./analyze_cancer_specific_p2g.R %s\n", cancer_type))
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("------------------------------------\n\n"))



# create output directories
dir_archr_output <- args$dir_archr_output # output/archr_output
dir_archr <- dirname(dir_archr_output) # output
dir_archr_rds <- sprintf("%s/rds", dir_archr)
if (dirname(args$dir_log) == args$dir_output) {
	dir_log <- args$dir_log
} else {
	dir_log <- sprintf("%s/%s", dir_output, basename(args$dir_log))
}
args$dir_log <- dir_log
dir_pdf <- sprintf("%s/pdf", dir_output)
dir_png <- sprintf("%s/png", dir_output)
dir_rds <- sprintf("%s/rds", dir_output)
args$dir_rds <- dir_rds
dir_tmp <- sprintf("%s/tmp", dir_output)
dir_tsv <- sprintf("%s/tsv", dir_output)
dir_xlsx <- sprintf("%s/xlsx", dir_output)

dir.create(dir_log, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_png, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rds, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tmp, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tsv, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_xlsx, showWarnings = FALSE, recursive = TRUE)




# set seed for reproducibility
set.seed(51)




# for debug: args <- list(); args$n_cores <- 20; dir_log <- "log_cancer_specific_p2g"; cancer_type <- "male-bc"; dimsToUse <- 1:30; reducedDims <- "Harmony"; df_output <- data.frame(); args$atac_clusters_cancer <- "0,2,6,10,12"; filename_archrproj_obj_final <- sprintf("%s/%s_final_archr_proj_archrgs.rds", dir_archr_rds, cancer_type);

suppressPackageStartupMessages(library(ChIPpeakAnno))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(stats))
suppressPackageStartupMessages(library(stringr))

switch(args$archr_package,
	"ArchR"={
		suppressPackageStartupMessages(library(ArchR))
	},
	"ArchR.dpeerlab"={
		suppressPackageStartupMessages(library(ArchR.dpeerlab))
	},
	{
		stop(sprintf("unsupported package: %s", args$archr_package))
	}
) # switch

# for parallel computing
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(parallel))

plan("multicore", workers = 8)
options(future.globals.maxSize = 5 * 1024^3) # 5GB






addArchRThreads(threads = args$n_cores)
addArchRGenome("hg38")




source("./r/archr/archr_v1.0.1_modified/archr_v1.0.1_modified.R")
source("./r/archr/archr_v0.9.5_modified/plotPeak2GeneHeatmap.distal.R")

source("./r/jupyter_message.R")
source("./r/utilities_for_sc_analyses.R")
source("./r/utilities_for_sc-atac-seq.R")





dimsToUse <- 1:args$max_dimstouse



atac_clusters_cancer <- c()
if (nchar(args$atac_clusters_cancer) > 0) {

	# atac_clusters_cancer.446B7L <- c(3,4,16)
	# atac_clusters_cancer.4CC61L <- c(0,2) # removed 19 since median of cnvs is 0
	# atac_clusters_cancer <- c(atac_clusters_cancer.446B7L, atac_clusters_cancer.4CC61L)

	atac_clusters_cancer <- strsplit(args$atac_clusters_cancer, ",")[[1]]
} # if
















cat(sprintf("\n------------------------------------\n"))
cat(sprintf("parameters\n"))

cat(sprintf("\tn_cores=%d\n", args$n_cores))
cat(sprintf("\tdir_output=%s\n", args$dir_output))
cat(sprintf("\tdir_log=%s\n", args$dir_log))
cat(sprintf("\tdir_seurat_obj=%s\n", args$dir_seurat_obj))
cat(sprintf("\tdir_archr_output=%s\n", args$dir_archr_output))
cat(sprintf("\tdir_rds=%s\n", args$dir_rds))
cat(sprintf("\tdir_tmp=%s\n", dir_tmp))
cat(sprintf("\tdir_xlsx=%s\n", dir_xlsx))
cat(sprintf("\tcancer_type=%s\n", args$cancer_type))

args <- update_args_cancer_type_standard(args)
cat(sprintf("\tcancer_type_standard=%s\n", args$cancer_type_standard))



# for seurat
cat(sprintf("\tmax_dimstouse=%d\n", args$max_dimstouse))

# for archr
cat(sprintf("\tarchr_package=%s\n", args$archr_package))

# arguments for cancer specific
cat(sprintf("\tclusters_cancer=%s\n", paste(atac_clusters_cancer, collapse=", ")))
cat(sprintf("\tf_exclude_encode_peaks=%s\n", args$f_exclude_encode_peaks))
cat(sprintf("\tf_exclude_cluster_highly_overlapped_with_normal_peaks=%s\n", args$f_exclude_cluster_highly_overlapped_with_normal_peaks))
cat(sprintf("\tf_exclude_cluster_low_nfrags=%s\n", args$f_exclude_cluster_low_nfrags))



cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("read peaks\n\n"))


gr_normal.peaks <- load_normal_peaks(args)


if (args$f_exclude_encode_peaks) {

	# read Encode peaks
	gr_encode <- load_encode_peaks(args)

} # if












###################################################################
# read archr proj obj
filename_archrproj_obj_final <- sprintf("%s/%s_archrproj_obj_final.rds", dir_archr_rds, cancer_type)
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("read %s\n", filename_archrproj_obj_final))
proj.archr <- readRDS(filename_archrproj_obj_final)




args$f_exclude_cluster_types <-  (args$f_exclude_cluster_epi_unassigned || args$f_exclude_cluster_highly_overlapped_with_normal_peaks || args$f_exclude_cluster_low_nfrags)

if (args$f_exclude_cluster_types) {

	# cluster.type_exclude
	cluster.type_exclude <- determine_cluster.type_exclude_with_archr_peakset(proj.archr, args)


	# subset ArchRProject
	proj.archr <- subset_archrproject(proj.archr, args, cluster.type_exclude)


} # if










# update reducedDims
df_output <- data.frame()

col_cluster_types <- "cluster.type"
reducedDims <- "LSI_ATAC"
f_multiome <- FALSE
if (args$harmony_theta[1] >= 0) {

  # harmony
  col_cluster_types <- "cluster.type.harmony" 
  reducedDims <- "Harmony"

} else if ("LSI_Combined" %in% names(proj.archr)) {

	# use LSI_Combined without harmony
	reducedDims <- "LSI_Combined"
	f_multiome <- TRUE

} else {

	# same as default
	reducedDims <- "LSI_ATAC"

} # if




cat(sprintf("\tcol_cluster_types=%s\n", col_cluster_types))
cat(sprintf("\treducedDims=%s\n", reducedDims))
df_output["reducedDims", "info"] <- reducedDims





if (length(atac_clusters_cancer) == 0) {

	cluster_types <- unique(proj.archr$predictedGroup_ArchR)
	f.epi <- grepl(pattern_tumor_epi, cluster_types)
	if (any(f.epi)) {
		cat(sprintf("\tdetermine atac_clusters_cancer with pattern_tumor_epi\n"))
	} else {
		cat(sprintf("\tdetermine atac_clusters_cancer with pattern_epi\n"))
		f.epi <- grepl(pattern_epi, cluster_types)
	} # if
	atac_clusters_cancer <- gsub("-.*$", "", cluster_types[f.epi])
	cat(sprintf("\t\tclusters_cancer=%s\n", paste(atac_clusters_cancer, collapse=", ")))

} # if

pattern_filename_atac_clusters_cancer <- paste(sprintf("X%s\\.", atac_clusters_cancer), collapse='|') # e.g. pattern_atac_clusters_cancer <- "X0\\.|X1\\.|X3\\."
pattern_atac_clusters_cancer <- paste(sprintf("^%s-", atac_clusters_cancer), collapse='|') # e.g. pattern_atac_clusters_cancer <- "0-|1-|3-"



# number of samples for each group
dt <- summarize_meta_data(proj.archr, "predictedGroup_ArchR", "Sample")




















# read cancer/epi peaks
filenames_epi <- list.files(path=sprintf("%s/PeakCalls", dir_archr_output), pattern="^X.*rds", all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

if (length(atac_clusters_cancer) > 0) {
	cat(sprintf("\tload cancer peaks with pattern_filename_atac_clusters_cancer\n"))
	f.epi <- grepl(pattern_filename_atac_clusters_cancer, filenames_epi)
} else {
	f.epi <- grepl(pattern_tumor_epi, filenames_epi)
	if (any(f.epi)) {
		cat(sprintf("\tload cancer peaks with pattern_tumor_epi\n"))
	} else {
		cat(sprintf("\tload epithelial peaks with pattern_epi\n"))
		f.epi <- grepl(pattern_epi, filenames_epi)
	}
} # if
filenames_epi <- filenames_epi[f.epi]




# peakcalls.epi.peakName
if (length(filenames_epi) > 0) {
	peakcalls.epi.peakName <- c()
	for (filename_epi in filenames_epi) {
		cat(sprintf("\tread %s\n", filename_epi))
		gr_peakcalls.epi <- readRDS(filename_epi)
		peakcalls.epi.peakName <- c(peakcalls.epi.peakName, sprintf("%s:%d-%d", seqnames(gr_peakcalls.epi), start(gr_peakcalls.epi), end(gr_peakcalls.epi)))
	} # for

	cat(sprintf("\t# of peakcalls.epi: %d\n", length(peakcalls.epi.peakName)))
} else {
	stop(sprintf("\tthere is no cancer/epi peaks in %s/PeakCalls.", dir_archr_output))
} # if









###################################################################
# permutation tests

if (args$f_addpeak2genelinks_rawpval) {
  soruce("./r/archr/archr_v0.9.5_modified/Archr_Peak_RawPval.R")
}

if (args$f_permutation_tests || args$f_permutation_tests_firstquart) {

  soruce("./r/archr/archr_v0.9.5_modified/Archr_Peak_RawPval.R")

  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("addPeak2GeneLinks_RawPval\n"))

  # Add p2g links (no restrictions on FDR, Correlation, Variance cutoff) with raw pvalue
  proj.archr <- addPeak2GeneLinks_RawPval(
    ArchRProj = proj.archr ,
    reducedDims = reducedDims, # default="IterativeLSI"
    useMatrix = "GeneIntegrationMatrix_ArchR",
    dimsToUse = dimsToUse,
    scaleDims = NULL,
    corCutOff = 0.75,
    cellsToUse = NULL,
    k = 100,
    knnIteration = 500,
    overlapCutoff = 0.8,
    maxDist = 250000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    predictionCutoff = 0.5, # default=0.4
    seed = 1,
    threads = max(floor(getArchRThreads() / 2), 1),
    verbose = TRUE,
    logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log, useLogs=TRUE)
  ) # addPeak2GeneLinks



  p2geneDF <- metadata(proj.archr@peakSet)$Peak2GeneLinks
  p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]

  peakset <- proj.archr@peakSet
  names(peakset) <- NULL
  df_peaks_from_peakset <- as.data.frame(peakset)
  p2geneDF$peakType <- df_peaks_from_peakset[p2geneDF$idxATAC, "peakType"]
  
  p2g.df.obs <- as.data.frame(p2geneDF)
  p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]

  first.quart <- summary(p2g.df.obs$RawPval)[2]

  # p2g_permutation_tests
  if (args$f_permutation_tests) {

	source("./r/perform_p2g_permutation_tests.R")

  }



  # p2g_permutation_tests-firstquart
  if (args$f_permutation_tests_firstquart) {

	source("./r/perform_p2g_permutation_tests-firstquart.R")

  }



} # if























#####################################################################
# Add p2g links (no restrictions on FDR, Correlation, Variance cutoff) with raw pvalue
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("addPeak2GeneLinks\n"))

if (args$f_addpeak2genelinks_rawpval) {

  proj.archr <- addPeak2GeneLinks_RawPval(
    ArchRProj = proj.archr,
    reducedDims = reducedDims, # default="IterativeLSI"
    useMatrix = "GeneIntegrationMatrix_ArchR",
    dimsToUse = dimsToUse,
    scaleDims = NULL,
    corCutOff = 0.75,
    cellsToUse = NULL,
    k = 100,
    knnIteration = 500,
    overlapCutoff = 0.8,
    maxDist = 250000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    predictionCutoff = 0.5, # default=0.4
    seed = 1,
    threads = max(floor(getArchRThreads() / 2), 1),
    verbose = TRUE,
    logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log, useLogs=TRUE)
  ) # addPeak2GeneLinks


} else {

  proj.archr <- addPeak2GeneLinks(
    ArchRProj = proj.archr,
    reducedDims = reducedDims, # default="IterativeLSI"
    useMatrix = "GeneIntegrationMatrix_ArchR",
    dimsToUse = dimsToUse,
    scaleDims = NULL,
    corCutOff = 0.75,
    cellsToUse = NULL,
    k = 100,
    knnIteration = 500,
    overlapCutoff = 0.8,
    maxDist = 250000,
    scaleTo = 10^4,
    log2Norm = TRUE,
    predictionCutoff = 0.5, # default=0.4
    addEmpiricalPval = TRUE, # added by H. Kim
    seed = 1,
    threads = max(floor(getArchRThreads() / 2), 1),
    verbose = TRUE,
    logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log, useLogs=TRUE)
  ) # addPeak2GeneLinks

} # if





if (args$f_exclude_cluster_types) {

	filename_archproj_obj_rds <- sprintf("%s/%s_archrproj_obj_p2gs.rds", dir_rds, cancer_type)
	saveRDS(proj.archr, filename_archproj_obj_rds)

} # if





store.prop <- numeric(0)
# All/Peak2GeneLinks/{seATAC-Group-KNN.rds,  seRNA-Group-KNN.rds}
test <- readRDS(sprintf("%s/Peak2GeneLinks/seATAC-Group-KNN.rds", dir_archr_output))
test <- test@metadata$KNNList@listData
for ( i in 1:length(test) ) {
  
  test[[i]] <- gsub("\\#.*","",test[[i]])
  num <- max(table(test[[i]]))
  store.prop[i] <- num/100

} # for

#saveRDS(store.prop, sprintf("%s/store_knn_proportions.rds", dir_rds))

pdf(sprintf("%s/hist_patient_purity_per_cell_aggregate.pdf", dir_pdf), width = 5,height = 3.5)
hist(store.prop,main="Distribtion of patient purity per cell aggregate")
dev.off()




# p2geneDF
p2geneDF <- metadata(proj.archr@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
# colnames of p2geneDF: idxATAC, idxRNA Correlation, Pval, FDR, VarQATAC, VarQRNA, EmpPval, EmpFDR, geneName, peakName



# peakset
peakset <- proj.archr@peakSet
peakset.peakName <- sprintf("%s:%d-%d", seqnames(peakset), start(peakset), end(peakset))

vec_cluster_type_peaks <- names(peakset) # e.g. c("0-Epithelial cells'", ...)
print(table(vec_cluster_type_peaks))
# > length(vec_cluster_type_peaks) [1] 175986



# df_peaks_from_peakset
names(peakset) <- NULL
df_peaks_from_peakset <- as.data.frame(peakset)
n_peaks <- nrow(df_peaks_from_peakset)



# update p2geneDF peakType
p2geneDF$peakType <- df_peaks_from_peakset[p2geneDF$idxATAC, "peakType"]

cat(sprintf("\tn_peaks: %d\n", n_peaks))
df_output["n_peaks", "info"] <- n_peaks

p2g.df.obs <- as.data.frame(p2geneDF)
p2g.df.obs <- p2g.df.obs[complete.cases(p2g.df.obs),]

n_p2g.obs <- nrow(p2g.df.obs)
cat(sprintf("\tn_p2g.obs: %d\n", n_p2g.obs))
df_output["n_p2g.obs", "info"] <- n_p2g.obs

n_peaks.p2g.obs <- length(unique(p2g.df.obs$peakName))
cat(sprintf("\tn_peaks.p2g.obs: %d\n", n_peaks.p2g.obs))
df_output["n_peaks.p2g.obs", "info"] <- n_peaks.p2g.obs



# save
saveRDS(p2g.df.obs, sprintf("%s/all_p2g_observed.rds", dir_rds))





pdf(sprintf("%s/hist_p2g_correlation.pdf", dir_pdf), width = 5, height = 3.5)
hist(p2g.df.obs$Correlation, col = "lightblue", main = paste0("Histogram of\n", nrow(p2g.df.obs), " P2G correlations"), xlab = "Correlation")
dev.off()



if (args$f_addpeak2genelinks_rawpval) {
  pdf(sprintf("%s/hist_p2g_pval.pdf", dir_pdf), width = 5, height = 3.5)
  hist(p2g.df.obs$RawPVal, col="lightblue", main = paste0("Histogram of\n", nrow(p2g.df.obs), " P2G p-values"), xlab = "p-value")
  abline(v=0.01, col = "red")
  dev.off()
} else {
  pdf(sprintf("%s/hist_p2g_fdr.pdf", dir_pdf), width = 5, height = 3.5)
  hist(p2g.df.obs$FDR, col="lightblue", main = paste0("Histogram of\n", nrow(p2g.df.obs), " P2G FDR"), xlab = "FDR")
  abline(v=0.01, col = "red")
  dev.off()
} # if








####################################################################
# filter




cat(sprintf("\n------------------------------------\n"))
cat(sprintf("p2g.df.hist\n"))
if (args$f_addpeak2genelinks_rawpval) {
	p2g.df.hist <- dplyr::filter(p2g.df.obs,
		 RawPVal <= 1e-12 &
		 Correlation >= 0.45)
} else {
	p2g.df.hist <- dplyr::filter(p2g.df.obs,
		 FDR <= 1e-12 &
		 Correlation >= 0.45)
} # if

pdf(sprintf("%s/hist_genes_per_peak.pdf", dir_pdf), width = 5,height = 3.5)
hist(table(p2g.df.hist$idxRNA), main="Distribution of genes per peaks")
dev.off()

pdf(sprintf("%s/hist_peaks_per_gene.pdf", dir_pdf), width = 5,height = 3.5)
hist(table(p2g.df.hist$idxATAC), main="Distribution of peaks per gene")
dev.off()














#########################################################
# loop for each peak type

peakTypes <- c("Promoter", "Distal")
#peakTypes <- c("Promoter", "Intronic", "Distal", "Exonic")

for (peakType_ in peakTypes) {

  peaktype <- tolower(peakType_)
  if (peaktype == "distal") peaktype <- "enhancer"

  cat(sprintf("\n\n------------------------------------\n"))
  cat(sprintf("%s\n", peakType_))
  cat(sprintf("------------------------------------\n\n"))

  ### subset to postive correlation P2Gs
  # p2g.df.sub: data.frame for p2g distal links with high correlation
  cat(sprintf("p2g.df.sub\n"))
  if (args$f_addpeak2genelinks_rawpval) {
    p2g.df.sub <- dplyr::filter(p2g.df.obs,
		RawPVal <= 1e-12 & 
		Correlation >= 0.45 &
		peakType == peakType_)
  } else {
    p2g.df.sub <- dplyr::filter(p2g.df.obs,
		FDR <= 1e-12 & 
		Correlation >= 0.45 &
		peakType == peakType_)
  } # if


  n_p2g.sub <- nrow(p2g.df.sub)
  var <- sprintf("n_p2g.%s", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_p2g.sub))
  df_output[var, "info"] <- n_p2g.sub

  n_peaks.p2g.sub <- length(unique(p2g.df.sub$peakName))
  var <- sprintf("n_peaks.p2g.%s", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_peaks.p2g.sub))
  df_output[var, "info"] <- n_peaks.p2g.sub




  ### plot peak2 gene heatmap
  # p2g.df.sub.plot
  p2g.df.sub$idx <- paste0(p2g.df.sub$idxATAC, "-", p2g.df.sub$idxRNA)
  p2g.df.sub.plot <- p2g.df.sub


  # heatmap_peak2gene before excluding epi. unassgined
  cat(sprintf("\tplotPeak2GeneHeatmap.distal with seed_kmeans=%g\n", args$seed_kmeans))
  test <- plotPeak2GeneHeatmap.distal(proj.archr,
		peaks = p2g.df.sub.plot,
		groupBy = "predictedGroup_ArchR",
		k = length(levels(factor(proj.archr$predictedGroup_ArchR))),
		corCutOff = .45,
		varCutOffATAC = 0,
		varCutOffRNA = 0,
		FDRCutOff = 1,
		returnMatrices = FALSE,
		nPlot =100000,
		seed = args$seed_kmeans
	) # plotPeak2GeneHeatmap.distal

  filename_heatmap <- sprintf("%s/heatmap_peak2gene_legend_%s.pdf", dir_pdf, peaktype)
  cat(sprintf("\tfilename_heatmap: %s\n", filename_heatmap))
  pdf(filename_heatmap, width = 8, height = 10)
  draw(test, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()



  # list_atac_rna_p2g
  cat(sprintf("\tplotPeak2GeneHeatmap.distal for matrices\n"))
  list_atac_rna_p2g <- plotPeak2GeneHeatmap.distal(proj.archr,
			peaks = p2g.df.sub.plot,
			groupBy = "predictedGroup_ArchR",
			k = length(levels(factor(proj.archr$predictedGroup_ArchR))),
			corCutOff = .45,
			varCutOffATAC = 0,
			varCutOffRNA = 0,
			FDRCutOff = 1,
			returnMatrices = TRUE,
			nPlot = 100000,
			seed = args$seed_kmeans
	) # plotPeak2GeneHeatmap.distal


  # Save P2G peaknames and Kmeans cluster for the genes of interest (specific to cancer cells)

  p2g.df.sub.plot$kmeans <- list_atac_rna_p2g$RNA$kmeansId
  saveRDS(p2g.df.sub.plot, sprintf("%s/p2g.df.sub.plot_%s.rds", dir_rds, peaktype))

  pdf(sprintf("%s/hist_genes_per_%s_peak.pdf", dir_pdf, peaktype), width = 5,height = 3.5)
  hist(table(p2g.df.sub.plot$idxRNA),main="Distribution of genes per distal peaks")
  dev.off()

  pdf(sprintf("%s/hist_%s_peaks_per_gene.pdf", dir_pdf, peaktype), width = 5,height = 3.5)
  hist(table(p2g.df.sub.plot$idxATAC),main="Distribution of distal peaks per gene")
  dev.off()






  # automatic selection of store.kemans developed by H. Kim
  cat(sprintf("\n\t------------------------------------\n"))
  cat(sprintf("\tdetermine store.kmeans\n\n"))



  
  kmeans_clusters <- sort(unique(list_atac_rna_p2g$ATAC$kmeansId))

  # determine store.kmeans, cluster.types_with_cancer_cells
  store.kmeans <- c()
  cluster.types_with_cancer_cells <- c()
  for (k in kmeans_clusters) {

    cat(sprintf("\tkmeans cluster=%d\n", k))
    idx_row <- which(list_atac_rna_p2g$ATAC$kmeansId == k)

    vec_score <- colSums(list_atac_rna_p2g$ATAC$matrix[idx_row, ])
    th_score <- quantile(vec_score, probs=c(0.9), na.rm=TRUE)
    idx_col_knngroup <- which(vec_score > th_score)
    
    tb <- table(list_atac_rna_p2g$ATAC$colData[idx_col_knngroup, "groupBy"])
    tb_sort <- sort(tb, decreasing=T)
    vec_cluster_type <- names(tb_sort)

    # f_cancer
    # pattern_atac_clusters_cancer was assigned by atac_clusters_cancer. when length(atac_clusters_cancer)=0, pattern_atac_clusters_cancer was determined by pattern_tumor_epi or pattern_epi (see above).
    f_cancer <- grepl(pattern_atac_clusters_cancer, vec_cluster_type)

    num_epi <- sum(tb_sort[f_cancer])
    num_others <- sum(tb_sort[!f_cancer])
    if (num_epi > num_others) {
	# number of cancer cells is larger than number of other cells
	store.kmeans <- c(store.kmeans, k)  
	cluster.types_with_cancer_cells <- union(cluster.types_with_cancer_cells, vec_cluster_type[f_cancer]) # H. Kim added
	print(tb_sort) 
    }    

  } # for

  cat(sprintf("\tstore.kmeans=%s\n", paste(store.kmeans, collapse=", ")))
  cat(sprintf("\tcluster.types_with_cancer_cells=%s\n", paste(cluster.types_with_cancer_cells, collapse=", ")))








  # p2g.df.sub.plot.cancer.enriched
  if (length(store.kmeans) == 0) {

	method_to_select_peaks <- "peakcalls.epi"
	#method_to_select_peaks <- "peakset.names"
	switch(method_to_select_peaks,
		"peakcalls.epi"={
			cat(sprintf("\tselect peaks with peakcalls.epi\n"))
			p2g.df.sub.plot.cancer.enriched <- p2g.df.sub.plot[p2g.df.sub.plot$peakName %in% peakcalls.epi.peakName,]
		},
		"peakset.names"={
			# select peaks with names(peaks)
			cat(sprintf("\tselect peaks with names(peaks)\n"))
			f.epi <- grepl(pattern_tumor_epi, vec_cluster_type_peaks)
			if (!any(f.epi)) {
				f.epi <- grepl(pattern_epi, vec_cluster_type_peaks)
			}
			p2g.df.sub.plot.cancer.enriched <- p2g.df.sub.plot[p2g.df.sub.plot$peakName %in% peakset.peakName[f.epi],]
		},
		{}
	) # switch

  } else {

	# extract P2Gs for kmeans clusters of interest
	# p2g.df.sub.plot.cancer.enriched
	p2g.df.sub.plot.cancer.enriched <- p2g.df.sub.plot[p2g.df.sub.plot$kmeans %in% store.kmeans,]

  } # if


  n_p2g.cancer.enriched <- nrow(p2g.df.sub.plot.cancer.enriched)
  var <- sprintf("n_p2g.%s.cancer.enriched", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_p2g.cancer.enriched))
  df_output[var, "info"] <- n_p2g.cancer.enriched

  n_peaks.p2g.cancer.enriched <- length(unique(p2g.df.sub.plot.cancer.enriched$peakName))
  var <- sprintf("n_peaks.p2g.%s.cancer.enriched", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_peaks.p2g.cancer.enriched))
  df_output[var, "info"] <- n_peaks.p2g.cancer.enriched


  # save
  saveRDS(p2g.df.sub.plot.cancer.enriched, sprintf("%s/cancer_enriched_%s_p2g_table.rds", dir_rds, peaktype))

  p2g <- GRanges(p2g.df.sub.plot.cancer.enriched$peakName)






















  # ol, widths, tot

  switch(args$cancer_type_standard,

	"bc"={
		# breast cancer
		widths.normal <- end(gr_normal.peaks) - start(gr_normal.peaks)
		if (args$f_exclude_encode_peaks) {
			widths.encode <- end(gr_encode) - start(gr_encode)
			widths <- c(widths.encode, widths.normal)
			ol <- findOverlapsOfPeaks(unique(p2g), gr_encode, unique(gr_normal.peaks), minoverlap = 1, connectedPeaks = "min")
		} else {
			widths <- widths.normal
			ol <- findOverlapsOfPeaks(unique(p2g), unique(gr_normal.peaks), minoverlap = 1, connectedPeaks = "min")
		} # if

	},

	"oc"={

		# ovarian cancer
		widths.2 <- end(gr_normal.peaks) - start(gr_normal.peaks)

		#oe.peaks <- readRDS("Ovarian_Epithelial_Cell_line_Peaks.rds")
		oe.peaks <- readRDS("./reference/normal_cell_lines/Ovarian_Epithelial_Cell_line_Peaks.rds")
		gr_oe.peaks <- GRanges(oe.peaks)
		widths.3 <- end(gr_oe.peaks) - start(gr_oe.peaks)

		if (args$f_exclude_encode_peaks) {
			widths.encode <- end(gr_encode) - start(gr_encode)
			widths <- c(widths.encode, widths.2, widths.3)
			ol <- findOverlapsOfPeaks(unique(p2g), gr_encode, unique(gr_oe.peaks), unique(gr_normal.peaks), minoverlap = 1, connectedPeaks = "min")
		} else {
			widths <- c(widths.2, widths.3)
			ol <- findOverlapsOfPeaks(unique(p2g), unique(gr_oe.peaks), unique(gr_normal.peaks), minoverlap = 1, connectedPeaks = "min")
		} # if

	},

	{
		cat(sprintf("H23K27ac narrow peaks are not available for %s normal cell line\n", cancer_type))
		widths.encode <- end(gr_encode) - start(gr_encode)
		widths <- widths.encode
		ol <- findOverlapsOfPeaks(unique(p2g), gr_encode, minoverlap = 1, connectedPeaks = "min")
	}

  ) # switch
  tot <- 3.3e+9*0.98/mean(widths)






  # save overlaps_of_peaks for manual drawing of venn diagram.
  saveRDS(ol, sprintf("%s/find_overlaps_of_%s_peaks_output_overlappingpeaks_obj.rds", dir_rds, peaktype))
  #saveRDS(ol$overlappingPeaks, sprintf("%s/overlapping_%s_peaks.rds", dir_rds, peaktype))



  # begin of modification by H. Kim
  #names <- names(ol$overlappingPeaks)[4:6]
  names <- names(ol$overlappingPeaks)
  print(names)
  # [1] "gr_encode///unique.gr_normal.peaks."  "unique.p2g.///unique.gr_normal.peaks." [3] "unique.p2g.///gr_encode"

  # select names start with unique.p2g
  f_p2g <- grepl("^unique.p2g", names)
  names <- names[f_p2g]

  print(names)
  # [1] "unique.p2g.///unique.gr_normal.peaks." "unique.p2g.///unique.oe.peaks." [3] "unique.p2g.///gr_encode" # for ov with --exclude_encode_peaks
  # [1] "unique.p2g.///unique.gr_normal.peaks." "unique.p2g.///gr_encode" # for male-bc with --exclude_encode_peaks
  # [1] "unique.p2g.///unique.gr_normal.peaks." # for male-bc

  #total <- c(ol$overlappingPeaks[[names[1]]]$overlapFeature, ol$overlappingPeaks[[names[2]]]$overlapFeature, ol$overlappingPeaks[[names[3]]]$overlapFeature)
  total <- c()
  for (n_name in 1:length(names)) {
	total <- c(total, ol$overlappingPeaks[[names[n_name]]]$overlapFeature)
  }
  # end of modification

  # names(table(total)): {"includeFeature", "inside", "overlapEnd", "overlapStart"}
  pdf(sprintf("%s/piechart_overlapping_%s_peaks.pdf", dir_pdf, peaktype))
  pie1(table(total))
  dev.off()






  ### all.overlap

  # begin of modification by H. Kim
  #encode.overlap <- ol$overlappingPeaks[[names[3]]]
  #levels(factor(encode.overlap$overlapFeature))
  #colnames(encode.overlap)[8:12] <- c("seqnames2","start2","end2","width2","strand2")
  #encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "includeFeature",encode.overlap$width2,"fill")
  #encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "inside",encode.overlap$width,encode.overlap$overlap)
  #encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapEnd",encode.overlap$end2 - encode.overlap$start,encode.overlap$overlap)
  #encode.overlap$overlap <- ifelse(encode.overlap$overlapFeature == "overlapStart",encode.overlap$end - encode.overlap$start2,encode.overlap$overlap)
  #levels(factor(encode.overlap$overlap))
  #encode.overlap$overlap <- as.numeric(encode.overlap$overlap)
  #hist(encode.overlap$overlap)

  modify_overlapping_peaks <- function(ov) {

	#levels(factor(ov$overlapFeature))
	colnames(ov)[8:12] <- c("seqnames2", "start2", "end2", "width2", "strand2")
	ov$overlap <- ifelse(ov$overlapFeature == "includeFeature", ov$width2,"fill")
	ov$overlap <- ifelse(ov$overlapFeature == "inside", ov$width, ov$overlap)
	ov$overlap <- ifelse(ov$overlapFeature == "overlapEnd", ov$end2 - ov$start, ov$overlap)
	ov$overlap <- ifelse(ov$overlapFeature == "overlapStart", ov$end - ov$start2,ov$overlap)
	#print(levels(factor(ov$overlap)))
	ov$overlap <- as.numeric(ov$overlap)
	#print(hist(ov$overlap))

        ov

  } # modify_overlapping_peaks


  all.overlap <- NULL
  for (n_name in 1:length(names)) {

	print(names[n_name])
	ov <- ol$overlappingPeaks[[names[n_name]]]
	ov <- modify_overlapping_peaks(ov)
	if (nrow(ov) > 3) {
		print(head(ov[1:3,]))
	}
	ov <- ov[, c(1:12, ncol(ov), grep("overlapFeature",colnames(ov)))]

	all.overlap <- rbind(all.overlap, ov)
  } # for
  # end of modification


  pdf(sprintf("%s/hist_all_overlaps_of_%s_peaks.pdf", dir_pdf, peaktype))
  hist(all.overlap$overlap)
  dev.off()


  types <- c("overlapEnd","overlapStart")
  all.overlap <- all.overlap[all.overlap$overlapFeature %in% types,]


  pdf(sprintf("%s/hist_partial_overlaps_of_%s_peaks.pdf", dir_pdf, peaktype))
  hist(all.overlap$overlap)
  dev.off()












  ### venn diagram

  pdf(sprintf("%s/venn_chippeakanno_overlaps_of_%s_peaks.pdf", dir_pdf, peaktype))
  makeVennDiagram(ol, totalTest = tot, connectedPeaks = "min")
  dev.off()


  venn <- makeVennDiagram(ol, totalTest = tot, connectedPeaks = "min")
  saveRDS(venn, sprintf("%s/venn_chippeakanno_overlaps_of_%s_peaks.rds", dir_rds, peaktype))













  ### find cancer specific peak to gene links 

  # > ol$venn_cnt
  #     unique.p2g.  gr_encode  unique.gr_normal.peaks. Counts
  #[1,]           0          0                        0      0
  #[2,]           0          0                        1   4773
  #[3,]           0          1                        0 790739
  #[4,]           0          1                        1  45700
  #[5,]           1          0                        0   1000
  #[6,]           1          0                        1     11
  #[7,]           1          1                        0   3448
  #[8,]           1          1                        1   1152

  # names(ol$peaklist): {"unique.gr_normal.peaks.", "gr_encode", "gr_encode///unique.gr_normal.peaks.", "unique.p2g.", "unique.p2g.///unique.gr_normal.peaks.", "unique.p2g.///gr_encode", "unique.p2g.///gr_encode///unique.gr_normal.peaks."}

  # cancer specific peaks after excluding encode peaks and other peaks observed in normal cell lines.
  peakNames.p2g.only <- paste0(ol$peaklist$unique.p2g.@seqnames, ":", ol$peaklist$unique.p2g.@ranges)

  # cancer specific distal peaks
  # p2g.df.sub.plot.cancer_specific
  f_select_peaks <- (p2g.df.sub.plot.cancer.enriched$peakName %in% peakNames.p2g.only)

  p2g.df.sub.plot.cancer_specific <- p2g.df.sub.plot.cancer.enriched[f_select_peaks,]

  n_p2g.cancer_specific <- nrow(p2g.df.sub.plot.cancer_specific)
  var <- sprintf("n_p2g.%s.cancer_specific", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_p2g.cancer_specific))
  df_output[var, "info"] <- n_p2g.cancer_specific

  n_peaks.p2g.cancer_specific <- length(unique(p2g.df.sub.plot.cancer_specific$peakName))
  var <- sprintf("n_peaks.p2g.%s.cancer_specific", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_peaks.p2g.cancer_specific))
  df_output[var, "info"] <- n_peaks.p2g.cancer_specific

  # save
  #saveRDS(p2g.df.sub.plot.cancer_specific, sprintf("%s/cancer_specific_%s_p2g_table.rds", dir_rds, peaktype))











  ### find normal specific peak to gene links 

  # normal epithelial specific peaks
  # unique.gr_normal.peaks.
  peakNames.gr_normal.only <- paste0(ol$peaklist$unique.gr_normal.peaks.@seqnames, ":", ol$peaklist$unique.gr_normal.peaks.@ranges)

  # p2g.df.sub.plot.normal_specific
  f_select_peaks <- (p2g.df.sub.plot.cancer.enriched$peakName %in% peakNames.gr_normal.only)

  p2g.df.sub.plot.normal_specific <- p2g.df.sub.plot.cancer.enriched[f_select_peaks,]

  n_p2g.normal_specific <- nrow(p2g.df.sub.plot.normal_specific)
  var <- sprintf("n_p2g.%s.normal_specific", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_p2g.normal_specific))
  df_output[var, "info"] <- n_p2g.normal_specific

  n_peaks.p2g.normal_specific <- length(unique(p2g.df.sub.plot.normal_specific$peakName))
  var <- sprintf("n_peaks.p2g.%s.normal_specific", peaktype)
  cat(sprintf("\t%s: %d\n", var, n_peaks.p2g.normal_specific))
  df_output[var, "info"] <- n_peaks.p2g.normal_specific

  # save
  saveRDS(p2g.df.sub.plot.normal_specific, sprintf("%s/normal_specific_%s_p2g_table.rds", dir_rds, peaktype))





















  # Plot P2G heatmap for cancer specific distal elements
  #################################################################




  # begin of modification by H. Kim
  #test <- plotPeak2GeneHeatmap.distal(proj.archr, peaks=p2g.df.sub.plot.cancer.enriched,groupBy = "cluster.new",k=length(levels(factor(proj.archr$cluster.new))), corCutOff = .45,varCutOffATAC = 0,varCutOffRNA = 0,FDRCutOff = 1,nPlot =100000,palGroup=cols, palATAC = paletteContinuous("solarExtra"), palRNA = paletteContinuous("solarExtra"))
  cat(sprintf("\tplotPeak2GeneHeatmap.distal\n"))
  test <- plotPeak2GeneHeatmap.distal(proj.archr,
		peaks = p2g.df.sub.plot.cancer_specific,
		groupBy = "predictedGroup_ArchR",
		k = length(levels(factor(proj.archr$predictedGroup_ArchR))),
		corCutOff = .45,
		varCutOffATAC = 0,
		varCutOffRNA = 0,
		FDRCutOff = 1,
		nPlot =100000
	) # plotPeak2GeneHeatmap.distal

  filename_heatmap <- sprintf("%s/heatmap_peak2gene_legend_%s-cancerspecific.pdf", dir_pdf, peaktype)
  cat(sprintf("\tfilename_heatmap: %s\n", filename_heatmap))
  pdf(filename_heatmap, width = 8, height = 10)
  draw(test, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
  dev.off()
  # end of modification



  cat(sprintf("\tplotPeak2GeneHeatmap.distal for matrices\n"))
  list_atac_rna_p2g <- plotPeak2GeneHeatmap.distal(proj.archr,
		peaks = p2g.df.sub.plot.cancer_specific,
		groupBy = "predictedGroup_ArchR",
		k = length(levels(factor(proj.archr$predictedGroup_ArchR))),
                corCutOff = .45,
		varCutOffATAC = 0,
		varCutOffRNA = 0,
		FDRCutOff = 1,
		returnMatrices = TRUE, # 
		nPlot = 100000
	) # plotPeak2GeneHeatmap.distal

  # save
  p2g.df.sub.plot.cancer_specific$kmeans <- list_atac_rna_p2g$RNA$kmeansId
  saveRDS(p2g.df.sub.plot.cancer_specific, sprintf("%s/cancer_specific_%s_p2g_table.rds", dir_rds, peaktype))

  # write.table
  write.table(p2g.df.sub.plot.cancer_specific[,c("peakName", "geneName", "Correlation", "FDR")], sprintf("%s/cancer_specific_%s_p2g_table.tsv", dir_tsv, peaktype), row.names = FALSE, col.names = TRUE, sep = "\t", quote = F)

  # plot
  pdf(sprintf("%s/hist_genes_per_%s_peak-cancerspecific.pdf", dir_pdf, peaktype), width = 5,height = 3.5)
  hist(table(p2g.df.sub.plot.cancer_specific$idxRNA), main="Distribution of genes per distal peaks")
dev.off()


  pdf(sprintf("%s/hist_%s_peaks_per_gene-cancerspecific.pdf", dir_pdf, peaktype), width = 5,height = 3.5)
  hist(table(p2g.df.sub.plot.cancer_specific$idxATAC), main="Distribution of distal peaks per gene")
  dev.off()










  # begin of addition by H. Kim




  source("./r/perform_p2g_overlap_analysis.R")
  source("./r/motif_enrichment_analysis.R")






  #####################################################################
  #  Motif and Feature Enrichment


  # reference:
  # https://www.archrproject.com/bookdown
  # /home/hkim77/francolab/scRNA_scATAC_Breast_Cancer_Project_2021/ER+_Plus_TrueNormal_Patients_scATAC-PeakCalling_chromVAR.R
  # /home/hkim77/francolab/scRNA_scATAC_Breast_Cancer_Project_2021/ER+_Plus_TrueNormal_Patients_scATAC_scRNA-Plotting.R

  cat(sprintf("\n\n------------------------------------\n"))
  cat(sprintf("Motif enrichment\n"))
  cat(sprintf("\n------------------------------------\n\n"))

  peakset <- proj.archr@peakSet
  peakset.peakName <- sprintf("%s:%d-%d", seqnames(peakset), start(peakset), end(peakset))






  # markerPeaks.all
  markerPeaks.all <- readRDS(sprintf("%s/%s_markerpeaks_proj.archr_normal-vs-cancer.rds", dir_archr_rds, cancer_type))
  # > dim(markerPeaks.all) [1] 163345      1

  df_rowdata <- rowData(markerPeaks.all)
  peakName.markerpeaks.all <- sprintf("%s:%d-%d", df_rowdata$seqnames, df_rowdata$start, df_rowdata$end)
  
  if (length(cluster.type_exclude) > 0) {
	f <- (peakName.markerpeaks.all %in% peakset.peakName)
	markerPeaks.all	 <- markerPeaks.all[f,]
	df_rowdata <- rowData(markerPeaks.all)
	peakName.markerpeaks.all <- sprintf("%s:%d-%d", df_rowdata$seqnames, df_rowdata$start, df_rowdata$end)
  } # if





  # peakName.for_comparison
  # target: enhancer peaks in normal vs. tumor cells
  # p2g.df.sub: pval <= 1e-12, cor >= 0.45, distal 
  #    enhancer peaks in tumor cells (epi)
  #       peaks in kmeans clusters where epithelial cells showed high reads in scATAC-seq and scRNA-seq
  #       exclude known peaks overlapped with peaks in normal epithelial cells
  #    enhancer peaks in normal epithelial cells
  #    enhancer peaks in normal non-epi cells

  # p2g.df.sub.plot.cancer.enriched: enhancer peaks in store.kmeans
  # p2g.df.sub.plot.cancer_specific: p2g.df.sub.plot.cancer.enriched without enhancer peaks in normal epithelial cells

  # peakName.non_cancer_specific: peaks removed (i.e. enhancer peaks overlapped with peaks in normal epithelial cells)
  peakName.non_cancer_specific <- setdiff(p2g.df.sub.plot.cancer.enriched$peakName, p2g.df.sub.plot.cancer_specific$peakName)


  str_motif_set <- "cisbp"
  list_out <- add_archr_anno_enrichemnt(proj.archr, p2g.df.sub.plot, markerPeaks.all, peakName.non_cancer_specific, "normal-vs-cancer", "non-cancer-specific", str_motif_set, df_output, args)
  df_output <- list_out$df_output





  # peakName.for_comparison: 
  #   enhancer peaks in tumor cells (epi) without peakName.non_cancer_specific
  #   enhancer peaks in normal epithelial cells
  #   enhancer peaks in normal non-epi cells  
  peakName.for_comparison <- setdiff(p2g.df.sub.plot$peakName, peakName.non_cancer_specific)


  list_out <- add_archr_anno_enrichemnt(proj.archr, p2g.df.sub.plot, markerPeaks.all, peakName.for_comparison, "normal-vs-cancer", "cancer-specific", str_motif_set, df_output, args)
  df_output <- list_out$df_output










  ###################################################################
  # ChromVAR deviations enrichment













  ###################################################################
  # 15.4 Identification of Positive TF-Regulators
  # https://www.archrproject.com/bookdown/identification-of-positive-tf-regulators.html

  vec_matrices <- getAvailableMatrices(proj.archr)
  if (any(vec_matrices == "MotifMatrix")) {

	# https://www.archrproject.com/reference/getGroupSE.html
	seGroupMotif <- getGroupSE(ArchRProj = proj.archr,
			 useMatrix = "MotifMatrix",
			 groupBy = "predictedGroup_ArchR",
			 divideN = TRUE,
			 scaleTo = NULL,
			 threads = getArchRThreads(),
			 verbose = FALSE,
			 logFile = createLogFile("getGroupSE", logDir=dir_log, useLogs=TRUE)
		) # getGroupSE

  

	f_cancer <- (colnames(seGroupMotif) %in% cluster.types_with_cancer_cells)

	seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
	rowData(seZ)$maxDelta <- Matrix::rowMeans(assay(seZ)[,f_cancer,drop=F], na.rm=T) - Matrix::rowMeans(assay(seZ)[,!f_cancer,drop=F], na.rm=T)
   
	# https://www.archrproject.com/reference/correlateMatrices.html
	corGIM_MM <- correlateMatrices(
		ArchRProj = proj.archr,
		useMatrix1 = "GeneIntegrationMatrix_ArchR",
		useMatrix2 = "MotifMatrix",
		reducedDims = reducedDims, # default="IterativeLSI"
		dimsToUse = dimsToUse,
		scaleDims = NULL,
		corCutOff = 0.75, # A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the 'corCutOff', it will be excluded from analysis.
		k = 100,
		knnIteration = 500,
		overlapCutoff = 0.8,
		seed = 1,
		threads = getArchRThreads(),
		verbose = TRUE,
		logFile = createLogFile("correlateMatrices", logDir=dir_log, useLogs=TRUE)
	) # correlateMatrices

	corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

	corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
	corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]

	# overlapping peaks
	corGIM_MM$peakName <- sprintf("%s:%s", corGIM_MM$GeneIntegrationMatrix_ArchR_seqnames, ifelse(corGIM_MM$GeneIntegrationMatrix_ArchR_strand == 1, sprintf("%d-%d:+", corGIM_MM$GeneIntegrationMatrix_ArchR_start, corGIM_MM$GeneIntegrationMatrix_ArchR_end), sprintf("%d-%d:-", corGIM_MM$GeneIntegrationMatrix_ArchR_end, corGIM_MM$GeneIntegrationMatrix_ArchR_start)))
  
	f.over <- GRanges(corGIM_MM$peakName) %over% GRanges(peakName.for_comparison)
	corGIM_MM <- corGIM_MM[f.over,]


	if (nrow(corGIM_MM) > 0) {

		# positive TFs in cancer clusters
		corGIM_MM$TFRegulator <- "NO"
		corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75, na.rm=TRUE))] <- "YES"

		cat(sprintf("positive TF regulators in cancer clusters\n"))
		print(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))

		ggPos <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
			geom_point() + 
			theme_ArchR() +
			geom_vline(xintercept = 0, lty = "dashed") + 
			scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
			xlab("Correlation To Gene Expression") +
			ylab("Max TF Motif Delta") +
			scale_y_continuous(
				expand = c(0,0), 
				limits = c(min(corGIM_MM$maxDelta)*1.05, max(corGIM_MM$maxDelta)*1.05)
			)

		ggsave(sprintf("%s/scatterplot_%s_positive-tf-regulators.cancer_clusters.%s.pdf", dir_pdf, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggPos)
		#ggsave(sprintf("%s/scatterplot_%s_positive-tf-regulators.cancer_clusters.%s.png", dir_png, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggPos)


		# negative TFs in cancer clusters
		corGIM_MM$TFRegulator <- "NO"
		corGIM_MM$TFRegulator[which(corGIM_MM$cor < -0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75, na.rm=TRUE))] <- "YES"

		cat(sprintf("negative TF regulators in cancer clusters\n"))
		print(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))

		ggNeg <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
			geom_point() + 
			theme_ArchR() +
			geom_vline(xintercept = 0, lty = "dashed") + 
			scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
			xlab("Correlation To Gene Expression") +
			ylab("Max TF Motif Delta") +
			scale_y_continuous(
				expand = c(0,0), 
				limits = c(min(corGIM_MM$maxDelta)*1.05, max(corGIM_MM$maxDelta)*1.05)
			)

		ggsave(sprintf("%s/scatterplot_%s_negative-tf-regulators.cancer_clusters.%s.pdf", dir_pdf, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggNeg)
		#ggsave(sprintf("%s/scatterplot_%s_negative-tf-regulators.cancer_clusters.%s.png", dir_png, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggNeg)
  

		# positive TFs in normal clusters
		corGIM_MM$TFRegulator <- "NO"
		corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta < quantile(corGIM_MM$maxDelta, 0.25, na.rm=TRUE))] <- "YES"

		cat(sprintf("positive TF regulators in normal clusters\n"))
		print(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))

		ggPos <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
			geom_point() + 
			theme_ArchR() +
			geom_vline(xintercept = 0, lty = "dashed") + 
			scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
			xlab("Correlation To Gene Expression") +
			ylab("Max TF Motif Delta") +
			scale_y_continuous(
				expand = c(0,0), 
				limits = c(min(corGIM_MM$maxDelta)*1.05, max(corGIM_MM$maxDelta)*1.05)
			)

		ggsave(sprintf("%s/scatterplot_%s_positive-tf-regulators.normal_clusters.%s.pdf", dir_pdf, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggPos)
		#ggsave(sprintf("%s/scatterplot_%s_positive-tf-regulators.normal_clusters.%s.png", dir_png, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggPos)


		# negatve TFs in normal clusters
		corGIM_MM$TFRegulator <- "NO"
		corGIM_MM$TFRegulator[which(corGIM_MM$cor < -0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta < quantile(corGIM_MM$maxDelta, 0.25, na.rm=TRUE))] <- "YES"

		cat(sprintf("negative TF regulators in normal clusters\n"))
		print(sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1]))

		ggNeg <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
			geom_point() + 
			theme_ArchR() +
			geom_vline(xintercept = 0, lty = "dashed") + 
			scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
			xlab("Correlation To Gene Expression") +
			ylab("Max TF Motif Delta") +
			scale_y_continuous(
				expand = c(0,0), 
				limits = c(min(corGIM_MM$maxDelta)*1.05, max(corGIM_MM$maxDelta)*1.05)
			)

		ggsave(sprintf("%s/scatterplot_%s_negative-tf-regulators.normal_clusters.%s.pdf", dir_pdf, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggNeg)
		#ggsave(sprintf("%s/scatterplot_%s_negative-tf-regulators.normal_clusters.%s.png", dir_png, str_motif_set, peaktype), width = 5, height = 5.5, plot=ggNeg)

	} # if (nrow(corGIM_MM) > 0)

  } # if MotifMatrix is available











} # for peakType_

















# write xlsx
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\twrite xlsx\n"))
file_name_xlsx <- sprintf("%s/%s_sc-atac-seq_cancer_specific_p2g_summary.xlsx", dir_xlsx, cancer_type)
openxlsx::write.xlsx(df_output, file_name_xlsx, row.names = TRUE, col.names = TRUE)

print(df_output)



# post-process
unlink(dir_tmp)
unlink("Rplots.pdf")
unlink(sprintf("%s/*", dir_log))
unlink("VennDiagram.*.log")



# message
cat(sprintf("------------------------------------\n"))
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("completed successfully~\n"))
cat(sprintf("------------------------------------\n\n"))


# end of addition
# end of script






