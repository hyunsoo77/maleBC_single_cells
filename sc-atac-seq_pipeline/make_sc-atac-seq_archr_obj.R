#!/usr/bin/env Rscript
#
# make_sc-atac-seq_archr_obj.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Apr.
#
# description:
#   1) scATAC-seq data processing
#	creating an ArchR project
#	dimensionality reduction
#	clustering
#	single-cell embeddings
#   2) Gene scores and marker genes
#		markersGS.archr, useMatrix = "ArchRGeneScore", groupBy = "ATAC_clusters"
#		markersGS.archr.pred, useMatrix = "ArchRGeneScore", groupBy = "predictedGroup_ArchR"
#   3) Pseudo-bulk replicates and calling peaks, identifying marker peaks
#	ATAC_clusters:
#		addGroupCoverages
#		addReproduciblePeakSet
#		addPeakMatrix
#		addBgdPeaks
#		markerPeaks, useMatrix = "PeakMatrix", groupBy = "ATAC_clusters"
#
#	predictedGroup_ArchR:
#		addGroupCoverages
#		addReproduciblePeakSet
#		addPeakMatrix
#		addBgdPeaks
#		markerPeaks.archr, useMatrix = "PeakMatrix",  groupBy = "predictedGroup_ArchR"
#
#   4) Motif and feature encrichment
#   5) ChromVAR deviations enrichment
#   6) Footprinting
#   7) Integative analysis
#	Co-accessiblity 
#	Peak2GeneLinkage
#
# requirement:
#   ArchR v1.0.1
#	or modified packages of ArchR v1.0.1
#		1) ArchR.dpeerlab
#
# input:
#   cancer_type: (e.g. male-bc, ovar)
#
# options:
#   --n_log:  0 < n_log < 10: light logging, 10 < n_log < 100: heavy logging, n_log > 100: reserved for future usage
#   --n_debug:  0 < n_debug < 10: light debugging, 10 < n_debug < 100: heavy debugging, n_debug > 100: saveArchR for debug after createArrowFiles()
#   --n_cores: (default=20) the number of cores for parallel computing
#   --dir_count: defualt="../count"
#   --dir_output: (default="./output")
#   --dir_log: (default="./output/log")
#
#   --dir_seurat_obj: (default="")
#   --f_use_scrnaseq_harmony_cluster: suse sc-rna-seq harmony cluster
#   --keep_most_common_cell_type_for_each_cluster: (default=FALSE)
#   --max_dimstouse: (default=30)
#   --seurat_resolution: (default=0.8)  https://www.archrproject.com/bookdown/clustering-using-seurats-findclusters-function.html
#   --cluster.types_with_cancer_cells: default="" clusters with cancer cells (mostly epithelial with CNV), e.g. 2,3,4,11,12,16,17,27,30,37,0-Fibroblast,27-Fibroblast,21-Ciliated,16-Fibroblast
#
#   --archr_package: {["ArchR"], "ArchR.dpeerlab"}
#   --f_load_archr_obj: load ArchRProj (default=FALSE)
#
#   --min_tss: (default=4)
#   --min_frags: (default=1000)
#   --th_log10_tssenrichment: (default="") --th_log10_tssenrichment sample1=0.9,sample2=0.9
#   --th_uncertainty_tssenrichment: (default="") --th_uncertainty_tssenrichment sample1=0.05,sample2=0.05
#   --th_log10_nfrags: (default="") --th_log10_nfrags sample1=3.2,sample2=3.2
#   --th_uncertainty_nfrags: (default="") --th_uncertainty_nfrags sample1=0.05,sample2=0.05
#   --colorlim: (default=2)
#
#   --seed.harmony: seed for computing harmony
#   --harmony_vars: string separated with comma (e.g. donot,technology) {dataset, species, tumor.type, donor, technology, condition, treatment}
#   --harmony_theta: theta for computing harmony
#   --harmony_lambda: lambda for computing harmony, default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.
#
#   --umap_n_neighbors: (default=30)
#   --umap_min_dist: (default=0.3)
#   --umap_metric: (default="cosine")
#   --umap_metric_doubletscores: (default="euclidean")
#
#   --max_peaks_per_group: (default=150,000) A numeric threshold for the maximum peaks to retain per group from groupBy in the union reproducible peak set.
#   --min_cells_to_call_peaks: (default=50) The minimum allowable number of unique cells that was used to create the coverage files on which peaks are called. This is important to allow for exclusion of pseudo-bulk replicates derived from very low cell numbers.
#   --macs2_cutoff: (default=1e-7) The numeric significance cutOff for the testing method indicated by method.
#
#   --make_tfsee_input: make tfsee input (obsoleted)
#   --f_save_tmp_rds: (default=FALSE)
#
# output:
#   ./output/archr_output/
#	Annotations
#	ArrowFiles
#	Background-Peaks.proj.atac.rds
#	Background-Peaks.proj.archr.rds
#	Embeddings
#	GroupCoverages
#		ATAC_clusters
#		predictedGroup_ArchR
#	LSI_ATAC
#	Peak2GeneLinks
#	PeakCalls
#		InsertionBeds
#		ReplicateCalls
#	Plots
#	RNAIntegration
#		GeneIntegrationMatrix_ArchR
#	Save-ArchR-Project.rds
#	
#   ./output/rds/
#	${cancer_type}_archrproj_obj_final.rds
#	${cancer_type}_markergenes_overlap_rna_and_atac_predictedgroup.rds
#	${cancer_type}_markerpeaks_proj.archr_normal-vs-cancer.rds
#	${cancer_type}_markerpeaks_proj.archr_predictedgroup.rds
#   ./output/tmp: tmp files
#
#
#
#
# usage:
# for male-ba
# ./make_sc-atac-seq_archr_obj.R --min_tss 0 --min_frags 500  male-bc
# ./make_sc-atac-seq_archr_obj.R --f_load_archr_obj TRUE  male-bc
#
# for ovar
# ./make_sc-atac-seq_archr_obj.R --min_tss 0 --min_frags 500 --max_dimstouse 50 --umap_n_neighbors 30 --umap_min_dist 0.3 --umap_metric cosine --umap_metric_doubletscores cosine  ovar 
# ./make_sc-atac-seq_archr_obj.R --f_load_archr_obj TRUE --max_dimstouse 50 --umap_n_neighbors 30 --umap_min_dist 0.3 --umap_metric cosine  ovar
#
#
#
#
# debug:
# library(ArchR); proj <- loadArchRProject(path = "./All", force = FALSE, showLogo = TRUE)
# library(Seurat); dir_seurat_obj <- "/datastore/nextgenout5/share/labs/francolab/hyunsoo.kim/sc-rna-seq/male-bc/run-20211030/rds"; cancer_type <- "male-bc"; rna <- readRDS(sprintf("%s/%s_sc-rna-seq_merged_seurat_obj.rds", dir_seurat_obj, cancer_type));
# args <- list(); args$harmony_vars <- ""; args$harmony_theta <- "0";
#
#
#
# reference:
#
# 0. turorials
#	https://www.archrproject.com/bookdown
# 1. papers
#	https://www.nature.com/articles/s41588-021-00790-6
# 2. codes
#	https://github.com/RegnerM2015/scENDO_scOVAR_2020
#





suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser(description="hs script",python_cmd="python")
# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--n_log", default=1, type="integer", help="n_log")
parser$add_argument("--n_debug", default=0, type="integer", help="n_debug")
parser$add_argument("-l", "--n_log_level", default=0, type="integer", help="n_log_level")
parser$add_argument("-nc", "--n_cores", default=20, type="double", help="the number of cores for parallel computing")


# arguments for input/output
parser$add_argument("-dc", "--dir_count", default="../count", type="character", help="direcotry for directories of samples")
parser$add_argument("-do", "--dir_output", default="./output", type="character", help="direcotry for output")
parser$add_argument("-dl", "--dir_log", default="./output/log", type="character", help="direcotry for log")
#parser$add_argument("-ff", "--figure_format", default="pdf", type="character", help="figure format") # coded to print figures with pdf format (pdf and png whenever it is possible since png is helpful for lighter PPT slides).


# arguments for seurat
parser$add_argument("-dr", "--dir_seurat_obj", default="", type="character", help="direcotry for seurat objects")
parser$add_argument("-fushc", "--f_use_scrnaseq_harmony_cluster", action="store_true", default=FALSE, help="use sc-rna-seq harmony cluster")
parser$add_argument("-kmcct", "--keep_most_common_cell_type_for_each_cluster", dest="f_keep_most_common_cell_type_for_each_cluster", action="store_true", default=FALSE, help="keep most common cell type for each cluster")
parser$add_argument("-res", "--seurat_resolution", default=0.8, type="double", help="Seurat FindClusters parameter of resolution")
parser$add_argument("-ct", "--cluster.types_with_cancer_cells", default="", type="character", help="clusters with cancer cells (mostly epithelial with CNV)")


# arguments for archr
parser$add_argument("-ar", "--archr_package", default="ArchR", type="character", help="ArchR package name")
parser$add_argument("--f_load_archr_obj", default=FALSE, type="logical", help="f_load_archr_obj")


# arguments for quality control
parser$add_argument("-mintss", "--min_tss", default=4, type="integer", help="minTSS")
parser$add_argument("-minfrags", "--min_frags", default=1000, type="integer", help="minFrags")
parser$add_argument("-tlt", "--th_log10_tssenrichment", default="", type="character", help="th_log10_tssenrichment")
parser$add_argument("-tut", "--th_uncertainty_tssenrichment", default="", type="character", help="th_uncertainty_tssenrichment")
parser$add_argument("-tln", "--th_log10_nfrags", default="", type="character", help="th_log10_nfrags")
parser$add_argument("-tun", "--th_uncertainty_nfrags", default="", type="character", help="th_uncertainty_nfrags")
parser$add_argument("-clim", "--colorlim", default=2.0, type="double", help="density color limit")

# arguments for LSI
parser$add_argument("-lsiir", "--lsi_iter_rna", default=2, type="double", help="the number of LSI iterations to perform for multiome scRNA")
parser$add_argument("-lsivfr", "--lsi_varfeatures_rna", default=2500, type="double", help="the number of N variable features to use for LSI. The top N features will be used based on the selectionMethod for multiome scRNA")
parser$add_argument("-lsiia", "--lsi_iter_atac", default=2, type="double", help="the number of LSI iterations to perform for scATAC")
parser$add_argument("-lsivfa", "--lsi_varfeatures_atac", default=25000, type="double", help="the number of N variable features to use for LSI. The top N features will be used based on the selectionMethod for scATAC")


# arguments for harmony
parser$add_argument("-sh", "--seed.harmony", default=51, type="double", help="seed for computing harmony")
parser$add_argument("-hv", "--harmony_vars", default="", type="character", help="harmony vars (e.g. 'donor,technology')")
parser$add_argument("-ht", "--harmony_theta", default="-1", type="character", help="theta for computing harmony")
parser$add_argument("-hl", "--harmony_lambda", default="-1", type="character", help="lambda for computing harmony, default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.")


# arguments for umap
parser$add_argument("-maxdim", "--max_dimstouse", default=30, type="integer", help="max dimsToUse")
parser$add_argument("-n", "--umap_n_neighbors", default=30, type="integer", help="UMAP number of neighbors")
parser$add_argument("-dist", "--umap_min_dist", default=0.3, type="double", help="UMAP min distance")
parser$add_argument("-metric", "--umap_metric", default="cosine", type="character", help="UMAP metrics")
parser$add_argument("-metricds", "--umap_metric_doubletscores", default="euclidean", type="character", help="UMAP metrics for addDoubletScores")


# arguments for archr/addReproduciblePeakSet
parser$add_argument("-mp", "--max_peaks_per_group", default=150000, type="integer", help="max peasks per group") # function default=150,000
parser$add_argument("-mc", "--min_cells_to_call_peaks", default=50, type="integer", help="min cells to call peaks") # function default=25
parser$add_argument("-mco", "--macs2_cutoff", default=1e-7, type="double", help="The numeric significance cutOff for the testing method indicated by method") # function default=0.1


# arguments for save
parser$add_argument("-mti", "--make_tfsee_input", dest="f_make_tfsee_input", action="store_true", default=FALSE, help="make TFSEE input")
parser$add_argument("-str", "--f_save_tmp_rds", action="store_true", default=FALSE, help="save temporary rds files")

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
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("------------------------------------\n\n"))



# create output directories
dir_archr_output <- sprintf("%s/archr_output", dir_output)
#dir_archr_output <- normalizePath(dir_archr_output)
args$dir_archr_output <- dir_archr_output
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
dir_xlsx <- sprintf("%s/xlsx", dir_output)

dir.create(dir_archr_output, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_log, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_png, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rds, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tmp, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_xlsx, showWarnings = FALSE, recursive = TRUE)


unlink(sprintf("%s/*", dir_log), recursive=T)
unlink(sprintf("%s/Embeddings", dir_archr_output), recursive=T)
unlink(sprintf("%s/GroupCoverages", dir_archr_output), recursive=T)
unlink(sprintf("%s/IterativeLSI", dir_archr_output), recursive=T)
unlink(sprintf("%s/Peak2GeneLinks", dir_archr_output), recursive=T)
unlink(sprintf("%s/PeakCalls", dir_archr_output), recursive=T)
unlink(sprintf("%s/Plots", dir_archr_output), recursive=T)
unlink(sprintf("%s/RNAIntegration", dir_archr_output), recursive=T)





# load subroutines
source("./r/jupyter_message.R")
source("./r/identify_cell_types.R")








# set seed for reproducibility
set.seed(51)







############################################################


# string
suppressPackageStartupMessages(library(stringr))

# plots
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(patchwork))

# heatmap
suppressPackageStartupMessages(library(ComplexHeatmap))

# clustering
suppressPackageStartupMessages(library(ConsensusClusterPlus))

# color map
suppressPackageStartupMessages(library(viridis))

# genome annotation
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(ensembldb))

# GSEA
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))

# sc-RNA-seq
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(SingleR))

# data structure
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))


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

plan("multicore", workers = args$n_cores)
options(future.globals.maxSize = (min(8, args$n_cores) * 1024^3)) # <= 8GB






# read modified ArchR functions
source("./r/archr/archr_v1.0.1_modified/archr_v1.0.1_modified.R")
source("./r/archr/archr_v0.9.5_modified/filterDoublets_modified.R")
source("./r/archr/archr_v1.0.1_modified/GgplotUtils.R")



#set.seed(3)
addArchRThreads(threads=args$n_cores) 
addArchRGenome("hg38")



encode.all <- read.delim("./reference/genome_annotation/GRCh38-cCREs.bed.gz", header=F)
colnames(encode.all)[1:3] <- c("seqnames", "start", "end")


# store patient metadata and colors:
# make patient sample metadata and color assignments 

sampleColors <- RColorBrewer::brewer.pal(11, "Paired")
sampleColors[11] <- "#8c8b8b"
pie(rep(1,11), col=sampleColors) 


# color patient tumors to resemble the cancer ribbon color 
sampleColors <- c(sampleColors[5],sampleColors[7],sampleColors[6],sampleColors[8],sampleColors[10],sampleColors[9],sampleColors[4],sampleColors[3],sampleColors[2],sampleColors[11],sampleColors[1])




f_load_archr_obj <- FALSE

dimsToUse <- 1:args$max_dimstouse
umap_n_neighbors <- args$umap_n_neighbors
umap_min_dist <- args$umap_min_dist
umap_metric <- args$umap_metric
umap_metric_doubletscores <- args$umap_metric_doubletscores

# fragments.tsv.gz
dirs <- list.dirs(path=args$dir_count, recursive = FALSE)
sample_names <- gsub("-ATAC", "", basename(dirs))
inputFiles <- c()
for (dir1 in dirs) {
  filename <- sprintf("%s/outs/fragments.tsv.gz", dir1)
  if (file.exists(filename)) {
	inputFiles <- c(inputFiles, filename)
	next
  } 
  filename <- sprintf("%s/fragments.tsv.gz", dir1)
  if (file.exists(filename)) {
	inputFiles <- c(inputFiles, filename)
	next
  } 
  filename <- sprintf("%s/outs/atac_fragments.tsv.gz", dir1)
  if (file.exists(filename)) {
	inputFiles <- c(inputFiles, filename)
	next
  } 
  filename <- sprintf("%s/atac_fragments.tsv.gz", dir1)
  if (file.exists(filename)) {
	inputFiles <- c(inputFiles, filename)
	next
  } 
} # for





cluster.types_with_cancer_cells <- c()
if (nchar(dir_seurat_obj) > 0) {

  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("read in matching scRNA-seq\n"))
  fname_seurat_obj_rds <- sprintf("%s/%s_sc-rna-seq_merged_seurat_obj.rds", dir_seurat_obj, cancer_type)
  cat(sprintf("\tread %s\n", fname_seurat_obj_rds))
  rna <- readRDS(fname_seurat_obj_rds)
  print(rna)

  sample_names <- unique(rna$Sample)

  str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)
  str_column_of_meta_data_harmony <- sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse=","))

  str_umap_reduction <- "umap"
  col_umap1 <- "UMAP_1"; col_umap2 <- "UMAP_2"
  col_cluster_types <- "cluster.type"

  if (DefaultAssay(rna) == "integrated") {

        # current RNA.snn.res came from assay="integrated"
	col_seurat_cluster <- str_column_of_meta_data_cluster
	col_cluster_types <- "cluster.type"

	# use assay="RNA"
	# assay=integrated does not work for addGeneIntegrationMatrix since it does not have slot="counts".
	DefaultAssay(rna) <- "RNA"

	if (args$f_use_scrnaseq_harmony_cluster) {
		# match sc-rna-seq harmony and archr harmony
		if ((args$harmony_theta[1] >= 0) && ("umapharmony" %in% names(rna))) {
			col_seurat_cluster <- str_column_of_meta_data_harmony
			str_umap_reduction <- "umapharmony"
			col_umap1 <- "umapharmony_1"; col_umap2 <- "umapharmony_2"
			col_cluster_types <- "cluster.type.harmony"
		}
	} # if

  } else {

	col_seurat_cluster <- str_column_of_meta_data_cluster
	col_cluster_types <- "cluster.type"
	
	if ((args$harmony_theta[1] >= 0) && ("umapharmony" %in% names(rna))) {
		col_seurat_cluster <- str_column_of_meta_data_harmony
		str_umap_reduction <- "umapharmony"
		col_umap1 <- "umapharmony_1"; col_umap2 <- "umapharmony_2"
		col_cluster_types <- "cluster.type.harmony"
	}

  } # if
  cat(sprintf("\tcol_seurat_cluster=%s\n", col_seurat_cluster))
  cat(sprintf("\tcol_cluster_types=%s\n", col_cluster_types))

  if (args$f_keep_most_common_cell_type_for_each_cluster) {

	# keep_most_common_cell_type_for_each_cluster
	rna <- keep_most_common_cell_type_for_each_cluster(rna, args, col_seurat_cluster=col_seurat_cluster, col_cluster_types=col_cluster_types, col_cell_types="cell.type", th_ratio=-1, min_ncell=-1, n_log=1)
	print(rna)

  } # if
  

  # load wilcox.markers
  filename_wilcox <- sprintf("%s/wilcox_degs/%s_%s_wilcox_degs.rds", dir_seurat_obj, cancer_type, col_cluster_types)
  cat(sprintf("\tread %s\n", filename_wilcox))
  Wilcox.markers <- readRDS(filename_wilcox)

  if (nchar(args$cluster.types_with_cancer_cells) > 0) {
	# name cancer cell populations of interest (mostly epithelial with CNV)
	fac_cluster.type <- factor(rna@meta.data[, col_cluster_types])
	levels_cluster.type <- levels(fac_cluster.type)
	items <- strsplit(args$cluster.types_with_cancer_cells, ", ")[[1]]
	for (item in items) {
		if (grepl("-", item)) {
			# 0-Fibroblast
			cluster.types_with_cancer_cells <- c(cluster.types_with_cancer_cells, item) 
		} else {
			# 0
			idx_level <- grep(sprintf("^%s-", item), levels_cluster.type)
			if (length(idx_level) > 0) {
				cluster.types_with_cancer_cells <- c(cluster.types_with_cancer_cells, levels_cluster.type[idx_level])
			}
		}
	} # for
  } # if

} # if




# parameters

# gene score matrix
str_cutOff_markersGS_for_heatmap = "FDR <= 0.01 & Log2FC >= 1.25"

# peak matrix
str_cutOff_markerPeaks_for_heatmap <- "FDR <= 1e-3 & Log2FC >= 3"





















cat(sprintf("\n------------------------------------\n"))
cat(sprintf("parameters\n"))

cat(sprintf("\tn_cores=%d\n", args$n_cores))
cat(sprintf("\tdir_count=%s\n", args$dir_count))
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
cat(sprintf("\tf_load_archr_obj=%s\n", args$f_load_archr_obj))

# for harmony
cat(sprintf("\tseed.harmony=%g\n", args$seed.harmony))
args$harmony_vars <- strsplit(args$harmony_vars, ", ")[[1]]
cat(sprintf("\tharmony_vars=%s\n", paste(args$harmony_vars, collapse=", ")))

args$harmony_theta <- as.numeric(strsplit(args$harmony_theta, ", ")[[1]])
cat(sprintf("\tharmony_theta=%s\n", paste(args$harmony_theta, collapse=", ")))

args$harmony_lambda <- as.numeric(strsplit(args$harmony_lambda, ", ")[[1]])
if (args$harmony_lambda[1] < 0) {
	# user did not provide with harmony_lambda
	args$harmony_lambda <- rep(1.0, length(args$harmony_theta))
}
cat(sprintf("\tharmony_lambda=%s\n", paste(args$harmony_lambda, collapse=", ")))

# for umap
cat(sprintf("\tumap_n_neighbors=%d\n", umap_n_neighbors))
cat(sprintf("\tumap_min_dist=%.1f\n", umap_min_dist))
cat(sprintf("\tumap_metric=%s\n", umap_metric))

# arguments for archr/addReproduciblePeakSet
cat(sprintf("\tmax_peaks_per_group=%d\n", args$max_peaks_per_group))
cat(sprintf("\tmin_cells_to_call_peaks=%d\n", args$min_cells_to_call_peaks))
cat(sprintf("\tmacs2_cutoff=%g\n", args$macs2_cutoff))











###########################################################
# PART 1: scATAC-seq data processing
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 1: scATAC-seq data processing\n"))
cat(sprintf("------------------------------------\n\n"))
###########################################################






if (args$f_load_archr_obj) {

  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("loadArchRProject\n"))
  proj <- loadArchRProject(path = dir_archr_output, force = FALSE, showLogo = TRUE)
  ArrowFiles <- getArrowFiles(ArchRProj = proj)

} else {

  # Create Arrow and ArchR project
  ##########################################################
  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("createArrowFiles\n"))
  cat(sprintf("\tminTSS=%d\n", args$min_tss))
  cat(sprintf("\tminFrags=%d\n", args$min_frags))
  cat(sprintf("\tinputFiles\n"))
  cat(sprintf("\t\t%s\n", inputFiles), sep="")

  # https://www.archrproject.com/reference/createArrowFiles.html
  ArrowFiles <- createArrowFiles(
    inputFiles = inputFiles,
    sampleNames = sample_names,
    outputNames = sprintf("%s/%s", dir_output, sample_names),
    #filterFrags = 1000, # default=1000 for ArchR v0.9.5, M. Regner: filterFrags = 0
    #filterTSS = 4, # default=4 for ArchR v0.9.5, M. Regner: filterTSS = 0 Dont set this too high because you can always increase later
    validBarcodes = NULL,
    geneAnnotation = getGeneAnnotation(),
    genomeAnnotation = getGenomeAnnotation(),
    minTSS = args$min_tss, # default=4, ArchR v0.9.5 does not have the parameter of minTSS, but it had filterTSS = 4, set minTSS = 0 when you want consistency with M. Regner: filterTSS = 0
    minFrags = args$min_frags, # default=1000, ArchR v0.9.5 minFrags=500, ArchR v1.0.1 minFrags=1000
    maxFrags = 1e+05, # default=1e+05, ArchR v0.9.5 maxFrags=100000
    QCDir = sprintf("%s/qc", args$dir_output),
    nucLength = 147, # The length in basepairs that wraps around a nucleosome. This number is used for identifying fragments as sub-nucleosome-spanning, mono-nucleosome-spanning, or multi-nucleosome-spanning.
    promoterRegion = c(2000, 100),
    TSSParams = list(),
    excludeChr = c("chrM", "chrY"),
    nChunk = 5,
    bcTag = "qname",
    gsubExpression = NULL,
    bamFlag = NULL,
    offsetPlus = 4,
    offsetMinus = -5,
    addTileMat = TRUE,
    TileMatParams = list(),
    addGeneScoreMat = FALSE,
    GeneScoreMatParams = list(),
    force = FALSE, # A boolean value indicating whether to force ArrowFiles to be overwritten if they already exist.
    threads = getArchRThreads(),
    parallelParam = NULL,
    subThreading = FALSE, # default=TRUE, When user do not pass "subThreading = F" to createArrowFiles(), arrow files can be created correctly, but later function calls which want to read or write arrow files will lead to this "HDF5. File accessibilty. Unable to open file." error.  https://github.com/GreenleafLab/ArchR/issues/248
    verbose = TRUE,
    cleanTmp = TRUE,
    logFile = createLogFile("createArrows", logDir=dir_log),
    filterFrags = NULL,
    filterTSS = NULL
  ) # createArrowFiles




  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("addDoubletScores\n"))
  cat(sprintf("\tumap_metric_doubletscores=%s\n", umap_metric_doubletscores))

  # addDoubletScores
  doubScores <- addDoubletScores(
    input = ArrowFiles,
    useMatrix = "TileMatrix",
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    nTrials = 5,
    dimsToUse = dimsToUse,
    LSIMethod = 1,
    scaleDims = FALSE, # A boolean that indicates whether to z-score the reduced dimensions for each cell during the LSI method performed for doublet determination. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs.
    corCutOff = 0.75,
    knnMethod = "UMAP",
    UMAPParams = list(n_neighbors=umap_n_neighbors, min_dist=umap_min_dist, metric=umap_metric_doubletscores, verbose =FALSE),
    LSIParams = list(outlierQuantiles = NULL, filterBias = FALSE),
    #outDir = getOutputDirectory(input), #  If null this returns "QualityControl" directory.
    outDir = sprintf("%s/qc", args$dir_output),
    threads = getArchRThreads(),
    force = FALSE, # If the UMAP projection is not accurate (when R < 0.8 for the reprojection of the training data - this occurs when you have a very homogenous population of cells), setting force=FALSE will return -1 for all doubletScores and doubletEnrichments. If you would like to override this (not recommended!), you can bypass this warning by setting force=TRUE.
    parallelParam = NULL,
    verbose = TRUE,
    logFile = createLogFile("addDoubletScores", logDir=dir_log)
  ) # addDoubletScores


  # ArchRProject
  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("ArchRProject\n"))
  # https://www.archrproject.com/reference/ArchRProject.html
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = dir_archr_output,
    copyArrows = TRUE, # A boolean value indicating whether ArrowFiles should be copied into outputDirectory. This is recommened so that you maintain an unaltered copy for later usage.
    geneAnnotation = getGeneAnnotation(),
    genomeAnnotation = getGenomeAnnotation(),
    showLogo = TRUE,
    threads = getArchRThreads()
  ) # ArchRProject


  if (args$n_debug > 100) {
    cat(sprintf("\tsaveArchRProject for debug\n"))
    saveArchRProject( ArchRProj = proj,
	outputDirectory = dir_archr_output,
	overwrite = TRUE,
	load = FALSE,
	dropCells = FALSE,
	logFile = createLogFile("saveArchRProject", logDir=dir_log), threads = getArchRThreads()
    ) # saveArchRProject
  } # if
  

} # if

























#######################################################
# Determine df_tss and df_depth for each sample
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("determine df_tss and df_depth for each sample\n"))

suppressPackageStartupMessages(library(mclust))

for (sample in sample_names) {

  proj.i <- proj[proj$Sample == sample]
  
  cat(sprintf("\n%s\n", sample))
  cat(sprintf("\t%s length(proj.i$nFrags)=%d\n", sample, length(proj.i$nFrags)))
  cat(sprintf("\t%s length(proj.i$TSSEnrichment)=%d\n", sample, length(proj.i$TSSEnrichment)))


  ### GMM for TSS per cell

 
  # modification of proj.i$TSS.cluster and TSS.cluster.uncertainty
  th_log10_tssenrichment <- NULL
  switch(sample,
	"3571DL"={ th_log10_tssenrichment <- 0.80 },
	"38FE7L"={ th_log10_tssenrichment <- 0.80 },
	"446B7L"={ th_log10_tssenrichment <- 0.90 },
	"4CC61L"={ th_log10_tssenrichment <- 0.90 },
	"52BC3L_multiome"={ th_log10_tssenrichment <- 0.90 },
	"5572CL_multiome"={ th_log10_tssenrichment <- 0.90 },
	{ 
		items <- strsplit(args$th_log10_tssenrichment, "[=,]")[[1]]
		idx <- match(sample, items)
		if (!is.na(idx)) {
			th_log10_tssenrichment <- as.numeric(items[idx+1])
		}
		
		if (is.null(th_log10_tssenrichment)) {
			TSS.clust <- tryCatch({
				TSS.clust <- Mclust(log10(proj.i$TSSEnrichment+1), G = 2)
        		}, error = function(e) {
               			print(paste('#error-handler-code'))
				cat(sprintf('%s',e))
				return(NULL)
			}, finally = {
			})

			if (!is.null(TSS.clust)) {
				proj.i$TSS.cluster <- TSS.clust$classification
				proj.i$TSS.cluster.uncertainty <- TSS.clust$uncertainty
				# no modificaiton
				th_log10_tssenrichment <- NULL
				f.cluster2 <- (proj.i$TSS.cluster == 2)
				log10_tssenrichment <- log10(proj.i$TSSEnrichment+1)
				cat(sprintf("\t%s proj.i$TSS.cluster was determined by Mclust(), min(log10_tssenrichment) in cluster2: %.2f\n", sample, min(log10_tssenrichment[f.cluster2], na.rm=T)))
			} else {
				stop(sprintf("no th_log10_tssenrichment for %s: try to use --th_log10_tssenrichment %s=0.90", sample, sample))
			}
		} # if
	}
  ) # switch

  if (!is.null(th_log10_tssenrichment)) {
  	cat(sprintf("\t%s proj.i$TSS.cluster was determined by log10(proj.i$TSSEnrichment+1) >= %.2f\n", sample, th_log10_tssenrichment))
  	proj.i$TSS.cluster <- ifelse(log10(proj.i$TSSEnrichment+1) >= th_log10_tssenrichment, "2", "1")
  	proj.i$TSS.cluster.uncertainty <- rep(NA, nrow(proj.i@cellColData))
  }

  # plot tss.png 
  gg <- ggPoint(
	x = log10(proj.i$nFrags),
	y = log10(proj.i$TSSEnrichment+1),
	color = as.character(proj.i$TSS.cluster),
	discrete = T,
	xlabel = TeX(r'($log_{10}$ (unique fragments))'),
	ylabel = TeX(r'($log_{10}$ (TSS Enrichment+1))'),
	rastr = TRUE
  ) + 
  ggtitle(paste0("GMM classification:\n", sample, " TSS Enrichment"))

  ggsave(sprintf("%s/%s_%s_atac_tss.pdf", dir_pdf, cancer_type, sample), width = 4, height = 4, plot=gg)
  ggsave(sprintf("%s/%s_%s_atac_tss.png", dir_png, cancer_type, sample), width = 4, height = 4, plot=gg)
  
  
  df.TSS <- data.frame(proj.i$cellNames, proj.i$TSS.cluster, proj.i$TSS.cluster.uncertainty, proj.i$TSSEnrichment)

  cat(sprintf("\t%s df.TSS was filered by proj.i.TSS.cluster=2\n", sample))
  df.TSS <- dplyr::filter(df.TSS, proj.i.TSS.cluster == "2")

  switch(sample,
	"sample"={ th_uncertainty <- NULL },
	{
		items <- strsplit(args$th_uncertainty_tssenrichment, "[=,]")[[1]]
		idx <- match(sample, items)
		if (!is.na(idx)) {
			th_uncertainty <- as.numeric(items[idx+1])
		} else {
			th_uncertainty <- 0.05
		} # if
	}
  ) # switch
  
  if (all(is.na(proj.i$TSS.cluster.uncertainty))) {
	th_uncertainty <- NULL
  }

  if (!is.null(th_uncertainty)) {
  	#df.TSS <- dplyr::filter(df.TSS, proj.i.TSS.cluster.uncertainty <= 0.05)
        cat(sprintf("\t%s df.TSS was filtered with TSS.cluster.uncertainty <= %.2f\n", sample, th_uncertainty))
  	df_tmp <- dplyr::filter(df.TSS, proj.i.TSS.cluster.uncertainty <= th_uncertainty)
	if (nrow(df_tmp) > 0) {
		df.TSS <- df_tmp
	} else {
		th_uncertainty <- quantile(proj.i$TSS.cluster.uncertainty, prob=0.20)
		cat(sprintf("\tfilter by proj.i.TSS.cluster.uncertainty <= %.2f\n", th_uncertainty))
		df.TSS <- dplyr::filter(df.TSS, proj.i.TSS.cluster.uncertainty <= th_uncertainty)
		print(head(df.TSS))
	}
  }

  # saveRDS
  filename <- sprintf("%s/df_tss_%s.rds", dir_rds, sample)
  cat(sprintf("\tsave RDS: %s\n", filename))
  #saveRDS(df.TSS,paste0("df_TSS_", sample, ".rds"))
  saveRDS(df.TSS, filename)
  





  ### GMM for fragments per cell

  # modifcation of proj.i$depth.cluster and depth.cluster.uncertainty
  th_log10_nfrags <- NULL
  switch(sample,
	"52BC3L_multiome"={ th_log10_nfrags <- 3.2 },
        {
		# no modification
		th_log10_nfrags <- NULL
		items <- strsplit(args$th_log10_nfrags, "[=,]")[[1]]
		idx <- match(sample, items)
		if (!is.na(idx)) {
			th_log10_nfrags <- as.numeric(items[idx+1])
		}
		
		if (is.null(th_log10_nfrags)) {
			depth.clust <- tryCatch({
  				depth.clust <- Mclust(log10(proj.i$nFrags), G = 2)
        		}, error = function(e) {
               			print(paste('#error-handler-code'))
				cat(sprintf('%s',e))
				return(NULL)
			}, finally = {
			})

			if (!is.null(depth.clust)) {
				proj.i$depth.cluster <- depth.clust$classification
				proj.i$depth.cluster.uncertainty <- depth.clust$uncertainty
				# no modificaiton
				th_log10_nfrags <- NULL
				f.cluster2 <- (proj.i$depth.cluster == 2)
				log10_nfrags <- log10(proj.i$nFrags)
				cat(sprintf("\t%s proj.i$depth.cluster was determined by Mclust(), min(log10_nfrags) in cluster2: %.2f\n", sample, min(log10_nfrags[f.cluster2], na.rm=T)))
			} else {
				stop(sprintf("no th_log10_nfrags for %s: try to use --th_log10_nfrags %s=3.2", sample, sample))
			}
		} # if
        }
  ) # switch

  if (!is.null(th_log10_nfrags)) {
  	cat(sprintf("\t%s proj.i$depth.cluster was determined by log10(proj.i$nFrags) >= %.2f\n", sample, th_log10_nfrags))
	proj.i$depth.cluster <- ifelse(log10(proj.i$nFrags) >= th_log10_nfrags, "2", "1")
	proj.i$depth.cluster.uncertainty <- rep(NA,nrow(proj.i@cellColData))
  }

  # plot depth.png
  gg <- ggPoint(
    x = log10(proj.i$nFrags),
    y = log10(proj.i$TSSEnrichment+1),
    color = as.character(proj.i$depth.cluster),
    xlabel = TeX(r'($log_{10}$ (unique fragments))'),
    ylabel = TeX(r'($log_{10}$ (TSS Enrichment+1))'),
    rastr = TRUE
  ) + 
  ggtitle(paste0("GMM classification:\n", sample, " Fragments"))

  ggsave(sprintf("%s/%s_%s_atac_depth.pdf", dir_pdf, cancer_type, sample), width = 4, height = 4, plot=gg)
  ggsave(sprintf("%s/%s_%s_atac_depth.png", dir_png, cancer_type, sample), width = 4, height = 4, plot=gg)

  df.depth <- data.frame(proj.i$cellNames, proj.i$depth.cluster, proj.i$depth.cluster.uncertainty, proj.i$nFrags)

  cat(sprintf("\t%s df.depth was filered by proj.i.depth.cluster=2\n", sample))
  df.depth <- dplyr::filter(df.depth, proj.i.depth.cluster == "2")
  switch(sample,
	"446B7L"={ th_uncertainty <- 0.10 },
	"4CC61L"={ th_uncertainty <- 0.10 },
	{
		items <- strsplit(args$th_uncertainty_nfrags, "[=,]")[[1]]
		idx <- match(sample, items)
		if (!is.na(idx)) {
			th_uncertainty <- as.numeric(items[idx+1])
		} else {
			th_uncertainty <- 0.05
		} # if
	}
  ) # switch

  if (all(is.na(proj.i$depth.cluster.uncertainty))) {
	th_uncertainty <- NULL
  }

  if (!is.null(th_uncertainty)) {
    #df.depth <- dplyr::filter(df.depth, proj.i.depth.cluster.uncertainty <= 0.05)
    cat(sprintf("\t%s df.depth was filtered with depth.cluster.uncertainty <= %.2f\n", sample, th_uncertainty))
    df_tmp <- dplyr::filter(df.depth, proj.i.depth.cluster.uncertainty <= th_uncertainty)
    if (nrow(df_tmp) > 0) {
	df.depth <- df_tmp
    } else {
	th_uncertainty <- quantile(proj.i$depth.cluster.uncertainty, prob=0.20)
	cat(sprintf("\tfilter by proj.i.depth.cluster.uncertainty <= %.2f\n", th_uncertainty))
	df.depth <- dplyr::filter(df.depth, proj.i.depth.cluster.uncertainty <= th_uncertainty)
	print(head(df.depth))
    }
  } # if



  # saveRDS
  filename <- sprintf("%s/df_depth_%s.rds", dir_rds, sample)
  cat(sprintf("\tsave RDS: %s\n", filename))
  #saveRDS(df.depth,paste0("df_depth_", sample, ".rds"))
  saveRDS(df.depth, filename)
  
  f_color_density <- TRUE
  x <- log10(proj.i$nFrags)
  y <- log10(proj.i$TSSEnrichment+1)
  h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y))
  if (any(h <= 0)) {
    # stop("bandwidths must be strictly positive")  https://rdrr.io/cran/MASS/src/R/kde2d.R
    f_color_density <- FALSE
  } # if

  gg <- ggPoint(
	x = log10(proj.i$nFrags),
	y = log10(proj.i$TSSEnrichment+1),
	colorDensity = f_color_density,
	#colorLimits = c(0, 2.25),
	colorLimits = c(0, args$colorlim+0.25),
	continuousSet = "sambaNight",
	xlabel = TeX(r'($log_{10}$ (unique fragments))'),
	ylabel = TeX(r'($log_{10}$ (TSS Enrichment+1))'),
	rastr = TRUE
    ) +
    geom_hline(yintercept = log10(min(df.TSS$proj.i.TSSEnrichment)+1), linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)), linetype = "dashed")+
    ggtitle(paste0("QC thresholds:\n", sample))

  ggsave(sprintf("%s/%s_%s_atac_qc.pdf", dir_pdf, cancer_type, sample), width = 4, height = 4, plot=gg)
  ggsave(sprintf("%s/%s_%s_atac_qc.png", dir_png, cancer_type, sample), width = 4, height = 4, plot=gg)
  
  gg <- ggPoint(
	x = log10(proj.i$nFrags),
	y = log10(proj.i$TSSEnrichment+1),
	color = proj.i$DoubletEnrichment,
	discrete = F,
	continuousSet = "sambaNight",
	xlabel = TeX(r'($log_{10}$ (unique fragments))'),
	ylabel = TeX(r'($log_{10}$ (TSS Enrichment+1))'),
	rastr = TRUE
  ) +
    geom_hline(yintercept = min(log10(df.TSS$proj.i.TSSEnrichment+1)), linetype = "dashed")+
    geom_vline(xintercept = min(log10(df.depth$proj.i.nFrags)), linetype = "dashed")+
    ggtitle(paste0("Doublet Enrichment:\n", sample))

  ggsave(sprintf("%s/%s_%s_atac_doublets.pdf", dir_pdf, cancer_type, sample),width = 4,height = 4, plot=gg)
  ggsave(sprintf("%s/%s_%s_atac_doublets.png", dir_png, cancer_type, sample),width = 4,height = 4, plot=gg)
  
} # for









##########################################################
# Filter out low quality cells, and remove doublets
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("Filter out low quality cells, and remove doublets\n"))

#list.depth <- list.files(pattern = "^df_depth")
list.depth <- list.files(path=dir_rds, pattern = "^df_depth", full.names = TRUE, recursive = FALSE)

df.depth <-  data.frame(cellNames=character(),
                        cluster=character(),
                        cluster.uncertainty=character(),
                        nFrags = character())

for (fname_rds in list.depth){
  cat(sprintf("\tread %s\n", fname_rds)) 
  df <- readRDS(fname_rds)
  colnames(df) <- c("cellNames", "cluster", "cluster.uncertainty", "nFrags")
  df.depth <- rbind(df.depth,df)
} # for

list.TSS <- list.files(path=dir_rds, pattern = "^df_tss", full.names = TRUE, recursive = FALSE)

df.TSS <-  data.frame(cellNames=character(),
                      cluster=character(),
                      cluster.uncertainty=character(),
                      TSSEnrichment = character())

for (fname_rds in list.TSS){
  cat(sprintf("\tread %s\n", fname_rds)) 
  df <- readRDS(fname_rds)
  colnames(df) <- c("cellNames", "cluster", "cluster.uncertainty", "TSSEnrichment")
  df.TSS <- rbind(df.TSS,df)
} # for


colnames(df.TSS) <- c("cellNames", "TSS.cluster", "TSS.cluster.uncertainty", "TSSEnrichment")
colnames(df.depth) <- c("cellNames", "depth.cluster", "depth.cluster.uncertainty", "nFrags")




# cellsPass
cellsPass <- intersect(df.TSS$cellNames, df.depth$cellNames)


# cellsFail
cellsFail <- proj$cellNames[!(proj$cellNames %in% cellsPass)]





# proj.filter
# screen for high quality barcodes (remove non cellular barcodes)
proj.filter <- proj[proj$cellNames %in% cellsPass]





df_output <- data.frame()
for (sample in sample_names) {

  proj.i <- proj[proj$Sample == sample]
  meta.data <- as.data.frame(getCellColData(proj.i))
  df_output["StartingNumCells", sample] <- nCells(proj.i)
  proj.i <- proj.filter[proj.filter$Sample == sample]
  df_output["PostThresholdingNumCells", sample] <- nCells(proj.i)

} # for








cat(sprintf("\tfilterDoublets\n"))
proj.filter <- filterDoublets(proj.filter,
		cutEnrich = 1,  # The minimum numeric cutoff for DoubletEnrichment. This number is equivalent to the number of simulated doublets identified as a nearest neighbor to the cell divided by the expected number given a random uniform distribution. https://www.archrproject.com/reference/filterDoublets.html
		cutScore = -Inf,
		filterRatio = 1 # The maximum ratio of predicted doublets to filter based on the number of pass-filter cells. For example, if there are 5000 cells, the maximum would be filterRatio * 5000^2 / (100000) (which simplifies to filterRatio * 5000 * 0.05). This filterRatio allows you to apply a consistent filter across multiple different samples that may have different percentages of doublets because they were run with different cell loading concentrations. The higher the filterRatio, the greater the number of cells potentially removed as doublets.  https://www.archrproject.com/reference/filterDoublets.html
) # filterDoublets





for (sample in sample_names) {
  proj.i <- proj.filter[proj.filter$Sample == sample]
  df_output["PostFilterDoubletsNumCells", sample] <- nCells(proj.i)
}

cat(sprintf("\tdf_output\n"))
print(df_output)




cat(sprintf("\tplotFragmentSizes\n"))
gg <- plotFragmentSizes(proj.filter,
		logFile = createLogFile("plotFragmentSizes", logDir=dir_log))+
	ggtitle("Fragment Size Histogram")
#gg <- plotFragmentSizes(proj.filter, threads=1)+ggtitle("Fragment Size Histogram")

ggsave(sprintf("%s/%s_atac_frags_hist.pdf", dir_pdf, cancer_type), width = 6, height = 4, plot=gg)
ggsave(sprintf("%s/%s_atac_frags_hist.png", dir_png, cancer_type), width = 6, height = 4, plot=gg)


cat(sprintf("\tplotTSSEnrichment\n"))
gg <- plotTSSEnrichment(proj.filter,
		logFile = createLogFile("plotTSSEnrichment", logDir=dir_log))+
	ggtitle("TSS Enrichment")

ggsave(sprintf("%s/%s_atac_tss.pdf", dir_pdf, cancer_type), width = 6, height = 4, plot=gg)
ggsave(sprintf("%s/%s_atac_tss.png", dir_png, cancer_type), width = 6, height = 4, plot=gg)




# save cell names
write.table(proj.filter$cellNames, file=sprintf("%s/qc/%s_cell_names_filtered.txt", dir_output, cancer_type), quote=FALSE, row.names=F, col.names=F)




if (!args$f_save_tmp_rds) {
	# remove temp files
	unlink(sprintf("%s/df_depth_*.rds", dir_rds))
	unlink(sprintf("%s/df_tss_*.rds", dir_rds))
} # if


















###################################################
# check multiome
cat(sprintf("\tcheck the overlap of cell names between scATAC and scRNA\n"))


# https://rdrr.io/github/GreenleafLab/ArchR/src/R/MatrixGeneExpression.R
allCells <- rownames(getCellColData(proj.filter))

if ("package:ArchR.dpeerlab" %in% (search())) {
	cellsInArrows <- unlist(lapply(ArrowFiles, ArchR.dpeerlab:::.availableCells), use.names=FALSE)
} else {
	cellsInArrows <- unlist(lapply(ArrowFiles, ArchR:::.availableCells), use.names=FALSE)
}

if (!is.null(allCells)) {
    cellsInArrows <- allCells
}

overlap_cell_names_between_atac_and_rna <- sum(cellsInArrows %in% colnames(rna)) / length(cellsInArrows)
#.logMessage("Overlap w/ scATAC = ", round(overlap_cell_names_between_atac_and_rna, 3), logFile = createLogFile("checkMultiome", logDir=dir_log), verbose = TRUE)
cat(sprintf("\t\toverlap of cell names between scATAC and scRNA: %.3f\n", overlap_cell_names_between_atac_and_rna))


f_multiome <- grepl("multiome", cancer_type)
#f_multiome <- (overlap_cell_names_between_atac_and_rna > 0.1)
cat(sprintf("\t\tf_multiome: %s\n", f_multiome))






if (f_multiome) {


  cat(sprintf("\t\toverlap of cell names between scATAC and scRNA for each sample\n"))
  for (sample in sample_names) {
	cells.rna.sample <- colnames(rna)[(rna$Sample == sample)]
	cells.atac.sample <- proj.filter$cellNames[(proj.filter$Sample == sample)]
	cells_overlapped <- cells.atac.sample[cells.atac.sample %in% cells.rna.sample]
	overlap_ratio_atac_and_rna.sample <- sum(cells.atac.sample %in% cells.rna.sample) / length(cells.atac.sample)
	cat(sprintf("\t\t\t%s: %g\n", sample, overlap_ratio_atac_and_rna.sample))
 	df_output["scRNA-seqNumCells", sample] <- length(cells_overlapped)
  } # for


  # use overlapped cells
  cat(sprintf("\tselect cells overlapped with scRNA-seq barcodes\n"))
  cells_overlapped <- proj.filter$cellNames[proj.filter$cellNames %in% colnames(rna)]
  proj.filter <- proj.filter[proj.filter$cellNames %in% cells_overlapped]
  cat(sprintf("\t\t# of cells: %d\n", nCells(proj.filter)))
  




  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("add scRNA\n"))

  # https://www.archrproject.com/reference/addGeneExpressionMatrix.html
  suppressPackageStartupMessages(library(SummarizedExperiment))
  counts_matrix <- GetAssayData(rna, assay="RNA", slot="counts")
  geneAnnotation <- getGeneAnnotation(proj.filter)
  idx <- match(rownames(counts_matrix), geneAnnotation$gene$symbol)
  f <- !is.na(idx); idx <- idx[f]
  se_rna <- SummarizedExperiment(assays=list(counts=counts_matrix[f,]), rowRanges=geneAnnotation$gene[idx], colData=rna@meta.data)

  proj.filter <- addGeneExpressionMatrix(
	input = proj.filter,
	seRNA = se_rna, # A a scRNA-seq SummarizedExperiment (gene x cell) to be integrated with the scATAC-seq data. Cell names from this object much match those of the cell names in the ArrowFiles/ArchRProject. We will add support shortly for Seurat Objects (see Seurat::as.SingleCellExperiment). The provided values MUST be in counts (integer), not log transformed.
	chromSizes = getChromSizes(proj.filter),
	excludeChr = c("chrM", "chrY"),
	scaleTo = 10000,
	verbose = TRUE,
	threads = getArchRThreads(),
	parallelParam = NULL,
	force = TRUE,
	logFile = createLogFile("addGeneExpressionMatrix",  logDir=dir_log)
  ) # addGeneExpressionMatrix


  # Filter Cells
  cat(sprintf("\tapply additional filter with TSSEnrichment > 6 & nFrags > 2500\n"))
  # TSSEnrichment log10(6) = 0.7781513, usually my choice was higher than this
  # nFrags log10(2500) = 3.39794, usually my choice was higher than this 
  proj.filter <- proj.filter[proj.filter$TSSEnrichment > 6 & proj.filter$nFrags > 2500 & !is.na(proj.filter$Gex_nUMI)]

  for (sample in sample_names) {
	cells.atac.sample <- proj.filter$cellNames[(proj.filter$Sample == sample)]
	# Filter1: proj.filter$TSSEnrichment > 6 & proj.filter$nFrags > 2500 & !is.na(proj.filter$Gex_nUMI)
 	df_output["addFilter1NumCells", sample] <- length(cells.atac.sample)
  } # for

  cat(sprintf("\t\t# of cells: %d\n", nCells(proj.filter)))

  





  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("Perform LSI reduction for RNA data only\n"))
  cat(sprintf("\taddIterativeLSI for LSI_RNA, lsi_iter_rna=%d, lsi_varfeatures_rna=%d\n", args$lsi_iter_rna, args$lsi_varfeatures_rna))

  # https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
  # https://www.archrproject.com/reference/addIterativeLSI.html
  proj.filter <- addIterativeLSI(
	ArchRProj = proj.filter, 
	useMatrix = "GeneExpressionMatrix", 
	name = "LSI_RNA",
	iterations = args$lsi_iter_rna, # default=2,  The number of LSI iterations to perform.
	#clusterParams = list( resolution = args$seurat_resolution, sampleCells = 10000, n.start = 10),
	clusterParams = list(resolution = c(2.0), sampleCells = 10000, maxClusters = 6, n.start = 10), # consistency with LSI_ATAC
	firstSelection = "variable", # default="top"
	depthCol = "Gex_nUMI", # default="nFrags"
	varFeatures = args$lsi_varfeatures_rna, # default=25000,  The number of N variable features to use for LSI. The top N features will be used based on the selectionMethod.
	dimsToUse = dimsToUse,
	LSIMethod = 2,  # default=2,  A number or string indicating the order of operations in the TF-IDF normalization. Possible values are: 1 or "tf-logidf", 2 or "log(tf-idf)", and 3 or "logtf-logidf".
	scaleDims = TRUE,
	corCutOff = 0.75,
	binarize = FALSE, # default=TRUE
	outlierQuantiles = c(0.02, 0.98),
	filterBias = TRUE,
	sampleCellsPre = 10000,
	projectCellsPre = FALSE,
	sampleCellsFinal = NULL,
	selectionMethod = "var",
	scaleTo = 10000,
	totalFeatures = 5e+05,
	filterQuantile = 0.995,
	excludeChr = c(),
	saveIterations = FALSE, # default=TRUE
	#UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE, fast_sgd = TRUE),
	UMAPParams = list(n_neighbors=umap_n_neighbors, min_dist=umap_min_dist, metric=umap_metric, verbose =FALSE, fast_sgd = TRUE),
	nPlot = 10000,
	#outDir = getOutputDirectory(ArchRProj),
	outDir = dir_archr_output,
	threads = getArchRThreads(),
	seed = 1,
	verbose = TRUE,
	force = TRUE, # default=FALSE
	logFile = createLogFile("addIterativeLSI_RNA", logDir=dir_log)
  ) # addIterativeLSI


} # if f_multiome








# Perform LSI reduction and clustering with ATAC data only
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("Perform LSI reduction and clustering with ATAC data only\n"))


####################################################
# Add LSI dimreduc
cat(sprintf("\taddIterativeLSI for LSI_ATAC, lsi_iter_atac=%d, lsi_varfeatures_atac=%d\n", args$lsi_iter_atac, args$lsi_varfeatures_atac))

# https://www.archrproject.com/reference/addIterativeLSI.html
proj.filter <- addIterativeLSI(
  ArchRProj = proj.filter,
  useMatrix = "TileMatrix",
  name = "LSI_ATAC",
  iterations = args$lsi_iter_atac, # default=2,  The number of LSI iterations to perform.
  #clusterParams = list( resolution = c(0.2), sampleCells = 10000, n.start = 10 ), # default for AachR v0.9.5
  clusterParams = list(resolution = c(2.0), sampleCells = 10000, maxClusters = 6, n.start = 10), # default for ArchR v1.0.1  	A list of Additional parameters to be passed to addClusters() for clustering within each iteration. These params can be constant across each iteration, or specified for each iteration individually. Thus each param must be of length == 1 or the total number of iterations - 1. PLEASE NOTE - We have updated these params to resolution=2 and maxClusters=6! To use previous settings use resolution=0.2 and maxClusters=NULL.
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = args$lsi_varfeatures_atac, # default=25000,  The number of N variable features to use for LSI. The top N features will be used based on the selectionMethod.
  dimsToUse = dimsToUse, # default=1:30
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75, # default=0.75
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 10000,
  projectCellsPre = FALSE,
  sampleCellsFinal = NULL,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 5e+05,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = FALSE, # default=TRUE
  #UMAPParams = list(n_neighbors=30, min_dist=0.3, metric="cosine", verbose =FALSE), # M. Regner for ArchR v0.9.5
  #UMAPParams = list(n_neighbors=40, min_dist=0.4, metric="cosine", verbose =FALSE, fast_sgd=TRUE), # default for ArchR v1.0.1
  UMAPParams = list(n_neighbors=umap_n_neighbors, min_dist=umap_min_dist, metric=umap_metric, verbose =FALSE, fast_sgd = TRUE),
  nPlot = 10000,
  outDir = dir_archr_output, # default=getOutputDirectory(proj.filter)
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  force = TRUE, # default=FALSE
  logFile = createLogFile("addIterativeLSI_ATAC", logDir=dir_log)
) # addIterativeLSI



reducedDims <- "LSI_ATAC"
embedding_name <- "UMAP"




if (f_multiome) {

	# https://www.archrproject.com/reference/addCombinedDims.html
	proj.filter <- addCombinedDims(proj.filter, 
			name =  "LSI_Combined",
			reducedDims = c("LSI_ATAC", "LSI_RNA"),
			dimWeights = NULL,
			dimsToUse = NULL,
			scaleDims = NULL,
			corCutOff = 0.75 # A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
		) # addCombinedDims


	reducedDims <- "LSI_Combined"
	embedding_name <- "UMAP"

} # if f_multiome





if (args$harmony_theta[1] >= 0) {

  # harmony
  harmony_vars <- args$harmony_vars
  if (length(harmony_vars) == 0) {
        harmony_vars <- "Sample"
  }

  harmony_vars <- harmony_vars[harmony_vars %in% colnames(rna@meta.data)]

  df_info <- as.data.frame(matrix("nonspecified", nrow=nCells(proj.filter), ncol=length(harmony_vars), dimnames=list(proj.filter$cellNames, harmony_vars)))
  for (sample in levels(factor(proj.filter$Sample))) {
    rna.sub <- rna[,rna$Sample == sample]
    idxcell <- BiocGenerics::which(proj.filter$Sample == sample)
    for (var in harmony_vars) {
	df_info[idxcell,var] <- rna.sub@meta.data[1,var]
    }
  } # for
  remove(list=c("rna.sub"))

  # https://www.archrproject.com/reference/addCellColData.html
  cat(sprintf("\taddCellColData\n"))
  for (var in harmony_vars) {
	proj.filter <- addCellColData(proj.filter,
			data=df_info[,var],
			name=var,
			cells=rownames(df_info),
			force=TRUE)	
  } # for

  # set.seed for RunHarmony()  https://github.com/immunogenomics/harmony/issues/13
  set.seed(args$seed.harmony)

  # https://github.com/GreenleafLab/ArchR/blob/master/R/Harmony.R
  # https://github.com/immunogenomics/harmony/blob/master/man/HarmonyMatrix.Rd
  cat(sprintf("\taddHarmony\n"))
  cat(sprintf("\t\tharmony_vars=%s\n", paste(harmony_vars, collapse=", ")))
  cat(sprintf("\t\tharmony_theta=%s\n", paste(args$harmony_theta, collapse=", ")))
  cat(sprintf("\t\tharmony_lambda=%s\n", paste(args$harmony_lambda, collapse=", ")))
  # https://www.archrproject.com/reference/addHarmony.html
  proj.filter <- addHarmony(
    ArchRProj = proj.filter,
    reducedDims = reducedDims, # The name of the reducedDims object (i.e. "IterativeLSI") to retrieve from the designated ArchRProject.
    dimsToUse = dimsToUse, # A vector containing the dimensions from the reducedDims object to use in clustering.
    scaleDims = NULL, # A boolean that indicates whether to z-score the reduced dimensions for each cell. This is useful forminimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. If set to NULL this will scale the dimensions based on the value of scaleDims when the reducedDims were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
    corCutOff = 0.75, # defalt=0.75, A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
    name = "Harmony",
    groupBy = harmony_vars,
    verbose = FALSE,
    force = TRUE,
    # begin of harmonyParams for harmony::HarmonyMatrix()
    theta = args$harmony_theta, # Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.
    lambda = args$harmony_lambda, # Ridge regression penalty parameter. Specify for each variable in vars_use. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.
    sigma = 0.1, # Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering.
    nclust = NULL, # Number of clusters in model. nclust=1 equivalent to simple linear regression.
    tau = 0, # Protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster.
    block.size = 0.05, # What proportion of cells to update during clustering. Between 0 to 1, default 0.05. Larger values may be faster but less accurate
    max.iter.harmony = 10,
    max.iter.cluster = 200, # default=200 in harmony::HarmonyMatrix()
    epsilon.cluster = 1e-05,
    epsilon.harmony = 1e-04
    # end of harmonyParams for harmony::HarmonyMatrix()
  ) # addHarmony

  reducedDims <- "Harmony"

} # if

cat(sprintf("\t\treducedDims=%s\n", reducedDims))







# addClusters

cat(sprintf("\taddClusters\n"))
cat(sprintf("\t\tseurat_resolution=%.1f\n", args$seurat_resolution))

# https://www.archrproject.com/reference/addClusters.html
proj.filter <- addClusters(
  input = proj.filter,
  reducedDims = reducedDims,
  name = "ATAC_clusters",
  sampleCells = NULL,
  seed = 1,
  method = "Seurat", # Supported methods are "Seurat" and "Scran".
  dimsToUse = dimsToUse,
  scaleDims = NULL, # A boolean value that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. If set to NULL this will scale the dimensions based on the value of scaleDims when the reducedDims were originally created during dimensionality reduction. This idea was introduced by Timothy Stuart.
  corCutOff = 0.75,
  knnAssign = 10,
  nOutlier = 5, # The minimum number of cells required for a group of cells to be called as a cluster. If a group of cells does not reach this threshold, then the cells will be considered outliers and assigned to nearby clusters.
  maxClusters = 25, # The maximum number of clusters to be called. If the number exceeds this the clusters are merged unbiasedly using hclust and cutree. This is useful for contraining the cluster calls to be reasonable if they are converging on large numbers. Useful in iterativeLSI as well for initial iteration. Default is set to 25.
  testBias = TRUE,
  filterBias = FALSE,
  biasClusters = 0.01,
  biasCol = "nFrags",
  biasVals = NULL,
  biasQuantiles = c(0.05, 0.95),
  biasEnrich = 10,
  biasProportion = 0.5,
  biasPval = 0.05,
  nPerm = 500,
  prefix = "C",
  ArchRProj = NULL,
  verbose = TRUE,
  tstart = NULL,
  force = TRUE,
  # for Seurat
  resolution = args$seurat_resolution,  # for Seurat FindClusters(rna, resolution=0.8) Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.  https://satijalab.org/seurat/reference/findclusters
  # for Scran
  # logFile
  logFile = createLogFile("addClusters", logDir=dir_log)
) # addClusters






# Add UMAP based on LSI dims
cat(sprintf("\taddUMAP\n"))
# https://www.archrproject.com/reference/addUMAP.html
proj.filter <- addUMAP(proj.filter,
  reducedDims = reducedDims,
  name = embedding_name,
  nNeighbors = umap_n_neighbors,
  minDist = umap_min_dist,
  metric = umap_metric, # default "cosine"
  dimsToUse = dimsToUse,
  scaleDims = NULL,
  corCutOff = 0.75,
  sampleCells = NULL,
  outlierQuantile = 0.9,
  saveModel = TRUE,
  verbose = TRUE,
  seed = 1,
  force = TRUE,
  threads =  getArchRThreads()
) # addUMAP







if (f_multiome) {

	cat(sprintf("\taddUMAP with LSI_RNA\n"))
	# https://www.archrproject.com/reference/addUMAP.html
	proj.filter <- addUMAP(proj.filter,
	  reducedDims = "LSI_RNA",
	  name = "UMAP_RNA",
	  nNeighbors=umap_n_neighbors,
	  minDist=umap_min_dist,
	  metric=umap_metric, # default "cosine"
	  dimsToUse=dimsToUse,
	  scaleDims = NULL,
	  corCutOff = 0.75,
	  sampleCells = NULL,
	  outlierQuantile = 0.9,
	  saveModel = TRUE,
	  verbose = TRUE,
	  seed = 1,
	  force = TRUE,
	  threads =  getArchRThreads()
	) # addUMAP


	cat(sprintf("\taddUMAP with LSI_ATAC\n"))
	# https://www.archrproject.com/reference/addUMAP.html
	proj.filter <- addUMAP(proj.filter,
	  reducedDims = "LSI_ATAC",
	  name = "UMAP_ATAC",
	  nNeighbors=umap_n_neighbors,
	  minDist=umap_min_dist,
	  metric=umap_metric, # default "cosine"
	  dimsToUse=dimsToUse,
	  scaleDims = NULL,
	  corCutOff = 0.75,
	  sampleCells = NULL,
	  outlierQuantile = 0.9,
	  saveModel = TRUE,
	  verbose = TRUE,
	  seed = 1,
	  force = TRUE,
	  threads =  getArchRThreads()
	) # addUMAP


} # if f_multiome



















######################################################################
# PART 2: Gene scores and marker genes
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 2: Gene scores and marker genes\n"))
cat(sprintf("------------------------------------\n\n"))
######################################################################







# Estimate gene activity in ATAC data and perform cluster type annotation: 
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("Estimate gene activity in ATAC data and perform cluster type annotation\n"))

# Add Gene activity matrix using ArchR model 
cat(sprintf("\taddGeneScoreMatrix\n"))

proj.filter <- addGeneScoreMatrix(proj.filter,
  genes = getGenes(proj.filter),
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "ArchRGeneScore",
  extendUpstream = c(1000, 1e+05),
  extendDownstream = c(1000, 1e+05),
  geneUpstream = 5000,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  tileSize = 500,
  ceiling = 4,
  geneScaleFactor = 5,
  scaleTo = 10000,
  excludeChr = c("chrY", "chrM"),
  blacklist = getBlacklist(proj.filter),
  threads = getArchRThreads(),
  parallelParam = NULL,
  #subThreading = TRUE,
  subThreading = FALSE,
  force = TRUE,
  logFile = createLogFile("addGeneScoreMatrix", logDir=dir_log)
)


#genescore.mat <- getMatrixFromProject(
#  ArchRProj = proj.filter,
#  useMatrix = "ArchRGeneScore",
#  useSeqnames = NULL,
#  verbose = TRUE,
#  binarize = FALSE,
#  threads = getArchRThreads(),
#  logFile = createLogFile("getMatrixFromProject", logDir=dir_log)
#)



# Add Gene activity matrix using Signac model 
# if ("package:ArchR.dpeerlab" %in% (search())) {
#	gene.ranges <- ArchR.dpeerlab::geneAnnoHg38
# } else {
#	gene.ranges <- ArchR::geneAnnoHg38
# }
# gene.ranges <- gene.ranges$genes
# genebodyandpromoter.coords <- Signac::Extend(x = gene.ranges, upstream = 2000, downstream = 0)
# genebodyandpromoter.coords$name <- genebodyandpromoter.coords$symbol
# proj.filter <- addFeatureMatrix(proj.filter,matrixName = "SignacGeneScore",features = genebodyandpromoter.coords,force = T)


# getAvailableMatrices
cat(sprintf("\tgetAvailableMatrices(proj.filter)\n"))
getAvailableMatrices(proj.filter)


if (args$f_save_tmp_rds) {
  # saveRDS
  filename <- sprintf("%s/%s_proj_lsi_and_genescores.rds", dir_rds, cancer_type)
  cat(sprintf("\tsave RDS: %s\n", filename))
  #saveRDS(proj.filter, "proj_LSI_AND_GeneScores.rds")
  saveRDS(proj.filter, filename)
} # if



# constrained integration to only align cells from the same patient tumor
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("Constrained integration to only align cells from the same patient tumor\n"))

groupList <- SimpleList()
for (sample in levels(factor(proj.filter$Sample))){
  
  rna.sub <- rna[,rna$Sample == sample]
  RNA.cells <- colnames(rna.sub)

  #idxSample <- BiocGenerics::which(proj.filter$Sample == sample)
  #ATAC.cells <- proj.filter$cellNames[idxSample]

  proj.meta <- as.data.frame(proj.filter@cellColData)
  proj.meta.sub <- proj.meta[proj.meta$Sample == sample,]
  ATAC.cells <- rownames(proj.meta.sub)

  groupList[[sample]] <- SimpleList(
    ATAC = ATAC.cells,
    RNA = RNA.cells
  )
} # for



####################################################
cat(sprintf("\taddGeneIntegrationMatrix\n"))



# Perform Seurat v3 label transfer between RNA/ATAC
# https://www.archrproject.com/reference/addGeneIntegrationMatrix.html



# proj.filter <- addGeneIntegrationMatrix(
#   ArchRProj = proj.filter,
#   useMatrix = "SignacGeneScore",
#   matrixName = "GeneIntegrationMatrix_Signac",
#   reducedDims = reducedDims,
#   seRNA = rna,
#   groupList = groupList,
#   addToArrow = T,
#   force= TRUE,
#   groupRNA = col_cluster_types,
#   nameCell = "predictedCell_Signac",
#   nameGroup = "predictedGroup_Signac",
#   nameScore = "predictedScore_Signac",
#   plotUMAP = F,
#   useImputation = F,
#   transferParams = list(dims = dimsToUse)
# )

#proj.filter <- addGeneIntegrationMatrix(
#  ArchRProj = proj.filter,
#  useMatrix = "ArchRGeneScore",
#  matrixName = "GeneIntegrationMatrix_ArchR",
#  reducedDims = reducedDims,
#  seRNA = rna,
#  groupList = groupList,
#  addToArrow = T,
#  force= TRUE,
#  groupRNA = col_cluster_types,
#  nameCell = "predictedCell_ArchR",
#  nameGroup = "predictedGroup_ArchR",
#  nameScore = "predictedScore_ArchR",
#  plotUMAP = F,
#  useImputation = F,
#  transferParams = list(dims = dismToUse)
#)

proj.filter <- addGeneIntegrationMatrix(
  ArchRProj = proj.filter,
  useMatrix = "ArchRGeneScore",
  matrixName = "GeneIntegrationMatrix_ArchR",
  reducedDims = reducedDims,
  seRNA = rna,
  groupATAC = NULL,
  groupRNA = col_cluster_types,
  groupList = groupList,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  dimsToUse = dimsToUse,
  scaleDims = NULL,
  corCutOff = 0.75,
  plotUMAP = FALSE,
  #UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
  UMAPParams = list(n_neighbors=umap_n_neighbors, min_dist=umap_min_dist, metric=umap_metric, verbose =FALSE),
  nGenes = 2000,
  useImputation = FALSE, # default=TRUE, A boolean value indicating whether to use imputation for creating the Gene Score Matrix prior to integration.
  reduction = "cca", # To assess the performance of each model, we used canonical correlation analysis to integrate scATAC-seq and scRNA-seq data from the same sample types and then compared the linked gene expression from scRNA-seq to the inferred gene scores from scATAC-seq  https://www.nature.com/articles/s41588-021-00790-6
  addToArrow = TRUE,
  scaleTo = 10000,
  genesUse = NULL,
  nameCell = "predictedCell_ArchR",
  nameGroup = "predictedGroup_ArchR",
  nameScore = "predictedScore_ArchR", # A column name to add to cellColData for the predicted scRNA-seq score in the specified ArchRProject. These scores represent the assignment accuracy of the group in the RNA cells. Lower scores represent ambiguous predictions and higher scores represent precise predictions.
  #transferParams = list(dims = dimsToUse), #  M. Regner for ArchR v0.9.5
  transferParams = list(dims = dimsToUse), # dims=NULL, Set of dimensions to use in the anchor weighting procedure. If NULL, the same dimensions that were used to find anchors will be used for weighting.  https://satijalab.org/seurat/reference/transferdata
  threads = getArchRThreads(),
  verbose = TRUE,
  force = TRUE,
  logFile = createLogFile("addGeneIntegrationMatrix", logDir=dir_log)
) # addGeneIntegrationMatrix









# number of samples for each group
dt <- summarize_meta_data(proj.filter, "predictedGroup_ArchR", "Sample")














# getAvailableMatrices
cat(sprintf("\tgetAvailableMatrices(proj.filter)\n"))
getAvailableMatrices(proj.filter)


if (args$f_save_tmp_rds) {
  # saveRDS
  filename <- sprintf("%s/%s_proj_lsi_genescores_annotations_int.rds", dir_rds, cancer_type)
  cat(sprintf("\tsave RDS: %s\n", filename))
  #saveRDS(proj.filter, "proj_LSI_GeneScores_Annotations_Int.rds")
  saveRDS(proj.filter, filename)
} # if









######################################################
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\tPlotting RNA/ATAC by sample, by cluster, by predicted label\n"))

colorBy <- "cellColData"

### make embedding highlighting by 1) Predicted group ArchR 2) Predicted group Signac 3) Sample 4) ATAC-only clusters 


# atac.archr.emb
atac.archr <- plotEmbedding(proj.filter,
  embedding = embedding_name,
  colorBy = colorBy,
  name = "predictedGroup_ArchR",
  log2Norm = NULL,
  imputeWeights = if (!grepl("coldata", tolower(colorBy[1]))) getImputeWeights(proj.filter),
  pal = NULL,
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding", logDir=dir_log)
)

atac.archr.emb <- as.data.frame(atac.archr$data)
atac.archr.emb$cluster.type.archr <- atac.archr.emb$color
atac.archr.emb$cluster.type.archr <- sub("-", ":", atac.archr.emb$cluster.type.archr)
atac.archr.emb$cluster.type.archr <- gsub(".*:", "", atac.archr.emb$cluster.type.archr)

#cat(sprintf("\thead(atac.archr.emb)\n"))
#head(atac.archr.emb)

atac.archr.emb$cluster.type.archr <- factor( atac.archr.emb$cluster.type.archr, levels = levels(as.factor(rna@meta.data[, col_cluster_types])) )

cat(sprintf("\thead(atac.archr.emb)\n"))
head(atac.archr.emb)





# atac.emb.sample
atac <- plotEmbedding(proj.filter,
  embedding = embedding_name,
  colorBy = colorBy,
  log2Norm = NULL,
  imputeWeights = if (!grepl("coldata", tolower(colorBy[1]))) getImputeWeights(proj.filter),
  pal = NULL,
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding", logDir=dir_log)
)

atac.emb.sample <- as.data.frame(atac$data)
atac.emb.sample$sample <- atac.emb.sample$color
atac.emb.sample$sample <- sub("-", ":", atac.emb.sample$sample )
atac.emb.sample$sample  <- gsub(".*:", "", atac.emb.sample$sample )

cat(sprintf("\thead(atac.emb.sample)\n"))
head(atac.emb.sample)





# atac.emb.cluster
atac <- plotEmbedding(proj.filter,
  embedding = embedding_name,
  colorBy = colorBy,
  name = "ATAC_clusters",
  log2Norm = NULL,
  imputeWeights = if (!grepl("coldata", tolower(colorBy[1]))) getImputeWeights(proj.filter),
  pal = NULL,
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  discreteSet = NULL,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding", logDir=dir_log)
)

atac.emb.cluster <- as.data.frame(atac$data)
atac.emb.cluster$sample <- atac.emb.cluster$color
atac.emb.cluster$sample <- sub("-", ":", atac.emb.cluster$sample )
atac.emb.cluster$sample  <- gsub(".*:", "", atac.emb.cluster$sample )

cat(sprintf("\thead(atac.emb.cluster)\n"))
head(atac.emb.cluster)

atac.emb.all <- cbind(atac.archr.emb[,c(1:2,4)],
                      atac.emb.sample[,4],
                      atac.emb.cluster[,4])

atac.emb.all$plain <- "Plain"

colnames(atac.emb.all) <- c(
	"UMAP1", "UMAP2", "Predicted.Group.ArchR",
        "Sample", "ATAC_clusters",
	"Blank")

cat(sprintf("\thead(atac.emb.all)\n"))
head(atac.emb.all)




gg <- ggplot(atac.emb.all,
  aes_string(x = "UMAP1", y="UMAP2", color = "Predicted.Group.ArchR"))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scATAC-seq: Predicted Group ArchR")+
  theme(plot.title = element_text(face = "bold"))+
  xlab(TeX(r'($UMAP_1$)'))+
  ylab(TeX(r'($UMAP_2$)'))+
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(sprintf("%s/%s_predictedgroup_archr_atac.pdf", dir_pdf, cancer_type), width = 8, height = 6, plot=gg)
ggsave(sprintf("%s/%s_predictedgroup_archr_atac.png", dir_png, cancer_type), width = 8, height = 6, plot=gg)




gg <- ggplot(atac.emb.all,
  aes_string(x = "UMAP1", y="UMAP2", color = "Sample"))+
  geom_point(size = .1)+
  theme_classic()+
  ggtitle(paste0("scATAC-seq: Sample"))+
  theme(plot.title = element_text(face = "bold"))+
  xlab(TeX(r'($UMAP_1$)'))+
  ylab(TeX(r'($UMAP_2$)'))+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)
  #guides(colour = guide_legend(override.aes = list(size=3)))+

ggsave(sprintf("%s/%s_sample_atac.pdf", dir_pdf, cancer_type), width = 8, height = 6, plot=gg)
ggsave(sprintf("%s/%s_sample_atac.png", dir_png, cancer_type), width = 8, height = 6, plot=gg)


  


prediction.scores <- data.frame(ArchR = proj.filter$predictedScore_ArchR)
var.list <- colnames(prediction.scores)

for (i in 1:length(var.list)) {
  
  gg <- ggplot(prediction.scores, aes_string(x = var.list[i]))+
    geom_histogram(binwidth = 0.025,fill="#000000", color="#e9ecef", alpha=0.9)+
    theme_classic()

  ggsave(sprintf("%s/%s_%s_atac.pdf", dir_pdf, cancer_type, var.list[i]), width = 8, height = 6, plot=gg)
  ggsave(sprintf("%s/%s_%s_atac.png", dir_png, cancer_type, var.list[i]), width = 8, height = 6, plot=gg)

} # for









# Plot matching scRNA-seq plots:
cat(sprintf("\tPlotting matching scRNA-seq plots\n"))
######################################################################
#rna.emb <- as.data.frame(rna@reductions$umap@cell.embeddings)
rna.emb <- as.data.frame(Embeddings(rna, reduction = str_umap_reduction))
rna.emb[, col_cluster_types] <- as.factor(rna@meta.data[, col_cluster_types])
rna.emb$sample <- rna$Sample

rna.emb[, col_cluster_types] <- factor(rna.emb[, col_cluster_types], levels = levels(atac.emb.all$Predicted.Group.ArchR))
rna.cell.plot <- ggplot(rna.emb,
  aes_string(x = col_umap1, y = col_umap2, color = col_cluster_types)) +
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scRNAseq")+
  theme(plot.title = element_text(face = "bold"))+
  xlab(TeX(r'($UMAP_1$)'))+
  ylab(TeX(r'($UMAP_2$)'))+
  theme(legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(sprintf("%s/%s_rna_all_labels.pdf", dir_pdf, cancer_type), width = 8, height = 6, plot=rna.cell.plot)
ggsave(sprintf("%s/%s_rna_all_labels.png", dir_png, cancer_type), width = 8, height = 6, plot=rna.cell.plot)
  


rna.emb$sample <- factor(rna.emb$sample, levels = levels(factor(atac.emb.all$Sample)))
rna.sample.plot <- ggplot(rna.emb,
  aes_string(x = col_umap1, y = col_umap2, color = "sample")) +
  geom_point(size = .1)+
  theme_classic()+
  ggtitle("scRNAseq")+
  theme(plot.title = element_text(face = "bold"))+
  xlab(TeX(r'($UMAP_1$)'))+
  ylab(TeX(r'($UMAP_2$)'))+
  theme(legend.key.size = unit(0.2, "cm"))+
  scale_color_manual(values = sampleColors)+
  guides(colour = guide_legend(override.aes = list(size=3)))

ggsave(sprintf("%s/%s_rna_all_samples.pdf", dir_pdf, cancer_type), width = 8, height = 6, plot=rna.sample.plot)
ggsave(sprintf("%s/%s_rna_all_samples.png", dir_png, cancer_type), width = 8, height = 6, plot=rna.sample.plot)


















# DEGs using ATAC labels with ArchRGeneScore matrix
cat(sprintf("\tgetMarkerFeatures ATAC_clusters\n"))
# markersGS.archr, useMatrix = "ArchRGeneScore", groupBy = "ATAC_clusters"
markersGS.archr <- getMarkerFeatures(
  ArchRProj = proj.filter,
  useMatrix = "ArchRGeneScore",
  groupBy = "ATAC_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures.ArchRGeneScore.ATAC_clusters", logDir=dir_log)
)


cat(sprintf("\tplotMarkerHeatmap\n"))
#heatmapGS.archr <- markerHeatmap(
heatmapGS.archr <- plotMarkerHeatmap(
  seMarker = markersGS.archr,
  #cutOff = "FDR <= 0.01 & Log2FC >= 0.5", # default
  cutOff = str_cutOff_markersGS_for_heatmap,
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2,2),
  #pal = NULL # default
  pal = viridis(n=256),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = NULL,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap.markersGS.archr", logDir=dir_log)
)


cat(sprintf("\tComplexHeatmap\n"))
ComplexHeatmap::draw(heatmapGS.archr, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.archr, name = "genescores-marker-heatmap_archr", width = 8, height = 32, ArchRProj = proj.filter, addDOC = FALSE) # output directory: All/Plots










# proj.filter.score_archr
# DEGs using predicted labels (removing small groups)
idxSample <- BiocGenerics::which(proj.filter$predictedScore_ArchR > 0.5)
cellsSample <- proj.filter$cellNames[idxSample]
proj.filter.score_archr <- proj.filter[cellsSample, ]

popular.groups <- summary(factor(proj.filter.score_archr$predictedGroup_ArchR))
popular.groups <- popular.groups[popular.groups > 10]
proj.filter.score_archr$Mode.Label <- ifelse(proj.filter.score_archr$predictedGroup_ArchR %in% names(popular.groups),TRUE,FALSE)

idxSample <- BiocGenerics::which(proj.filter.score_archr$Mode.Label == TRUE)
cellsSample <- proj.filter.score_archr$cellNames[idxSample]
proj.filter.score_archr <- proj.filter.score_archr[cellsSample, ]




# DEGs using predicted labels
cat(sprintf("\tgetMarkerFeatures predictedGroup_ArchR\n"))
# markersGS.archr.pred, useMatrix = "ArchRGeneScore", groupBy = "predictedGroup_ArchR"
markersGS.archr.pred <- getMarkerFeatures(
  ArchRProj = proj.filter.score_archr,
  useMatrix = "ArchRGeneScore",
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures.ArchRGeneScore.predictedGroup_ArchR", logDir=dir_log)
) # getMarkerFeatures

# save RDS
#saveRDS(markersGS.archr.pred, sprintf("%s/%s_markersgs_proj.filter.score_archr_predictedgroup_archr.rds", dir_rds, cancer_type))






cat(sprintf("\tplotMarkerHeatmap\n"))
#heatmapGS.archr.pred <- markerHeatmap(
heatmapGS.archr.pred <- plotMarkerHeatmap(
  seMarker = markersGS.archr.pred,
  #cutOff = "FDR <= 0.01 & Log2FC >= 0.5", # default
  cutOff = str_cutOff_markersGS_for_heatmap,
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2,2),
  #pal = NULL # default
  pal = viridis(n=256),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = NULL,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap.markersGS.archr.pred", logDir=dir_log)
)


cat(sprintf("\tComplexHeatmap\n"))
ComplexHeatmap::draw(heatmapGS.archr.pred, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.archr.pred, name = "genescores-marker-heatmap_archr_pred", width = 8, height = 32, ArchRProj = proj.filter.score_archr, addDOC = FALSE) # output directory: All/Plots















#################################################################
# 
# 
# # Signac
# #########################################################
# # DEGs using ATAC labels
# markersGS.signac <- getMarkerFeatures(
#   ArchRProj = proj.filter,
#   useMatrix = "SignacGeneScore",
#   groupBy = "ATAC_clusters",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# 
# heatmapGS.signac <- markerHeatmap(
#   seMarker = markersGS.signac,
#   cutOff = str_cutOff_markersGS_for_heatmap,
#   labelMarkers =NULL,
#   transpose = F,
#   pal =  viridis(n=256),
#   limits = c(-2,2)
# )
# ComplexHeatmap::draw(heatmapGS.signac, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapGS.signac , name = "GeneScores-Marker-Heatmap_signac", width = 8, height = 6, ArchRProj = proj.filter, addDOC = FALSE)
# 
# 
# 
# 
# proj.filter.score_signac
# # DEGs using predicted labels (removing small groups)
# idxSample <- BiocGenerics::which(proj.filter$predictedScore_Signac > 0.5)
# cellsSample <- proj.filter$cellNames[idxSample]
# proj.filter.score_signac <- proj.filter[cellsSample, ]
# 
# popular.groups <- summary(factor(proj.filter.score_signac$predictedGroup_Signac))
# popular.groups <- popular.groups[popular.groups > 10]
# proj.filter.score_signac$Mode.Label <- ifelse(proj.filter.score_signac$predictedGroup_Signac %in% names(popular.groups),TRUE,FALSE)
# 
# idxSample <- BiocGenerics::which(proj.filter.score_signac$Mode.Label == TRUE)
# cellsSample <- proj.filter.score_signac$cellNames[idxSample]
# proj.filter.score_signac <- proj.filter.score_signac[cellsSample, ]
# 
# # DEGs using predicted labels
# markersGS.signac.pred <- getMarkerFeatures(
#   ArchRProj = proj.filter.score_signac,
#   useMatrix = "SignacGeneScore",
#   groupBy = "predictedGroup_Signac",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# heatmapGS.signac.pred <- markerHeatmap(
#   seMarker = markersGS.signac.pred,
#   cutOff = str_cutOff_markersGS_for_heatmap,
#   labelMarkers =NULL,
#   transpose = F,
#   pal =  viridis(n=256),
#   limits = c(-2,2)
# )
# ComplexHeatmap::draw(heatmapGS.signac.pred, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# plotPDF(heatmapGS.signac.pred, name = "GeneScores-Marker-Heatmap_Signac_pred", width = 8, height = 6, ArchRProj = proj.filter.score_signac, addDOC = FALSE)












######################################################################
# PART 3: Pseudo-bulk replicates and calling peaks
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 3: Pseudo-bulk replicates and calling peaks\n"))
cat(sprintf("------------------------------------\n\n"))
######################################################################




############################
# ATAC clusters

# peak analysis with all ATAC cells
#proj.atac <- proj.filter

# peak analysis with selected cells
proj.atac <- proj.filter.score_archr
#proj.atac <- proj.filter.score_signac

cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\tATAC clusters\n"))

cat(sprintf("\taddGroupCoverages\n"))

if ("package:ArchR.dpeerlab" %in% (search())) {
  proj.atac <- addGroupCoverages(ArchRProj = proj.atac,
		groupBy = "ATAC_clusters",
		useLabels = TRUE,
		minCells = 40,
		maxCells = 500,
		maxFragments = 25 * 10^6,
		minReplicates = 2,
		maxReplicates = 5,
		sampleRatio = 0.8,
		kmerLength = 6,
		maxFragmentLength = 147,
		threads = getArchRThreads(),
		returnGroups = FALSE,
		parallelParam = NULL,
		force = TRUE, # force=T
		verbose = TRUE,
		logFile = createLogFile("addGroupCoverages.ATAC_clusters", logDir=dir_log)
)
} else {
  # https://www.archrproject.com/reference/addGroupCoverages.html
  proj.atac <- addGroupCoverages(ArchRProj = proj.atac,
		groupBy = "ATAC_clusters",
		useLabels = TRUE,
		minCells = 40,
		maxCells = 500,
		maxFragments = 25 * 10^6,
		minReplicates = 2,
		maxReplicates = 5,
		sampleRatio = 0.8,
		kmerLength = 6,
		threads = getArchRThreads(),
		returnGroups = FALSE,
		parallelParam = NULL,
		force = TRUE, # force=T
		verbose = TRUE,
		logFile = createLogFile("addGroupCoverages.ATAC_clusters", logDir=dir_log)
  )
} # if






peakMethod <- "Macs2"
pathToMacs2 <- findMacs2()

cat(sprintf("\taddReproduciblePeakSet\n"))
proj.atac <- addReproduciblePeakSet(
  ArchRProj = proj.atac,
  groupBy = "ATAC_clusters",
  peakMethod = peakMethod, # default="Macs2"
  reproducibility = "2", # default="2",  A string that indicates how peak reproducibility should be handled. This string is dynamic and can be a function of n where n is the number of samples being assessed. For example, reproducibility = "2" means at least 2 samples must have a peak call at this locus and reproducibility = "(n+1)/2" means that the majority of samples must have a peak call at this locus.
  peaksPerCell = 500, # default=500
  maxPeaks = args$max_peaks_per_group, # function default=150,000, script default=150,000, A numeric threshold for the maximum peaks to retain per group from groupBy in the union reproducible peak set.
  minCells = args$min_cells_to_call_peaks, # function default=25, script default=50, The minimum allowable number of unique cells that was used to create the coverage files on which peaks are called. This is important to allow for exclusion of pseudo-bulk replicates derived from very low cell numbers.
  excludeChr = c("chrM", "chrY"),
  #pathToMacs2 = if (tolower(peakMethod) == "macs2") findMacs2() else NULL,
  pathToMacs2 = pathToMacs2,
  genomeSize = NULL,
  shift = -75,
  extsize = 150, # The number of basepairs to extend the MACS2 fragment after shift has been applied. When combined with extsize this allows you to create proper fragments, centered at the Tn5 insertion site, for use with MACS2 (see MACS2 documentation).
  method = if (tolower(peakMethod) == "macs2") "q" else "p",
  cutOff = args$macs2_cutoff, # function default=0.1, script default=1e-7, The numeric significance cutOff for the testing method indicated by method (see MACS2 documentation).
  additionalParams = "--nomodel --nolambda",
  extendSummits = 250, # The number of basepairs to extend peak summits (in both directions) to obtain final fixed-width peaks. For example, extendSummits = 250 will create 501-bp fixed-width peaks from the 1-bp summits.
  promoterRegion = c(2000, 100), # default=c(2000, 100)
  genomeAnnotation = getGenomeAnnotation(proj.atac),
  geneAnnotation = getGeneAnnotation(proj.atac),
  plot = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addReproduciblePeakSet", logDir=dir_log)
) # addReproduciblePeakSet






cat(sprintf("\taddPeakMatrix\n"))
if ("package:ArchR.dpeerlab" %in% (search())) {
  proj.atac <- addPeakMatrix(proj.atac,
    ceiling = 4, # default=4, ceiling=10^9 for some pupurposes
    maxFragmentLength = 147, 
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    logFile = createLogFile("addPeakMatrix.proj", logDir=dir_log)
  )
} else {
  # https://www.archrproject.com/reference/addPeakMatrix.html
  proj.atac <- addPeakMatrix(proj.atac,
    ceiling = 4,
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    logFile = createLogFile("addPeakMatrix.proj", logDir=dir_log)
  )
} # if





cat(sprintf("\taddBgdPeaks\n"))
# https://www.archrproject.com/reference/addBgdPeaks.html
proj.atac <- addBgdPeaks(proj.atac,
	nIterations = 50,
	w = 0.1,
	binSize = 50,
	seed = 1,
	method = "chromVAR",
	outFile = file.path(dir_archr_output, "Background-Peaks.proj.atac.rds"),
	force = TRUE
  ) # addBgdPeaks



# markePeaks with PeakMatrix
cat(sprintf("\tgetMarkerFeatures\n"))
# markerPeaks, useMatrix = "PeakMatrix", groupBy = "ATAC_clusters"
markerPeaks <- getMarkerFeatures(
  ArchRProj = proj.atac,
  useMatrix = "PeakMatrix",
  groupBy = "ATAC_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures.PeakMatrix.ATAC_clusters", logDir=dir_log)
) # getMarkerFeatures





cat(sprintf("\tplotMarkerHeatmap\n"))
#heatmapPeaks<- markerHeatmap(
heatmapPeaks<- plotMarkerHeatmap(
  seMarker = markerPeaks,
  #cutOff = "FDR <= 0.01 & Log2FC >= 0.5", # default
  cutOff = str_cutOff_markerPeaks_for_heatmap,
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2,2),
  #pal = NULL # default
  pal = viridis(n=256),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = NULL,
  nLabel = 15,
  nPrint = 15, # If provided seMarker is from "GeneScoreMatrix" print the top n genes for each group based on how uniquely up-regulated the gene is.
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap.markerPeaks", logDir=dir_log)
)

cat(sprintf("\tComplexHeatmap\n"))
ComplexHeatmap::draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "markerpeaks_atac_clusters", width = 8, height = 32, ArchRProj = proj.atac, addDOC = FALSE) # output directory: All/Plots























############################
# ArchR predicted labels

cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\tArchR predicted labels\n"))

# proj.filter.score_archr
# DEGs using predicted labels (removing small groups)
idxSample <- BiocGenerics::which(proj.filter$predictedScore_ArchR >= 0.5)
cellsSample <- proj.filter$cellNames[idxSample]
proj.filter.score_archr <- proj.filter[cellsSample, ]

# select groups with # of cells > 10
popular.groups <- summary(factor(proj.filter.score_archr$predictedGroup_ArchR))
popular.groups <- popular.groups[popular.groups > 10]
proj.filter.score_archr$Mode.Label <- ifelse(proj.filter.score_archr$predictedGroup_ArchR %in% names(popular.groups), TRUE, FALSE)

idxSample <- BiocGenerics::which(proj.filter.score_archr$Mode.Label == TRUE)
cellsSample <- proj.filter.score_archr$cellNames[idxSample]
proj.filter.score_archr <- proj.filter.score_archr[cellsSample, ]






# # proj.filter.score_signac
# # DEGs using predicted labels (removing small groups)
# idxSample <- BiocGenerics::which(proj.filter$predictedScore_Signac > 0.5)
# cellsSample <- proj.filter$cellNames[idxSample]
# proj.filter.score_signac <- proj.filter[cellsSample, ]
#
# popular.groups <- summary(factor(proj.filter.score_signac$predictedGroup_Signac))
# popular.groups <- popular.groups[popular.groups > 10]
# proj.filter.score_signac$Mode.Label <- ifelse(proj.filter.score_signac$predictedGroup_Signac %in% names(popular.groups),TRUE,FALSE)
#
# idxSample <- BiocGenerics::which(proj.filter.score_signac$Mode.Label == TRUE)
# cellsSample <- proj.filter.score_signac$cellNames[idxSample]
# proj.filter.score_signac <- proj.filter.score_signac[cellsSample, ]
#





# peak analysis with selected cells
proj.archr <- proj.filter.score_archr
#proj.archr <- proj.filter.score_signac

for (sample in sample_names) {
	cells.atac.sample <- proj.archr$cellNames[(proj.archr$Sample == sample)]
	# Filter2: predictedScore_ArchR >= 0.5, select groups with # of cells > 10
 	df_output["PostFilter2NumCells", sample] <- length(cells.atac.sample)
} # for



# proj.archr
cat(sprintf("\taddGroupCoverages\n"))
if ("package:ArchR.dpeerlab" %in% (search())) {
  proj.archr <- addGroupCoverages(ArchRProj = proj.archr,
		groupBy = "predictedGroup_ArchR",
		useLabels = TRUE,
		minCells = 40,
		maxCells = 500,
		maxFragments = 25 * 10^6,
		minReplicates = 2,
		maxReplicates = 5,
		sampleRatio = 0.8,
		kmerLength = 6,
		maxFragmentLength = 147,
		threads = getArchRThreads(),
		returnGroups = FALSE,
		parallelParam = NULL,
		force = TRUE, # force=T
		verbose = TRUE,
		logFile = createLogFile("addGroupCoverages.predictedGroup_ArchR", logDir=dir_log)
  )
} else {
  # https://www.archrproject.com/reference/addGroupCoverages.html
  proj.archr <- addGroupCoverages(ArchRProj = proj.archr,
		groupBy = "predictedGroup_ArchR",
		useLabels = TRUE,
		minCells = 40,
		maxCells = 500,
		maxFragments = 25 * 10^6,
		minReplicates = 2,
		maxReplicates = 5,
		sampleRatio = 0.8,
		kmerLength = 6,
		threads = getArchRThreads(),
		returnGroups = FALSE,
		parallelParam = NULL,
		force = TRUE, # force=T
		verbose = TRUE,
		logFile = createLogFile("addGroupCoverages.predictedGroup_ArchR", logDir=dir_log)
  )
} # if







cat(sprintf("\taddReproduciblePeakSet\n"))
proj.archr <- addReproduciblePeakSet(
  ArchRProj = proj.archr,
  groupBy = "predictedGroup_ArchR",
  peakMethod = peakMethod,
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = args$max_peaks_per_group, # function default=150,000, script default=150,000
  minCells = args$min_cells_to_call_peaks, # function default=25, script default=50
  excludeChr = c("chrM", "chrY"),
  #pathToMacs2 = if (tolower(peakMethod) == "macs2") findMacs2() else NULL,
  pathToMacs2 = pathToMacs2,
  genomeSize = NULL,
  shift = -75,
  extsize = 150,
  method = if (tolower(peakMethod) == "macs2") "q" else "p",
  cutOff = args$macs2_cutoff, # function default=0.1, script default=1e-7
  additionalParams = "--nomodel --nolambda",
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  genomeAnnotation = getGenomeAnnotation(proj.archr),
  geneAnnotation = getGeneAnnotation(proj.archr),
  plot = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addReproduciblePeakSet", logDir=dir_log)
) # addReproduciblePeakSet








cat(sprintf("\taddPeakMatrix\n"))
if ("package:ArchR.dpeerlab" %in% (search())) {
  proj.archr <- addPeakMatrix(proj.archr,
    ceiling = 4, # default=4, ceiling=10^9 for some pupurposes
    maxFragmentLength = 147, 
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    logFile = createLogFile("addPeakMatrix.proj.archr", logDir=dir_log)
  )
} else {
  # https://www.archrproject.com/reference/addPeakMatrix.html
  proj.archr <- addPeakMatrix(proj.archr,
    ceiling = 4,
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    logFile = createLogFile("addPeakMatrix.proj.archr", logDir=dir_log)
  )
} # if



cat(sprintf("\taddBgdPeaks\n"))
# https://www.archrproject.com/reference/addBgdPeaks.html
proj.archr <- addBgdPeaks(proj.archr,
	nIterations = 50,
	w = 0.1,
	binSize = 50,
	seed = 1,
	method = "chromVAR",
	outFile = file.path(dir_archr_output, "Background-Peaks.proj.archr.rds"),
	force = TRUE
  ) # addBgdPeaks



cat(sprintf("\tgetMarkerFeatures\n"))
# markerPeaks.archr, useMatrix = "PeakMatrix",  groupBy = "predictedGroup_ArchR"
markerPeaks.archr <- getMarkerFeatures(
  ArchRProj = proj.archr,
  useMatrix = "PeakMatrix",
  groupBy = "predictedGroup_ArchR",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures.PeakMatrix.predictedGroup_ArchR", logDir=dir_log)
) # getMarkerFeatures

# save RDS
saveRDS(markerPeaks.archr, sprintf("%s/%s_markerpeaks_proj.archr_predictedgroup.rds", dir_rds, cancer_type))




cat(sprintf("\tplotMarkerHeatmap\n"))
#heatmapPeaks.archr <- markerHeatmap(
heatmapPeaks.archr <- plotMarkerHeatmap(
  seMarker = markerPeaks.archr,
  #cutOff = "FDR <= 0.01 & Log2FC >= 0.5", # default
  cutOff = str_cutOff_markerPeaks_for_heatmap,
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2,2),
  #pal = NULL # default
  pal = viridis(n=256),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = NULL,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap.markerPeaks.archr", logDir=dir_log)
)

cat(sprintf("\tComplexHeatmap\n"))
ComplexHeatmap::draw(heatmapPeaks.archr, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks.archr , name = "markerpeaks_archr_predicted_labels", width = 8, height = 32, ArchRProj =  proj.archr, addDOC = FALSE) # output directory: All/Plots












### peaks enriched in epithelial cells

if (length(cluster.types_with_cancer_cells) == 0) {
	# collect cluster.type with epithelial cells
	cluster.type.u <- unique(proj.archr$predictedGroup_ArchR)

	f.tumor.epi <- grepl(pattern_tumor_epi, cluster.type.u)
 	if (any(f.tumor.epi)) {
		cluster.types_with_cancer_cells <- cluster.type.u[f.tumor.epi]
	} else {
		f.epi <- grepl(pattern_epi, cluster.type.u)
		cluster.types_with_cancer_cells <- cluster.type.u[f.epi]
	}
} # if



cat(sprintf("\tcluster.types_with_cancer_cells: %s\n", paste(cluster.types_with_cancer_cells, collapse=", ")))


# proj.archr$Cancer.Group
proj.archr$Cancer.Group <- ifelse(proj.archr$predictedGroup_ArchR %in% cluster.types_with_cancer_cells, "Cancer", "Normal")


cat(sprintf("\tgetMarkerFeatures\n"))
# https://www.archrproject.com/reference/getMarkerFeatures.html
# markerPeaks.for_comparison, useMatrix = "PeakMatrix", groupBy = "Cancer.Group"
markerPeaks.for_comparison <- getMarkerFeatures(
  ArchRProj = proj.archr,
  useMatrix = "PeakMatrix",
  groupBy = "Cancer.Group",
  useGroups = "Cancer",
  bgdGroups = "Normal",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 500,
  scaleTo = 10^4,
  threads = getArchRThreads(),
  k = 100,
  bufferRatio = 0.8,
  binarize = FALSE,
  useSeqnames = NULL,
  verbose = TRUE,
  logFile = createLogFile("getMarkerFeatures.PeakMatrix.Cancer.Group", logDir=dir_log)
) # getMarkerFeatures

# save RDS
saveRDS(markerPeaks.for_comparison, sprintf("%s/%s_markerpeaks_proj.archr_normal-vs-cancer.rds", dir_rds, cancer_type))




f_plotMarkerPeaks.cancer <- FALSE
#Error in hclust(get_dist(submat, distance), method = method) :
#  NA/NaN/Inf in foreign function call (arg 10)

if (f_plotMarkerPeaks.cancer) {

  cat(sprintf("\tplotMarkerHeatmap\n"))
  heatmapPeaks.cancer <- plotMarkerHeatmap(
    seMarker = markerPeaks.for_comparison,
    #cutOff = "FDR <= 0.01 & Log2FC >= 0.5", # default
    cutOff = str_cutOff_markerPeaks_for_heatmap,
    log2Norm = TRUE,
    scaleTo = 10^4,
    scaleRows = TRUE,
    plotLog2FC = TRUE, # Must use plotLog2FC = TRUE when ncol(seMarker) <= 2
    limits = c(-2,2),
    #pal = NULL # default
    pal = viridis(n=256),
    binaryClusterRows = TRUE,
    clusterCols = TRUE,
    labelMarkers = NULL,
    nLabel = 15,
    nPrint = 15,
    labelRows = FALSE,
    returnMatrix = FALSE,
    transpose = FALSE,
    invert = FALSE,
    logFile = createLogFile("plotMarkerHeatmap.markerPeaks.for_comparison", logDir=dir_log)
  ) # heatmapPeaks.cancer

  cat(sprintf("\tComplexHeatmap\n"))
  ComplexHeatmap::draw(heatmapPeaks.cancer, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmapPeaks.cancer, name = "markers_cancer-enriched_peaks", width = 8, height = 32, ArchRProj =  proj.archr, addDOC = FALSE) # output directory: All/Plots

} # if (f_plotMarkerPeaks.cancer)




































# RNA heatmap
####################################################################
topN <- Wilcox.markers %>% group_by(cluster) %>% dplyr::filter(p_val_adj <= 0.01) %>% top_n(30, desc(avg_log2FC))  # H. Kim changed avg_logFC to avg_log2FC

# Seurat Heatmap (top N DEGs per cluster)
#rna.sub <- subset(rna,downsample = 300)
#DoHeatmap(rna.sub,features = top30$gene,draw.lines = F,size = 2,disp.min = -2,disp.max = 2)+scale_fill_viridis()

## All 2000 variable features
#DoHeatmap(rna.sub,features = VariableFeatures(rna),draw.lines = F,size = 2,disp.min = -2,disp.max = 2)+scale_fill_viridis()


# Downsample cells from each cluster
rna.sub <- subset(rna,downsample =300)
rna.sub <- NormalizeData(rna.sub)
rna.sub <- rna.sub[rownames(rna.sub) %in% topN$gene,]
rna.sub <- ScaleData(rna.sub,features = rownames(rna.sub))

mat <- rna.sub@assays$RNA@scale.data

cluster_anno <- rna.sub@meta.data[, col_cluster_types]

# gg <- ggplot_build(rna.cell.plot)
#
# data <- as.data.frame(gg$data)
# names <- data$group
# names <- names -1
# data$group <- names
# #data$group <- data$colour
#
#
# cols <- unique(data$colour)
# names <- as.numeric(unique(data$group))
# names(cols) <- names
#
# cols <- cols[order(factor(names(cols), levels=levels(cluster_anno)))]

col_fun = circlize::colorRamp2(c(-2, 0, 2),viridis(n = 3))

# https://www.rdocumentation.org/packages/ComplexHeatmap/versions/1.10.2/topics/Heatmap
heatmapRNA <- Heatmap(mat, name = "Expression",
        column_split = factor(cluster_anno),
        cluster_columns =T,
        show_column_dend = F,
        cluster_column_slices = T,
        column_title_gp = gpar(fontsize = 8),
        column_gap = unit(0.1, "mm"),
        cluster_rows = T,
        show_row_dend = FALSE,
        col = col_fun,
        column_title_rot = 90,
        show_column_names = F)

ComplexHeatmap::draw(heatmapRNA, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#plotPDF(heatmapRNA, name = "heatmap_rna", width = 8, height = 6) # output directory: ./Plots
plotPDF(heatmapRNA, name = "heatmap_rna", width = 8, height = 32, ArchRProj = proj.archr) # output directory: All/Plots








###############################################################
#source("./r/archr/archr_v0.9.5_modified/Archr_Peak_Null.R")
# Overlap RNA/ATAC DEG hits 1) ArchR 2) Signac


# markers from scATACseq
# https://www.archrproject.com/reference/getMarkers.html
markerList.atac.up <- getMarkers(markersGS.archr.pred,
			cutOff = "FDR <= 0.01 & Log2FC >= 1.0",
			n = NULL,
			returnGR = FALSE
	)



markerList.atac.dn <- getMarkers(markersGS.archr.pred,
			cutOff = "FDR <= 0.01 & Log2FC <= -1.0",
			n = NULL,
			returnGR = FALSE
	)


# markers from scRNAseq
#Wilcox.markers <- readRDS(sprintf("%s/wilcox_degs/%s_wilcox_degs.rds", dir_seurat_obj, cancer_type))
#Wilcox.markers$cluster <- str_replace(Wilcox.markers$cluster, "/", "_")

# make marker RNA list
markerList.rna.up <- list()
markerList.rna.dn <- list()
cluster.types <- levels(factor(rna@meta.data[, col_cluster_types]))
for ( cluster.type1 in cluster.types) {

   markerList.rna.up[[cluster.type1]] <- Wilcox.markers[Wilcox.markers$cluster == cluster.type1 &
	Wilcox.markers$p_val_adj <= 0.01 &
	Wilcox.markers$avg_log2FC >= 1.0,] # H Kim changed avg_logFC to avg_log2FC.

   #cluster.type1 <- gsub(" ", "_", cluster.type1)
   #saveRDS(rna.cluster.up[[cluster.type1]], sprintf("%s/%s_rna_cluster_markers_%s_up.rds", dir_rds, cancer_type, cluster.type1))

   markerList.rna.dn[[cluster.type1]] <- Wilcox.markers[Wilcox.markers$cluster == cluster.type1 &
	Wilcox.markers$p_val_adj <= 0.01 &
	Wilcox.markers$avg_log2FC <= -1.0,] # H Kim changed avg_logFC to avg_log2FC.
   #saveRDS(rna.cluster.dn[[cluster.type1]], sprintf("%s/%s_rna_cluster_markers_%s_dn.rds", dir_rds, cancer_type, cluster.type1))

} # for


gene.hits.up <- list()
gene.hits.dn <- list()
for (cluster.type1 in cluster.types) {

   atac.clust.up <- markerList.atac.up[[cluster.type1]]
   rna.clust.up <- markerList.rna.up[[cluster.type1]]
   gene.hits.up[[cluster.type1]] <- intersect(atac.clust.up$name, rna.clust.up$gene)

   atac.clust.dn <- markerList.atac.dn[[cluster.type1]]
   rna.clust.dn <- markerList.rna.dn[[cluster.type1]]
   gene.hits.dn[[cluster.type1]] <- intersect(atac.clust.dn$name, rna.clust.dn$gene)

} # for



# list_marker_info
list_marker_info <- list()
list_marker_info$markerList.rna.up <- markerList.rna.up
list_marker_info$markerList.rna.dn <- markerList.rna.dn
list_marker_info$markerList.atac.up <- markerList.atac.up
list_marker_info$markerList.atac.dn <- markerList.atac.dn
list_marker_info$gene.hits.up <- gene.hits.up
list_marker_info$gene.hits.dn <- gene.hits.dn

# save RDS
saveRDS(list_marker_info, sprintf("%s/%s_markergenes_overlap_rna_and_atac_predictedgroup.rds", dir_rds, cancer_type))





# update df_output
for (sample in sample_names) {

  proj.i <- proj.archr[proj.archr$Sample == sample]
  meta.data <- as.data.frame(getCellColData(proj.i))
  # Sample, TSSEnrichment, ReadsInTSS, ReadsInPromoter, ReadsInBlacklist, PromoterRatio, PassQC, NucleosomeRatio, nMultiFrags, nMonoFrags, nFrags, nDiFrags, DoubletScore, DoubletEnrichment, BlacklistRatio, ATAC_clusters, predictedCell_ArchR, predictedGroup_ArchR, predictedScore_ArchR, Mode.Label, ReadsInPeaks, FRIP, Cancer.Group
  # nCount_RNA, nFeature_RNA, percent.mt, Phase, RNA_snn_res.0.2, seurat_clusters, epi_krt_epcam, epi_normal_tumor, normal_epi_type, epi_type, SingleR.HTAPP_toolbox, SingleR.HPCA, SingleR.BED, celltype.cycling, PScore, cell.type, CNV.value, CNV.corr, CNV.Pos, CNV.type, CNV.genes, cluster.type, RNA_harmony_th.0, cluster.type.harmony
  df_output["TotalReadsInPeaks", sample] <- sum(meta.data$ReadsInPeaks)

} # for









################################################################
# PART 4: Motif and Feature Enrichment
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 4: Motif and Feature Enrichment\n"))
cat(sprintf("------------------------------------\n\n"))
################################################################

# reference:
# https://www.archrproject.com/bookdown
# /home/hkim77/francolab/scRNA_scATAC_Breast_Cancer_Project_2021/ER+_Plus_TrueNormal_Patients_scATAC-PeakCalling_chromVAR.R
# /home/hkim77/francolab/scRNA_scATAC_Breast_Cancer_Project_2021/ER+_Plus_TrueNormal_Patients_scATAC_scRNA-Plotting.R


str_motif_set <- "cisbp"
# https://www.archrproject.com/reference/addMotifAnnotations.html
proj.archr <- addMotifAnnotations(ArchRProj = proj.archr,
			motifSet = str_motif_set, # The motif set to be used for annotation. Options include: (i) "JASPAR2016", "JASPAR2018", "JASPAR2020" which gives the 2016, 2018 or 2020 version of JASPAR motifs or (ii) one of "cisbp", "encode", or "homer" which gives the corresponding motif sets from the chromVAR package.
			name = "Motif",
			logFile = createLogFile("addMotifAnnotations", logDir=dir_log)
		)




n_markerpeaks.for_comparison <- nrow(markerPeaks.for_comparison)
cat(sprintf("\tn_markerpeaks.for_comparison: %d\n", n_markerpeaks.for_comparison))
df_output["n_markerpeaks.for_comparison", "info"] <- n_markerpeaks.for_comparison




### up

cat(sprintf("\tpeakAnnoEnrichment for up-regulated motifs\n"))
motifsUp <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr,
    peakAnnotation = "Motif",
    #cutOff = "FDR < 0.1 & Log2FC > 0.5",
    cutOff = "FDR < 0.01 & Log2FC > 1.0",
    logFile = createLogFile("peakAnnoEnrichmentUp", logDir=dir_log)
  )


df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

n_motifs.up <- nrow(df)
cat(sprintf("\t\tn_motifs.up: %d\n", n_motifs.up))
df_output["n_motifs.up", "info"] <- n_motifs.up

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab(TeX(r'(-$log_{10} (P_{adj})$ Motif Enrichment)')) +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(sprintf("%s/scatterplot_%s_motif_up.cancer_clusters.pdf", dir_pdf, str_motif_set), width = 5, height = 5, plot=ggUp)
ggsave(sprintf("%s/scatterplot_%s_motif_up.cancer_clusters.png", dir_png, str_motif_set), width = 5, height = 5, plot=ggUp)

cutOff_mlog10Padj <- 5
if (length(which(df$mlog10Padj < cutOff_mlog10Padj)) > 0) {

  heatmapEM <- tryCatch({
		plotEnrichHeatmap(motifsUp, n = 20,
			cutOff = cutOff_mlog10Padj, # default=20
			transpose = TRUE,
			returnMatrix = FALSE,
			logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log))
	}, error = function(e) {
		print(paste('#error-handler-code'))
		cat(sprintf('%s',e))
		return(NULL)
	}, finally = {
	})

  if (!is.null(heatmapEM)) {
	# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
	plotPDF(heatmapEM, name = sprintf("heatmap_%s_motif_up.cancer_clusters", str_motif_set), width = 8, height = 5, ArchRProj = proj.archr, addDOC = FALSE)
  } # if

} # if



# 12.3 ArchR enrichment
str_collection <- "EncodeTFBS"
#str_collection <- "ATAC"
#str_collection <- "Codex"

# https://www.archrproject.com/reference/addArchRAnnotations.html
cat(sprintf("\taddArchRAnnotations for %s\n", str_collection))
proj.archr <- addArchRAnnotations(ArchRProj = proj.archr,
	collection = str_collection,
	force = TRUE,
	logFile = createLogFile("addArchRAnnotations", logDir=dir_log)
  )

cat(sprintf("\tpeakAnnoEnrichment for up-regulated motifs with %s\n", str_collection))
motifsUp <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr,
    peakAnnotation = str_collection,
    #cutOff = "FDR < 0.1 & Log2FC > 0.5",
    cutOff = "FDR < 0.01 & Log2FC > 1.0",
    logFile = createLogFile(sprintf("peakAnnoEnrichment_%s_Up", str_collection), logDir=dir_log)
  )


df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

n_motifs.up.EncodeTFBS <- nrow(df)
cat(sprintf("\t\tn_motifs.up.EncodeTFBS: %d\n", n_motifs.up.EncodeTFBS))
df_output["n_motifs.up.EncodeTFBS", "info"] <- n_motifs.up.EncodeTFBS

if (length(which(df$mlog10Padj < cutOff_mlog10Padj)) > 0) {

  heatmapEM <- tryCatch({
	plotEnrichHeatmap(motifsUp, n = 20,
		cutOff = cutOff_mlog10Padj, # default=20
		transpose = TRUE,
		returnMatrix = FALSE,
		logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log))
	}, error = function(e) {
		print(paste('#error-handler-code'))
		cat(sprintf('%s',e))
		return(NULL)
	}, finally = {
	})

  if (!is.null(heatmapEM)) {
	# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
	plotPDF(heatmapEM, name = sprintf("heatmap_%s_enrichemnt_cancer_clusters_up", str_collection), width = 8, height = 6, ArchRProj = proj.archr, addDOC = FALSE)
  } # if

} # if









### down


# https://www.archrproject.com/reference/peakAnnoEnrichment.html
# This function will perform hypergeometric enrichment of a given peak annotation within the defined marker peaks.
cat(sprintf("\tpeakAnnoEnrichment for down-regulated motifs\n"))
motifsDo <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr,
    peakAnnotation = "Motif",
    #cutOff = "FDR < 0.1 & Log2FC < -0.5",
    cutOff = "FDR < 0.01 & Log2FC < -1.0",
    logFile = createLogFile("peakAnnoEnrichment_Do", logDir=dir_log)
  )


df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

n_motifs.down <- nrow(df)
cat(sprintf("\t\tn_motifs.down: %d\n", n_motifs.down))
df_output["n_motifs.down", "info"] <- n_motifs.down


ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab(TeX(r'(-$log_{10} (P_{adj})$ Motif Enrichment)')) +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(sprintf("%s/scatterplot_%s_motif_down.cancer_clusters.pdf", dir_pdf, str_motif_set), width = 5, height = 5, plot=ggDo)
ggsave(sprintf("%s/scatterplot_%s_motif_down.cancer_clusters.png", dir_png, str_motif_set), width = 5, height = 5, plot=ggDo)

if (length(which(df$mlog10Padj < cutOff_mlog10Padj)) > 0) {

  heatmapEM <- tryCatch({
	plotEnrichHeatmap(motifsDo, n = 20,
		cutOff = cutOff_mlog10Padj, # default=20
		transpose = TRUE,
		returnMatrix = FALSE,
		logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log))
	}, error = function(e) {
		print(paste('#error-handler-code'))
		cat(sprintf('%s',e))
		return(NULL)
	}, finally = {
	})

  if (!is.null(heatmapEM)) {
	# ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
	plotPDF(heatmapEM, name = sprintf("heatmap_%s_motif_down.cancer_clusters", str_motif_set), width = 8, height = 5, ArchRProj = proj.archr, addDOC = FALSE)
  } # if

} # if





### up_down
plotPDF(ggUp, ggDo, name = sprintf("scatterplots_%s_motif_updown.cancer_clusters", str_motif_set), width = 5, height = 5, ArchRProj = proj.archr, addDOC = FALSE)









################################################################
# PART 5. ChromVAR deviations enrichment
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 5: ChromVAR deviations enrichment\n"))
cat(sprintf("------------------------------------\n\n"))
################################################################



# https://www.archrproject.com/reference/addBgdPeaks.html
#proj.archr <- addBgdPeaks(proj.archr)


# custom peak annotations
#EncodePeaks <- c(
#  Encode_K562_GATA1 = "https://www.encodeproject.org/files/ENCFF632NQI/@@download/ENCFF632NQI.bed.gz",
#  Encode_GM12878_CEBPB = "https://www.encodeproject.org/files/ENCFF761MGJ/@@download/ENCFF761MGJ.bed.gz",
#  Encode_K562_Ebf1 = "https://www.encodeproject.org/files/ENCFF868VSY/@@download/ENCFF868VSY.bed.gz",
#  Encode_K562_Pax5 = "https://www.encodeproject.org/files/ENCFF339KUO/@@download/ENCFF339KUO.bed.gz"
#)
#if("ChIP" %ni% names(proj.archr@peakAnnotation)){
#    proj.archr <- addPeakAnnotations(ArchRProj = proj.archr, regions = EncodePeaks, name = "ChIP")
#}



idxSample <- BiocGenerics::which(proj.archr$predictedGroup_ArchR %in% cluster.types_with_cancer_cells)
cellsSample <- proj.archr$cellNames[idxSample]
proj.archr.cancer <- proj.archr[cellsSample, ]


# We are now ready to compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. This function has an optional parameter called matrixName that allows us to define the name of the deviations matrix that will be stored in the Arrow files. If we do not provide a value to this parameter, as in the example below, this function creates a matrix name by adding the word “Matrix” to the name of the peakAnnotation. The example below creates a deviations matrix in each of our Arrow files called “MotifMatrix”.
str_peakAnnotation <- "Motif" # {"Motif", "EncodeTFBS", "ATAC"}

proj.archr.cancer <- addDeviationsMatrix(
  ArchRProj = proj.archr.cancer, 
  peakAnnotation = str_peakAnnotation,
  force = TRUE,
  logFile = createLogFile(sprintf("addDeviationsMatrix_%s", str_peakAnnotation), logDir=dir_log)
) # addDeviationsMatrix


#motif.mat <- getMatrixFromProject(
#  ArchRProj = proj.archr.cancer,
#  useMatrix = "MotifMatrix",
#  useSeqnames = NULL,
#  verbose = TRUE,
#  binarize = FALSE,
#  threads = getArchRThreads(),
#  logFile = createLogFile("getMatrixFromProject", logDir=dir_log)
#)


# str_peak_annotation_matrix <- "MotifMatrix"
str_peak_annotation_matrix <- sprintf("%sMatrix", str_peakAnnotation)
plotVarDev <- getVarDeviations(proj.archr.cancer, name = str_peak_annotation_matrix, plot = TRUE)

plotPDF(plotVarDev, name = sprintf("scatterplot_variable_%s_deviation_scores_cancer_clusters", str_peakAnnotation), width = 5, height = 5, ArchRProj = proj.archr.cancer, addDOC = FALSE)


# https://www.archrproject.com/bookdown/motif-deviations.html
#motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
#markerMotifs <- getFeatures(proj.archr.cancer, select = paste(motifs, collapse="|"), useMatrix = str_peak_annotation_matrix)












################################################################
# PART 6. Footprinting with ArchR
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 6: Footprinting with ArchR\n"))
cat(sprintf("------------------------------------\n\n"))
################################################################


#motifPositions <- getPositions(proj.archr.cancer)
#motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
#markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
#markerMotifs

# addGroupCoverages() has been done above.
#projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "Clusters2")
#seFoot <- getFootprints(
#  ArchRProj = proj.archr.cancer, 
#  positions = motifPositions[markerMotifs], 
#  groupBy = "predictedGroup_ArchR"
#)

#plotFootprints(
#  seFoot = seFoot,
#  ArchRProj = proj.archr.cancer, 
#  normMethod = "Divide", # {"Subtract", "Divide", "None"}
#  plotName = "Footprints-Divide-Bias",
#  addDOC = FALSE,
#  smoothWindow = 5
#)

#seTSS <- getFootprints(
#  ArchRProj = proj.archr.cancer, 
#  positions = GRangesList(TSS = getTSS(proj.archr)), 
#  groupBy = "predictedGroup_ArchR",
#  flank = 2000
#)

#plotFootprints(
#  seFoot = seTSS,
#  ArchRProj = proj.archr.cancer, 
#  normMethod = "None",
#  plotName = "TSS-No-Normalization",
#  addDOC = FALSE,
#  flank = 2000,
#  flankNorm = 100
#)










################################################################
# PART 7. Integrative Analysis with ArchR
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("PART 7: Integrative Analysis with ArchR\n"))
cat(sprintf("------------------------------------\n\n"))
################################################################




# Add Coaccessiblity and Peak2Gene links:
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("addCoAccessibility\n"))

proj.archr <- addCoAccessibility(
  ArchRProj = proj.archr,
  reducedDims = reducedDims,
  dimsToUse = dimsToUse,
  scaleDims = NULL,
  corCutOff = 0.75,
  k = 100,
  knnIteration = 500,
  overlapCutoff = 0.8,
  maxDist = 1e+05,
  scaleTo = 10^4,
  log2Norm = TRUE,
  seed = 1,
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("addCoAccessibility", logDir=dir_log)
)


# addPeak2GeneLinks
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("addPeak2GeneLinks\n"))

proj.archr <- addPeak2GeneLinks(
  ArchRProj = proj.archr ,
  reducedDims = reducedDims,
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
  addEmpiricalPval = TRUE, # added for v1.0.1
  seed = 1,
  threads = max(floor(getArchRThreads()/2), 1),
  verbose = TRUE,
  logFile = createLogFile("addPeak2GeneLinks", logDir=dir_log)
) # addPeak2GeneLinks








# update proj.archr@cellColData
md_rna <- rna@meta.data
md_atac <- proj.archr@cellColData
cols_rna <- colnames(md_rna)
cols_atac <- colnames(md_atac)
cols <- setdiff(cols_rna, cols_atac)
if (f_multiome) {

	cat(sprintf("update proj.archr@cellColData with multiome\n"))
	idx <- match(rownames(md_atac), rownames(md_rna))
	f <- !is.na(idx); idx <- idx[f]
	md_atac[f, cols] <- md_rna[idx, cols]
	# update proj.archr@cellColData
	proj.archr@cellColData <- md_atac	

	n_cells <- nrow(md_atac)
	n_cells_correctly_predicted_group <- length(which(md_atac$cluster.type == md_atac$predictedGroup_ArchR))
	cat(sprintf("\tcorrectly predictedGroup_ArchR: %d/%d=%g\n", n_cells_correctly_predicted_group, n_cells, n_cells_correctly_predicted_group/n_cells))

	n_cells_correctly_predicted_cells <- length(which(rownames(md_atac) == md_atac$predictedCell_ArchR))
	cat(sprintf("\tcorrectly predictedCell_ArchR: %d/%d=%g\n", n_cells_correctly_predicted_cells, n_cells, n_cells_correctly_predicted_cells/n_cells))


} else {

	cat(sprintf("update proj.archr@cellColData with predicted cells\n"))
	idx <- match(md_atac$predictedCell_ArchR, rownames(md_rna))
	f <- !is.na(idx); idx <- idx[f]
	md_atac[f, cols] <- md_rna[idx, cols]
	# update proj.archr@cellColData
	proj.archr@cellColData <- md_atac	

} # if f_multiome


# number of samples for each group
#dt <- summarize_meta_data(proj.archr, "predictedGroup_ArchR", "Sample")











################################################################
# save archproj obj


f_save_final_rds <- TRUE
if (f_save_final_rds) {
  # saveRDS
  filename <- sprintf("%s/%s_archrproj_obj_final.rds", dir_rds, cancer_type)
  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("save RDS: %s\n", filename))
  saveRDS(proj.archr, filename)
} # if










# # For every cluster, visualize the peaks coaccessible with marker genes
# 
# gg <- ggplot_build(rna.cell.plot)
# 
# data <- as.data.frame(gg$data)
# data$cluster.type <- rna.emb[, col_cluster_types]
# 
# 
# cols <- unique(data$colour)
# names(cols) <- unique(data$cluster.type)
# 
# for (i in names(gene.hits)){
#   genes <- gene.hits[[i]]
# 
#   p <- plotBrowserTrack(
#     ArchRProj = proj.archr,
#     groupBy = "predictedGroup_ArchR",
#     geneSymbol = genes,
#     upstream = 50000,
#     downstream = 50000,
#     loops = getCoAccessibility(
#       ArchRProj = proj.archr,
#       corCutOff = 0.70,  # ArchR v0.9.5 default=0.5
#       resolution = 1,    # ArchR v0.9.5 default=1000
#       returnLoops = TRUE
#     ),
#     features = makeGRangesFromDataFrame(encode.all),
#     pal = cols
#   )
#   plotPDF(plotList = p,
#           name = paste0("ArchR_CoAscblty_Hits",i, ".pdf"),
#           ArchRProj = proj.archr,
#           addDOC = T, width = 5, height = 5)
# 
# 
#   p <- plotBrowserTrack(
#     ArchRProj = proj.archr,
#     groupBy = "predictedGroup_ArchR",
#     geneSymbol = genes,
#     upstream = 50000,
#     downstream = 50000,
#     loops = getPeak2GeneLinks(
#       ArchRProj = proj.archr,
#       corCutOff = 0.70,   # ArchR v0.9.5 default=0.45
#       FDRCutOff = 1e-04,  # ArchR v0.9.5 default=1e-4
#       resolution = 1,     # ArchR v0.9.5 default=1000
#       returnLoops = TRUE
#     ),
#     features = makeGRangesFromDataFrame(encode.all),
#     pal = cols
# 
#   )
#   loops <- getPeak2GeneLinks(
#     ArchRProj = proj.archr,
#     corCutOff = 0.70,
#     FDRCutOff = 1e-04,
#     resolution = 1,
#     returnLoops = F
#   )
#   #saveRDS(loops, paste0("./rds", i, "_empirical_pval_archr.rds"))
#   saveRDS(loops, sprintf("%s/%s_empirical_pval_archr.rds", dir_rds, i))
# 
#   plotPDF(plotList = p,
#           name = paste0("ArchR_Peak2Gene_Hits_",i, ".pdf"),
#           ArchRProj = proj.archr,
#           addDOC = T, width = 6, height = 8)
# 
# }







if (args$f_make_tfsee_input) {
	
	source("./r/make_tfsee_input.R")
	dir_bam_subset <- sprintf("%s/bam_subset_archr", dir_output)
	out <- make_tfsee_input(proj.archr.cancer, dir_bam_subset)

} # if












cat(sprintf("------------------------------------\n"))
cat(sprintf("\tprint meta data\n\n"))

proj.meta <- as.data.frame(proj.archr@cellColData)
print(head(proj.meta))



if (args$n_debug < 100) {

  # save archr proj
  cat(sprintf("\n------------------------------------\n"))
  cat(sprintf("\tsaveArchRProject\n"))
  saveArchRProject( ArchRProj = proj.archr,
	outputDirectory = dir_archr_output,
	overwrite = TRUE,
	load = FALSE,
	dropCells = TRUE, # A boolean indicating whether to drop cells that are not in ArchRProject from corresponding Arrow Files.
	logFile = createLogFile("saveArchRProject", logDir=dir_log),
	threads = getArchRThreads()
  ) # saveArchRProject

} # if


# write xlsx
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\twrite xlsx\n"))
file_name_xlsx <- sprintf("%s/%s_sc-atac-seq_pipeline_summary.xlsx", dir_xlsx, cancer_type)
openxlsx::write.xlsx(df_output, file_name_xlsx, row.names = TRUE, col.names = TRUE)

print(df_output)



# post-process
unlink(dir_tmp)
unlink("Rplots.pdf")
# for reducing disk usage
unlink(sprintf("%s/*.arrow", dir_output))
unlink(sprintf("%s/GroupCoverages", dir_archr_output), recursive = TRUE)



# message
cat(sprintf("------------------------------------\n"))
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("completed successfully~\n"))
cat(sprintf("------------------------------------\n\n"))




# end of script


