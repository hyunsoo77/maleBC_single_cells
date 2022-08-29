#!/usr/bin/env Rscript
#
# make_sc-rna-seq_merged_seurat_obj.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Feb.
#
# options:
#   --n_cores: (default=8) the number of cores for parallel computing
#   --dir_otuput: (default="./output")
#   --dir_log: (default="./output/log")
#
#   --dir_seurat_obj: (default="./output/rds")
#   --filename_to_select_barcodes: (default="")  filename to select barcodes
#   --pattern_to_select_cell_types: (default="")  pattern to select cell types in seurat obj
#   --max_n_cells: (default=1e9) select n cells with higher number of genes expressed
#   --min_n_cells: (default=50) ignore samples with low number of cells
#   --type_parsing_rds_filename: (default="gsub") type for parsing rds file name
#   --cancer_type_for_parsing_rds_filename: (default="")
#   --samples_to_exclude: (default=""): samples to exclude e.g. sample1_to_exclude,sample2_to_exclude
#
#   --method_sctransform_vst: (default="") {"poisson", "glmGamPoi", ...}, use sctransform::vst() when method_sctransform_vst is not empty
#   --method_integration: {["cca"], "rpca", "rlsi", "none"}, use rpca for large sample sizes
#   --k.anchor: (default=5) increase the strength of alignment by increasing the k.anchor parameter  https://satijalab.org/seurat/articles/integration_rpca.html
#   --type_integration_anchor_features: {["all"], "2000"}, type of integration anchor features
#   --method_batch_effect_correction: (default="") pattern of {"cca", "rpca", "rlsi", "harmony", "none"} for "cca|rpca|rlsi": integrate data for assay="integrated" and slot="data", harmony: reduction=harmony for cell-based analysis
#   --batch_keys_for_reference: (default="") batch keys for reference, used in seurat object integration (e.g. sample1_good_quality,sample2_good_quality)
#   --k.weight: (default=100) number of neighbors to consider when weighting anchors
#
#   --vars.to.regress: (default="percent.mt") e.g. "nCount_RNA,nFeature_RNA,percent.mt,S.Score,G2M.Score"
#
#   --max_dimstouse: (default=30) dimsToUse <- 1:args$max_dimstouse
#   --seurat_resolution: (default: 0.8) Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.   https://satijalab.org/seurat/reference/findclusters
#   --seed.seurat_clustering: (default: 0) 
#
#   --seed.harmony: (default=51) seed for computing harmony
#   --type_parsing_rds_filename_for_dataset: type for parsing rds file name for dataset
#   --type_parsing_rds_filename_for_species: type for parsing rds file name for species
#   --type_parsing_rds_filename_for_tumor.type: type for parsing rds file name for tumor.type
#   --type_parsing_rds_filename_for_donor: type for parsing rds file name for donor
#   --type_parsing_rds_filename_for_technology: type for parsing rds file name for technology
#   --type_parsing_rds_filename_for_condition: type for parsing rds file name for condition
#   --type_parsing_rds_filename_for_treatment: type for parsing rds file name for treatment
#
#   --harmony_theta: (default=-1) theta for computing harmony, compute harmony when harmony_theta[1] >= 0
#   --harmony_lambda: lambda for computing harmony, default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction
#   --seed.symphony: (default=-1) seed for computing symphony, compute symphony when seed.symphony > 0
#
#   --type_signatures: {[""], "panglaodb", "panglaodb_plus_others"}, type of signatures to get biomarkers for RunPCA
#
#   --umap_n_neighbors: (defualt=30)
#   --umap_min_dist: (default=0.3)
#   --umap_metric: (default="cosine")
#
#   --method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", ["singler_blueprint_encode"]}
#
#   --f_save_jackstraw_png: (default=FALSE) flag to save jackstraw png
#   --f_save_cnv_pc1_boxplot_png: (default=FALSE) flag to save cnv_pc1_boxplot.png
#   --f_save_infercnv_rds: (default=FALSE) flag to save infercnv rds
#   --f_save_umap_figure: (default=FALSE) flag to save umap figure
#   --f_save_vlnplot_png: (default=FALSE) flag to save violinplot png
#   
#   --str_update_cell_types: (default="") update cell types in get_most_common_cell_type_for_each_cluster(), (e.g. "B-cells,Mast cells"), execute "rna <- identify_cell_types(rna, "panglaodb", args)" for "B-cells", execute "rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)" for using "Mast cells". this feature was turned off for more accurate cell type annotations.
#
#   --multik: logical to execute multik for clustering labels at optimial K
#   --multik_reps: (default=10) # of repeats of subsampling and consensusing clustering to determine optimal Ks
#   --no_diet_seurat: do not diet seurat obj, but store full seurat obj
#
# input:
#   cancer_type: (e.g. brca1-mut-bc, er+-bc, her2+-bc, male-bc, normal-breast, tnbc, tr-bc, matt-oc)
#
# output:
#   $dir_seurat_obj/${cancer_type}_sc-rna-seq_merged_seurat_obj.rds
#   xlsx/${cancer_type}_sc-rna-seq_pipeline_summary.xlsx
#
# usage:
# ./make_sc-rna-seq_merged_seurat_obj.R male-bc
#
# reference:
# this script is based on the origical script written by Matt Regner  https://github.com/RegnerM2015/scENDO_scOVAR_2020
# https://satijalab.org/
# https://github.com/satijalab/seurat-wrappers
#
#
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser(description="hs script",python_cmd="python")

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("--n_log", default=1, type="integer", help="n_log")
parser$add_argument("--n_debug", default=0, type="integer", help="n_debug")
parser$add_argument("-l","--n_log_level", default=0, type="integer", help="n_log_level")
parser$add_argument("-nc","--n_cores", default=8, type="integer", help="the number of cores for parallel computing")

# arguments for input/output
parser$add_argument("-do", "--dir_output", default="./output", type="character", help="direcotry for output")
parser$add_argument("-dl", "--dir_log", default="./output/log", type="character", help="direcotry for log")
#parser$add_argument("-ff", "--figure_format", default="pdf", type="character", help="figure format") # coded to print figures with pdf format (pdf and png whenever it is possible since png is helpful for lighter PPT slides).


# arguments for seurat
parser$add_argument("-dso","--dir_seurat_obj", default="./output/rds", type="character", help="direcotry for seurat objects")
parser$add_argument("-fb","--filename_to_select_barcodes", default="", type="character", help="filename to select barcodes")
parser$add_argument("-pct","--pattern_to_select_cell_types", default="", type="character", help="pattern to select cell types")
parser$add_argument("-maxc","--max_n_cells", default=1e9, type="double", help="select n cells for each sample with higher number of genes expressed")
parser$add_argument("-minc","--min_n_cells", default=50, type="double", help="ignore samples with low number of cells")
parser$add_argument("-tprdsfn","--type_parsing_rds_filename", default="gsub", type="character", help="type for parsing rds file name")
parser$add_argument("-ctprdsfn","--cancer_type_for_parsing_rds_filename", default="", type="character", help="cancer type for parsing rds file name")
parser$add_argument("-ste","--samples_to_exclude", default="", type="character", help="samples to exclude e.g. sample1_to_exclude,sample2_to_exclude")


# seurat/batch effect correction
parser$add_argument("-mv","--method_sctransform_vst", default="", type="character", help="method for sctransform::vst()")
parser$add_argument("-mi","--method_integration", default="cca", type="character", help="method for integration {[cca], rpca, rlsi, none}, use rpca for large sample sizes")
parser$add_argument("-ka","--k.anchor", default=5, type="integer", help="increase the strength of alignment by increasing the k.anchor parameter")
parser$add_argument("-tiaf","--type_integration_anchor_features", default="all", type="character", help="type of integration anchor features")
#parser$add_argument("-int","--integrated", dest="f_integrated", default=TRUE, action="store_true", help="integrate data for assay=\"integrated\" and slot=\"data\"")
#parser$add_argument("-nint","--no_integrated", dest="f_integrated", action="store_false", help="do not integrate data for assay=\"integrated\" and slot=\"data\"")
parser$add_argument("-mbr","--method_batch_effect_correction", default="", type="character", help="pattern of {\"cca\", \"rpca\", \"rlsi\", \"harmony\", \"none\"} for \"cca|rpca|rlsi\": integrate data for assay=\"integrated\" and slot=\"data\", harmony: reduction=\"harmony\" for cell-based analysis")
parser$add_argument("-bkr","--batch_keys_for_reference", default="", type="character", help="batch keys for reference, used in seurat object integration (e.g. sample1_good_quality,sample2_good_quality)")
parser$add_argument("-kw","--k.weight", default=100, type="integer", help="number of neighbors to consider when weighting anchors")


# seurat/scaling
parser$add_argument("-vr","--vars.to.regress", default="percent.mt", type="character", help="vars.to.regress for ScaleData()")

# seurat/pca
parser$add_argument("-maxdim","--max_dimstouse", default=30, type="integer", help="max dimsToUse")

# seurat/clustering
parser$add_argument("-res","--seurat_resolution", default=0.8, type="double", help="Seurat FindClusters parameter of resolution")
parser$add_argument("-ssc","--seed.seurat_clustering", default=0, type="double", help="Seurat FindClusters seed")



# arguments for harmony
parser$add_argument("-sh","--seed.harmony", default=51, type="double", help="seed for computing harmony")
parser$add_argument("-ds","--type_parsing_rds_filename_for_dataset", default="", type="character", help="type for parsing rds file name for dataset")
parser$add_argument("-sp","--type_parsing_rds_filename_for_species", default="", type="character", help="type for parsing rds file name for species")
parser$add_argument("-tt","--type_parsing_rds_filename_for_tumor.type", default="", type="character", help="type for parsing rds file name for tumor.type")
parser$add_argument("-donor","--type_parsing_rds_filename_for_donor", default="", type="character", help="type for parsing rds file name for donor")
parser$add_argument("-tech","--type_parsing_rds_filename_for_technology", default="", type="character", help="type for parsing rds file name for technology")
parser$add_argument("-cond","--type_parsing_rds_filename_for_condition", default="", type="character", help="type for parsing rds file name for condition")
parser$add_argument("-treat","--type_parsing_rds_filename_for_treatment", default="", type="character", help="type for parsing rds file name for treatment")


parser$add_argument("-ht","--harmony_theta", default="-1", type="character", help="theta for computing harmony")
parser$add_argument("-hl", "--harmony_lambda", default="-1", type="character", help="lambda for computing harmony, default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction.")


# arguments for symphony
parser$add_argument("-ss","--seed.symphony", default=-1, type="double", help="seed for computing symphony")

# argument for pca
parser$add_argument("-tsig","--type_signatures", default="", type="character", help="type of signatures ({[\"\"], \"panglaodb\", \"panglaodb_plus_others\"}) to get biomarkers for RunPCA")

# arguments for umap
parser$add_argument("-nn","--umap_n_neighbors", default=30, type="integer", help="UMAP number of neighbors")
parser$add_argument("-dist","--umap_min_dist", default=0.3, type="double", help="UMAP min distance")
parser$add_argument("-metric","--umap_metric", default="cosine", type="character", help="UMAP metrics")

# arguments for SingleR cell typing
parser$add_argument("-mict","--method_to_identify_cell_types", default="singler_blueprint_encode", type="character", help="method to identify cell types")

# arguments for infercnv
parser$add_argument("-sjp","--f_save_jackstraw_png", action="store_true", default=FALSE, help="save jackstraw png")
parser$add_argument("-scpbp","--f_save_cnv_pc1_boxplot_png", action="store_true", default=FALSE, help="save cnv_pc1_boxplot.png")
parser$add_argument("-sir","--f_save_infercnv_rds", action="store_true", default=FALSE, help="save InferCNV rds")

# arguments for plots
parser$add_argument("-sup","--f_save_umap_figure", action="store_true", default=FALSE, help="save umap figures")
parser$add_argument("-svp","--f_save_vlnplot_png", action="store_true", default=FALSE, help="save violinplot png")

# arguments for renaming cell types
parser$add_argument("-uct","--str_update_cell_types", default="", type="character", help="update cell types in get_most_common_cell_type_for_each_cluster() (e.g. \"B-cells,Mast cells\"")

# arguments for others
parser$add_argument("-cts","--cancer_type_standard", default="", type="character", help="standard cancer type")
parser$add_argument("-mk", "--multik", dest="f_multik", action="store_true", default=FALSE, help="execute multik for clustering labels at optimial K")
parser$add_argument("-mkr","--multik_reps", default=10, type="integer", help="# of repeats of subsampling and consensusing clustering to determine optimal Ks")
parser$add_argument("-cpdb", "--cellphonedb", dest="f_cellphonedb", action="store_true", default=FALSE, help="execute cellphonedb")
parser$add_argument("-nds", "--no_diet_seurat", dest="f_diet_seurat", action="store_false", default=TRUE, help="do not diet seurat object for reducing memory/disk usage")
parser$add_argument("-dsl", "--diet_seurat_level", default=1, type="double", help="diet seurat object level")



parser$add_argument("cancer_type", nargs=1, help="cancer_type")



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


dimsToUse <- 1:args$max_dimstouse
args$dimsToUse <- dimsToUse






# create output directory
if (dirname(args$dir_log) == args$dir_output) {
        dir_log <- args$dir_log
} else {
        dir_log <- sprintf("%s/%s", dir_output, args$dir_log)
}
args$dir_log <- dir_log
dir_pdf <- sprintf("%s/pdf", dir_output)
dir_png <- sprintf("%s/png", dir_output)
dir_rdata <- sprintf("%s/rdata", dir_output)
if (dirname(args$dir_seurat_obj) == args$dir_output) {
        dir_rds <- args$dir_seurat_obj
} else {
        dir_rds <- sprintf("%s/%s", dir_output, args$dir_seurat_obj)
}
args$dir_rds <- dir_rds
dir_xlsx <- sprintf("%s/xlsx", dir_output)

dir.create(dir_log, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_png, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rds, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_xlsx, showWarnings = FALSE, recursive = TRUE)





source("./r/utilities_for_sc_analyses.R")



cat(sprintf("\n------------------------------------\n"))
cat(sprintf("parameters\n"))

cat(sprintf("\tn_cores=%s\n", args$n_cores))
cat(sprintf("\tdir_output=%s\n", args$dir_output))
cat(sprintf("\tdir_log=%s\n", args$dir_log))
cat(sprintf("\tdir_seurat_obj=%s\n", args$dir_seurat_obj))
cat(sprintf("\tdir_rds=%s\n", args$dir_rds))
cat(sprintf("\tcancer_type=%s\n", args$cancer_type))

args <- update_args_cancer_type_standard(args)
cat(sprintf("\tcancer_type_standard=%s\n", args$cancer_type_standard))


cat(sprintf("\tfilename_to_select_barcodes=%s\n", args$filename_to_select_barcodes))
cat(sprintf("\tpattern_to_select_cell_types=%s\n", args$pattern_to_select_cell_types))
cat(sprintf("\tmax_n_cells=%d\n", args$max_n_cells))
cat(sprintf("\tmin_n_cells=%d\n", args$min_n_cells))
cat(sprintf("\ttype_parsing_rds_filename=%s\n", args$type_parsing_rds_filename))
cat(sprintf("\tcancer_type_for_parsing_rds_filename=%s\n", args$cancer_type_for_parsing_rds_filename))
cat(sprintf("\tsamples_to_exclude=%s\n", args$samples_to_exclude))


# for seurat
cat(sprintf("\tmethod_sctransform_vst=%s\n", args$method_sctransform_vst))
pattern_seurat_batch_effect_correction <- "cca|rpca|rlsi"
cat(sprintf("\tmethod_integration=%s\n", args$method_integration))
cat(sprintf("\tk.anchor=%d\n", args$k.anchor))
cat(sprintf("\ttype_integration_anchor_features=%s\n", args$type_integration_anchor_features))


# args$method_batch_effect_correction
if (nchar(args$method_batch_effect_correction) == 0) {
	args$method_batch_effect_correction <- args$method_integration
}
if (args$harmony_theta[1] < 0) {
	args$method_batch_effect_correction <- gsub("[,]*harmony|harmony[,]*", "", args$method_batch_effect_correction)
} else {
	if (!grepl("harmony", args$method_batch_effect_correction)) {
		args$method_batch_effect_correction <- paste(args$method_batch_effect_correction, "harmony", sep=",")
	}
} # if


cat(sprintf("\tmethod_batch_effect_correction=%s\n", args$method_batch_effect_correction))
cat(sprintf("\tbatch_keys_for_reference=%s\n", args$batch_keys_for_reference))
cat(sprintf("\tk.weight=%d\n", args$k.weight))






# seurat/scaling
if (nchar(args$vars.to.regress) > 0) {
  args$vars.to.regress <- strsplit(args$vars.to.regress, ",")[[1]]
  cat(sprintf("\tvars.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
} else {
  cat(sprintf("\tvars.to.regress=%s\n", args$vars.to.regress))
  args$vars.to.regress <- NULL
}



# seurat/pca
cat(sprintf("\tmax_dimstouse=%d\n", args$max_dimstouse))



# seurat/clustering
cat(sprintf("\tseurat_resolution=%g\n", args$seurat_resolution))
cat(sprintf("\tseed.seurat_clustering=%g\n", args$seed.seurat_clustering))
str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)
cat(sprintf("\tstr_column_of_meta_data_cluster=%s\n", str_column_of_meta_data_cluster))

if (args$f_multik) {
	str_column_of_meta_data_cluster_multik <- "RNA_multik"
	cat(sprintf("\tstr_column_of_meta_data_cluster_multik=%s\n", str_column_of_meta_data_cluster_multik))
}






# for harmony
cat(sprintf("\tseed.harmony=%g\n", args$seed.harmony))
cat(sprintf("\ttype_parsing_rds_filename_for_dataset=%s\n", args$type_parsing_rds_filename_for_dataset))
cat(sprintf("\ttype_parsing_rds_filename_for_species=%s\n", args$type_parsing_rds_filename_for_species))
cat(sprintf("\ttype_parsing_rds_filename_for_tumor.type=%s\n", args$type_parsing_rds_filename_for_tumor.type))
cat(sprintf("\ttype_parsing_rds_filename_for_donor=%s\n", args$type_parsing_rds_filename_for_donor))
cat(sprintf("\ttype_parsing_rds_filename_for_technology=%s\n", args$type_parsing_rds_filename_for_technology))
cat(sprintf("\ttype_parsing_rds_filename_for_condition=%s\n", args$type_parsing_rds_filename_for_condition))
cat(sprintf("\ttype_parsing_rds_filename_for_treatment=%s\n", args$type_parsing_rds_filename_for_treatment))



args$harmony_theta <- as.numeric(strsplit(args$harmony_theta, ",")[[1]])
cat(sprintf("\tharmony_theta=%s\n", paste(args$harmony_theta, collapse="_")))
str_column_of_meta_data_harmony <- sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse="_"))
cat(sprintf("\tstr_column_of_meta_data_harmony=%s\n", str_column_of_meta_data_harmony))

if (args$f_multik) {
	str_column_of_meta_data_harmony_multik <- sprintf("RNA_harmony_th.%s_multik", paste(args$harmony_theta, collapse="_"))
	cat(sprintf("\tstr_column_of_meta_data_harmony_multik=%s\n", str_column_of_meta_data_harmony_multik))
}


args$harmony_lambda <- as.numeric(strsplit(args$harmony_lambda, ", ")[[1]])
if (args$harmony_lambda[1] < 0) {
        # user did not provide with harmony_lambda
        args$harmony_lambda <- rep(1.0, length(args$harmony_theta))
}
cat(sprintf("\tharmony_lambda=%s\n", paste(args$harmony_lambda, collapse=", ")))


cat(sprintf("\tseed.symphony=%s\n", args$seed.symphony))

# for pca
cat(sprintf("\ttype_signatures=%s\n", args$type_signatures))

# for umap
cat(sprintf("\tumap_n_neighbors=%d\n", args$umap_n_neighbors))
cat(sprintf("\tumap_min_dist=%.1f\n", args$umap_min_dist))
cat(sprintf("\tumap_metric=%s\n", args$umap_metric))

# for cell typing
cat(sprintf("\tmethod_to_identify_cell_types=%s\n", args$method_to_identify_cell_types))
cat(sprintf("\tstr_update_cell_types=%s\n", args$str_update_cell_types))





# set seed for reproducibility
set.seed(51)


# load subroutines

source("./r/rowr.R")
source("./r/batch_correction_seurat.R")
source("./r/batch_correction_harmony.R")
source("./r/enrichment_analysis.R")
source("./r/find_markers.R")
source("./r/identify_cell_types.R")
source("./r/jupyter_message.R")
source("./r/merge_seurat_objects.R")
source("./r/plot_sc_clusters.R")
source("./r/stacked_violin.R")
source("./r/run_multik.R")
source("./r/find_clusters_after_merging.R")
#source("./r/update_rna_pc1_umap_cnv_for_merging.R")


# load libraries

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ConsensusClusterPlus))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(psych))


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(symphony))


source("./r/seurat/seurat_reproducible.R")




if (args$f_multik) {
  suppressPackageStartupMessages(library(MultiK))
} else {
  suppressPackageStartupMessages(library(future))
  plan("multicore", workers = args$n_cores)
  options(future.globals.maxSize = 100 * 1024^3) # 100GB
} # if





merged_sample_id <- sprintf("%s_merged", cancer_type)
















###########################################################
# Part 1: Merge individually processed scRNA-seq datasets
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 1: Merge individually processed scRNA-seq datasets\n"))
###########################################################

# list processed scRNA-seq Seurat objects made in the previous R scripts
filenames_rds <- list.files(path=args$dir_rds, pattern = "*_sample_seurat_obj.rds", full.names = T)



# merge_seurat_objects
list_out <- merge_seurat_objects(filenames_rds, args)
rna <- list_out$rna
args <- list_out$args
df_sample.meta <- list_out$df_sample.meta








##### normalization/variable features/scale data

cat(sprintf("------------------------------------\n"))


# args$list_cancer_type_specific_info 
args$list_cancer_type_specific_info <- get_cancer_type_specific_info(args)

# store mitochondrial percentage in object meta data
rna <- identify_cell_types(rna, "percent.mt", args)



if (nchar(args$method_sctransform_vst) == 0) {

  if (grepl(pattern_seurat_batch_effect_correction, args$method_batch_effect_correction)) {

    # batch correction
    cat(sprintf("\tBatch correction with Seurat\n\n"))
    rna <- batch_correction_with_seurat_lognormalize(rna, args, vars.to.regress = args$vars.to.regress)

    cat(sprintf("\tDrop0\n"))
    rna <- update_seurat_obj_after_batch_correction_seurat(rna)

    #cat(sprintf("\tFindVariableFeatures\n"))
    # Warning message: In FindVariableFeatures.Assay(object = assay.data, selection.method = selection.method,  : selection.method set to 'vst' but count slot is empty; will use data slot instead
    # We do not support the identification of variable features on integrated data. https://github.com/satijalab/seurat/issues/1528
    #rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

    # cell cycle after bactch correction
    rna <- identify_cell_types(rna, "cell_cycle", args)

    cat(sprintf("\tScaleData var.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
    #all.genes <- rownames(rna)
    #rna <- ScaleData(rna, features = all.genes)
    # regressing out percent.mt is slow when features = all.genes
    #rna <- ScaleData(rna, features = all.genes, vars.to.regress = "percent.mt", verbose = FALSE)

    # There is no need to regress out sample-specific attributes like percent.mito after integration, since batch correction takes care of this at the feature level.  https://github.com/satijalab/seurat/issues/3579
    #rna <- ScaleData(rna, verbose = FALSE)

    rna <- ScaleData(rna, vars.to.regress = args$vars.to.regress, verbose = FALSE)

  } else {

    # normalize data without batch correction
    cat(sprintf("\tNormalizeData\n"))
    rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

    cat(sprintf("\tFindVariableFeatures\n"))
    rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

    # cell cycle after normalization
    rna <- identify_cell_types(rna, "cell_cycle", args)

    cat(sprintf("\tScaleData var.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
    #all.genes <- rownames(rna)
    #rna <- ScaleData(rna, features = all.genes)
    # regressing out percent.mt is slow when features = all.genes.
    #rna <- ScaleData(rna, features = all.genes, vars.to.regress = "percent.mt", verbose = FALSE)
    rna <- ScaleData(rna, vars.to.regress = args$vars.to.regress, verbose = FALSE)

  } # if


} else {

  # Perform SCTransform since args$method_sctransform_vst is given.

  # normalize data without batch correction for cell cycle
  cat(sprintf("\tNormalizeData\n"))
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

  # cell cycle
  rna <- identify_cell_types(rna, "cell_cycle", args)


  if (grepl(pattern_seurat_batch_effect_correction, args$method_batch_effect_correction)) {

    # batch correction with SCTransform
    cat(sprintf("\tBatch correction with Seurat SCTransform\n\n"))
    rna <- batch_correction_with_seurat_sctransform(rna, args, vars.to.regress = args$vars.to.regress)

    cat(sprintf("\tDrop0\n"))
    rna <- update_seurat_obj_after_batch_correction_seurat(rna)

  } else {

    # normalize data without batch correction

    cat(sprintf("\tSCTransform\n"))
    # https://satijalab.org/seurat/reference/sctransform
    # https://satijalab.org/seurat/articles/sctransform_vignette.html
    # Apply sctransform normalization
    # * Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
    # * Transformed data will be available in the SCT assay, which is set as the default after running sctransform
    # * During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage

    # The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure. It can be invoked by specifying method="glmGamPoi".
    rna <- SCTransform(rna,
		 method = args$method_sctransform_vst, # default method for sctransform::vst() = "poisson".  https://rdrr.io/cran/sctransform/man/vst.htm
		 vars.to.regress = args$vars.to.regress,
		 verbose = FALSE)

  } # fi

} # if
 







cat(sprintf("------------------------------------\n"))
cat(sprintf("\tRunPCA\n\n"))

npcs <- 50
features_for_pca <- get_genes_from_signatures(args)
if (length(features_for_pca) > 0) {
	cat(sprintf("\t\t# of genes for PCA: %d\n", length(features_for_pca)))
	if ( length(features_for_pca) <= ncol(rna) ) {
		npcs <- min(npcs, length(features_for_pca)-1)
	} else {
		npcs <- min(npcs, ncol(rna)-1)
	}
} # if



# https://satijalab.org/seurat/reference/runpca
rna <- RunPCA(rna,
  assay = NULL,  # Name of Assay PCA is being run on, when=NULL, DefaultAssay(rna)
  features = features_for_pca, # default=NULL, Features to compute PCA on. If features=NULL, PCA will be run using the variable features for the Assay. Note that the features must be present in the scaled data. Any requested features that are not scaled or have 0 variance will be dropped, and the PCA will be run using the remaining features.
  npcs = npcs, # default=50
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = TRUE, # Print the top genes associated with high/low loadings for the PCs
  ndims.print = 1:5, # PCs to print genes for
  nfeatures.print = 30, # Number of genes to print for each PC
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42 # default=42
) # RunPCA

str_reduction <- "pca"









if (grepl("harmony", args$method_batch_effect_correction)) {

  # compute harmony with theta
  cat(sprintf("------------------------------------\n"))
  if (grepl(pattern_seurat_batch_effect_correction, args$method_batch_effect_correction)) {
  	cat(sprintf("\tCompute Harmony\n\n"))
  } else {
	args$method_batch_removal <- "harmony"
  	cat(sprintf("\tBatch correction with Harmony\n\n"))
  }

  rna <- run_harmony(rna, args, vars.to.regress = args$vars.to.regress)
  if ("harmony" %in% names(rna)) {
	str_reduction <- "harmony"
  } else {
  	cat(sprintf("\tskip bath correction with Harmony\n\n"))
  }
  
  if (args$seed.symphony > 0) {
	f_save <- save_symphony_reference_rds(rna, args)
  }

} # if



















# score cells for cell type specific gene signatures
#cat(sprintf("------------------------------------\n"))
#cat(sprintf("\tScore cells for cell type specific gene signaturesl\n"))

#rna <- identify_cell_types(rna, "panglaodb", args)
#rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)



rna <- identify_cell_types(rna, "celltype_cycling", args)









###########################################################
# Part 2: FindClusters
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 2: Find clusters\n"))
###########################################################


rna <- find_clusters_after_merging(rna, args)
#rna <- update_rna_pc1_umap_cnv_for_merging(rna, args)










###########################################################
# Part 3: Cell type annotation of clusters 
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 3: Cell type annotation of clusters\n"))
###########################################################



if (DefaultAssay(rna) == "integrated") {
	str_umap_reduction <- "umap"
	#Idents(rna) <- "RNA_snn_res.0.7"
	Idents(rna) <- str_column_of_meta_data_cluster
} else {
	str_umap_reduction <- "umap"
	Idents(rna) <- str_column_of_meta_data_cluster
	if (grepl("harmony", args$method_batch_effect_correction) && ("harmony" %in% names(rna))) {
		str_umap_reduction <- "umapharmony"
		Idents(rna) <- str_column_of_meta_data_harmony
	}
} # if
cat(sprintf("\tstr_umap_reduction=%s\n", str_umap_reduction))




if (args$f_save_umap_figure) {

  ### Visualize clusters and SingleR annotations 
  cat(sprintf("\tDimPlot\n"))

  #DimPlot(rna, reduction = str_umap_reduction, group.by = "Sample", label = T)+ggsave("Patient_UMAP.pdf",width = 6,height = 4)
  #DimPlot(rna, reduction = str_umap_reduction, group.by = str_column_of_meta_data_cluster, label = T)+ggsave("Cluster_UMAP.pdf",width = 6,height = 4)
  #DimPlot(rna, reduction = str_umap_reduction, group.by = "cell.type",label = T)+ggsave("SingleR_UMAP.pdf",width = 6,height = 4)
  #DimPlot(rna, reduction = str_umap_reduction, group.by = "SingleR.HPCA",label = T)+ggsave("HPCA_UMAP.pdf",width = 8,height = 4)
  #DimPlot(rna, reduction = str_umap_reduction, group.by = "SingleR.BED",label = T)+ggsave("BED_UMAP.pdf",width = 8,height = 4)

  gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = "Sample", label = T)
  #ggsave(sprintf("%s/umap_after_%s_%s_samples.pdf", dir_pdf, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
  ggsave(sprintf("%s/umap_after_%s_%s_samples.png", dir_png, str_reduction, cancer_type), width = 6, height = 4, plot=gg)

  gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = str_column_of_meta_data_cluster, label = T)
  #ggsave(sprintf("%s/umap_after_%s_%s_clusters.pdf", dir_pdf, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
  ggsave(sprintf("%s/umap_after_%s_%s_clusters.png", dir_png, str_reduction, cancer_type), width = 6, height = 4, plot=gg)

  if ("cell.type" %in% names(rna@meta.data)) {
    # merged
    gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = "cell.type", label = T)
    #ggsave(sprintf("%s/umap_after_%s_%s_singler.pdf", dir_pdf, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
    ggsave(sprintf("%s/umap_after_%s_%s_singler.png", dir_png, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
  } # if

  if ("SingleR.HTAPP_toolbox" %in% names(rna@meta.data)) {
    # 1) Slyper et al. Nat. Medicine 2020 scRNA-seq ovarian tumor
    gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = "SingleR.HTAPP_toolbox", label = T)
    #ggsave(sprintf("%s/umap_after_%s_%s_singler_htapp_toolbox.pdf", dir_pdf, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
    ggsave(sprintf("%s/umap_after_%s_%s_singler_htapp_toolbox.png", dir_png, str_reduction, cancer_type), width = 6, height = 4, plot=gg)
  } # if

  if ("SingleR.HPCA" %in% names(rna@meta.data)) {
    # 2) Human Primary Cell Atlas Data (microarray)
    gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = "SingleR.HPCA", label = T)
    #ggsave(sprintf("%s/umap_after_%s_%s_singler_hpca.pdf", dir_pdf, str_reduction, cancer_type), width = 8, height = 4, plot=gg)
    ggsave(sprintf("%s/umap_after_%s_%s_singler_hpca.png", dir_png, str_reduction, cancer_type), width = 8, height = 4, plot=gg)
  }

  if ("SingleR.BED" %in% names(rna@meta.data)) {
    # 3) BluePrint Encode (bulk RNA-seq) 
    gg <- DimPlot(rna, reduction = str_umap_reduction, group.by = "SingleR.BED", label = T)
    #ggsave(sprintf("%s/umap_after_%s_%s_singler_blueprint_encode.pdf", dir_pdf, str_reduction, cancer_type), width = 8,height = 4, plot=gg)
    ggsave(sprintf("%s/umap_after_%s_%s_singler_blueprint_encode.png", dir_png, str_reduction, cancer_type), width = 8,height = 4, plot=gg)
  }


} # if






# get_most_common_cell_type_for_each_cluster

rna$cluster.type <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_cluster)

if (str_column_of_meta_data_harmony %in% colnames(rna@meta.data)) {
	rna$cluster.type.harmony <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_harmony)
}

if (args$f_multik) {

	rna$cluster.type_multik <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_cluster_multik)

	if (str_column_of_meta_data_harmony_multik %in% colnames(rna@meta.data)) {
		rna$cluster.type.harmony_multik <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_harmony_multik)
	}

} # if




# col_seurat_cluster, col_cluster_types
if (DefaultAssay(rna) == "integrated") {
	col_seurat_cluster <- str_column_of_meta_data_cluster
	col_cluster_types <- "cluster.type"
} else {
	col_seurat_cluster <- str_column_of_meta_data_cluster
	col_cluster_types <- "cluster.type"
        if (grepl("harmony", args$method_batch_effect_correction) && ("harmony" %in% names(rna))) {
		col_seurat_cluster <- str_column_of_meta_data_harmony
		col_cluster_types <- "cluster.type.harmony"
        }
} # if
cat(sprintf("\tcol_seurat_cluster=%s\n", col_seurat_cluster))
cat(sprintf("\tcol_cluster_types=%s\n", col_cluster_types))




if (args$f_multik) {
	Idents(rna) <- col_cluster_types_multik
} else {
	Idents(rna) <- col_cluster_types
}












####################################################################

if (args$f_save_umap_figure) {

  #DimPlot(rna, reduction = str_umap_reduction, group.by = col_cluster_types, label = F)+ggsave("cell_type_UMAP.pdf",width = 10,height = 4)
  gg <- print_umap_cell_type(rna, str_umap_reduction = str_umap_reduction, col_cluster_types = col_cluster_types, width=10, height=7, filename=sprintf("%s/umap_after_%s_%s_cell_types.png", dir_png, str_reduction, cancer_type))

} # if




f_plot_cnv_value_corr <- TRUE
if (f_plot_cnv_value_corr) {

	df <- rna@meta.data
	f.cnv <- !is.na(df$CNV.corr) & df$CNV.Pos == "11"
	df <- df[f.cnv,,drop=F]
	f.epi <- grepl(pattern_tumor_epi, df[, col_cluster_types])
	df <- df[f.epi,,drop=F]
	if (nrow(df) > 0) {
        	gg <- ggplot(df, aes_string(x="CNV.value", y="CNV.corr", color=col_cluster_types)) +
			geom_point(size=2.0, alpha=0.5) +
			#scale_color_manual(values=c("cancer"="red", "unassigned"="grey", "normal"="blue")) +
			ggtitle(sprintf("inferCNV")) +
			xlab("CNV values") +
			ylab("kendall tau with top 5% CN altered cells") +
			theme_bw()
        	ggsave(sprintf("%s/scatterplot_infercnv_cna_vs_cor.pdf", dir_pdf), width = 6, height = 5, plot=gg)
		vecx <- layer_scales(gg)$x$range$range
		vecy <- layer_scales(gg)$y$range$range

		cluster_types_epi <- unique(df[, col_cluster_types])
		for (cluster_type_epi in cluster_types_epi) {
			idx <- which(df[, col_cluster_types] == cluster_type_epi)
			df1 <- df[idx,,drop=F]
        		gg <- ggplot(df1, aes(x=CNV.value, y=CNV.corr, color=Sample)) +
				geom_point(size=2.0, alpha=0.5) +
				#scale_color_manual(values=c("cancer"="red", "unassigned"="grey", "normal"="blue")) +
				ggtitle(cluster_type_epi) +
				xlim(vecx) +
				ylim(vecy) +
				xlab("CNV values") +
				ylab("kendall tau with top 5% CN altered cells") +
				theme_bw()

        		ggsave(sprintf("%s/scatterplot_infercnv_cna_vs_cor_%s.pdf", dir_pdf, gsub(" ", "_", tolower(cluster_type_epi))), width = 6, height = 5, plot=gg)
		} # for
	} # if

} # if











# perform DEGs analysis with cell type annotated clusters 
Wilcox.markers <- find_all_markers(rna, args, assay=NULL, col_cluster_types="cluster.type", f_save_rds=TRUE, dir_rds=dir_rds)

if ("cluster.type.harmony" %in% colnames(rna@meta.data)) {
  Wilcox.markers <- find_all_markers(rna, args, assay="RNA", col_cluster_types="cluster.type.harmony", f_save_rds=TRUE, dir_rds=dir_rds)
}






###########################################################
# get_diet_seurat_obj
rna <- get_diet_seurat_obj(rna, args, diet_level=1)




cat(sprintf("------------------------------------\n"))
cat(sprintf("\tprint meta data\n\n"))
# extract meta data
md <- rna@meta.data %>% as.data.table
print(head(md))


cat(sprintf("------------------------------------\n"))
cat(sprintf("\tcount the number of cells per unique combinations of samples and %s\n\n", col_seurat_cluster))

dt <- summarize_meta_data(rna, "Sample", col_seurat_cluster)



cat(sprintf("------------------------------------\n"))
cat(sprintf("\tcount the number of cells per unique combinations of cell types and %s\n\n", col_seurat_cluster))

dt <- summarize_meta_data(rna, "cell.type", col_seurat_cluster)




cat(sprintf("------------------------------------\n"))
cat(sprintf("\tcount epithelial cell numbers with %s\n\n", "cell.type"))
idx <- grep(pattern_epi, rna@meta.data[, "cell.type"])
table(rna@meta.data[idx, "cell.type"])




cat(sprintf("------------------------------------\n"))
cat(sprintf("\tnumber of samples for each %s\n\n", col_cluster_types))

dt <- summarize_meta_data(rna, col_cluster_types, "Sample")











###########################################################
# save seurat object 
cat(sprintf("------------------------------------\n"))
cat(sprintf("\tsave seurat object\n"))
path_seurat_obj <- sprintf("%s/%s_sc-rna-seq_merged_seurat_obj.rds", dir_rds, cancer_type)
cat(sprintf("\tsaveRDS(rna, '%s')\n", path_seurat_obj))
saveRDS(rna, path_seurat_obj)























###########################################################
cat(sprintf("------------------------------------\n"))
cat(sprintf("\twrite.xlsx\n"))


df_merged.meta <- data.frame(TotalCells = length(colnames(rna)),
                          NumClusters = length(levels(as.factor(Idents(rna)))),
                          stringsAsFactors = FALSE)

if ("cluster.type.harmony" %in% colnames(rna@meta.data)) {
  df_merged.meta$NumClustersHarmony <- length(levels(as.factor(rna$cluster.type.harmony)))
}

if ("cluster.type_multik" %in% colnames(rna@meta.data)) {
  df_merged.meta$NumClustersMultiK <- length(levels(as.factor(rna$cluster.type_multik)))
}

if ("cluster.type.harmony_multik" %in% colnames(rna@meta.data)) {
  df_merged.meta$NumClustersHarmonyMultiK <- length(levels(as.factor(rna$cluster.type.harmony_multik)))
}


df_merged.meta <- as.data.frame(t(df_merged.meta))
colnames(df_merged.meta) <- merged_sample_id

wb <- createWorkbook()
n_sheet <- 0

n_sheet <- n_sheet + 1
addWorksheet(wb, sheetName = "sample", gridLines = TRUE)
writeDataTable(wb, sheet = n_sheet, x = df_sample.meta, colNames = TRUE, rowNames = TRUE)

n_sheet <- n_sheet + 1
addWorksheet(wb, sheetName = "merged", gridLines = TRUE)
writeDataTable(wb, sheet = n_sheet, x = df_merged.meta, colNames = TRUE, rowNames = TRUE)
#setColWidths(wb, sheet = n_sheet, cols = 1, widths = 50)

filename_xlsx <- sprintf("%s/%s_sc-rna-seq_pipeline_summary.xlsx", dir_xlsx, cancer_type)
#openxlsx::write.xlsx(df_merged.meta, sprintf("%s/%s_sc-rna-seq_pipeline_summary.xlsx", dir_xlsx, cancer_type), row.names = T, col.names = TRUE)
saveWorkbook(wb, filename_xlsx, overwrite = TRUE)

print(df_merged.meta)







if (args$f_cellphonedb) {
	# cellphonedb
	cat(sprintf("------------------------------------\n"))
	cat(sprintf("\tcellphonedb\n\n"))

	source("r/run_cellphonedb.R")
	out <- run_cellphonedb(rna, args, "cell.type", cancer_type, max_n_cells_for_cellphonedb = 3000)
} # if






# diet seurat object for deeper levels
if (args$diet_seurat_level > 1) {

	cat(sprintf("\n\tdiet_seurat_level=%d\n", args$diet_seurat_level))
	rna <- get_diet_seurat_obj(rna, args)

	dir_rds_diet <- sprintf("%s/diet", dir_rds)
	dir.create(dir_rds_diet, showWarnings = FALSE, recursive = TRUE)
	path_seurat_obj <- sprintf("%s/%s_sc-rna-seq_sample_seurat_obj.rds", dir_rds_diet, sample_id)
	cat(sprintf("\tsaveRDS(rna, '%s')\n", path_seurat_obj))
	saveRDS(rna, path_seurat_obj)

} # if



# post-process
unlink("Rplots.pdf")


cat(sprintf("\n------------------------------------\n"))
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("completed successfully~\n"))
cat(sprintf("------------------------------------\n\n"))





