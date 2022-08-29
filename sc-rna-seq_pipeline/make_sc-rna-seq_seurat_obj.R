#!/usr/bin/env Rscript
#
# make_sc-rna-seq_seurat_obj.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, May
#
# options:
#   --n_cores: (default=8) the number of cores for parallel computing
#   --dir_count: (default="../count")
#   --gene.column: (default=2) Specify which column of genes.tsv or features.tsv to use for gene names
#   --cell.column: (default=1) Specify which column of barcodes.tsv to use for cell names
#   --dir_otuput: (default="./output")
#   --dir_log: (default="./output/log")
#   --dir_seurat_obj: (default="./output/rds")
#
#   --type_qc: {"arguments", ["mad_arguments"]}, type of quality control
#   --min_ncount_rna: (default=-1) minimal # of RNA counts
#   --min_nfeature_rna: (default=-1) minimal $ of detected genes
#   --th_percent.mt: (default=25) # less than 25% mitochondrial counts
#
#   --n_mad_log_counts: (default=2) # remove outliers < -2mad
#   --n_mad_log_features: (default=2) # remove outliers < -2mad
#   --n_mad_log1p_mito: (default=2) # remove outliers > 2mad
#
#   --vars.to.regress: (default="percent.mt") e.g. "nCount_RNA,nFeature_RNA,percent.mt,S.Score,G2M.Score"
#   --method_sctransform_vst: (default="") {"poisson", "glmGamPoi", ...}, use sctransform::vst() when method_sctransform_vst is not empty
#   --max_dimstouse: (default=30) dimsToUse <- 1:args$max_dimstouse
#   --seurat_resolution: (default: 0.8) Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.   https://satijalab.org/seurat/reference/findclusters
#   --seed.seurat_clustering: (default: 0) 
#
#   --umap_n_neighbors: (defualt=30)
#   --umap_min_dist: (default=0.3)
#   --umap_metric: (default="cosine")
#
#   --method_to_identify_cell_types: {"singler_htapp_toolbox", "singler_hpca_cellidx", ["singler_blueprint_encode"]}
#   --method_to_identify_reference_clusters: {"none", "immune_plasma", "nonepi", "nonepi_krt_epcam", "nonepi_singler", "immune_singler", ["immune_plus_endo_singler"], "normal-like_singler", "normal_sample_epi", "normal-like"}, "none" to skip infercnv
#   --method_to_update_cell_types: {[epithelial_cell_types], normal_epithelial_cell_types, cancer_epithelial_cell_types, cancer_normal_epithelial_cell_types}
#   --method_to_identify_subtypes: {"brca_cell_atlas_scsubtyper", "brca_pam50", "brca_pam50_genefu", ["none"]}
#   --type_centroids: {["4centroids"], "5centroids"}
#   --prefix_reference_clusters: (default="auto")
#
#   --doubletdecon.rhop: (default=1.1) try rhop=0.5 to remove the error message of "no locations are finite".  https://github.com/EDePasquale/DoubletDecon/issues/26
#   --doublet.rate: (default=-1) the doublet formation rates based on the number of cells recovered. For the number of cells recovered, you can use the estimated number of cells from the Cell Ranger output.  https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
#   --max_n_cells: (default=Inf) select n cells with higher number of genes expressed
#   --f_save_doublet_png: (default=FALSE) flag to save double png
#   --f_save_doublet_rds: (default=FALSE) flag to save double rds
#
#   --f_save_jackstraw_png: (default=FALSE) flag to save jackstraw png
#   --f_save_cnv_pc1_boxplot_png: (default=FALSE) flag to save cnv_pc1_boxplot.png
#   --f_infercnv_observations_only_epithelial: (default=TRUE)
#   --method_barcodes_reference: {"nonepi_krt_epcam", "normal_epithelial_cells", ["nonepi_singler"]}, method to select reference cells
#   --method_barcodes_observation: {"epi_krt_epcam", "tumor_epithelial_cells", ["epi_singler"]}, additional filter to select observation barcodes for infercnv
#   --min_num_cells_per_cluster: (default=10)
#   --type_infercnv_argset: {["vignettes"], "wiki", "matt", "default", "none"}, "none" to skip infercnv
#   --infercnv_cutoff: (default=0.1) infercnv cutoff, use 1 for smart-seq, 0.1 for 10x-genomics
#   --method_to_update_seurat_obj_with_infercnv: {["peroulab"], "cna_cor_per_cluster", "default"}
#   --infercnv_pos_notpos: classify tumor and non-tumor with infercnv, otherwise classify tumor, non-tumor, and unassigned (default).
#
#   --method_to_determine_th_cna_value_corr: {["ng_2021_and_thresholding", "ng_2021", "fixed"} method_to_determine_th_cna_value_cor, use --th_cna_value and --th_cna_corr when --method_to_determine_th_cna_value_corr "fixed"
#   --th_cna_value: (default=0.05) threshold value to determine the first character of CNV.Pos (1x: CNV_value > th_cna_value)
#   --th_cna_corr: (default=0.4) threshold value to determine the second character of CNV.Pos (x1: cor.estimate > th_cna_corr)
#
#   --bayes_max_pnormal: (default=0.5) By default this threshold is set to 0.5, given this any CNV region that has a posterior probability of being of a normal state greater than 0.5 is relabeled as "normal" and no longer considered an identified CNV region.  https://github-wiki-see.page/m/broadinstitute/infercnv/wiki/Infercnv-i6-HMM-type
#   --min_prob_normal: (default=0.05) call tumor cells if probability of state 3: 1x: neutral < min_prob_normal
#
#   --f_save_infercnv_rds: (default=FALSE) flag to save infercnv rds
#
#   --cancer_type_standard: standard symbols of cancer type (e.g. "bc" for  brca1-mut-bc, er+-bc, her2+-bc, male-bc, normal-breast, tnbc, tr-bc)
#   --cancer_subtype: cancer subtype (e.g. idc, ilc, mbc, nm, dcis, ln_mets), NM: nonmalignant, DCIS: ductual carinoma in situ, IDC: invasiv ductal carcinoma, LN mets: lymph node metastases
#   --no_diet_seurat: do not diet seurat obj, but store full seurat obj
#   --diet_seurat_level: (default=1)
#
#
# input:
#   cancer_type: (e.g. brca1-mut-bc, er+-bc, her2+-bc, male-bc, normal-breast, tnbc, tr-bc, matt-oc)
#   sample_id: (e.g. 3BAE2L, 3E5CFL)
#
# usage:
# ./make_sc-rna-seq_seurat_obj.R male-bc 446B7L
# ./make_sc-rna-seq_seurat_obj.R male-bc 4CC61L
# ./make_sc-rna-seq_seurat_obj.R matt-oc 3BAE2L
# ./make_sc-rna-seq_seurat_obj.R matt-oc 3E5CFL
#
# output:
# ./$dir_seurat_obj/${sample_id}_sc-rna-seq_seurat_obj.rds
#
# reference:
# this script is based on the origical script written by Matt Regner  https://github.com/RegnerM2015/scENDO_scOVAR_2020
# https://satijalab.org/
#
# 0. tutorials
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# https://ucdavis-bioinformatics-training.github.io/2019-single-cell-RNA-sequencing-Workshop-UCD_UCSF/scrnaseq_analysis/scRNA_Workshop-PART5.html
# https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_integ_clustering_analysis.html
#
# 1. gene signatures
# for immune, reference/ESTIMATE_signatures.csv  https://www.frontiersin.org/articles/10.3389/fonc.2019.01212/full
# for other cell types, PanglaoDB_markers_27_Mar_2020.tsv.gz  https://panglaodb.se/, https://zenodo.org/record/4438615#.YUnuDP5KiX3
#
# 2. copy number estination
# 1) InferCNV  https://bioconductor.org/packages/devel/bioc/vignettes/infercnv/inst/doc/inferCNV.html
#
# 3. cell type estimation with SingleR
# 1) A single-cell and single-nucleus RNA-Seq toolbox for fresh and frozen human tumors, Slyper et al., Nat. Medicine, 2020  https://www.nature.com/articles/s41591-020-0844-1, Finally, we provide a website that displays a comprehensive analysis summary for each sample tested (https://tumor-toolbox.broadinstitute.org).
# 2) Human Primary Cell Atlas Data (microarray)
# 3) BluePrint Encode (bulk RNA-seq)
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
parser$add_argument("-dc", "--dir_count", default="../count", type="character", help="direcotry for directories of samples")
parser$add_argument("-gc", "--gene.column", default=2, type="integer", help="gene column")
parser$add_argument("-cc", "--cell.column", default=1, type="integer", help="cell column")
parser$add_argument("-do", "--dir_output", default="./output", type="character", help="direcotry for output")
parser$add_argument("-dl", "--dir_log", default="./output/log", type="character", help="direcotry for log")
#parser$add_argument("-ff", "--figure_format", default="pdf", type="character", help="figure format") # coded to print figures with pdf format (pdf and png whenever it is possible since png is helpful for lighter PPT slides).


# arguments for seurat
parser$add_argument("-ds","--dir_seurat_obj", default="./output/rds", type="character", help="direcotry for seurat objects")

# seurat/filtering cells
parser$add_argument("-tqc","--type_qc", default="mad_arguments", type="character", help="type of quality control")
parser$add_argument("-minc","--min_ncount_rna", default=-1, type="double", help="minimum of nCount_RNA")
parser$add_argument("-minf","--min_nfeature_rna", default=-1, type="double", help="minimum of nFeature_RNA")
parser$add_argument("-thpmt","--th_percent.mt", default=25, type="double", help="threshold for mitochrondrial count percentage")

parser$add_argument("-madc","--n_mad_log_counts", default=2, type="double", help="remove outliers <-2mad")
parser$add_argument("-madf","--n_mad_log_features", default=2, type="double", help="remove outliers <-2mad")
parser$add_argument("-madm","--n_mad_log1p_mito", default=2, type="double", help="remove outliers >2mad")

# seurat/normalization
parser$add_argument("-vr","--vars.to.regress", default="percent.mt", type="character", help="vars.to.regress for ScaleData()")
parser$add_argument("-mv","--method_sctransform_vst", default="", type="character", help="method for sctransform::vst()")

# seurat/pca
parser$add_argument("-maxdim","--max_dimstouse", default=30, type="integer", help="max dimsToUse")

# seurat/clustering
parser$add_argument("-res","--seurat_resolution", default=0.8, type="double", help="Seurat FindClusters parameter of resolution")
parser$add_argument("-ssc","--seed.seurat_clustering", default=0, type="double", help="Seurat FindClusters seed")

# arguments for umap
parser$add_argument("-nn","--umap_n_neighbors", default=30, type="integer", help="UMAP number of neighbors")
parser$add_argument("-dist","--umap_min_dist", default=0.3, type="double", help="UMAP min distance")
parser$add_argument("-metric","--umap_metric", default="cosine", type="character", help="UMAP metrics")

# arguments for the identification of reference clusters
parser$add_argument("-mict","--method_to_identify_cell_types", default="singler_blueprint_encode", type="character", help="method to identify cell types")
parser$add_argument("-i","--method_to_identify_reference_clusters", default="immune_plus_endo_singler", type="character", help="method to identify reference clusters")
parser$add_argument("-r","--prefix_reference_clusters", default="auto", type="character", help="prefix for reference clusters")


# arguments for cell type identification
parser$add_argument("-uct","--method_to_update_cell_types", default="epithelial_cell_types", type="character", help="method to update cell types")
parser$add_argument("-mist","--method_to_identify_subtypes", default="none", type="character", help="method to identify subtypes")
parser$add_argument("-tc","--type_centroids", default="4centroids", type="character", help="type of centroids")



# arguments for doublet removal
parser$add_argument("-mdr","--method_to_remove_doublets", default="DoubletFinder", type="character", help="method to remove doublts")

# arguments for doubletdecon
parser$add_argument("-dds","--doubletdecon.species", default="hsa", type="character", help="doubletdecon species")
parser$add_argument("-rhop","--doubletdecon.rhop", default=1.1, type="double", help="doubletdecon rhop") # take the trial and error approach, such as trying rhop values of 0.5, 1, and 1.5. If any of these work then I will usually tweak them up or down until I get the software to generate a cluster merging heatmap that shows the merging I expect (based on my knowledge of the cell types in the data — more on that in a protocol that is currently in progress!). It would be best to kill the run after the heatmap pops up to save you time. Again, this is why the DoubletDeconUI for Windows would be tremendously helpful!  https://github.com/EDePasquale/DoubletDecon/issues/26

# arguments for doubletfinder
parser$add_argument("-d","--doublet.rate", default=-1, type="double", help="doublet rate") # the doublet formation rates based on the number of cells recovered. For the number of cells recovered, you can use the estimated number of cells from the Cell Ranger output.  https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
parser$add_argument("-nmc", "--max_n_cells", default=1e9, type="double", help="select n cells with higher number of genes expressed")
parser$add_argument("-sdp", "--f_save_doublet_png", default=FALSE, type="logical", help="save doublet png")
parser$add_argument("-sdr", "--f_save_doublet_rds", default=FALSE, type="logical", help="save doublet rds")







# arguments for infercnv
parser$add_argument("-sjp","--f_save_jackstraw_png", default=FALSE, type="logical", help="save jackstraw png")
parser$add_argument("-scpbp","--f_save_cnv_pc1_boxplot_png", default=FALSE, type="logical", help="save cnv_pc1_boxplot.png")
parser$add_argument("-foe","--f_infercnv_observations_only_epithelial", default=TRUE, type="logical", help="InferCNV observations for only epithelial cells")
parser$add_argument("-mr","--method_barcodes_reference", default="nonepi_singler", type="character", help="method to select reference cells") # {"nonepi_krt_epcam", "normal_epithelial_cells", ["nonepi_singler"]}
parser$add_argument("-mo","--method_barcodes_observation", default="epi_singler", type="character", help="method to select observation cells") # {"epi_krt_epcam", "tumor_epithelial_cells", ["epi_singler"]}

parser$add_argument("-m","--min_num_cells_per_cluster", default=10, type="integer", help="minimum # of cells to define a cluster for InferCNV")
parser$add_argument("-t","--type_infercnv_argset", default="vignettes", type="character", help="InferCNV argument set")
parser$add_argument("-c","--infercnv_cutoff", default=0.1, type="double", help="InferCNV cutoff, use 1 for smart-seq, 0.1 for 10x-genomics")
parser$add_argument("-muso","--method_to_update_seurat_obj_with_infercnv", default="peroulab", type="character", help="method to update seurat obj with infercnv")
parser$add_argument("-ipnp","--infercnv_pos_notpos", dest="f_infercnv_pos_notpos", action="store_true", default=FALSE, help="classify tumor and non-tumor with infercnv, otherwise classify tumor, non-tumor, and unassigned (default)")




# for infercnv_wo_hmm
parser$add_argument("-mtdth","--method_to_determine_th_cna_value_corr", default="ng_2021_and_thresholding", type="character", help="method_to determine th_cna_value and th_cna_corr")
parser$add_argument("-tcv","--th_cna_value", default=0.05, type="double", help="threshold value to determine the first character of CNV.Pos (1x: CNV_value > th_cna_value)")
parser$add_argument("-tcc","--th_cna_corr", default=0.4, type="double", help="threshold value to determine the second character of CNV.Pos (x1: cor.estimate > th_cna_corr)")




# for infercnv_hmm
parser$add_argument("-b","--bayes_max_pnormal", default=0.5, type="double", help="maximum P(Normal) allowed for a CNV prediction according to BayesNet. (de-fault=0.5, note zero turns it off)sim_methodmethod for calibrating CNV levels in the i6 HMM (default: ’meanvar’)sim_foregrounddon’t use... for debugging, developer option.") # CNV regions identified by the HMM are filtered out if the CNV region's posterior probability of being normal exceeds a specified threshold. This combats possibility of miss identified CNV's by removing CNV's that are most likely to be normal and not a true CNV events. By default this threshold is set to 0.5, given this any CNV region that has a posterior probability of being of a normal state greater than 0.5 is relabeled as "normal" and no longer considered an identified CNV region.  https://github-wiki-see.page/m/broadinstitute/infercnv/wiki/Infercnv-i6-HMM-type
parser$add_argument("-n","--min_prob_normal", default=0.05, type="double", help="call tumor cells if probability of state 3: 1x: neutral < min_prob_normal")




# for infercnv amp/del definition
parser$add_argument("-mamp", "--min_amp", default=1.20, type="double", help="infercnv modified expression minimal value to call copy number amplification")
parser$add_argument("-mdel", "--max_del", default=0.80, type="double", help="infercnv modified expression maximal value to call copy number deletion")




# for infercnv save rds
parser$add_argument("-sir","--f_save_infercnv_rds", default=FALSE, type="logical", help="save InferCNV rds")




# arguments for others
parser$add_argument("-cts","--cancer_type_standard", default="", type="character", help="standard cancer type")
parser$add_argument("-cst","--cancer_subtype", default="", type="character", help="cancer subtype (e.g. idc, ilc, mbc, nm, dcis, ln_mets)")
parser$add_argument("-cpdb", "--cellphonedb", dest="f_cellphonedb", action="store_true", default=FALSE, help="execute cellphonedb")
parser$add_argument("-nds", "--no_diet_seurat", dest="f_diet_seurat", action="store_false", default=TRUE, help="do not diet seurat object for reducing memory/disk usage")
parser$add_argument("-dsl", "--diet_seurat_level", default=1, type="double", help="diet seurat object level")

parser$add_argument("cancer_type", nargs=1, help="cancer type")
parser$add_argument("sample_id", nargs=1, help="sample_id")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()


dir_output <- args$dir_output
#figure_format <- args$figure_format
dir_seurat_obj <- args$dir_seurat_obj
cancer_type <- args$cancer_type


cat(sprintf("------------------------------------\n"))
#cat(sprintf("./make_sc-rna-seq_seurat_obj.R %s %s\n", args$cancer_type, args$sample_id))
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
dir_tsv <- sprintf("%s/tsv", dir_output)
dir_xlsx <- sprintf("%s/xlsx", dir_output)

dir.create(dir_log, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_pdf, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_png, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_rds, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tsv, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_xlsx, showWarnings = FALSE, recursive = TRUE)





source("./r/utilities_for_sc_analyses.R")



cat(sprintf("\n------------------------------------\n"))
cat(sprintf("parameters\n"))


cat(sprintf("\tn_cores=%s\n", args$n_cores))
cat(sprintf("\tdir_count=%s\n", args$dir_count))
cat(sprintf("\tgene.column=%d\n", args$gene.column))
cat(sprintf("\tcell.column=%d\n", args$cell.column))
cat(sprintf("\tdir_output=%s\n", args$dir_output))
cat(sprintf("\tdir_log=%s\n", args$dir_log))
cat(sprintf("\tdir_seurat_obj=%s\n", args$dir_seurat_obj))
cat(sprintf("\tdir_rds=%s\n", args$dir_rds))
cat(sprintf("\tcancer_type=%s\n", args$cancer_type))

args <- update_args_cancer_type_standard(args)
cat(sprintf("\tcancer_type_standard=%s\n", args$cancer_type_standard))
cat(sprintf("\tcancer_subtype=%s\n", args$cancer_subtype))


# for filtering cells


# for seurat

if (nchar(args$vars.to.regress) > 0) {
	args$vars.to.regress <- strsplit(args$vars.to.regress, ",")[[1]]
	cat(sprintf("\tvars.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
} else {
  cat(sprintf("\tvars.to.regress=%s\n", args$vars.to.regress))
  args$vars.to.regress <- NULL
}

cat(sprintf("\tmethod_sctransform_vst=%s\n", args$method_sctransform_vst))
cat(sprintf("\tmax_dimstouse=%d\n", args$max_dimstouse))

cat(sprintf("\tseurat_resolution=%.1f\n", args$seurat_resolution))
cat(sprintf("\tseed.seurat_clustering=%g\n", args$seed.seurat_clustering))


# for umap
cat(sprintf("\tumap_n_neighbors=%d\n", args$umap_n_neighbors))
cat(sprintf("\tumap_min_dist=%.1f\n", args$umap_min_dist))
cat(sprintf("\tumap_metric=%s\n", args$umap_metric))



# for the identification of reference clusters
cat(sprintf("\tmethod_to_identify_cell_types=%s\n", args$method_to_identify_cell_types))
cat(sprintf("\tmethod_to_identify_reference_clusters=%s\n", args$method_to_identify_reference_clusters))
if (args$prefix_reference_clusters == "auto") {
	switch(args$method_to_identify_reference_clusters,
		"immune_plasma"={ args$prefix_reference_clusters <- "immune." },
		"nonepi"={ args$prefix_reference_clusters <- "non-epi." },
		"nonepi_krt_epcam"={ args$prefix_reference_clusters <- "non-epi." },
		"nonepi_singler"={ args$prefix_reference_clusters <- "non-epi." },
		"immune_singler"={ args$prefix_reference_clusters <- "immune." },
		"immune_plus_endo_singler"={ args$prefix_reference_clusters <- "immendo." },
		"normal-like_singler"={ args$prefix_reference_clusters <- "normal-like." },
		{
			# default way to determine reference.clusters to bypass the reference
			# i.e. immune_plus_endo_singler	
			args$prefix_reference_clusters <- "immendo."
		}
	) # switch
} # if
cat(sprintf("\tprefix_reference_clusters=%s\n", args$prefix_reference_clusters))
cat(sprintf("\tmethod_to_update_cell_types=%s\n", args$method_to_update_cell_types))
cat(sprintf("\tmethod_to_identify_subtypes=%s\n", args$method_to_identify_subtypes))
cat(sprintf("\ttype_centroids=%s\n", args$type_centroids))

# doublet removal
cat(sprintf("\tdoubletdecon.rhop=%.1f\n\n", args$doubletdecon.rhop))
cat(sprintf("\tmax_n_cells=%d\n", args$max_n_cells))

# infercnv
cat(sprintf("\tf_infercnv_observations_only_epithelial=%s\n", args$f_infercnv_observations_only_epithelial))
cat(sprintf("\tmethod_barcodes_reference=%s\n", args$method_barcodes_reference))
cat(sprintf("\tmethod_barcodes_observation=%s\n", args$method_barcodes_observation))
cat(sprintf("\tmin_num_cells_per_cluster=%d\n", args$min_num_cells_per_cluster))
cat(sprintf("\ttype_infercnv_argset=%s\n", args$type_infercnv_argset))
cat(sprintf("\tinfercnv_cutoff=%s\n", args$infercnv_cutoff))
cat(sprintf("\tmethod_to_update_seurat_obj_with_infercnv=%s\n", args$method_to_update_seurat_obj_with_infercnv))
cat(sprintf("\tf_infercnv_pos_notpos=%s\n", args$f_infercnv_pos_notpos))

# for infercnv_wo_hmm
cat(sprintf("\tmethod_to_determine_th_cna_value_corr=%s\n", args$method_to_determine_th_cna_value_corr))
cat(sprintf("\tth_cna_value=%g\n", args$th_cna_value))
cat(sprintf("\tth_cna_corr=%g\n", args$th_cna_corr))

# for infercnv_hmm
cat(sprintf("\tbayes_max_pnormal=%.2f\n", args$bayes_max_pnormal))
cat(sprintf("\tmin_prob_normal=%.2f\n", args$min_prob_normal))




# set seed for reproducibility
set.seed(51)



# load subroutines

source("./r/rowr.R")
source("./r/convert_symbol.R")
source("./r/enrichment_analysis.R")
source("./r/find_markers.R")
source("./r/identify_cell_types.R")
source("./r/identify_reference_clusters.R")
source("./r/jupyter_message.R")
source("./r/merge_seurat_objects.R")
source("./r/monocle3_garnett_utilities.R")
source("./r/markercount_utilities.R")
source("./r/run_infercnv.R")


source("./r/find_clusters_before_doublet_removal.R")
source("./r/find_clusters_after_doublet_removal.R")
#source("./r/update_rna_pc1_umap_cnv_before_doublet_removal.R")
#source("./r/update_rna_pc1_umap_cnv_after_doublet_removal.R")







# for general utils

#suppressPackageStartupMessages(library(base))  # R.oo::attach(), R.oo::detach(), R.oo::load), R.oo::save(), R.utils::cat(), R.utils::commandArgs(), R.utils::getOption(), R.utils::inherits(), R.utils::isOpen(), R.utils::nullfile(), R.utils::parse(), R.utils::warnings(), spam::backsolve(), spam::forwardsolve()
#suppressPackageStartupMessages(library(methods))  # R.oo::getClasses(), R.oo::getMethods()
#suppressPackageStartupMessages(library(R.methodsS3))  # R.oo::throw()

suppressPackageStartupMessages(library(R.utils))  # spam::cleanup()
#suppressPackageStartupMessages(library(tidyr)) # S4Vectors::expand(), R.utils::extract()
# You have loaded plyr after dplyr - this is likely to cause problems.
# If you need functions from both plyr and dplyr, please load plyr first, then dplyr: library(plyr); library(dplyr)
suppressPackageStartupMessages(library(plyr)) 
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(utils))  # R.utils::timestamp()



# for data structure
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(data.table))



# for multivariate analysis
suppressPackageStartupMessages(library(stats)) # pcaMethods::loadings(), gplots::lowess()
#suppressPackageStartupMessages(library(stats4)) # spam::mle()
suppressPackageStartupMessages(library(psych)) # pcaMethods::pca(), fields::describe()
suppressPackageStartupMessages(library(ICIKendallTau)) 



# for display
#suppressPackageStartupMessages(library(IRdisplay)) # spam::display()



# for graph 
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ComplexHeatmap)) # R.utils::draw()
suppressPackageStartupMessages(library(ConsensusClusterPlus))
suppressPackageStartupMessages(library(scales)) # viridis::viridis_pal()
suppressPackageStartupMessages(library(ggplot2)) # limSolve::resolution()



# for enrichment analysis
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(clusterProfiler))



# for single cell analysis
suppressPackageStartupMessages(library(SingleCellExperiment)) # R.oo::trim(), GenomicRanges::trim(), IRanges::trim() 
suppressPackageStartupMessages(library(scater)) # DeconRNASeq::multiplot()
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(infercnv))
suppressPackageStartupMessages(library(DoubletDecon))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(SingleR))


source("./r/seurat/seurat_reproducible.R")



# for parallel computing
suppressPackageStartupMessages(library(future))
#suppressPackageStartupMessages(library(purrr)) # Signac::reduce(), scales::discard(), data.table::transpose(), plyr::compact(), clusterProfiler::simplify(), GenomicRanges::reduce(), IRanges::reduce()
#suppressPackageStartupMessages(library(furrr))

plan("multicore", workers = args$n_cores)
options(future.globals.maxSize = (min(8, args$n_cores)*1024^3)) # <= 8GB










###########################################################
# Part 1: scRNA-seq processing before doublet removal
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 1: scRNA-seq processing before doublet removal\n"))
cat(sprintf("------------------------------------\n"))
###########################################################



# Load the RNA dataset
if (grepl(".txt.gz", args$sample_id)) {

	# read gene x cell count matrix
	filename <- sprintf("%s/%s", args$dir_count, args$sample_id)
	cat(sprintf("\tread.table(%s)\n", filename))
	df <- read.table(filename, sep="\t", header=T, row.names=1, check.names=FALSE)
	counts <- as(as.matrix(df), "dgCMatrix")

	sample_id <- gsub("_Expression_Matrix.txt.gz", "", args$sample_id)
	sample_id <- sprintf("%s_%s", args$cancer_type, sample_id)
	

} else {

	sample_id <- sprintf("%s_%s", args$cancer_type, args$sample_id)
	dir_sample <- args$sample_id

	#counts <- Read10X_h5(filename=sprintf("../count/%s/outs/filtered_feature_bc_matrix.h5", dir_sample))
	file_name_h5 <- sprintf("%s/%s/outs/filtered_feature_bc_matrix.h5", args$dir_count, dir_sample)
	if (file.exists(file_name_h5)) {
		cat(sprintf("\tRead10X_h5(%s)\n", file_name_h5))
		counts <- Read10X_h5(file_name_h5)
	} else {
		dir_sample1 <- sprintf("%s/%s", args$dir_count, dir_sample)
		cat(sprintf("\tRead10X(%s)\n", dir_sample1))
		counts <- Read10X(dir_sample1, gene.column = args$gene.column, cell.column = args$cell.column)
	}

} # if


if (!is.null(names(counts))) {
	if ("Gene Expression" %in% names(counts)) {
		# for multiome kit
		counts <- counts[["Gene Expression"]]
	}
}

# Initialize the Seurat object with the raw (non-normalized data).
rna <- CreateSeuratObject(
	counts = counts, # unnormalize data such as raw counts or TPMs
	project = sample_id,
	assay = "RNA",
	min.cells = 3, # genes not present in at least 3 cells are removed
	min.features = 0,
	names.field = 1,
	names.delim = "_",
	meta.data = NULL)


# print information of seurat obj
#
# e.g.)
# An object of class Seurat
# 17704 features across 6468 samples within 1 assay
# Active assay: RNA (17704 features, 0 variable features)

log_obj(rna, tab=2)

n_cells <- ncol(rna)

if (args$doublet.rate < 0) {
  # automatic assignment of doublet.rate with # of cells
  # https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x.php
  #       Multiplet rate (%)      # of Cell Loaded        # of Cell Recovered
  #       0.40%   800     500
  #       0.80%   1,600   1,000
  #       1.60%   3,200   2,000
  #       2.30%   4,800   3,000
  #       3.10%   6,400   4,000
  #       3.90%   8,000   5,000
  #       4.60%   9,600   6,000
  #       5.40%   11,200  7,000
  #       6.10%   12,800  8,000
  #       6.90%   14,400  9,000
  #       7.60%   16,000  10,000
  if (n_cells > 9500) { args$doublet.rate <- 0.076 }
  else if (n_cells > 8500) { args$doublet.rate <- 0.069 } 
  else if (n_cells > 7500) { args$doublet.rate <- 0.061 } 
  else if (n_cells > 6500) { args$doublet.rate <- 0.054 } 
  else if (n_cells > 5500) { args$doublet.rate <- 0.046 } 
  else if (n_cells > 4500) { args$doublet.rate <- 0.039 } 
  else if (n_cells > 3500) { args$doublet.rate <- 0.031 } 
  else if (n_cells > 2500) { args$doublet.rate <- 0.023 } 
  else if (n_cells > 1500) { args$doublet.rate <- 0.016 } 
  else if (n_cells > 750) { args$doublet.rate <- 0.008 } 
  else { args$doublet.rate <- 0.004 }
} # if

cat(sprintf("\tn_cells: %d\n", n_cells))
cat(sprintf("\tdoublet.rate=%.3f\n\n", args$doublet.rate))



PreQCNumCells <- length(colnames(rna))
cat(sprintf("\t# of cells at pre-QC: %d\n", PreQCNumCells))





# args$list_cancer_type_specific_info 
args$list_cancer_type_specific_info <- get_cancer_type_specific_info(args)

# store mitochondrial percentage in object meta data
rna <- identify_cell_types(rna, "percent.mt", args)




# print stats
cat(sprintf("\tnCount_RNA distribution:\n"))
log_obj(round(quantile(rna$nCount_RNA, probs=seq(0.1,0.9,by=0.1)),1), tab=2)
cat(sprintf("\tnFeature_RNA distribution:\n"))
log_obj(round(quantile(rna$nFeature_RNA, probs=seq(0.1,0.9,by=0.1)),1), tab=2)
cat(sprintf("\tpercent.mt distribution:\n"))
log_obj(round(quantile(rna$percent.mt, probs=seq(0.1,0.9,by=0.1)),1), tab=2)




# filtering cells
# QC metrics: nCount_RNA, nFeature_RNA, and percent mitochrondial counts


# arguments: use min_ncount_rna, min_nfeature_rna, th_percent.mt
#args$type_qc <- "arguments"

# mad_and_th_percent.mt: use n_mad_log_counts, n_mad_log_features, th_percent.mt
# useful to include more cells with th_percent.mt=25 (e.g. GSE176078)
# --type_qc mad_arguments --n_mad_log1p_mito -1 --th_percent.mt 25

# mad_arguments: use n_mad_log_counts, n_mad_log_features, n_mad_log1p_mito, min_ncount_rna, min_nfeature_rna, th_percent.mt
# useful to include high quality cells (e.g. GSE161529)
# --type_qc mad_arguments --n_mad_log_counts 2 --n_mad_log_features 2 --n_mad_log1p_mito 2 --min_ncount_rna -1 --min_nfeature_rna -1 --th_percent.mt 25



cat(sprintf("\n\ttype_qc=%s\n", args$type_qc))
switch(args$type_qc,

	"arguments"={

		cat(sprintf("\tmin_ncount_rna=%d\n", args$min_ncount_rna))
		cat(sprintf("\tmin_nfeature_rna=%d\n", args$min_nfeature_rna))
		cat(sprintf("\tth_percent.mt=%d\n", args$th_percent.mt))


		f <- (rna$nCount_RNA < args$min_ncount_rna)
		TooLow_nCount_RNA <- length(which(f))
		cat(sprintf("\t\t# of cells with nCount_RNA < %g: %d\n", args$min_ncount_rna, TooLow_nCount_RNA))

		f <- (rna$nFeature_RNA < args$min_nfeature_rna)
		TooLow_nFeature_RNA <- length(which(f))
		cat(sprintf("\t\t# of cells with nFeature_RNA < %g: %d\n", args$min_nfeature_rna, TooLow_nFeature_RNA))

		f <- (rna$percent.mt > args$th_percent.mt)
		TooHigh_percent.mt <- length(which(f))
		cat(sprintf("\t\t# of cells with percent.mt > %g: %d\n", args$th_percent.mt, TooHigh_percent.mt))

		rna <- subset(rna, subset = (nCount_RNA >= args$min_ncount_rna &
               		 nFeature_RNA >= args$min_nfeature_rna & 
              		 percent.mt <= args$th_percent.mt))

	},

	"mad_arguments"={

		# outliers are >2 MADs
		cat(sprintf("\tn_mad_log_counts=%d\n", args$n_mad_log_counts))
		cat(sprintf("\tn_mad_log_features=%d\n", args$n_mad_log_features))
		cat(sprintf("\tn_mad_log1p_mito=%d\n", args$n_mad_log1p_mito))

		# https://github.com/davismcc/archive-scater/blob/master/R/qc.R
		f_nCount_RNA_outlier_mad <- isOutlier(log(rna$nCount_RNA), nmads = args$n_mad_log_counts, type = "lower", log = FALSE)
		f_nFeature_RNA_outlier_mad <- isOutlier(log(rna$nFeature_RNA), nmads = args$n_mad_log_features, type = "lower", log = FALSE)
		f_percent_mt_outlier_mad <- isOutlier(log1p(rna$percent.mt), nmads = args$n_mad_log1p_mito, type = "higher", log = FALSE)


		# ncount
		TooLow_nCount_RNA <- NA
		f_cells_exclude_ncount <- rep(FALSE, n_cells)
		if (args$n_mad_log_counts > 0) {
			f_cells_exclude_ncount <- f_nCount_RNA_outlier_mad
			TooLow_nCount_RNA <- length(which(f_cells_exclude_ncount))
			cat(sprintf("\t\t# of outliers log(nCount_RNA) < %d*MAD: %d\n", args$n_mad_log_counts, TooLow_nCount_RNA))
			list_out <- get_med_mad(log(rna$nCount_RNA), nmads = args$n_mad_log_counts, type = "lower", log = FALSE, n_log = 1)
			if (any(f_cells_exclude_ncount)) {
				cat(sprintf("\t\t\toutlier max(nCount_RNA)=%d\n", max(rna$nCount_RNA[f_cells_exclude_ncount])))
			}
		}

		if (args$min_ncount_rna > 0) {
			f <- (rna$nCount_RNA < args$min_ncount_rna)
			TooLow_nCount_RNA_ <- length(which(f))
			cat(sprintf("\t\t# of cells with nCount_RNA < %g: %d\n", args$min_ncount_rna, TooLow_nCount_RNA_))
			if (TooLow_nCount_RNA_ > 0) {
				f_cells_exclude_ncount[f] <- TRUE
			}

			TooLow_nCount_RNA <- length(which(f_cells_exclude_ncount))
			cat(sprintf("\t\t# of cells with outliers or nCount_RNA < %g: %d\n", args$min_ncount_rna, TooLow_nCount_RNA))
			if (any(f_cells_exclude_ncount)) {
				cat(sprintf("\t\t\tmax(nCount_RNA in removed cells)=%.1f\n", max(rna$nCount_RNA[f_cells_exclude_ncount])))
			}
		}

		# nfeature
		TooLow_nFeature_RNA <- NA
		f_cells_exclude_nfeature <- rep(FALSE, n_cells)
		if (args$n_mad_log_features > 0) {
			f_cells_exclude_nfeature <- f_nFeature_RNA_outlier_mad
			TooLow_nFeature_RNA <- length(which(f_cells_exclude_nfeature))
			cat(sprintf("\t\t# of outliers log(nFeature_RNA) < %d*MAD: %d\n", args$n_mad_log_features, TooLow_nFeature_RNA))
			list_out <- get_med_mad(log(rna$nFeature_RNA), nmads = args$n_mad_log_features, type = "lower", log = FALSE, n_log = 1)
			if (any(f_cells_exclude_nfeature)) {
				cat(sprintf("\t\t\toutlier max(nFeature_RNA)=%d\n", max(rna$nFeature_RNA[f_cells_exclude_nfeature])))
			}
		}

		if (args$min_nfeature_rna > 0) {
			f <- (rna$nFeature_RNA < args$min_nfeature_rna)
			TooLow_nFeature_RNA_ <- length(which(f))
			cat(sprintf("\t\t# of cells with nFeature_RNA < %g: %d\n", args$min_nfeature_rna, TooLow_nFeature_RNA_))
			if (TooLow_nFeature_RNA_ > 0) {
				f_cells_exclude_nfeature[f] <- TRUE
			}

			TooLow_nFeature_RNA <- length(which(f_cells_exclude_nfeature))
			cat(sprintf("\t\t# of cells with outliers or nFeature_RNA < %g: %d\n", args$min_nfeature_rna, TooLow_nFeature_RNA))
			if (any(f_cells_exclude_nfeature)) {
				cat(sprintf("\t\t\tmax(nFeature_RNA in removed cells)=%.1f\n", max(rna$nFeature_RNA[f_cells_exclude_nfeature])))
			}
		}

		# percent_mt
		TooHigh_percent.mt <- NA	
		f_cells_exclude_percent.mt <- rep(FALSE, n_cells)
		if (args$n_mad_log1p_mito > 0) {
			f_cells_exclude_percent.mt <- f_percent_mt_outlier_mad
			TooHigh_percent.mt <- length(which(f_cells_exclude_percent.mt))
			cat(sprintf("\t\t# of outliers log1p(percent.mt) > %d*MAD: %d\n", args$n_mad_log1p_mito, TooHigh_percent.mt))
			list_out <- get_med_mad(log1p(rna$percent.mt), nmads = args$n_mad_log1p_mito, type = "higher", log = FALSE, n_log = 1)
			if (any(f_cells_exclude_percent.mt)) {
				cat(sprintf("\t\t\toutlier min(percent.mt)=%.1f\n", min(rna$percent.mt[f_cells_exclude_percent.mt])))
			}
		}

		if (args$th_percent.mt > 0) {
			f <- (rna$percent.mt > args$th_percent.mt)
			TooHigh_percent.mt_ <- length(which(f))
			cat(sprintf("\t\t# of cells with percent.mt > %g: %d\n", args$th_percent.mt, TooHigh_percent.mt_))
			if (TooHigh_percent.mt_ > 0) {
				f_cells_exclude_percent.mt[f] <- TRUE
			}

			TooHigh_percent.mt <- length(which(f_cells_exclude_percent.mt))
			cat(sprintf("\t\t# of cells with outliers or percent.mt > %g: %d\n", args$th_percent.mt, TooHigh_percent.mt))
			if (any(f_cells_exclude_percent.mt)) {
				cat(sprintf("\t\t\tmin(percent.mt in removed cells)=%.1f\n", min(rna$percent.mt[f_cells_exclude_percent.mt])))
			}
		}

		rna$f_cells_exclude_ncount <- f_cells_exclude_ncount
		rna$f_cells_exclude_nfeature <- f_cells_exclude_nfeature
		rna$f_cells_exclude_percent.mt <- f_cells_exclude_percent.mt

		rna <- subset(rna, subset = (f_cells_exclude_ncount == FALSE) &
               		 (f_cells_exclude_nfeature == FALSE) & 
              		 (f_cells_exclude_percent.mt == FALSE))
	},

	{}
) # switch


PostQCNumCells <- ncol(rna)
cat(sprintf("\t# of cells after post QC: %d\n", PostQCNumCells))




 



if (nchar(args$method_sctransform_vst) == 0) {

  # normalization
  cat(sprintf("\tNormalizeData\n"))
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

  # feature selection
  cat(sprintf("\tFindVariableFeatures\n"))
  rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

  # cell cycle after normalization
  rna <- identify_cell_types(rna, "cell_cycle", args)

  # scaling
  cat(sprintf("\tScaleData var.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
  #all.genes <- rownames(rna)
  #rna <- ScaleData(rna, features = all.genes)
  # regressing out percent.mt is slow when features = all.genes
  #rna <- ScaleData(rna, features = all.genes, vars.to.regress="percent.mt", verbose = FALSE)
  rna <- ScaleData(rna, vars.to.regress = args$vars.to.regress, verbose = FALSE)

} else {

  # normalization
  cat(sprintf("\tNormalizeData\n"))
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

  # cell cycle after normalization
  rna <- identify_cell_types(rna, "cell_cycle", args)

  # https://satijalab.org/seurat/reference/sctransform
  # https://satijalab.org/seurat/articles/sctransform_vignette.html

  # Apply sctransform normalization
  # * Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
  # * Transformed data will be available in the SCT assay, which is set as the default after running sctransform
  # * During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage

  # The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure. It can be invoked by specifying method="glmGamPoi".
  cat(sprintf("\tSCTransform\n"))
  rna <- SCTransform(rna, method = args$method_sctransform_vst, vars.to.regress = args$vars.to.regress, verbose = FALSE)

} # if



# RunPCA
cat(sprintf("\tRunPCA\n"))
# https://satijalab.org/seurat/reference/runpca
rna <- RunPCA(rna,
  assay = NULL,  # Name of Assay PCA is being run on, when=NULL, DefaultAssay(rna)
  features = NULL, # default=NULL, Features to compute PCA on. If features=NULL, PCA will be run using the variable features for the Assay. Note that the features must be present in the scaled data. Any requested features that are not scaled or have 0 variance will be dropped, and the PCA will be run using the remaining features.
  npcs = 50, # default=50
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = FALSE, # Print the top genes associated with high/low loadings for the PCs
  ndims.print = 1:5, # PCs to print genes for
  nfeatures.print = 30, # Number of genes to print for each PC
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42 # default=42
) # RunPCA


log_obj(rna, tab=2)




# Score cells for cell type specific gene signatures
cat(sprintf("------------------------------------\n"))
cat(sprintf("\tScore cells for cell type specific gene signatures before doublet removal\n"))

if (grepl("singler", args$method_to_identify_reference_clusters)) {
	# immune_singler, immune_plus_endo_singler, normal-like_singler
	rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)
} else {
	#rna <- identify_cell_types(rna, "panglaodb", args)
	rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)
}







cat(sprintf("------------------------------------\n"))
cat(sprintf("\tFind clusters before doublet removal\n"))

rna <- find_clusters_before_doublet_removal(rna, args)
#rna <- update_rna_pc1_umap_cnv_before_doublet_removal(rna, args)








###########################################################
# Part 2: Doublet detection and removal
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 2: Doublet detection and removal\n"))
###########################################################


#' # Doublet detection: 1) DoubletDecon 2) DoubletFinder 3) Take intersection of calls

### Doublet Decon
# Add if else statement to regress out nCount RNA if needed before running DoubletDecon
idents.length <- length(levels(Idents(rna)))
for (i in levels(Idents(rna))){
  i <- as.numeric(i)
  levels(Idents(rna))[i] <- (i - 1)
}




decon.doublets <- c()
if (grepl("DoubletDecon", args$method_to_remove_doublets)) {

  # DoubletDecon is senstitive to slot=data
  cat(sprintf("------------------------------------\n"))
  cat(sprintf("\tDoubletDecon\n"))
  cat(sprintf("\tDoubletDecon\n"), file=stderr())

  cat(sprintf("\t\tImproved_Seurat_Pre_Process\n"))
  # https://github.com/EDePasquale/DoubletDecon/blob/master/R/Improved_Seurat_Pre_Process.R
  newFiles <- Improved_Seurat_Pre_Process(rna, num_genes=50, write_files=F, data_type="counts")

  vec_rhop <- unique(c(args$doubletdecon.rhop, seq(1.0, 0.5, by=-0.1)))

  for (rhop in vec_rhop) {

    cat(sprintf("\t\trhop=%g\n", rhop))
    results <- tryCatch({

	  out <- capture.output(
	  out <- capture.output(
	    results <- Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile,
		groupsFile=newFiles$newGroupsFile,
		filename=sample_id,
		#location="./DoubletDecon",
		location=sprintf("%s/%s_doubletdecon_", dir_log, sample_id),
		fullDataFile=NULL,
		removeCC=FALSE,
		species=args$doubletdecon.species, # Species as scientific species name, KEGG ID, three letter species abbreviation, or NCBI ID. Default is "mmu". https://bioinformatics.stackexchange.com/questions/4937/are-these-standard-species-abbreviations-and-how-to-look-up-others
		rhop=rhop, # default=1.1
		write=F,
		PMF=TRUE,
		useFull=FALSE,
		heatmap=FALSE,
		centroids=TRUE,
		num_doubs=100,
		only50=FALSE,
		min_uniq=4,
		nCores=args$n_cores)
	  ,append = FALSE, type = "output", split = FALSE)
	  ,append = FALSE, type = "message", split = FALSE)

	  results

	}, error = function(e) {
		cat(sprintf("\t\t#error-handler-code"))
		# slurm_er+bc_GSM4909299_ER-MH0114-T3.err
		# Unable to perform mcl function for blacklist clustering, please try a different rhop.
		if (!grepl("please try a different rhop|no locations are finite", e)) {
    			cat(sprintf("\t\t%s",e))
		}
		return(NULL)
	}, finally = {
    }) # tryCatch
  
    if (!is.null(results)) {
	break
    }
  } # for


  if (is.null(results)) {
  	cat(sprintf("\t\tdecon.doublets failure\n"))
	decon.doublets <- c()
  } else {
	decon.doublets <- rownames(results$Final_doublets_groups)
	decon.doublets <- gsub("\\.","-",decon.doublets)
  } # if

  n_decon.doublets <- length(decon.doublets)
  cat(sprintf("\t\t# of decon.doublets: %d\n", n_decon.doublets))


} # if












### DoubletFinder

DF.doublets <- c()
if (grepl("DoubletFinder", args$method_to_remove_doublets)) {


  cat(sprintf("------------------------------------\n"))
  cat(sprintf("\tDoubletFinder\n"))
  cat(sprintf("\tDoubletFinder\n"), file=stderr())


  # Add modfieid Parameter sweep function to regress out nCOunt RNA if needed
  # Add if else statement to regress out nCount RNA if needed before running DoubletFinder

  cat(sprintf("\t\tparamSweep_v3()\n"))
  cat(sprintf("\t\tparamSweep_v3()\n"), file=stderr())
  # Estimate Doublets
  # pK Identification (no ground-truth)
  # https://github.com/chris-mcginnis-ucsf/DoubletFinder/blob/master/R/paramSweep_v3.R
  sweep.res.list <- tryCatch({

		out <- capture.output(
		out <- capture.output(
		  sweep.res.list <- paramSweep_v3(rna, PCs = dimsToUse, sct = FALSE, num.cores=args$n_cores)
	  	,append = FALSE, type = "output", split = FALSE)
	  	,append = FALSE, type = "message", split = FALSE)

		sweep.res.list

	}, error = function(e) {
  		cat(sprintf("\t\tDF.doublets failure\n"))
		cat(sprintf("\t\t\t#error-handler-code"))
		# slurm_er+bc_GSM4909299_ER-MH0114-T3.err
		# Error in names(sweep.res.list) <- name.vec : 'names' attribute [186] must be the same length as the vector [157]
		# Calls: paramSweep_v3
		# In addition: Warning message: In mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3,  : scheduled core 1 did not deliver a result, all values of the job will be affected
    		cat(sprintf("\t\t\t%s",e))
		return(NULL)
	}, finally = {
  }) # tryCatch


  if (is.null(sweep.res.list)) {

	cat(sprintf("------------------------------------\n"))
	cat(sprintf("\tUse doublets called by DoubletDecon\n"))

	doublets <- decon.doublets

	pK.1 <- NA

  } else {

	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

	bcmvn<- find.pK(sweep.stats)

	pK.1 <- as.numeric(unfactor(dplyr::arrange(bcmvn, desc(BCmetric))$pK[1]))

	nExp_poi <- round(args$doublet.rate * length(colnames(counts)))  # Assuming 4.6% doublet formation rate - tailor for your dataset

	cat(sprintf("\t\tdoubletFinder_v3()\n"))
	cat(sprintf("\t\tdoubletFinder_v3()\n"), file=stderr())

	# Run doubletFinder with pk.1
	rna <- doubletFinder_v3(rna, PCs = dimsToUse, pN = 0.25, pK = pK.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

	doublet.column <- paste0("DF.classifications_0.25_", pK.1, "_", nExp_poi)
	doublet.calls <- rna[[doublet.column]]
	colnames(doublet.calls) <- "Call"

	rna.dub <- dplyr::filter(doublet.calls, Call == "Doublet")
	rna.singlet <- dplyr::filter(doublet.calls, Call == "Singlet")

	DF.doublets <- rownames(rna.dub)

  } # if

  n_DF.doublets <- length(DF.doublets)
  cat(sprintf("\t\t# of DF.doublets: %d\n", n_DF.doublets))

} # if










### Intersect doublet calls

cat(sprintf("------------------------------------\n"))

n_decon.doublets <- length(decon.doublets)
n_DF.doublets <- length(DF.doublets)

if ((n_decon.doublets > 0) && (n_DF.doublets > 0)) {
	cat(sprintf("\tIntersect doublet calls\n"))
	doublets <- intersect(decon.doublets, DF.doublets)
} else if ((n_decon.doublets > 0) && (n_DF.doublets == 0)) {
	cat(sprintf("\tUse doublets called by DoubletDecon\n"))
	doublets <- decon.doublets
} else if ((n_decon.doublets == 0) && (n_DF.doublets > 0)) {
	cat(sprintf("\tUse doublets called by DoubletFinder.\n"))
	doublets <- DF.doublets
} else {
	doublets <- c()
}

if (length(doublets) > 0) {
	rna@meta.data$Doublet.Call <- ifelse(rownames(rna@meta.data) %in% doublets, "TRUE", "FALSE")
}
cat(sprintf("\t# of doublets: %d\n", length(doublets)))




if (args$f_save_doublet_png) {
  #FeatureScatter(rna,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "Doublet.Call")+ggsave("Doublet_calls_FeatureScatter.png")
  gg <- FeatureScatter(rna,feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "Doublet.Call")
  ggsave(sprintf("%s/%s_doublet_calls_featurescatter.pdf", dir_pdf, sample_id), plot=gg)
  #ggsave(sprintf("%s/%s_doublet_calls_featurescatter.png", dir_png, sample_id), plot=gg)

  #DimPlot(rna,group.by = "Doublet.Call")+ggsave("Doublet_calls_UMAP.png")
  gg <- DimPlot(rna,group.by = "Doublet.Call")
  ggsave(sprintf("%s/%s_doublet_calls_umap.pdf", dir_pdf, sample_id), plot=gg)
  #ggsave(sprintf("%s/%s_doublet_calls_umap.png", dir_png, sample_id), plot=gg)
} # if


if (args$f_save_doublet_rds) {

  path_rds <- sprintf("%s/%s_intersection_doubletdecon_doubletfinder_doublets.rds", dir_rds, sample_id)
  cat(sprintf("\tsaveRDS(rna, '%s')\n", path_rds))
  saveRDS(doublets, path_rds)

  # save doublet information
  save(decon.doublets, DF.doublets, doublets, file=sprintf("%s/%s_doubletdecon_doubletfinder_doublets.rdata", dir_rdata, sample_id))

} # if


# cells after filtering
cells <- colnames(rna)
































###########################################################
# Part 3: scRNA-seq processing after doublet removal
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 3: scRNA-seq processing after doublet removal\n"))
###########################################################

# Load the RNA dataset

if (grepl(".txt.gz", args$sample_id)) {

	# read gene x cell count matrix
	filename <- sprintf("%s/%s", args$dir_count, args$sample_id)
	df <- read.table(filename, sep="\t", header=T, row.names=1, check.names=FALSE)
	counts.init <- as(as.matrix(df), "dgCMatrix")

	args$sample_id <- gsub("_Expression_Matrix.txt.gz", "", args$sample_id)
	sample_id <- sprintf("%s_%s", args$cancer_type, args$sample_id)

} else {

	#counts.init <- Read10X_h5(filename= "./filtered_feature_bc_matrix.h5")
	#counts.init <- Read10X_h5(filename=sprintf("../count/%s/outs/filtered_feature_bc_matrix.h5", dir_sample))
	file_name_h5 <- sprintf("%s/%s/outs/filtered_feature_bc_matrix.h5", args$dir_count, dir_sample)
	if (file.exists(file_name_h5)) {
		counts.init <- Read10X_h5(file_name_h5)
	} else {
		dir_sample1 <- sprintf("%s/%s", args$dir_count, dir_sample)
		counts.init <- Read10X(dir_sample1, gene.column = args$gene.column, cell.column = args$cell.column)
	}

} # if





# initialize the Seurat object with the raw (non-normalized data).
cat(sprintf("------------------------------------\n"))
cat(sprintf("\tCreateSeuratObject\n"))

if (!is.null(names(counts.init))) {
	if ("Gene Expression" %in% names(counts.init)) {
		# for multiome kit
		counts.init <- counts.init[["Gene Expression"]]
	}
}

rna <- CreateSeuratObject(
	counts = counts.init, # unnormalize data such as raw counts or TPMs
	project = sample_id,
	assay = "RNA",
	min.cells = 3, # genes not present in at least 3 cells are removed
	min.features = 0,
	names.field = 1,
	names.delim = "_",
	meta.data = NULL)

log_obj(rna, tab=2)




# remove doublets
cat(sprintf("\tremove doublets\n"))
rna <- rna[,cells]
rna <- rna[,!(colnames(rna) %in% doublets)]


n_cells_after_doublet_removal <- length(colnames(rna))
cat(sprintf("\t# of cells after doublet removal: %d\n", n_cells_after_doublet_removal))

if (n_cells_after_doublet_removal > args$max_n_cells) {

	col_cell_types <- get_column_name_for_cell_types(rna, args)
        idx_cells <- select_cells_with_high_counts(rna, col_cell_types, max_n_cells = args$max_n_cells, min_n_cells_per_cell_type = 20, assay = NULL, slot = "counts")
	rna <- rna[, idx_cells]

	n_cells <- length(colnames(rna))
	cat(sprintf("\t# of cells after applying max_n_cells: %d\n", n_cells))

} # if





# store mitochondrial percentage in object meta data
rna <- identify_cell_types(rna, "percent.mt", args)




if (nchar(args$method_sctransform_vst) == 0) {

  # normalization
  cat(sprintf("\tNormalizeData\n"))
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

  # feature selection
  cat(sprintf("\tFindVariableFeatures\n"))
  rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)

  # cell cycle after normalization
  rna <- identify_cell_types(rna, "cell_cycle", args)

  # scaling
  cat(sprintf("\tScaleData var.to.regress=%s\n", paste(args$vars.to.regress, collapse=", ")))
  #all.genes <- rownames(rna)
  #rna <- ScaleData(rna, features = all.genes)
  # regressing out percent.mt is slow when features = all.genes
  #rna <- ScaleData(rna, features = all.genes, vars.to.regress="percent.mt", verbose=FALSE)
  rna <- ScaleData(rna, vars.to.regress = args$vars.to.regress, verbose = FALSE)

} else {

  # normalization
  cat(sprintf("\tNormalizeData\n"))
  rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000)

  # cell cycle after normalization
  rna <- identify_cell_types(rna, "cell_cycle", args)

  # https://satijalab.org/seurat/articles/sctransform_vignette.html
  cat(sprintf("\tSCTransform\n"))
  rna <- SCTransform(rna, method = args$method_sctransform_vst, vars.to.regress =  args$vars.to.regress, verbose = FALSE)

} # if


# RunPCA
# https://satijalab.org/seurat/reference/runpca
cat(sprintf("\tRunPCA\n"))
rna <- RunPCA(rna,
  assay = NULL, # Name of Assay PCA is being run on, when=NULL, DefaultAssay(rna)
  features = NULL, # default=NULL, Features to compute PCA on. If features=NULL, PCA will be run using the variable features for the Assay. Note that the features must be present in the scaled data. Any requested features that are not scaled or have 0 variance will be dropped, and the PCA will be run using the remaining features.
  npcs = 50, # default=50
  rev.pca = FALSE,
  weight.by.var = TRUE,
  verbose = FALSE, # Print the top genes associated with high/low loadings for the PCs
  ndims.print = 1:5, # PCs to print genes for
  nfeatures.print = 30, # Number of genes to print for each PC
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 42 # default=42
) # RunPCA


log_obj(rna, tab=2)







# Score cells for cell type specific gene signatures
cat(sprintf("------------------------------------\n"))
cat(sprintf("\tScore cells for cell type specific gene signatures after doublet removal\n"))
if (grepl("singler", args$method_to_identify_reference_clusters)) {
	# immune_singler, immune_plus_endo_singler, normal-like_singler
	rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)
} else {
	rna <- identify_cell_types(rna, "panglaodb", args)
	rna <- identify_cell_types(rna, "cell_type_specific_gene_signatures", args)
}





cat(sprintf("------------------------------------\n"))
cat(sprintf("\tFind clusters after doublet removal\n"))

rna <- find_clusters_after_doublet_removal(rna, args)
#rna <- update_rna_pc1_umap_cnv_after_doublet_removal(rna, args)





# set idents from a value in object metadata
#Idents(rna)<- "RNA_snn_res.0.7"
str_column_of_meta_data_cluster <- sprintf("RNA_snn_res.%g", args$seurat_resolution)
Idents(rna) <- str_column_of_meta_data_cluster


# perform DEGs analysis with cell type annotated clusters 
Wilcox.markers <- find_all_markers(rna, args, assay=NULL, col_cluster_types=str_column_of_meta_data_cluster, sample_id=sample_id, f_save_rds=TRUE, dir_rds=dir_rds)









###########################################################
# Part 4: SingleR cell typing 
cat(sprintf("\n\n------------------------------------\n"))
cat(sprintf("Part 4: SingleR cell typing\n"))
cat(sprintf("------------------------------------\n"))
###########################################################

# SingleR labeling of celltypes
# SingleR annotation

# Set reference paths:
# Slyper et al. GSM4186985 (sample ID: HTAPP-624-SMP-3212)
# Processed using default Seurat parameters 


if (args$method_to_identify_cell_types != "singler_htapp_toolbox") {
  # 1) Slyper et al. Nat. Medicine 2020 scRNA-seq ovarian tumor 
  #ref.data.cancer_type <- readRDS("./HTAPP_ovarian_reference.rds")
  rna <- identify_cell_types(rna, "singler_htapp_toolbox", args)
}

if (args$method_to_identify_cell_types != "singler_hpca_cellidx") {
  # 2) Human Primary Cell Atlas Data (microarray)
  #ref.data.HPCA <- readRDS("./HPCA_celldex.rds")
  rna <- identify_cell_types(rna, "singler_hpca_cellidx", args)
}

if (args$method_to_identify_cell_types != "singler_blueprint_encode") {
  # 3) BluePrint Encode (bulk RNA-seq)
  #ref.data.BED <- readRDS("./BluePrintEncode_celldex.rds")
  rna <- identify_cell_types(rna, "singler_blueprint_encode", args)
}


# update cell types
rna <- update_cell_types_after_infercnv(rna, args)










###########################################################


# get_diet_seurat_obj
rna <- get_diet_seurat_obj(rna, args, diet_level=1)





cat(sprintf("------------------------------------\n"))
cat(sprintf("\tprint meta data\n\n"))
# extract meta data
md <- rna@meta.data %>% as.data.table
print(head(md))




cat(sprintf("------------------------------------\n"))
cat(sprintf("\tcount the number of cells per unique combinations of cell types and %s\n\n", str_column_of_meta_data_cluster))

col_cell_types <- get_column_name_for_cell_types(rna, args)
dt <- summarize_meta_data(rna, col_cell_types, str_column_of_meta_data_cluster)




cat(sprintf("------------------------------------\n"))
cat(sprintf("\tcount epithelial cell numbers with %s\n", col_cell_types))
idx <- grep(pattern_epi, rna@meta.data[, col_cell_types])
table(rna@meta.data[idx, col_cell_types])


# gene expression for subtypes
cat(sprintf("\tinspect gene expression for subtypes\n"))
switch(args$method_to_update_cell_types,
        "normal_epithelial_cell_types"={
                # all normal samples
		# https://www.nature.com/articles/s41467-018-04334-1
		genes_inspect <- c("KRT14", "KRT18", "LTF")
	},
        "epithelial_cell_types"={
                # set tumor/normal epithelial cells with CNV.Pos
		genes_inspect <- c("KRT14", "KRT18", "LTF")
	},
        "cancer_epithelial_cell_types"={
                # all tumor samples
		# https://breast-cancer-research.biomedcentral.com/articles/10.1186/bcr2635
		genes_inspect <- c("KRT14", "ERBB2", "KRT18")
	},
        "cancer_normal_epithelial_cell_types"={
                # normal + tumor samples
		genes_inspect <- c("KRT14", "ERBB2", "KRT18", "LTF")
	},
	{ 
		genes_inspect <- c()
	}
) # switch

if (length(genes_inspect) > 0) {
	for (gene in genes_inspect) {
		df <- inspect_gene_expr_distributions(rna, gene, n_log=1)
	} # for
} # if












###########################################################
# save seurat object 

cat(sprintf("------------------------------------\n"))
cat(sprintf("\tsave seurat object\n"))

path_seurat_obj <- sprintf("%s/%s_sc-rna-seq_sample_seurat_obj.rds", dir_rds, sample_id)
cat(sprintf("\tsaveRDS(rna, '%s')\n", path_seurat_obj))
saveRDS(rna, path_seurat_obj)




















###########################################################
cat(sprintf("------------------------------------\n"))
cat(sprintf("\twrite.xlsx\n"))


output.meta <- data.frame(

			# starting cells
			StartingNumCells = length(colnames(counts.init)),

			# QC
			min_nCount_RNA = args$min_ncount_rna,
                        min_nFeature_RNA = args$min_nfeature_rna,
                        th_percent.mt = args$th_percent.mt,

                        nMADLogCounts = args$n_mad_log_counts,
                        nMADLogFeatures = args$n_mad_log_features,
                        nMADLog1pMito = args$n_mad_log1p_mito,

			TooLow_nCount_RNA = TooLow_nCount_RNA,
			TooLow_nFeature_RNA = TooLow_nFeature_RNA,
			TooHigh_percent.mt = TooHigh_percent.mt,

                        PostQCNumCells = PostQCNumCells,

			# doublet removal
                        ExpectedDoubletFraction = args$doublet.rate,
                        ObservedDoubletFraction = length(doublets)/length(colnames(counts.init)),
                        PostDoubletNumCells = length(colnames(rna)),
                        #DoubletFinderpK = pK.1,

			# number of clusters
                        NumClusters = length(levels(Idents(rna))),

			# stats
                        SumCounts = sum(rna$nCount_RNA),
                        MinNumCounts = min(rna$nCount_RNA),
                        MaxNumCounts = max(rna$nCount_RNA),
                        MedianNumbCounts = median(rna$nCount_RNA),

                        MinNumFeats = min(rna$nFeature_RNA),
                        MaxNumFeats = max(rna$nFeature_RNA),
                        MedianNumbFeats = median(rna$nFeature_RNA),

                        stringsAsFactors = FALSE

		) # data.frame


# remove rows
cols <- colnames(output.meta)
switch(args$type_qc,
        "arguments"={
                f <- !grepl("nMAD", cols)
		cols <- cols[f]
        },
        "mad_arguments"={
		f <- rep(TRUE, length(cols))
		if (args$min_ncount_rna < 0) {
                	f[cols == "min_nCount_RNA"] <- FALSE
		}
		if (args$min_nfeature_rna < 0) {
                	f[cols == "min_nFeature_RNA"] <- FALSE
		}
		if (args$th_percent.mt < 0) {
                	f[cols == "th_percent.mt"] <- FALSE
		}
		if (args$n_mad_log_counts < 0) {
                	f[cols == "nMADLogCounts"] <- FALSE
		}
		if (args$n_mad_log_features < 0) {
                	f[cols == "nMADLogFeatures"] <- FALSE
		}
		if (args$n_mad_log1p_mito < 0) {
                	f[cols == "nMADLog1pMito"] <- FALSE
		}
		cols <- cols[f]
	},
	{}
) # switch
output.meta <- output.meta[,cols]






output <- as.data.frame(t(output.meta))
colnames(output) <- sample_id
file_name_xlsx <- sprintf("%s/%s_sc-rna-seq_pipeline_summary.xlsx", dir_xlsx, sample_id)
openxlsx::write.xlsx(output, file_name_xlsx, row.names = T, col.names = TRUE)

print(output.meta)





if (args$f_cellphonedb) {
	# cellphonedb
	cat(sprintf("------------------------------------\n"))
	cat(sprintf("\tcellphonedb\n\n"))

	source("r/run_cellphonedb.R")
	out <- run_cellphonedb(rna, args, col_cell_types, sample_id, max_n_cells_for_cellphonedb = 2000)
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
#cat(sprintf("./make_sc-rna-seq_seurat_obj.R %s %s\n", args$cancer_type, args$sample_id))
script_name <- sub(".*=", "", commandArgs()[4])
args_ <- commandArgs(trailingOnly = TRUE)
cat(sprintf("%s %s\n", script_name, paste(args_, collapse = " ")))
cat(sprintf("completed successfully~\n"))
cat(sprintf("------------------------------------\n\n"))





