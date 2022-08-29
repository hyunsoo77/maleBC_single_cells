#!/usr/bin/env Rscript
#
# write_tile_matrices.R
# author: H. Kim
# date created: 2022, May
# date last modified: 2022, May
#
# usage:
# ./r/write_tile_matrices.R --dir_count ../count_male-bc --dir_output output_male-bc male-bc
#
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


parser$add_argument("cancer_type", nargs=1, help="cancer type")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

dir_output <- args$dir_output
cancer_type <- args$cancer_type


suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(HDF5Array))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(stringr))



dir_archr_output <- sprintf("%s/archr_output", dir_output)
#dir_archr_output <- normalizePath(dir_archr_output)
args$dir_archr_output <- dir_archr_output
dir_geo <- sprintf("%s/geo", dir_output)
dir_geo_processed_data <- sprintf("%s/processed_data", dir_geo)
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

# create geo submission directory
#dir.create(dir_geo, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_geo_processed_data, showWarnings = FALSE, recursive = TRUE)






filenames_arrow <- list.filenames_arrow(path=sprintf("%s/ArrowFiles", dir_archr_output), dir_pattern = ".arrow", full.names=TRUE)
filenames_rds <- c()

# for loop
for (filename_arrow in filenames_arrow){

  cat(sprintf("%s\n", filename_arrow))

  # https://www.archrproject.com/reference/getMatrixFromArrow.html
  se_tile_matrix <- getMatrixFromArrow(ArrowFile = filename_arrow,
	useMatrix = "TileMatrix",
	useSeqnames = NULL,
	cellNames = NULL,
	ArchRProj = NULL,
	verbose = FALSE,
 	binarize = TRUE, # defualt=FALSE, A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
	logFile = createLogFile("getMatrixFromProject", logDir=dir_log)
   ) # getMatrixFromArrow

  sample <- gsub(".arrow", "", basename(filename_arrow))
  filename_rds <- sprintf("%s_sc-atac-seq_tile_matrix.rds", sample)
  filenames_rds <- c(filenames_rds, filename_rds)
  path_rds <- sprintf("%s/%s", dir_geo_processed_data, filename_rds)
  saveRDS(se_tile_matrix, file = filesname_rds)

} # for



# write xlsx
cat(sprintf("\n------------------------------------\n"))
cat(sprintf("\twrite xlsx\n"))
df_rds_files <- data.frame(files=filenames_rds)
file_name_xlsx <- sprintf("%s/processed_files_sc-atac-seq_tile_matrix.xlsx", dir_geo)
openxlsx::write.xlsx(df_rds_files, file_name_xlsx, row.names = TRUE, col.names = TRUE)





