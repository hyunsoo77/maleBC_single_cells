# 
# utilities_for_sc-atac-seq.R
# date craeted: 2021, Oct.
# date last modified: 2021, Dec.
#
# content:
# load_normal_peaks()
# load_encode_peaks()
# determine_cluster.type_exclude_with_archr_peakset()
# subset_archrproject()
#
# reference:
# 




# load_normal_peaks
load_normal_peaks <- function(args) {

  switch(args$cancer_type_standard,

	"bc"={
		# breast cancer
		filename_rds <- "./reference/normal_cell_lines_breast/HMEC_H3K27ac_peaks.rds"
		cat(sprintf("\tread %s\n", filename_rds))
		gr_normal.peaks <- readRDS(filename_rds)
	},

	"oc"={
		# ovarian cancer
		filename_rds <- "./reference/normal_cell_lines/Fallopian_Tube_Cell_line_Peaks.rds"
		cat(sprintf("\tread %s\n", filename_rds))
		normal.peaks <- readRDS(filenames_rds)
		gr_normal.peaks <- GRanges(normal.peaks)
	},

	{
		cat(sprintf("H23K27ac narrow peaks are not available for %s normal cell line\n", cancer_type))
	}

  ) # switch

  gr_normal.peaks

} # load_normal_peaks








# load_encode_peaks
load_encode_peaks <- function(args) {

  # read Encode peaks
  filename_bed <- "./reference/genome_annotation/GRCh38-cCREs.bed.gz"
  cat(sprintf("\tread %s\n", filename_bed))
  encode.all <- read.delim(filename_bed, header=F)
  colnames(encode.all)[1:3] <- c("seqnames","start","end")
  gr_encode <- GRanges(encode.all)

  gr_encode

} # load_encode_peaks








# determine_cluster.type_exclude_with_archr_peakset
determine_cluster.type_exclude_with_archr_peakset <- function(proj.archr, args) {

  df_meta.data <- as.data.frame(getCellColData(proj.archr))
  cluster_types <- naturalsort::naturalsort(unique(df_meta.data$predictedGroup_ArchR))


  # filtering peakset
  cat(sprintf("\n\tselect peaks\n"))

  # number of samples for each group
  dt <- summarize_meta_data(proj.archr, "predictedGroup_ArchR", "Sample", tab=2)

  peakset <- proj.archr@peakSet

  cluster_types_in_peakset <- naturalsort::naturalsort(unique(names(peakset)))
  cat(sprintf("\t\tcluster types in peakset: %s\n", paste(cluster_types_in_peakset, collapse=", ")))
  cluster.type_exclude <- c()

  cluster.type_exclude <- setdiff(cluster_types, cluster_types_in_peakset)
  cat(sprintf("\t\tcluster types not in peakset: %s\n", paste(cluster.type_exclude, collapse=", ")))

  if (args$f_exclude_cluster_epi_unassigned) {

	cat(sprintf("\t\tcheck epi. unassigned\n"))
	for (cluster in naturalsort::naturalsort(unique(names(peakset)))) {
		if (grepl("Epi. Unassigned", cluster)) {
  			cluster.type_exclude <- union(cluster.type_exclude, cluster)
			idx <- which(names(peakset) == cluster)
			cat(sprintf("\t\t\t%s: %d\n", cluster, length(idx)))
		}
	} # for

  } # if


  if (args$f_exclude_cluster_highly_overlapped_with_normal_peaks) {
	# check overlap with normal peaks
	cat(sprintf("\t\tcheck overlap with normal peaks\n"))
	for (cluster in naturalsort::naturalsort(unique(names(peakset)))) {
		idx <- which(names(peakset) == cluster)
		gr_cluster <- peakset[idx]
		idx_ov <- which(gr_cluster %over% gr_normal.peaks)
		ratio_overlapped_with_normal_peaks <- length(idx_ov)/length(gr_cluster)
		str_tail <- ""
		if (ratio_overlapped_with_normal_peaks > args$max_ratio_overlapped_with_normal_peaks) {
			# male-bc: 9-Epithelial cells: 4789/6166=0.776679
  			cluster.type_exclude <- union(cluster.type_exclude, cluster)
			str_tail <- sprintf(" > %g", args$max_ratio_overlapped_with_normal_peaks)
		}
		cat(sprintf("\t\t\t%s: %d/%d=%g%s\n", cluster,
			length(idx_ov),
			length(gr_cluster),
			ratio_overlapped_with_normal_peaks, str_tail))
	} # for
  } # if


  if (args$f_exclude_cluster_low_nfrags) {
	# check average nFrags
	cat(sprintf("\t\tcheck mean of nFrags\n"))
	if (args$min_nfrags < 0) {
		# automatically determine min_nfrags
		#args$min_nfrags <- median(df_meta.data$nFrags) - 0.3*stats::mad(df_meta.data$nFrags)
		args$min_nfrags <- mean(df_meta.data$nFrags) - 0.5*stats::sd(df_meta.data$nFrags)
		#args$min_nfrags <- exp(mean(log(df_meta.data$nFrags)) - 0.4*stats::sd(log(df_meta.data$nFrags)))
	}

	for (cluster in naturalsort::naturalsort(unique(df_meta.data$predictedGroup_ArchR))) {
		idx <- which(df_meta.data$predictedGroup_ArchR == cluster)
		mean_nfrags <- mean(df_meta.data[idx, "nFrags"], na.rm=T)
		str_tail <- ""
		if (mean_nfrags < args$min_nfrags) {
			# male-bc: 2-Epithelial cells mean_nfrags=4872.894
			# male-bc: 10-Epithelial cells mean_nfrags=4099.144	
  			cluster.type_exclude <- union(cluster.type_exclude, cluster)
			str_tail <- sprintf(" < %g", args$min_nfrags)
		}
		cat(sprintf("\t\t\t%s: mean_nFrag=%g%s\n", cluster,
			mean_nfrags, str_tail))
	} # for
  } # if

  if (length(cluster.type_exclude) > 0) {
	# update peakset
	f.exclude <- names(peakset) %in% cluster.type_exclude
	peakset <- peakset[!f.exclude]
  } # if

  cluster_types_in_peakset <- naturalsort::naturalsort(unique(names(peakset)))
  cat(sprintf("\t\tcluster types in peakset after peak selection: %s\n", paste(cluster_types_in_peakset, collapse=", ")))

  cluster.type_exclude <- naturalsort::naturalsort(cluster.type_exclude)

  cluster.type_exclude


} # determine_cluster.type_exclude_with_archr_peakset











# subset_archrproject
subset_archrproject <- function(proj.archr, args, cluster.type_exclude) {


  cat(sprintf("\n\tsubset_archrproject\n"))
  df_meta.data <- as.data.frame(getCellColData(proj.archr))
  
  cat(sprintf("\t\t# of cells: %d\n", nrow(df_meta.data)))
  cat(sprintf("\t\tcluster.type_exclude: %s\n", paste(cluster.type_exclude, collapse=", ")))
  f <- (!df_meta.data$predictedGroup_ArchR %in% cluster.type_exclude)
  cell_names <- rownames(df_meta.data)[f]

  dir_output_subset <- sprintf("%s/archr_output", args$dir_output)

  if (args$f_subset_archrproject_force_to_update) {

	cat(sprintf("\t\trm -rf %s\n", dir_output_subset))
	unlink(dir_output_subset, recursive = TRUE, force = FALSE)

	# https://www.archrproject.com/reference/subsetArchRProject.html
	cat(sprintf("\t\tsubsetArchRProject\n"))
	proj.archr <- subsetArchRProject(
		ArchRProj = proj.archr,
		cells = cell_names, # A vector of cells to subset ArchRProject by. Alternatively can provide a subset ArchRProject.
		outputDirectory = dir_output_subset,
		dropCells = TRUE, # A boolean indicating whether to drop cells that are not in ArchRProject from corresponding Arrow Files.
		logFile = createLogFile("subsetArchRProject", logDir=args$dir_log), # The path to a file to be used for logging ArchR output.
		threads = getArchRThreads(),
		force = args$f_subset_archrproject_force_to_update # If output directory exists overwrite.
  	) # subsetArchRProject

"Copying ArchRProject to new outputDirectory : /datastore/nextgenout5/share/labs/francolab/hyunsoo.kim/sc-atac-seq/male-bc/run-20220523/output_p2g_male-bc/archr_output
Copying Arrow Files...
Getting ImputeWeights
No imputeWeights found, returning NULL
Copying Other Files...
Copying Other Files (1 of 9): Annotations
Copying Other Files (2 of 9): Background-Peaks.proj.archr.rds
Copying Other Files (3 of 9): Background-Peaks.proj.atac.rds
Copying Other Files (4 of 9): Embeddings
Copying Other Files (5 of 9): LSI_ATAC
Copying Other Files (6 of 9): Peak2GeneLinks
Copying Other Files (7 of 9): PeakCalls
Copying Other Files (8 of 9): Plots
Copying Other Files (9 of 9): RNAIntegration
Saving ArchRProject...
Loading ArchRProject...
Successfully loaded ArchRProject!"

	#Error in subsetArchRProject(ArchRProj = proj.archr, cells = cell_names,  : outputDirectory exists! Please set force = TRUE to overwrite existing dir ectory!

  } else {

	filename_rds <- sprintf("%s/Save-ArchR-Project.rds", dir_output_subset)
	cat(sprintf("\t\treadRDS('%s')\n", filename_rds))
	proj.archr <- readRDS(filename_rds)

  } # if

  cat(sprintf("\t\t# of cells: %d\n\n", nCells(proj.archr)))


  proj.archr


} # subset_archrproject












