#
# identify_reference_clusters.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2021, Dec.
#
# content:
# identify_reference_clusters()
# update_rna_idents_with_reference_clusters()
#
# comment:
# called by make_sc-rna-seq_seurat_obj.R
#
#

suppressPackageStartupMessages(library(psych))





# identify_reference_clusters
# input:
#   args$method_to_identify_reference_clusters {["immune_plasma"], "nonepi", "nonepi_krt_epcam"}
#
# usage:
# immune.clusters <- identify_reference_clusters(rna, "matt")
identify_reference_clusters <- function(rna, args) {

  cat(sprintf("\t\t[identify_reference_clusters()]\n"))

  clusters <- levels(Idents(rna))

  switch(args$method_to_identify_reference_clusters,

	"immune_plasma"={

	  # identify reference clusters by relative enrichment of ESTIMATE immune signature
	  test <- VlnPlot(rna, features = "immune.2")
	  data <- psych::describeBy(test$data$immune.2, test$data$ident, mat = TRUE)
	  data.immune <- dplyr::filter(data, median > 0.1)

	  test <- VlnPlot(rna, features = "plasma.7")
	  data <- psych::describeBy(test$data$plasma.7, test$data$ident, mat = TRUE)
	  data.plasma <- dplyr::filter(data, median > 0.1)

	  immune.clusters <- intersect(data.immune$group1, clusters)
	  plasma.clusters <- intersect(data.plasma$group1, clusters)

	  reference.clusters <- unique(append(immune.clusters, plasma.clusters))
	},

	"nonepi"={

	  # refernce.clusters may contain immune, plasma cells
	  # but, reference.clusters cannot have normal epithelial cells
	  test <- VlnPlot(rna, features = "epithelial.5")
	  data <- psych::describeBy(test$data$epithelial.5, test$data$ident, mat = TRUE)
	  #data.epi <- dplyr::filter(data, median > 0.1)
	  data.nonepi <- dplyr::filter(data, median < 0)

	  reference.clusters <- intersect(data.nonepi$group1, clusters)

	},

	"nonepi_krt_epcam"={

		# basal marker: KRT5, KRT14, KRT17
		# luminal marker: KRT8/18, KRT19, EPCAM
		rna <- suppressWarnings(suppressMessages( AddModuleScore(rna, features=list(c("KRT14", "KRT5", "AKR1B1", "CAV1", "GYPC", "MYL9", "MYLK", "PPP1R14A"), c("KRT18", "KRT19", "CD24", "AZGP1", "C15orf48", "DEFB1", "PDZK1IP1", "RHOV")), nbin=2, ctrl=1, assay=NULL, name=c("epi_basal.", "epi_luminal."), search=TRUE) ))

		vln.df <- VlnPlot(rna, features = "epi_basal.1")
		data <- psych::describeBy(vln.df$data$epi_basal.1, vln.df$data$ident, mat = TRUE, quant=c(.25, .75))
		#data.epi_basal <- dplyr::filter(data, Q0.25 > 0)
		data.nonepi_basal <- dplyr::filter(data, Q0.25 <= 0)

		vln.df <- VlnPlot(rna, features = "epi_luminal.2")
		data <- psych::describeBy(vln.df$data$epi_luminal.2, vln.df$data$ident, mat = TRUE, quant=c(.25, .75))
		#data.epi_luminal <- dplyr::filter(data, Q0.25 > 0)
		data.nonepi_luminal <- dplyr::filter(data, Q0.25 <= 0)

		clusters.nonepi <- intersect(data.nonepi_basal$group1, data.nonepi_luminal$group1)
		reference.clusters <- intersect(clusters.nonepi, clusters)

	},	

	"nonepi_singler"={

	  # refernce.clusters may contain immune, plasma cells
	  # but, reference.clusters cannot have normal epithelial cells
	  vec_cell_types <- get_vec_cell_types(rna, args)
	  reference.clusters <- c()
	  for (cluster in clusters) {
		idx <- which(Idents(rna) == cluster)
		tb <- sort(table(vec_cell_types[idx]), decreasing=T)
		if (!grepl(pattern_epi, names(tb)[1])) {
	  		reference.clusters <- c(reference.clusters, cluster)
		}	
	  }
	},

	"immune_singler"={

	  # refernce.clusters may contain immune, plasma cells
	  # but, reference.clusters cannot have normal epithelial cells
	  vec_cell_types <- get_vec_cell_types(rna, args)
	  reference.clusters <- c()
	  for (cluster in clusters) {
		idx <- which(Idents(rna) == cluster)
		tb <- sort(table(vec_cell_types[idx]), decreasing=T)
		if (grepl(pattern_immune_cells, names(tb)[1])) {
	  		reference.clusters <- c(reference.clusters, cluster)
		}	
	  }

 	},

	"immune_plus_endo_singler"={

	  # refernce.clusters may contain immune, plasma cells
	  # but, reference.clusters cannot have normal epithelial cells
	  vec_cell_types <- get_vec_cell_types(rna, args)
	  reference.clusters <- c()
	  for (cluster in clusters) {
		idx <- which(Idents(rna) == cluster)
		tb <- sort(table(vec_cell_types[idx]), decreasing=T)
		if (grepl(pattern_immune_plus_endothelial_cells, names(tb)[1])) {
	  		reference.clusters <- c(reference.clusters, cluster)
		} else {
			# Epi, Fibro
		}
	  }

 	},

	"normal-like_singler"={

	  # assume normal-like contains most normal epithelial cells.
	  vec_cell_types <- get_vec_cell_types(rna, args)
	  reference.clusters <- c()
	  for (cluster in clusters) {
		idx <- which(Idents(rna) == cluster)
		tb <- sort(table(vec_cell_types[idx]), decreasing=T)
		if (grepl("Normal-like", names(tb)[1])) {
	  		reference.clusters <- c(reference.clusters, cluster)
		}	
	  }

	},

	{
		# default way to determine reference.clusters to bypass the reference
		# i.e. immune_singler
	  	# refernce.clusters may contain immune, plasma cells
		# but, reference.clusters cannot have normal epithelial cells
		vec_cell_types <- get_vec_cell_types(rna, args)
		reference.clusters <- c()
		for (cluster in clusters) {
			idx <- which(Idents(rna) == cluster)
			tb <- sort(table(vec_cell_types[idx]), decreasing=T)
			if (grepl(pattern_immune_plus_endothelial_cells, names(tb)[1])) {
		  		reference.clusters <- c(reference.clusters, cluster)
			}	
		}

	}
  ) # switch

  cat(sprintf("\t\treference.clusters: %d/%d, %s\n", length(reference.clusters), length(clusters), paste(reference.clusters, collapse=", ")))
  idx_ref <- which(Idents(rna) %in% reference.clusters)
  cat(sprintf("\t\tn_reference_cells: %d\n", length(idx_ref)))
  
  observation.clusters <- setdiff(clusters, reference.clusters)
  cat(sprintf("\t\tobservation.clusters: %d/%d, %s\n", length(observation.clusters), length(clusters), paste(observation.clusters, collapse=", ")))
  idx_obs <- which(Idents(rna) %in% observation.clusters)
  cat(sprintf("\t\tn_observation_cells: %d\n", length(idx_obs)))

  reference.clusters

} # identify_reference_clusters









# update_rna_idents_with_reference_clusters
#
# input:
#   colname_group: {predoublet.idents, postdoublet.idents}
# output:
#   0 --> immune.0
#
# usage:
# reference.clusters <- identify_reference_clusters(rna, args)
# rna <- update_rna_idents_with_reference_clusters(rna, reference.clusters, "predoublet.idents", args)
#
update_rna_idents_with_reference_clusters <- function(rna, reference.clusters, colname_group, args) {

    for (i in 1:length(reference.clusters)) {

      j <- which(levels(Idents(rna)) == reference.clusters[i])
      #levels(Idents(rna))[j] <- paste0("immune.", reference.clusters[i])
      levels(Idents(rna))[j] <- paste0(args$prefix_reference_clusters, reference.clusters[i])

    } # for

    # update meta.data
    rna@meta.data[, colname_group] <- Idents(rna)

    rna

} # update_rna_idents_with_reference_clusters




