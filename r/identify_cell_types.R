#
# identify_cell_types.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2021, Dec.
#
# content:
# identify_cell_types()
# update_cluster.ids()
# get_most_common_cell_type_for_each_cluster()
# expand_epi_cells()
# update_cell_types()
# identify_cell_types_with_multiple_methods()
#
# comment:
# called by make_sc-rna-seq_seurat_obj.R
#
# debug:
# library(Seurat); source("r/identify_cell_types.R"); args <- c(); args$cancer_type <- "brca-mut-bc"; rna <- identify_cell_types(rna, "normal_epithelial_cells", args)
#
#

suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(GSA))
suppressPackageStartupMessages(library(SummarizedExperiment))






# for subroutines
source("./r/utilities_for_sc_analyses.R")





# cell type specific gene signatures

#mast_cells_3genes <- c("TPSB2", "TPSAB1", "KIT") # Three genes are small number of genes, so it is possible not to detect mast cells when the number of feautres of Seurat obj is less than 21,000.
# KIT is also expressed in normal breast epithelium.  https://aacrjournals.org/clincancerres/article/10/1/178/165808/KIT-CD117-Positive-Breast-Cancers-Are-Infrequent

mast_cells_11genes <- c("TPSB2", "TPSAB1", "TPSD1", "TESPA1", "RGS13", "SLC18A2", "CPA3", "MS4A2", "HPGDS", "ADCYAP1", "HDC") # TPSB2, TPSAB1, TPSD1 and CPA3 encode the well-known mast cell proteases. MS4A2, HPGDS and HDC were also reported as genes encoding well-defined mast cell markers.  In the 69 asthma patients, the mast cell-specific 11-gene signature was significantly correlated with mast cell numbers in biopsies (Pearson, r = 0.42, P < .001, FDR = 0.0084).  https://onlinelibrary.wiley.com/doi/10.1111/cea.13732

# mast_cells_6genes <- c("TPSB2", "TPSAB1", "RGS13", "CPA3", "MS4A2", "HPGDS") from er+bc-53FBAL dataset.






# PanglaoDB
#
tsv=gzfile("reference/cell_type_signatures/PanglaoDB_markers_27_Mar_2020.tsv.gz")
panglaodb <- read.csv(tsv,header=T,sep = "\t")
panglaodb <- dplyr::filter(panglaodb,species == "Hs" | species == "Mm Hs")# Human subset
panglaodb <- dplyr::filter(panglaodb,organ == "Connective tissue" |
                             organ == "Epithelium" |
                             organ == "Immune system" |
                             organ == "Reproductive"|
                             organ == "Vasculature" |
                             organ == "Smooth muscle"
)
panglaodb <- split(as.character(panglaodb$official.gene.symbol), panglaodb$cell.type)



# ESTIMATE signatures
# Supplementary Table 1. Gene signatures for stromal and immune score estimation.  https://www.frontiersin.org/articles/10.3389/fonc.2019.01212/full#SM1
ESTIMATE.signatures <- "reference/cell_type_signatures/ESTIMATE_signatures.csv"



# breast epithelial cell markers
#load("reference/list_markers_for_breast_normal_epi_cells_hkim.rda")







# get_cancer_type_specific_info
#
# output:
#   list_out:
#     filename_rds_ref_htapp:
#     df_sig:
#     cols:
#     subtypes:
#     df_cna: copy number alteration
#
#   update global variables:
#     pattern_normal_cells
#     pattern_tumor_cells
#
# usage:
# list_cancer_type_specific_info <- get_cancer_type_specific_info(args)
get_cancer_type_specific_info <- function(args) {

  list_out <- list()

  type_normal_epi_type <- "unknown"
  filename_rds_ref_htapp <- NULL
  filename_rds_ref_major <- NULL
  filename_garnett_marker_file <- NULL
  filename_garnett_classifier_prefix <- NULL
  filename_rds_centroids <- NULL
  df_sig <- NULL
  cols <- NULL
  subtypes <- NULL
  df_cna <- NULL

  switch(args$cancer_type_standard,
	"bc"={
		# brca1-mut-bc, er+bc, her2+bc, male-bc, normal-breast, tr-bc (tamoxifen resistance breast cancer)

		#type_normal_epi_type <- "peroulab"
		type_normal_epi_type <- "swarbricklab"
		#type_normal_epi_type <- "kessenbrocklab"

 		# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7220853/
		# https://tumor-toolbox.broadinstitute.org/
		filename_rds_ref_htapp <- "reference/sc-rna-seq_references/GSM4186971_HTAPP-254-SMP-571_mbc_sc-rna-seq_reference_seurat_obj.rds"
		# A single-cell and spatially resolved atlas of human breast cancers  https://www.nature.com/articles/s41588-021-00911-1
		filename_rds_ref_major <- "reference/sc-rna-seq_references/brca-mini-atlas_sc-rna-seq_reference_seurat_obj.rds"
		#filename_garnett_marker_file <- "reference/garnett/breast_cancer_TME_markers_v4_XCell_revised_modified.txt"
		filename_garnett_marker_file <- "reference/garnett/brca-mini-atlas_garnett_marker_file_major_hkim.txt"
		filename_garnett_classifier_prefix <- "reference/garnett/brca-mini-atlas_garnett_classifier"
		# Supplementary Table 4.  scSubtype gene lists  Gene lists used to define the single-cell scSubtype molecular subtype classifier, one for each scSubtype (Basal_SC, Her2E_SC, LumA_SC and LumB_SC).
		if (!is.null(args$method_to_identify_subtypes)) {
			switch(args$method_to_identify_subtypes,
				"none"={
				},
				"brca_cell_atlas_scsubtyper"={
					df_sig <- read.csv("reference/scsubtyper/NatGen_Supplementary_table_S4.csv", header = T)
			 		cols <- c("Basal_SC", "Her2E_SC", "LumA_SC", "LumB_SC")
					subtypes <- c("Basal", "Her2E", "LumA", "LumB")
				},
				"scbcsubtype"={
					# swarbricklab_breast_cancer_epi
					filename_rds_centroids <- sprintf("reference/centroids_for_singler/scbcsubtype_%s.rds", args$type_centroids)
					if (file.exists(filename_rds_centroids)) {
						se_centroid <- readRDS(filename_rds_centroids)
						df_sig <- assay(se_centroid)	
						# cols <- c("Basal", "CLow", "Her2E", "LumA", "LumB", "Normal-like", "Normal")
						cols <- colnames(df_sig)
					}
				},
				"brca_pam50_tcga_centroids"={
					filename_rds_centroids <- "reference/centroids_for_singler/tcga_brca_molecular_subtype_centroids.rds"
					se_centroid <- readRDS(filename_rds_centroids)
					df_sig <- assay(se_centroid)	
					# cols <- c("Basal", "CLow", "Her2E", "LumA", "LumB", "Normal-like", "Normal")
					cols <- colnames(df_sig)
					subtypes <- cols
				},
				{
					if (grepl("pam50", args$method_to_identify_subtypes)) {
						# brca_pam50_cor, brca_pam50_genefu
						suppressPackageStartupMessages(library(genefu))
						data(pam50)
						df_sig <- as.data.frame(pam50$centroids)
						cols <- c("Basal", "Her2E", "LumA", "LumB", "Normal-like")
						subtypes <- cols
					} # if
	 			}
       			) # switch
     		} # if

     		df_cna <- read.table("reference/copy_number_alterations/tcga_brca_cna_hkim.tsv", sep="\t", header=T, row.names=1)
	},

	"oc"={

		# matt-oc
		# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7220853/
		# https://tumor-toolbox.broadinstitute.org/
		# Slyper et al. GSM4186985 (sample ID: HTAPP-624-SMP-3212)
		filename_rds_ref_htapp <- "reference/sc-rna-seq_references/GSM4186985_HTAPP-624-SMP-3212_oc_sc-rna-seq_reference_seurat_obj.rds"

	},

	"atll"={

		# adult T-cell leukemia/lymphoma (ATLL)
		pattern_normal_cells <<- "(^|-)NK cells|(^|-)Monocytes|(^|-)DC|(^|-)Macrophages|(^|-)Mast cells|(^|-)B-cells"
		pattern_tumor_cells <<- "(^|-)T-cell|(^|-)CD4\\+ T-cells|(^|-)CD8\\+ T-cells|(^|-)Treg-cells"

	},

	{

		stop(sprintf("no reference data defined for %s", args$cancer_type))

	}

  ) # switch

  list_out$type_normal_epi_type <- type_normal_epi_type
  list_out$filename_rds_ref_htapp <- filename_rds_ref_htapp
  list_out$filename_rds_ref_major <- filename_rds_ref_major
  list_out$filename_garnett_marker_file <- filename_garnett_marker_file
  list_out$filename_garnett_classifier_prefix <- filename_garnett_classifier_prefix
  list_out$filename_rds_centroids <- filename_rds_centroids
  list_out$df_sig <- df_sig
  list_out$cols <- cols
  list_out$subtypes <- subtypes
  list_out$df_cna <- df_cna

  list_out

} # get_cancer_type_specific_info














# identify_cell_types
#
# usage:
# identify_cell_types(rna, "singler_blueprint_encode")
identify_cell_types <- function(rna, method, args, n_log=1) {

  if (n_log > 0) {
    cat(sprintf("\tidentify_cell_types(rna, method=%s)\n", method))
  }

  #list_cancer_type_specific_info <- get_cancer_type_specific_info(args)
  list_cancer_type_specific_info <- args$list_cancer_type_specific_info

  immune.stromal <- read.csv(ESTIMATE.signatures, header = F)

  # from ESTIMATE.signatures
  stromal <- immune.stromal$V1[1:141]
  immune <- immune.stromal$V1[142:282]

  # from PanglaoDB_markers_27_Mar_2020.tsv.gz
  fibroblast <- panglaodb$Fibroblasts
  endothelial <- panglaodb$`Endothelial cells`
  epithelial <- panglaodb$`Epithelial cells`
  smooth <- panglaodb$`Smooth muscle cells`
  plasma <- panglaodb$`Plasma cells`
  b_cells <- panglaodb$`B cells`

  # other cell types in PanglaoDB
  mast_cells <- panglaodb$`Mast cells`
  macrophages <- panglaodb$Macrophages
  dendritic_cells <- panglaodb$`Dendritic cells`
  t_cells <- panglaodb$`T cells`
  nk_cells <- panglaodb$`NK cells`

  # store default assay for SingleR
  default_assay <- DefaultAssay(rna)

  if ((DefaultAssay(rna) == "integrated") && !is.null(args$type_integration_anchor_features)) {
	if (args$type_integration_anchor_features != "all") {
		# when type_integration_anchor_features="2000", the number of genes in assay="integreated" may not be enough for calling cell types, so use assay of RNA.
		DefaultAssay(rna) <- "RNA"
	}
  } # if

  rna.sce <- as.SingleCellExperiment(rna)

  switch(method,

	"percent.mt"={

		if (n_log > 0) {
			cat(sprintf("\t\tPercentageFeatureSet\n"))
		}
		# https://satijalab.org/seurat/reference/percentagefeatureset
		# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
		#rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
		rna <- PercentageFeatureSet(rna, pattern = "^MT-", col.name="percent.mt")

	},

	"cell_cycle"={

		# https://satijalab.org/seurat/archive/v3.0/cell_cycle_vignette.html
		# modified gene names: FAM64A --> PIMREG, HN1 --> JPT1
		exp.mat <- read.table(file = "./reference/cell_cycle/nestorawa_for_cell_cycle_expression_matrix_modified.txt", header = TRUE, as.is = TRUE, row.names = 1)

		# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can segregate this list into markers of G2/M phase and markers of S phase
		s.genes <- Seurat::cc.genes$s.genes
		g2m.genes <- Seurat::cc.genes$g2m.genes
		
		# modified gene names to standard names.
		idx <- match("MLF1IP", s.genes)
		if (length(idx) == 1) {
			s.genes[idx] <- "CENPU"
		}

		idx <- match(c("FAM64A", "HN1"), g2m.genes)
		if (length(idx) == 2) {
			g2m.genes[idx] <- c("PIMREG", "JPT1")
		}

		# https://satijalab.org/seurat/reference/cellcyclescoring
		# A Seurat object with the following columns added to object meta data: S.Score, G2M.Score, and Phase.
		if (n_log > 0) {
			cat(sprintf("\t\tCellCycleScoring\n"))
		}
		rna <- CellCycleScoring(rna, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
		# Insufficient data values to produce 24 bins? check if normalized before calling CellCycleScoring()  https://github.com/satijalab/seurat/issues/1181

	},


	"cell_type_specific_gene_signatures"={

		# https://satijalab.org/seurat/reference/addmodulescore
		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}

		features <- list(mast_cells_11genes)
		name <- c("mast_cells.")
		idx_mast_genes <- which(mast_cells_11genes %in% rownames(rna.sce))
		if (length(idx_mast_genes) < 2) {
			return(rna)
		}
		rna <- suppressWarnings( suppressMessages( AddModuleScore(rna,
			features = features,
			pool = NULL,
			nbin = 24, # Number of bins of aggregate expression levels for all analyzed features
			ctrl = 100, # Number of control features selected from the same bin per analyzed feature
			k = FALSE,
			assay = NULL,
			name = name,
			seed=1,
			search = TRUE # Search for symbol synonyms for features in features that don't match features in object? Searches the HGNC's gene names database; see UpdateSymbolList for more details
		) ) ) # AddModuleScore

		

  		# usage:
		# mast_cells.1 from other sources

	},

	"panglaodb"={
		# https://satijalab.org/seurat/reference/addmodulescore
		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}

		features <- list(stromal, immune, fibroblast, endothelial, epithelial, smooth, plasma, b_cells)
		name <- c("stromal.","immune.","fibroblast.","endothelial.", "epithelial.","smooth.","plasma.","b_cells.")
		rna <- suppressWarnings( suppressMessages( AddModuleScore(rna,
			features = features,
			pool = NULL,
			nbin = 24, # Number of bins of aggregate expression levels for all analyzed features
			ctrl = 100, # Number of control features selected from the same bin per analyzed feature
			k = FALSE,
			assay = NULL,
			name = name,
			seed=1,
			search = TRUE # Search for symbol synonyms for features in features that don't match features in object? Searches the HGNC's gene names database; see UpdateSymbolList for more details
		) ) ) # AddModuleScore

		

  		# usage:
		# stromal.1, immune.2 from ESTIMATE.signatures
  		# fibroblast.3, endothelial.4, epithelial.5, smooth.6, plasma.7, b_cells.8, mast_cells.9 from PanglaoDB_markers_27_Mar_2020.tsv.gz
		# warning: epithelial.5 does not have good discrimination power. use epi_krt_epcam

	},


	"panglaodb_others"={

		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}
		rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features = list(mast_cells, macrophages, dendritic_cells, t_cells, nk_cells), name = c("Mast.", "Macrophage.", "DC.", "T.", "NK."), search = TRUE) ) )

		# usage: "Mast.1","Macrophage.2", "DC.3","T.4","NK.5"
	},





	"singler_htapp_toolbox"={
		filename_rds_ref_htapp <- list_cancer_type_specific_info$filename_rds_ref_htapp
		if (!is.null(filename_rds_ref_htapp)) {

			cat(sprintf("\t\tread %s\n", filename_rds_ref_htapp))
			ref.data.cancer_type <- readRDS(filename_rds_ref_htapp)

			# MBC:
			# B cell: 447
			# Epithelial cell: 93
			# Macrophage: 280
			# NK cell: 191
			# T cell: 1711

			#log_obj(table(ref.data.cancer_type$annotate), tab=2)

			# 1) Slyper et al. Nat. Medicine 2020 scRNA-seq ovarian tumor
			ref.data.cancer_type <- as.SingleCellExperiment(ref.data.cancer_type)

			# Use de.method wilcox with scRNA-seq reference b/c the reference data is more sparse
			#predictions.ovar.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=ref.data.cancer_type, labels=ref.data.cancer_type$annotate,de.method = "wilcox")
			predictions.cancer_type.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=ref.data.cancer_type, labels=ref.data.cancer_type$annotate, de.method = "wilcox")

			rna$SingleR.HTAPP_toolbox <- predictions.cancer_type.sc$pruned.labels
			#log_obj(table(rna$SingleR.HTAPP_toolbox), tab=2)

		} # if
	},

	"singler_hpca_cellidx"={
		ref.data.HPCA <- readRDS("reference/celldex/HPCA_celldex.rds")
		predictions.HPCA.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=ref.data.HPCA, labels=ref.data.HPCA$label.main)
		rna$SingleR.HPCA <- predictions.HPCA.sc$pruned.labels
		
	},

	"singler_blueprint_encode"={
		ref.data.BED <- readRDS("reference/celldex/BluePrintEncode_celldex.rds")
		predictions.BED.sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=ref.data.BED, labels=ref.data.BED$label.main)
		rna$SingleR.BED <- predictions.BED.sc$pruned.labels
	},



	"singler_swarbricklab_breast_normal_epi"={
		se_centroid <- readRDS("reference/centroids_for_singler/swarbricklab_breast_normal_epi_3centroids.rds")
		predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
		rna$breast_normal_epi_type <- predictions_sc$pruned.labels
	},

	"singler_kessenbrocklab_breast_normal_epi"={
		se_centroid <- readRDS("reference/centroids_for_singler/kessenbrocklab_breast_normal_epi_centroids.rds")
		predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
		rna$breast_normal_epi_type <- predictions_sc$pruned.labels
	},

	"singler_scbcsubtype"={
		# singler_swarbricklab_breast_cancer_epi
		se_centroid <- readRDS(list_cancer_type_specific_info$filename_rds_centroids)
		predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
		rna$breast_cancer_epi_type <- predictions_sc$pruned.labels
	},

	"biomarkers"={

		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}
		rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features=list(c("ESR1"), c("PGR"), c("ERBB2")), nbin=24, ctrl=100, assay=NULL, name=c("esr1.", "pgr.", "erbb2."), search=T) ) )
		rna$biomarker <- ""
		rna$biomarker[rna$esr1.1 > 1.0] <- "ESR1"
		rna$biomarker[rna$pgr.2 > 1.0] <- "PGR"
		rna$biomarker[rna$erbb2.3 > 1.0] <- "ERBB2"
	},

	"basal_luminal_epithelial_cells"={
		# Keratins as markers that distinguish normal and tumor-derived mammary epithelial cells. It is shown here that normal cells produce keratins K5, K6, K7, K14, and K17, whereas tumor cells produce mainly keratins K8, K18, and K19. In immortalized cells, which are preneoplastic or partially transformed, the levels of K5 mRNA and protein are lower than in normal cells, whereas the amount of K18 is increased. Thus, K5 is an important marker in the tumorigenic process, distinguishing normal from tumor cells, and decreased K5 expression correlates with tumorigenic progression.  https://www.pnas.org/content/87/6/2319.long
		# KRT5, KRT6A, KRT14, KRT17: basal epithelial markers, silent in ER+ breast cancer
		# KRT8, KRT18: luminal epithelial markers
		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}
		rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features=list(c("KRT5", "KRT6A", "KRT7", "KRT14", "KRT17"), c("KRT8", "KRT18", "KRT19")), nbin=24, ctrl=100, assay=NULL, name=c("basal_epi.", "luminal_epi."), search=TRUE) ) )

		rna$epi_basal_luminal <- ""
		rna$epi_basal_luminal[rna$basal_epi.1 > 1.0] <- "basal_epi"
		rna$epi_basal_luminal[rna$luminal_epi.2 > 1.0] <- "luminal_epi"
	},

	"epi_krt_epcam"={

		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}
		# basal marker: KRT5, KRT14, KRT17
		# luminal marker: KRT8/18, KRT19, EPCAM
		rna <- suppressWarnings(suppressMessages( AddModuleScore(rna, features=list(c("KRT14", "KRT5", "AKR1B1", "CAV1", "GYPC", "MYL9", "MYLK", "PPP1R14A"), c("KRT18", "KRT19", "CD24", "AZGP1", "C15orf48", "DEFB1", "PDZK1IP1", "RHOV")), nbin=2, ctrl=1, assay=NULL, name=c("epi_basal.", "epi_luminal."), search=TRUE) ))

		rna$epi_krt_epcam <- "nonepi"
		rna$epi_krt_epcam[rna$epi_basal.1 > 1.0] <- "epi"
		rna$epi_krt_epcam[rna$epi_luminal.2 > 1.0] <- "epi"

	},

	"normal_epithelial_cells"={

		switch(list_cancer_type_specific_info$type_normal_epi_type,

			"peroulab"={
				fname_gmt <- "reference/gmt/epithelial_signatures_sc_from_peroulab.gmt"
				source("./r/enrichment_analysis.R") # for GSA.read.gmt_modified
				gset_obj <- GSA.read.gmt_modified(fname_gmt)
				names <- paste0(tolower(gset_obj$geneset.names), ".")
				#names <- list("Luminal_Progenitor_", "Basal_", "Mature_Luminal_", "Epithelial_")

				if (n_log > 0) {
					cat(sprintf("\t\tAddModuleScore for %s\n", method))
					cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
				}
				rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features=gset_obj$genesets, nbin=24, ctrl=100, assay=NULL, name=paste0(tolower(gset_obj$geneset.names), "."), search=T) ) )
				# usage: normal_luminalprog.1, normal_basal.2, normal_matureluminal.3, epithelialcell_signature.4
			},

			"swarbricklab"={
				filename_centroid_rds <- "reference/centroids_for_singler/swarbricklab_breast_normal_epi_3centroids.rds"
				cat(sprintf("\t\treadRDS('%s')\n", filename_centroid_rds))
				se_centroid <- readRDS(filename_centroid_rds)
				predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
				rna$breast_normal_epi_type <- predictions_sc$pruned.labels

				# usage: rna$breast_normal_epi_type
			},

			"kessenbrocklab"={
				filename_centroid_rds <- "reference/centroids_for_singler/kessenbrocklab_breast_normal_epi_5centroids.rds"
				cat(sprintf("\t\treadRDS('%s')\n", filename_centroid_rds))
				se_centroid <- readRDS(filename_centroid_rds)
				predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
				rna$breast_normal_epi_type <- predictions_sc$pruned.labels

				# usage: rna$breast_normal_epi_type

				# https://www.nature.com/articles/s41467-018-04334-1
				list_markers <- list("BEp_MaSCs"=c("TCF4", "ZEB1"))
				rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features=list_markers, nbin=24, ctrl=100, assay=NULL, name=paste0(tolower(names(list_markers)), "."), search=T) ) )
				# usage: bep_mascs.1, bep_myo.2, lep_secretory.3
				f.bep <- grepl("BEp", rna$breast_normal_epi_type)
				rna$breast_normal_epi_type[f.bep & (rna$bep_mascs.1 > 1.0)] <- "BEp_MaSCs"
			},

			"kessenbrocklab_test"={
				filename_centroid_rds <- "reference/centroids_for_singler/kessenbrocklab_breast_normal_epi_3centroids.rds" 
				cat(sprintf("\t\treadRDS('%s')\n", filename_centroid_rds))
				se_centroid <- readRDS(filename_centroid_rds)
				predictions_sc <- SingleR(test=rna.sce, assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)
				rna$breast_normal_epi_type <- predictions_sc$pruned.labels
				# https://www.nature.com/articles/s41467-018-04334-1
				list_markers <- list("BEp_MaSCs"=c("TCF4", "ZEB1"), "BEp_myo"=c("TAGLN", "MYLK", "KRT6A", "ID1", "ACTG2", "ACTA2", "IGFBP4", "CALML3", "CAPS", "SPARC"), "LEp_secretory"=c("S100A8", "LTF", "FDCSP", "SERPINB4", "OVOS2", "AKR1C3", "AKR1C2", "AKR1C1", "FBLN5", "OLFM4"))
				rna <- suppressWarnings( suppressMessages( AddModuleScore(rna, features=list_markers, nbin=24, ctrl=100, assay=NULL, name=paste0(tolower(names(list_markers)), "."), search=T) ) )
				# usage: bep_mascs.1, bep_myo.2, lep_secretory.3
				f.bep <- grepl("BEp", rna$breast_normal_epi_type)
				f.lep_prog <- grepl("LEp_prog", rna$breast_normal_epi_type)
				rna$breast_normal_epi_type[f.bep & (rna$bep_mascs.1 > 1.0)] <- "BEp_MaSCs"
				rna$breast_normal_epi_type[f.bep & (rna$bep_myo.2 > 1.0)] <- "BEp_myo"
				rna$breast_normal_epi_type[f.lep_prog & (rna$lep_secretory.3 > 1.0)] <- "LEp_secretory"
				f.lep <- grepl("LEp$", rna$breast_normal_epi_type)
				rna$breast_normal_epi_type[f.lep] <- "LEp_hormone"
			},
			{}
		) # switch

	},

	"cancer_epithelial_cells"={

		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}

		cols <- list_cancer_type_specific_info$cols
		list_genes <- list()
		for (col in cols) {
		  list_genes[[col]] <- select_genes_with_seurat_obj(rna, df_sig[,col])	
		}

		# scale data
		if (n_log > 0) {
			cat(sprintf("\t\tScaleData for %s\n", method))
			cat(sprintf("\t\tScaleData for %s\n", method), file=stderr())
		}

		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			idx_cells <- 1:ncol(rna)
		} # if

		# scaling data
		# results depend on the collection of samples due to scaling
		# https://satijalab.org/seurat/reference/scaledata
		rna.tmp <- ScaleData(rna[,idx_cells],
				features = unlist(list_genes),
				assay = NULL,
				vars.to.regress = NULL,
				split.by = NULL,
				do.scale = TRUE,
				do.center = TRUE,
				verbose = FALSE)

		source("./r/seurat/seurat_utilities.R") # for AddModuleScore_modified
		if (n_log > 0) {
			cat(sprintf("\t\tAddModuleScore for %s\n", method))
			cat(sprintf("\t\tAddModuleScore for %s\n", method), file=stderr())
		}

		rna.tmp <- suppressWarnings( suppressMessages( 
			AddModuleScore_modified(rna.tmp,
				features = list_genes,
				pool = NULL,
				nbin = 24,
				ctrl = 10, # default=100
				assay = NULL,
				slot = "scale.data",
				name = paste0(cols, "."),
				seed = 1,
				search = T) 
			) )

		cols <- paste0(cols, ".", 1:length(cols))
		rna@meta.data[idx_cell, cols] <- rna.tmp@meta.data[,cols] 
		# usage: Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4

		remove(list = c("rna.tmp"))

	},


	"celltype_garnett"={

		
		#for (type_classifier in c("major", "minor", "subset")) {
		for (type_classifier in c("major")) {
			filename_garnett_classifier <- sprintf("%s_%s.rds", list_cancer_type_specific_info$filename_garnett_classifier_prefix, type_classifier)
			if (!file.exists(filename_garnett_classifier)) {
				next
			}

			#cds <- convert_seurat2monocle(rna, method_create_cds = "as.cell_data_set")
			cds <- convert_seurat2monocle(rna, method_create_cds = "new_cell_data_set")

 			cat(sprintf("\t\treadRDS(%s)\n", filename_garnett_classifier))
			garnett_classifier <- readRDS(filename_garnett_classifier)

			set.seed(42)
 			cat(sprintf("\t\tclassify_cells\n"))
			cds <- classify_cells(cds, garnett_classifier,
				db = org.Hs.eg.db,
				cluster_extend = TRUE,
				cds_gene_id_type = "SYMBOL")

			#rna <- AddMetaData(object = rna, metadata = pData(cds)[, c("cell_type", "cluster_ext_type")], col.name = c(sprintf("celltype_garnett.%s", type_classifier), sprintf("celltype_gatnett.%s.ext", type_classifier)))
			rna <- AddMetaData(object = rna, metadata = pData(cds)[, "cell_type"], col.name = sprintf("celltype_garnett.%s", type_classifier))

		} # for

	},


	"celltype_td"={

		filename_rds_ref_major <- list_cancer_type_specific_info$filename_rds_ref_major
		if (!is.null(filename_rds_ref_major) && file.exists(filename_rds_ref_major)) {


			# reference based cell type data tansfer
 			cat(sprintf("\t\treadRDS(%s)\n", filename_rds_ref_major))
			rna.reference <- readRDS(filename_rds_ref_major)

			# features
			method_to_determine_features_for_dimension_reduction <- "ref.vf"
			switch(method_to_determine_features_for_dimension_reduction,
				"null"={
					# default=NULL, Features to use for dimensional reduction. If not specified, set as variable features of the reference object which are also present in the query.
					# rna.reference, variable features: could be large (e.g. 10781)
					# rna: variable features: 2000
					features <- NULL
				},
				"ref.vf"={
					rna.tmp <- FindVariableFeatures(rna.reference, nfeatures=2000)
					features <- intersect(VariableFeatures(rna.tmp), rownames(rna))
					# rna.reference, variable features: 2000
					# rna: row.names: may be > 15,000
				},
				"ref.vf_query.vf"={
					rna.tmp <- FindVariableFeatures(rna.reference, nfeatures=2000)
					features <- intersect(VariableFeatures(rna.tmp), VariableFeatures(rna))
					# rna.reference, variable features: 2000
					# rna: variable features: 2000
				},
				"query.vf"={
					features <- intersect(rownames(rna.reference), VariableFeatures(rna))
					# rna.reference, row.names: may be > 27,000
					# rna: variable features: 2000
				},
				{}
			) # switch

			# https://learn.gencore.bio.nyu.edu/seurat-integration-and-label-transfer/
 			cat(sprintf("\t\tFindTransferAnchors with n_features=%d\n", length(features)))
			# https://satijalab.org/seurat/reference/findtransferanchors
			anchors <- FindTransferAnchors(reference = rna.reference,
					query = rna,
					normalization.method = "LogNormalize",
					recompute.residuals = TRUE,
					reference.assay = NULL,
					reference.neighbors = NULL,
					query.assay = NULL,
					reduction = "pcaproject", # Dimensional reduction to perform when finding anchors. Options are: pcaproject: Project the PCA from the reference onto the query. We recommend using PCA when reference and query datasets are from scRNA-seq
					reference.reduction = NULL, # Name of dimensional reduction to use from the reference if running the pcaproject workflow. Optionally enables reuse of precomputed reference dimensional reduction. If NULL (default), use a PCA computed on the reference object.
					project.query = FALSE,
					features = features, # default=NULL, Features to use for dimensional reduction. If not specified, set as variable features of the reference object which are also present in the query.
					scale = TRUE,
					npcs = args$max_dimstouse, # default=30
					l2.norm = TRUE,
					dims = args$dimsToUse, # default=1:30
					k.anchor = 5,
					k.filter = 200,
					k.score = 30,
					max.features = 200,
					nn.method = "annoy",
					n.trees = 50,
					eps = 0,
					approx.pca = TRUE,
					mapping.score.k = NULL,
					verbose = TRUE
				) # FindTransferAnchors

			cols <- colnames(rna.reference@meta.data)
			if ("celltype_major" %in% cols) {
				# major
 				cat(sprintf("\t\tTransferData with celltype_major\n"))
				predictions <- TransferData(anchorset = anchors, refdata = rna.reference$celltype_major, dims = args$dimsToUse)
				rna <- AddMetaData(object = rna, metadata = predictions[, c("predicted.id", "prediction.score.max")], col.name = c("celltype_td.major", "celltype_td.major.score.max"))
			}
	
			if ("celltype_minor" %in% cols) {
				# minor
 				cat(sprintf("\t\tTransferData with celltype_minor\n"))
				predictions <- TransferData(anchorset = anchors, refdata = rna.reference$celltype_minor, dims = args$dimsToUse)
				rna <- AddMetaData(object = rna, metadata = predictions[, c("predicted.id", "prediction.score.max")], col.name = c("celltype_td.minor", "celltype_td.minor.score.max"))
			}

			if ("celltype_subset" %in% cols) {
				# subset
 				cat(sprintf("\t\tTransferData with celltype_subset\n"))
				predictions <- TransferData(anchorset = anchors, refdata = rna.reference$celltype_subset, dims = args$dimsToUse)
				rna <- AddMetaData(object = rna, metadata = predictions[, c("predicted.id", "prediction.score.max")], col.name = c("celltype_td.subset", "celltype_td.subset.score.max"))
			}

		} # if

	},


	"celltype_markercount"={

 		cat(sprintf("\t\tclassify_cells\n"))
		df_res <- classify_cells_with_markercount(rna, args, list_cancer_type_specific_info$filename_garnett_marker_file)

		rna <- AddMetaData(object = rna, metadata = df_res[, "cell_type_pred"], col.name = "celltype_markercount")

	},


	"celltype_markercount_ref"={
		
		# not implemented yet

	},



	"celltype_cycling"={

		# Averaged normalized expression of 11 genes34 (BIRC5, CCNB1, CDC20, NUF2, CEP55, NDC80, MKI67, PTTG1, RRM2, TYMS and UBE2C), independent of the SCSubtype gene lists, was used to compute the proliferation score.  https://www.nature.com/articles/s41588-021-00911-1
		rna[["PScore"]] <- NA
		rna[["celltype.cycling"]] <- ""

		df_sig <- data.frame(proliferation=c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C")) # https://www.nature.com/articles/s41588-021-00911-1
		df_score <- calculate_sig_scores(rna, args, df_sig, "epi")
		if (!is.null(df_score)) {
			barcodes <- colnames(df_score)
			pscore_epi <- as.numeric(df_score[1,])
			rna@meta.data[barcodes, "PScore"] <- pscore_epi
			rna@meta.data[barcodes[which(pscore_epi > 0)], "celltype.cycling"] <- "cycling"
		}

		df_score <- calculate_sig_scores(rna, args, df_sig, "nonepi")
		if (!is.null(df_score)) {
			barcodes <- colnames(df_score)
			pscore_nonepi <- as.numeric(df_score[1,])
			rna@meta.data[barcodes, "PScore"] <- pscore_nonepi
			rna@meta.data[barcodes[which(pscore_nonepi > 0)], "celltype.cycling"] <- "cycling"
		}

	},

	


	"brca_cell_atlas_scsubtyper"={

		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}
		cols <- list_cancer_type_specific_info$cols

		temp_allgenes <- c()
		for (col in cols) {
			temp_allgenes <- c(temp_allgenes, as.vector(df_sig[,col]))
		}
		temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])
		f.ok <- (temp_allgenes %in% rownames(rna))

		barcodes <- colnames(rna)
		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\t\tusing CNV positive cells\n"))
			barcodes <- barcodes[grep("11", rna$CNV.Pos)]
		} # if
 		cat(sprintf("\t\t# of barcodes: %d\n", length(barcodes)))



		# reference: https://www.nature.com/articles/s41588-021-00911-1
		# reading in the training datasets
		dir_train <- "./reference/scsubtyper/training_dataset"
		Her2combinedNov2019 <- readRDS(sprintf("%s/Her2combinedNov2019.rds", dir_train)) # seurat obj 18713x5279, assays: RNA (default), integrated (2000), dimension reduction: pca, umap
		LumAcombinedNov2020 <- readRDS(sprintf("%s/LumAcombinedNov2020.rds", dir_train)) # seurat obj 17619x5771, assays: RNA (default), integrated (2000), dimension reduction: pca, umap
		LumBcombinedNov2019 <- readRDS(sprintf("%s/LumBcombinedNov2019.rds", dir_train)) # seurat obj 16507x2483, assays: RNA (default), integrated (2000), dimension reduction: pca, umap
		BasalTrainingdataset2020 <- readRDS(sprintf("%s/BasalTrainingdataset2020.rds", dir_train)) # seurat obj 20226x3000, assays: RNA (default), integrated (2000), dimension reduction: none

		# rna.tmp
		rna.tmp <- rna[, barcodes]
		rna.tmp <- RenameCells(rna.tmp, new.names=paste0("scsubtyper_", colnames(rna.tmp)))
		barcodes.tmp <- colnames(rna.tmp)
		rna.tmp$orig.ident <- args$sample_id
		rna.tmp <- subset(rna.tmp, subset = nFeature_RNA > 200)

		type_merge <- "merge"
		switch(type_merge,
			"merge"={
 				cat(sprintf("\t\tmerge\n"))
				tumors.combined <- merge(Her2combinedNov2019, y=c(LumAcombinedNov2020, LumBcombinedNov2019, BasalTrainingdataset2020, rna.tmp))
			},
			"integrate"={
 				cat(sprintf("\t\tNormalizeData\n"))
				rna.tmp <- NormalizeData(rna.tmp, verbose = FALSE)
 				cat(sprintf("\t\tFindVariableFeatures\n"))
				rna.tmp <- FindVariableFeatures(rna.tmp, selection.method = "vst", nfeatures = 3000)

				# integrating the testing with the training datasets
		 		cat(sprintf("\t\tFindIntegrationAnchors\n"))
				tumors.anchors <- FindIntegrationAnchors(object.list = list(BasalTrainingdataset2020, LumAcombinedNov2020, Her2combinedNov2019, LumBcombinedNov2019, rna.tmp), dims = 1:50)

		 		cat(sprintf("\t\tIntegrateData\n"))
				tumors.combined <- IntegrateData(anchorset = tumors.anchors, dims = 1:50)
			},
			{}
		) # switch
		Idents(object = tumors.combined) <- "orig.ident"


		# performing normalizing and scaling the whole dataset
		DefaultAssay(tumors.combined) <- "RNA" # Aatish's code
		#DefaultAssay(tumors.combined) <- "integrated" # test
		assay.tmp <- DefaultAssay(tumors.combined)
		tumors.combined <- PercentageFeatureSet(tumors.combined, pattern = "^MT-", col.name="percent.mt")
		switch(assay.tmp,
			"RNA"={
 				cat(sprintf("\t\tNormalizeData assay=%s\n", assay.tmp))
				tumors.combined <- NormalizeData(object = tumors.combined)
 				cat(sprintf("\t\tFindVariableFeatures assay=%s\n", assay.tmp))
				tumors.combined <- FindVariableFeatures(object = tumors.combined, selection.method = "vst", nfeatures = 3000)
			},
			{}
		) # switch




		# scatterplot
		f_plot_scatterplot <- FALSE
		if (f_plot_scatterplot) {
			plot1 <- FeatureScatter(object = tumors.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
			plot2 <- FeatureScatter(object = tumors.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
			#CombinePlots(plots = list(plot1, plot2))
			#patchwork::wrap_plots(list(plot1, plot2), ncol = 2)
			ggsave(sprintf("%s/scatterplot_scsubtyper_ncount_vs_persent.mt_%s_%s.pdf", dir_pdf, args$cancer_type, args$sample_id), width = 7, height = 4, plot=plot1)
			ggsave(sprintf("%s/scatterplot_scsubtyper_ncount_vs_nfeature_%s_%s.pdf", dir_pdf, args$cancer_type, args$sample_id), width = 7, height = 4, plot=plot2)
		} # if

		all.genes <- rownames(x = tumors.combined)

		method_to_scale <- "no.regress"
		switch(method_to_scale,
			"no.regress"={
 				cat(sprintf("\t\tScaleData assay=%s\n", assay.tmp))
				tumors.combined <- ScaleData(object = tumors.combined, features = all.genes)
			},
			"vars.to.regress"={
				# regressing out percent.mt is slow when features = all.genes
 				cat(sprintf("\t\tScaleData assay=%s, var.to.regress=%s\n", assay.tmp, paste(args$vars.to.regress, collapse=", ")))
				tumors.combined <- ScaleData(object = tumors.combined, features = all.genes, vars.to.regress = args$vars.to.regress, verbose = FALSE)
			},
			{}
		) # switch

		# signature score calculations
		#tocalc <- as.data.frame(tumors.combined@assays$RNA@scale.data)
		tocalc <- GetAssayData(object = tumors.combined,
				 assay = assay.tmp,
				 slot = "scale.data")

		outdat <- matrix(0,
				nrow = ncol(df_sig),
				ncol = ncol(tocalc),
				dimnames = list(colnames(df_sig),
				colnames(tocalc)))

		#for(i in 1:nrow(sigdat)) {
		for(i in 1:ncol(df_sig)) {
			#sigdat[i,!is.na(sigdat[i,])]->module
			#row <- as.character(unlist(module))
		  	genes <- as.character(df_sig[,i])
			genes <- unique(genes[genes != ""])
			idx_genes <- which(rownames(tocalc) %in% genes)
			temp <- apply(tocalc[idx_genes,], 2,
				function(x){
					mean(as.numeric(x), na.rm=TRUE)
				} )
			outdat[i,] <- as.numeric(temp)
		} # for

		final <- outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
		final <- as.data.frame(final)
		is.num <- sapply(final, is.numeric)
		final[is.num] <- lapply(final[is.num], round, 4)
		finalm <- as.matrix(final)
		center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
			get_average <- function(v) sum(v * row.w)/sum(row.w)
			average <- apply(x, 2, get_average)
			sweep(x, 2, average)
		}
		finalmt <- as.data.frame(t(finalm))
		finalm.sweep.t <- center_sweep(finalmt)
		finalnames <- colnames(finalm.sweep.t)[max.col(finalm.sweep.t, ties.method="first")]

		# writing out the scores and calls
		#write.table(finalm.sweep.t, sprintf("%s/scsubtyper_%s_%s_scores.tsv", dir_tsv, args$cancer_type, args$sample_id), sep="\t")
		#write.table(Finalnames, sprintf("%s/scsubtyper_%s_%s_calls.tsv", dir_pdf, args$cancer_type, args$sample_id), sep="\t")

		# update scores
		finalm.sweep.t$SCSubtypeCall <- finalnames
		cols <- paste0(cols, ".", 1:length(cols))
		colnames(finalm.sweep.t) <- c(cols, "SCSubtypeCall")

		# display stats
		if ("SingleR.HTAPP_toolbox" %in% colnames(rna.tmp@meta.data)) {
 			cat(sprintf("\t\tSingleR.HTAPP_toolbox with Epi.\n"))
			idx <- grep("Epi", rna.tmp@meta.data[barcodes.tmp, "SingleR.HTAPP_toolbox"])
			log_obj(table(finalm.sweep.t[barcodes.tmp[idx], "SCSubtypeCall"]), tab=2)
		}

		if ("SingleR.HPCA" %in% colnames(rna.tmp@meta.data)) {
 			cat(sprintf("\t\tSingleR.HPCA with Epi.\n"))
			idx <- grep("Epi", rna.tmp@meta.data[barcodes.tmp, "SingleR.HPCA"])
			log_obj(table(finalm.sweep.t[barcodes.tmp[idx], "SCSubtypeCall"]), tab=2)
		}

		if ("SingleR.BED" %in% colnames(rna.tmp@meta.data)) {
 			cat(sprintf("\t\tSingleR.BED with Epi.\n"))
			idx <- grep("Epi", rna.tmp@meta.data[barcodes.tmp, "SingleR.BED"])
			log_obj(table(finalm.sweep.t[barcodes.tmp[idx], "SCSubtypeCall"]), tab=2)
		}

		# update rna@meta.data
		rna@meta.data[barcodes, cols] <- finalm.sweep.t[barcodes.tmp, cols] 
		# usage: Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4

		remove(list = c("final", "finalmt", "finalm.sweep.t"))

	},

	"brca_cell_atlas_scsubtyper_with_signatures"={

		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}
		cols <- list_cancer_type_specific_info$cols

		temp_allgenes <- c()
		for (col in cols) {
			temp_allgenes <- c(temp_allgenes, as.vector(df_sig[,col]))
		}
		temp_allgenes <- unique(temp_allgenes[!temp_allgenes == ""])
		f.ok <- (temp_allgenes %in% rownames(rna))

		# scale data
		if (n_log > 0) {
			cat(sprintf("\t\tScaleData for %s\n", method))
			cat(sprintf("\t\tScaleData for %s\n", method), file=stderr())
		}

		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			idx_cells <- 1:ncol(rna)
		} # if

		# scaling data
		# results depend on the collection of samples due to scaling
		# https://satijalab.org/seurat/reference/scaledata
		rna.tmp <- ScaleData(rna[,idx_cells],
				features = temp_allgenes[f.ok],
				assay = NULL,
				vars.to.regress = NULL,
				split.by = NULL,
				do.scale = TRUE,
				do.center = TRUE,
				verbose = FALSE)

		#tocalc <- as.data.frame(rna.tmp@assays$RNA@scale.data)
		tocalc <- GetAssayData(object = rna.tmp,
				 assay = NULL,
				 slot = "scale.data")
		remove(list = c("rna.tmp"))

		#calculate mean scsubtype scores
		outdat <- matrix(0,
                 	nrow=ncol(df_sig),
                 	ncol=ncol(tocalc),
                 	dimnames=list(colnames(df_sig),
                               colnames(tocalc)))

		#for(i in 1:nrow(sigdat)) {
		for (i in 1:ncol(df_sig)) {
  
		  # sigdat[i,!is.na(sigdat[i,])]->module
		  # row <- as.character(unlist(module))
		  genes <- as.character(df_sig[,i])
		  genes <- unique(genes[genes != ""])
		  idx_genes <- which(rownames(tocalc) %in% genes)
		  temp <- apply( tocalc[idx_genes,], 2,
			function(x) {
			  mean(as.numeric(x), na.rm=TRUE)
			} )
		  outdat[i,] <- as.numeric(temp)

		} # for

		final <- outdat[which(rowSums(outdat,na.rm=TRUE)!=0),]
		final <- as.data.frame(final)
		is.num <- sapply(final, is.numeric)
		final[is.num] <- lapply(final[is.num], round, 4)
		finalm <- as.matrix(final)

		# scaling scores function before calling the highest call
		center_sweep <- function(x, row.w = rep(1, nrow(x))/nrow(x)) {
		  get_average <- function(v) sum(v * row.w)/sum(row.w)
		  average <- apply(x, 2, get_average)
		  sweep(x, 2, average)
		}

		# obtaining the highest call
		finalmt <- as.data.frame(t(finalm))
		finalm.sweep.t <- center_sweep(finalmt)
		finalnames <- colnames(finalm.sweep.t)[max.col(finalm.sweep.t, ties.method="first")]

		# update scores
		finalm.sweep.t$SCSubtypeCall <- finalnames
		cols <- paste0(cols, ".", 1:length(cols))
		colnames(finalm.sweep.t) <- c(cols, "SCSubtypeCall")
		rna@meta.data[idx_cells, cols] <- finalm.sweep.t[,cols] 
		# usage: Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4

		remove(list = c("final", "finalmt", "finalm.sweep.t"))

	},


	"scbcsubtype"={

		# swarbricklab_breast_cancer_epi
		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			#idx_cells <- 1:ncol(rna)
			vec_cell_types <- get_vec_cell_types(rna, args)
			idx_cells <- grep(pattern_epi, vec_cell_types)
		} # if

		med_nfeatures <- median(rna$nFeature_RNA)
                mtx <- GetAssayData(object = rna, assay = "RNA", slot = "data")
		
		cat(sprintf("\t\treadRDS('%s')\n", list_cancer_type_specific_info$filename_rds_centroids))
		se_centroid <- readRDS(list_cancer_type_specific_info$filename_rds_centroids)
		mtx_centroid <- assay(se_centroid)
		df_rowdata <- rowData(se_centroid)
		syms_centroid <- intersect(rownames(mtx_centroid), rownames(mtx))
		mtx_centroid <- mtx_centroid[syms_centroid,]
		syms_centroid_dge <- syms_centroid
		mtx_centroid_dge <- mtx_centroid
		if (nrow(df_rowdata) > 0) {
			if ("f_dge" %in% colnames(df_rowdata)) {
				syms_centroid_dge <- intersect(syms_centroid, rownames(df_rowdata)[df_rowdata$f_dge])
				mtx_centroid_dge <- mtx_centroid[syms_centroid_dge,]
			}
		}

		# mtx_expr: logcount for syms_centroid and tumor cells
		mtx_expr <- mtx[syms_centroid, idx_cells]
		subtypes <- colnames(mtx_centroid)
		subtypes_syms_centroid <- subtypes[max.col(mtx_centroid, ties.method="first")]
		subtypes_syms_centroid_dge <- subtypes[max.col(mtx_centroid_dge, ties.method="first")]
		nv_subtypes_syms_centroid <- subtypes_syms_centroid
		names(nv_subtypes_syms_centroid) <- rownames(mtx_centroid)
		
		suppressPackageStartupMessages(library(ICIKendallTau))
		mtx.corr <- matrix(0, ncol(mtx_expr), length(subtypes), dimnames=list(colnames(mtx_expr), subtypes))
		mtx.nnzero <- matrix(0, ncol(mtx_expr), length(subtypes), dimnames=list(colnames(mtx_expr), subtypes))
		for (i in seq(subtypes)) {
			cat(sprintf("\t\t%s: compute Kendall's tau between %s centroid and logcount\n", subtypes[i], subtypes[i]))
			vec_centroid <- mtx_centroid[,i]
			#syms_centroid_subtype <- syms_centroid[(subtypes_syms_centroid == subtypes[i])]
			syms_centroid_subtype <- syms_centroid_dge[(subtypes_syms_centroid_dge == subtypes[i])]
			mtx.nnzero[,i] <- diff(mtx_expr[syms_centroid_subtype,]@p) # nnz per column
			for (j in 1:ncol(mtx_expr)) {
				vec_sample <- mtx_expr[,j]
				if (length(unique(vec_sample)) == 1) next	
				nv_cor <- ici_kt(vec_centroid, vec_sample)
				mtx.corr[j, i] <- nv_cor[["tau"]]
				# nv.cor[["pvalue"]]
			}
		} # for


		# initializing rna$breast_cancer_epi_type for calling subtypes
		rna$breast_cancer_epi_type <- NA

		list_rules <- metadata(se_centroid)[["list_rules"]]
		if (!is.null(list_rules)) {

			for (subtype in names(list_rules)) {

				list_info <- list_rules[[subtype]]
				if (!"min_expr" %in% names(list_info)) next

				nv_min_expr <- list_info[["min_expr"]]
				if (length(nv_min_expr) < 1) next

				# special treatment
				if (med_nfeatures < 1250) {
					# normalized log counts can be small due to lower # of genes per cells.
					next
				}

				if (subtype == "Basal") {
					if (grepl("mbc", args$cancer_subtype)) {
						# metastatic breast cancer cells usually contain basal cells.	
						next
					}
				} # if

				# confident exclusion of subtype candidates by rules
				idx_subtype <- max.col(mtx.corr, ties.method="first")
				subtypes_called.1st <- subtypes[idx_subtype]
				vec_corr.1st <- mtx.corr[cbind(1:nrow(mtx.corr), idx_subtype)]
				for (gene in names(nv_min_expr)) {
					# some genes need to be highly expresed for calling a subtype, e.g. Her2E ERBB2 > 1.0
					min_v <- nv_min_expr[gene]
					genes <- strsplit(gene, "\\|")[[1]]
					genes <- genes[genes %in% rownames(mtx_expr)]
					if (length(genes) < 1) {
						vec_expr <- rep(0, ncol(mtx_expr))
					} else if (length(genes) == 1) {
						vec_expr <- mtx_expr[genes,,drop=T]
					} else {
						vec_expr <- sparseMatrixStats::colMaxs(mtx_expr[genes,,drop=F], na.rm=T)
					}
					idx_low_expr <- which(vec_expr < min_v)
					idx_low_expr_in_subtype <- which((vec_expr < min_v) & (subtypes_called.1st == subtype))
 					cat(sprintf("\t\tdo not call %s when %s < %g: %d", subtype, gene, min_v, length(idx_low_expr_in_subtype)))

					type_reset <- "check_next_choice"
					switch(type_reset,
						"strict"={
							mtx.corr[idx_low_expr, subtype] <- -100
 							cat(sprintf(", reset: %d\n", length(intersect(idx_low_expr_in_subtype, idx_low_expr))))
						},
						"check_next_choice"={
							# check next choice
							mtx.corr.2nd <- mtx.corr
							mtx.corr.2nd[idx_low_expr, subtype] <- -100
							idx_subtype.2nd <- max.col(mtx.corr.2nd, ties.method="first")
							subtypes_called.2nd <- subtypes[idx_subtype.2nd]
							vec_corr.2nd <- mtx.corr.2nd[cbind(1:nrow(mtx.corr), idx_subtype.2nd)]
							# 1st choice was not very confident, corr(2nd choice) was not very low and it was close to corr(1st choice) --> it's okay to choose 2nd choice..
							# CID45171, CID44991 Her2E: many cells ERBB2 < 1, ignore Her2E calls when corr < 0.30 and the next choice would be still reasonable.
							f <- (vec_corr.1st[idx_low_expr] < 0.30) & (vec_corr.2nd[idx_low_expr] > 0.10) & (vec_corr.1st[idx_low_expr] - vec_corr.2nd[idx_low_expr] < 0.10)
							if (min_v > 0) {
								# exclude zero expression
								# CD4523 (MBC, TNBC): many cells ERBB2 == 0, ignore Her2E calls since the similar expression pattern might be risen by other routes than ERBB2 amp, e.g., other RTK-MAPK signals. 
								f <- f | (vec_expr <= 0)
							}
							# ERBB2 can be 0 or very low when nCount_RNA <= 1000 or nFeature_RNA <= 500.
							f <- f & (rna@meta.data[idx_cells, "nCount_RNA"] > 1000)
							f <- f & (rna@meta.data[idx_cells, "nFeature_RNA"] > 500)
							mtx.corr[idx_low_expr[f], subtype] <- -100
 							cat(sprintf(", reset: %d\n", length(intersect(idx_low_expr_in_subtype, idx_low_expr[f]))))
						},
						{}
					) # switch
				} # for gene
			} # for subtype

			for (subtype in names(list_rules)) {
				list_info <- list_rules[[subtype]]
				if (!"normal_epi_type" %in% names(list_info) || (!"normal_epi_type" %in% colnames(rna@meta.data))) next
				if (length(list_info[["normal_epi_type"]]) < 1) next

				# confident exclusion of subtype candidates by rules
				idx_subtype <- max.col(mtx.corr, ties.method="first")
				subtypes_called.1st <- subtypes[idx_subtype]
				vec_corr.1st <- mtx.corr[cbind(1:nrow(mtx.corr), idx_subtype)]

				# some subtypes were related with normal epithelial types, e.g. Basal subtype comes from BEp or LEp_prog.
				pattern_normal_epi_type <- paste(list_info[["normal_epi_type"]], collapse="|")
				idx_not_normal_epi_type <- which(!grepl(pattern_normal_epi_type, rna@meta.data[idx_cells, "normal_epi_type"]))
				idx_not_normal_epi_type_in_subtype <- which(!grepl(pattern_normal_epi_type, rna@meta.data[idx_cells, "normal_epi_type"]) & (subtypes_called.1st == "Basal"))
 				cat(sprintf("\t\tdo not call %s when epithelical cells are not %s: %d", subtype, pattern_normal_epi_type, length(idx_not_normal_epi_type_in_subtype)))

				type_reset <- "check_next_choice"
				switch(type_reset,
					"strict"={
						mtx.corr[idx_not_normal_epi_type, subtype] <- -100
 						cat(sprintf(", reset: %d\n", length(intersect(idx_not_normal_epi_type_in_subtype, idx_not_normal_epi_type))))
					},
					"check_next_choice"={
						# check next choice
						mtx.corr.2nd <- mtx.corr
						mtx.corr.2nd[idx_not_normal_epi_type, subtype] <- -100
						idx_subtype.2nd <- max.col(mtx.corr.2nd, ties.method="first")
						subtypes_called.2nd <- subtypes[idx_subtype.2nd]
						vec_corr.2nd <- mtx.corr.2nd[cbind(1:nrow(mtx.corr), idx_subtype.2nd)]
						# 1st choice was not very confident, corr(2nd choice) was not very low and it was close to corr(1st choice) --> it's okay to choose 2nd choice..
						# CID4513 Basal: many LEp
						f <- (vec_corr.1st[idx_not_normal_epi_type] < 0.30) & (vec_corr.2nd[idx_not_normal_epi_type] > 0.10) & (vec_corr.1st[idx_not_normal_epi_type] - vec_corr.2nd[idx_not_normal_epi_type] < 0.10)
						mtx.corr[idx_not_normal_epi_type[f], subtype] <- -100
 						cat(sprintf(", reset: %d\n", length(intersect(idx_not_normal_epi_type_in_subtype, idx_not_normal_epi_type[f]))))
					},
					{}
				) # switch
			} # for subtype
		} # if



		# call subtypes
		idx_subtype <- max.col(mtx.corr, ties.method="first")
		rna@meta.data[idx_cells, "breast_cancer_epi_type"] <- subtypes[idx_subtype]
		vec_corr <- mtx.corr[cbind(1:nrow(mtx.corr), idx_subtype)]
		vec_confidence_check <- rep(TRUE, length(vec_corr))

		# subsequential classifcation
		list_biomarkers <- metadata(se_centroid)[["list_biomarkers"]]
		if (!is.null(list_biomarkers)) {
			for (subtype in names(list_biomarkers)) {
				list_info <- list_biomarkers[[subtype]]
				idx_cells_subtype <- idx_cells
				if ("subtype of" %in% names(list_info)) {
					x_subtype <- which(rna@meta.data[idx_cells, "breast_cancer_epi_type"] == list_info[["subtype of"]])
					idx_cells_subtype <- idx_cells_subtype[x_subtype]
				}
				syms_expr <- list_info[["expressed"]]	
				syms_expr <- syms_expr[syms_expr %in% rownames(mtx)]
				syms_not_expr <- list_info[["not expressed"]]	
				syms_not_expr <- syms_not_expr[syms_not_expr %in% rownames(mtx)]
				if (length(syms_expr) == 0) {
					min_expr <- rep(0, length(idx_cells_subtype))
					mean_expr <- rep(0, length(idx_cells_subtype))
				} else if  (length(syms_expr) == 1) {
					min_expr <- min(mtx[syms_expr, idx_cells_subtype], na.rm=T)
					mean_expr <- mean(mtx[syms_expr, idx_cells_subtype], na.rm=T)
				} else {
					min_expr <- sparseMatrixStats::colMins(mtx[syms_expr, idx_cells_subtype], na.rm=T)
					mean_expr <- colMeans(mtx[syms_expr, idx_cells_subtype], na.rm=T)
				}

				if (length(syms_not_expr) == 0) {
					max_not_expr <- rep(0, length(idx_cells_subtype))
					mean_not_expr <- rep(0, length(idx_cells_subtype))
				} else if (length(syms_not_expr) == 1) {
					max_not_expr <- max(mtx[syms_not_expr, idx_cells_subtype], na.rm=T)
					mean_not_expr <- mean(mtx[syms_not_expr, idx_cells_subtype], na.rm=T)
				} else {
					max_not_expr <- sparseMatrixStats::colMaxs(mtx[syms_not_expr, idx_cells_subtype], na.rm=T)
					mean_not_expr <- colMeans(mtx[syms_not_expr, idx_cells_subtype], na.rm=T)
				}
				x <- which((min_expr > 1.0) & (max_not_expr < 0.1))
 				cat(sprintf("\t\tcall %s when min(%s) > 1.0 and max(%s) < 0.1: %d\n", subtype, paste(syms_expr, collapse=", "), paste(syms_not_expr, collapse=", "), length(idx_cells_subtype[x])))
				#x <- which((mean_expr > 1.0) & (mean_not_expr < 0.1))
 				#cat(sprintf("\t\tcall %s when mean(%s) > 1.0 and mean(%s) < 0.1: %d\n", subtype, paste(syms_expr, collapse=", "), paste(syms_not_expr, collapse=", "), length(idx_cells_subtype[x])))
				rna@meta.data[idx_cells_subtype[x], "breast_cancer_epi_type"] <- subtype
				idx <- match(idx_cells_subtype[x], idx_cells)
				f <- !is.na(idx); idx <- idx[f]
				vec_confidence_check[idx] <- FALSE
			} # for
		} # if


		# consider confidence
		rna$corr <- NA
		rna$corr[idx_cells] <- vec_corr
		df <- inspect_gene_expr_distributions(rna, "corr", col_cell_types="breast_cancer_epi_type", n_log=1)
		method_to_exclude_subtype_calls <- "none"
		#method_to_exclude_subtype_calls <- "corr.magnitude"
		switch(method_to_exclude_subtype_calls,
			"none"={
				# do not exclude subtype calls
			},
			"corr.magnitude"={
				# the effect of corr.magnitude is minimal for good dataset.
				vec_corr <- mtx.corr[cbind(1:nrow(mtx.corr), idx_subtype)]
				f_exclude <- (vec_corr < 0.1)
				f_exclude <- f_exclude & vec_confidence_check
 				cat(sprintf("\t\tdo not call subtype when corr < 0.1: %d\n", length(which(f_exclude))))
				rna@meta.data[idx_cells[f_exclude], "breast_cancer_epi_type"] <- "Epi. Tumor"
			},
			"ncount_nfeature"={
				f_exclude <- (rna@meta.data[idx_cells, "nCount_RNA"] < 1000)
				f_exclude <- f_exclude | (rna@meta.data[idx_cells, "nFeature_RNA"] < 500)
				f_exclude <- f_exclude & vec_confidence_check
 				cat(sprintf("\t\tdo not call subtype when nCount_RNA < 1000 or nFeature_RNA < 500: %d\n", length(which(f_exclude))))
				rna@meta.data[idx_cells[f_exclude], "breast_cancer_epi_type"] <- "Epi. Tumor"
			},
			"nnzero"={
				vec_confidence <- mtx.nnzero[cbind(1:ncol(mtx_expr), idx_subtype)]
				#print(head(vec_confidence,100))
				f_exclude <- (vec_confidence < 2)
				f_exclude <- f_exclude & vec_confidence_check
 				cat(sprintf("\t\tdo not call subtype when nnzero < 2: %d\n", length(which(f_exclude))))
				rna@meta.data[idx_cells[f_exclude], "breast_cancer_epi_type"] <- "Epi. Tumor"
			},
			{ }
		) # switch


		# display information
		df <- inspect_gene_expr_distributions(rna, "nCount_RNA", col_cell_types="breast_cancer_epi_type", n_log=1)
		df <- inspect_gene_expr_distributions(rna, "nFeature_RNA", col_cell_types="breast_cancer_epi_type", n_log=1)
 		cat(sprintf("\t\tnonexisting centroid genes: %s\n", paste(setdiff(rownames(mtx_centroid), syms_centroid), collapse=", ")))
		for (subtype in subtypes) {
			f_cells_subtype <- (rna@meta.data[idx_cells, "breast_cancer_epi_type"] == subtype)
			if (length(which(f_cells_subtype)) == 0) next
			#syms_centroid_subtype <- syms_centroid[subtypes_syms_centroid == subtype]
			syms_centroid_subtype <- syms_centroid_dge[subtypes_syms_centroid_dge == subtype]
			if (length(syms_centroid_subtype) < 2) {
				f_syms_centroid_more_expr <- mean(mtx_expr[syms_centroid_subtype, f_cells_subtype, drop=F]) > mean(mtx_expr[syms_centroid_subtype, !f_cells_subtype, drop=F])
			} else {
				f_syms_centroid_more_expr <- rowMeans(mtx_expr[syms_centroid_subtype, f_cells_subtype, drop=F]) > rowMeans(mtx_expr[syms_centroid_subtype, !f_cells_subtype, drop=F])
			}
			if (any(!f_cells_subtype)) {
				cat(sprintf("\t\t%s genes higher expressed than other subtypes: %s\n", subtype, paste(syms_centroid_subtype[f_syms_centroid_more_expr], collapse=", ")))
			}
		} # for

		# usage: rna$breast_cancer_epi_type
	},


	"brca_pam50_tcga_centroids"={

		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			#idx_cells <- 1:ncol(rna)
			vec_cell_types <- get_vec_cell_types(rna, args)
			idx_cells <- grepl(pattern_epi, vec_cell_types)
		} # if

		filename_rds_centroids <- "reference/centroids_for_singler/tcga_brca_molecular_subtype_centroids.rds"
		cat(sprintf("\t\treadRDS('%s')\n", filename_rds_centroids))
		se_centroid <- readRDS(filename_rds_centroids)
		se_centroid.sc <- SingleR(test=rna.sce[,idx_cells], assay.type.test="logcounts", assay.type.ref="logcounts", ref=se_centroid, labels=se_centroid$label.main)

		rna$brac_pam50 <- NA
		rna@meta.data[idx_cells, "brca_pam50"] <- se_centroid.sc$pruned.labels

		# usage: rna$brca_pam50
	},

	"brca_pam50_cor"={

		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			#idx_cells <- 1:ncol(rna)
			vec_cell_types <- get_vec_cell_types(rna, args)
			idx_cells <- grep(pattern_epi, vec_cell_types)
		} # if

		mtx <- GetAssayData(object = rna[,idx_cells],
				 assay = NULL,
				 slot = "data")

		suppressPackageStartupMessages(library(genefu))
		data(pam50)

		# df_pam50
		df_pam50 <- as.data.frame(pam50$centroids)
		cols <- c("Basal", "Her2E", "LumA", "LumB", "Normal-like")
		colnames(df_pam50) <- cols

		# df_map
		# columns: probe, probe.centroids, EntrezGene.ID
		df_map <- pam50$centroids.map
                suppressPackageStartupMessages(library(AnnotationDbi))
                suppressPackageStartupMessages(library(org.Hs.eg.db))
                suppressMessages(suppressWarnings(
                        map.id.alias <- AnnotationDbi::select(org.Hs.eg.db, keys=as.character(df_map$EntrezGene.ID), columns=c("SYMBOL"), keytype="ENTREZID") ))
		idx <- match(df_map$EntrezGene.ID, map.id.alias$ENTREZID)
		f <- !is.na(idx); idx <- idx[f]
		df_map[f,"probe"] <- map.id.alias[idx,"SYMBOL"]
		rownames(df_map) <- df_map$probe

		idx_probes <- match(df_map$probe, rownames(mtx))
		f_probes <- !is.na(idx_probes); idx_probes <- idx_probes[f_probes]
		# update rownames of df_pam50
		rownames(df_pam50) <- df_map[,"probe"]

     		cols <- list_cancer_type_specific_info$cols
		cols <- paste0(cols, ".", 1:length(cols))
		rna@meta.data[idx_cells, cols] <- NA

		method_subtyping <- "cor"
		switch(method_subtyping,
			"cor"={
				# using all tumor cells
				mtx_cor <- cor((as.matrix(mtx[idx_probes,])), (as.matrix(df_pam50)), method="kendal") # pearson, spearman, kendal
				# update scores
				rna@meta.data[idx_cells, cols] <- mtx_cor

			},
			{}
		) # switch

		# usage: Basal.1, Her2E.2, LumA.3, LumB.4, Normal-like.5
	},

	"brca_pam50_genefu"={

		if ("CNV.Pos" %in% colnames(rna@meta.data)) {
 			cat(sprintf("\tusing CNV positive cells\n"))
			idx_cells <- grep("11", rna$CNV.Pos)
		} else {
			idx_cells <- 1:ncol(rna)
		} # if

		mtx <- GetAssayData(object = rna[,idx_cells],
				 assay = NULL,
				 slot = "data")

		suppressPackageStartupMessages(library(Matrix))
		suppressPackageStartupMessages(library(genefu))
		data(pam50)

		# df_map
		# columns: probe, probe.centroids, EntrezGene.ID
		df_map <- pam50$centroids.map
                suppressPackageStartupMessages(library(AnnotationDbi))
                suppressPackageStartupMessages(library(org.Hs.eg.db))
                suppressMessages(suppressWarnings(
                        map.id.alias <- AnnotationDbi::select(org.Hs.eg.db, keys=as.character(df_map$EntrezGene.ID), columns=c("SYMBOL"), keytype="ENTREZID") ))
		idx <- match(df_map$EntrezGene.ID, map.id.alias$ENTREZID)
		f <- !is.na(idx); idx <- idx[f]
		df_map[f,"probe"] <- map.id.alias[idx,"SYMBOL"]
		rownames(df_map) <- df_map$probe

		idx_probes <- match(df_map$probe, rownames(mtx))
		f_probes <- !is.na(idx_probes); idx_probes <- idx_probes[f_probes]

     		cols <- list_cancer_type_specific_info$cols
		cols <- paste0(cols, ".", 1:length(cols))
		rna@meta.data[idx_cells, cols] <- NA

		method_subtyping <- "pam50"
		switch(method_subtyping,

			"pam50"={
				# using all tumor cells
				# https://github.com/bhklab/genefu/blob/master/R/molecular.subtyping.R
				PAM50Preds <- intrinsic.cluster.predict(
					sbt.model = pam50,
					data = t(mtx[idx_probes, ]),
					annot = df_map[f_probes,], # Matrix of annotations with at least one column named "EntrezGene.ID" (for ssp, scm, AIMS, and claudinLow models) or "Gene.Symbol" (for the intClust model), dimnames being properly defined.
					do.mapping = TRUE, # TRUE if the mapping through Entrez Gene ids must be performed (in case of ambiguities, the most variant probe is kept for each gene), FALSE otherwise.
					verbose = FALSE
				) # intrinsic.cluster.predict

				# update scores
				rna@meta.data[idx_cells, cols] <- PAM50Preds$subtype.proba 
			},

			"pam50_chunk"={
				for (x in bit::chunk(from=1, to=ncol(mtx), by=1000)) {
					idx_cells_chunk <- x[1]:x[2]
					suppressMessages(suppressWarnings(
						PAM50Preds <- intrinsic.cluster.predict(sbt.model = pam50, data = t(mtx[idx_probes, idx_cells_chunk]), annot = df_map[f_probes,], do.mapping=TRUE)
					))

					# update scores
					rna@meta.data[idx_cells[idx_cells_chunk], cols] <- PAM50Preds$subtype.proba 

				} # for
			},

			{}
		) # switch

		# usage: Basal.1, Her2E.2, LumA.3, LumB.4, Normal-like.5

	},

	{}
  ) # swtich


  # restore default assay
  DefaultAssay(rna) <- default_assay

  rna

} # identify_cell_types






















# update_cluster.ids
# determine mot common cell type with enrichment score
#
# comment:
# called by get_most_common_cell_type_for_each_cluster()
update_cluster.ids <- function(rna, cluster.ids, cell_type) {

  switch(cell_type,
    "B-cells"={
      # assess B cell enrichment to potentially rename clusters
      # b_cells.8: see identify_cell_types <- function(rna, "panglaodb", args)
      col_feature <- "b_cells.8"
      # rename cluster if median enrichment score is greater than 0.225
      th_median_enrichment_score <- 0.225
      col_cell_type <- "b.cell"
     },
    "Mast cells"={
      # Mast cells were not defined by HTAPP, HPCA, or BED.
      # assess mast cell enrichment to potentially rename clusters
      # mast_cells.1 signature: see identify_cell_types <- function(rna, "cell_type_specific_gene_signatures", args)
      col_feature <- "mast_cells.1"
      if (!col_feature %in% colnames(rna@meta.data)) {
		return(cluster.ids)
      }
      # rename cluster if median enrichment score is greater than 0.225
      th_median_enrichment_score <- 0.225
      col_cell_type <- "mast.cell"
    },
    {}
  ) # switch
    
  vln.df <- VlnPlot(rna, features = col_feature)
  df_summary <- describeBy(vln.df$data[, col_feature], vln.df$data$ident, mat = TRUE)
  df_summary <- dplyr::filter(df_summary, median > th_median_enrichment_score)

  if( nrow(df_summary) > 0 ) {

      # assign
      # rna@meta.data[,col_cell_type] <- ifelse(fac_cluster %in% as.character(df_summary$group1), TRUE, FALSE)
      for (i in 1:nrow(df_summary)) {
        cluster.ids[[df_summary$group1[i]]] <- cell_type
      } # for

  } # if

  cluster.ids

} # update_cluster.ids












# get_most_common_cell_type_for_each_cluster
# assign most common cell type for each cluster
#
# input:
#    rna: seurat obj
#    args:
#      str_update_cell_types: e.g. "B-cells,Mast cells"
#    col_seurat_cluster: (e.g. sprintf("RNA_snn_res.%g", args$seurat_resolution) or sprintf("RNA_harmony_th.%s", paste(args$harmony_theta, collapse="_"))
#    col_cell_types: (default="SingleR") column of cell types
#      switch(args$method_to_identify_cell_types,
#		"singler_htapp_toolbox"={ obj[[i]]$SingleR <- obj[[i]]$SingleR.HTAPP_toolbox },
#		"singler_hpca_celldex"={ obj[[i]]$SingleR <- obj[[i]]$SingleR.HPCA },
#		"singler_blueprint_encode"={ obj[[i]]$SingleR <- obj[[i]]$SingleR.BED }
#      )
#      see make_sc-rna-seq_merged_seurat_obj.R
#
# output:
#    cell.type: string vector
#
# comment:
# call by make_sc-rna-seq_merged_seurat_obj.R
#
# usage:
# rna$cell.type <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_cluster); Idents(rna) <- rna$cell.type
# rna$cell.type.harmony <- get_most_common_cell_type_for_each_cluster(rna, args, str_column_of_meta_data_harmony)
get_most_common_cell_type_for_each_cluster <- function(rna, args, col_seurat_cluster, col_cell_types="cell.type") {

  cat(sprintf("\tget_most_common_cell_type_for_each_cluster(%s)\n", col_seurat_cluster))

  Idents(rna) <- col_seurat_cluster

  fac_cluster <- rna[[col_seurat_cluster]][,1]
  names(fac_cluster) <- rownames(rna[[col_seurat_cluster]])

  #cells <- rna@meta.data %>% dplyr::group_by_at(col_seurat_cluster) %>% dplyr::count(cell.type)
  cells <- rna@meta.data %>% dplyr::group_by_at(col_seurat_cluster) %>% dplyr::count(across(all_of(col_cell_types)))
  cluster.ids <- rep("fill", length(levels(Idents(rna))))
  names(cluster.ids) <- levels(Idents(rna))

  for ( i in factor(cells[[col_seurat_cluster]]) ) {
    cells.sub <- cells %>% dplyr::filter(get(col_seurat_cluster) == i) %>% dplyr::arrange(desc(n))
    #cluster.ids[[i]] <- cells.sub$cell.type[1]
    cluster.ids[[i]] <- cells.sub[[1, col_cell_types]]
  } # for


  if (!is.null(args$str_update_cell_types) && (nchar(args$str_update_cell_types) > 0)) {
    cell_types <- strsplit(args$str_update_cell_types, ",")[[1]]
    for (cell_type in cell_types) {
	# determine mot common cell type with enrichment score (e.g. B-cells, Mast cells)
	# alternative approach by H. Kim:
	# B-cells: SingleR.BED has reference centroid vector for B-cells
	# Mast cells: try to annotate mast cells with mast cell specific gene signature, but do not touch epithelial cells: make_sc-rna-seq_seurat_obj.R: identify_cell_types_with_multiple_methods(): rna <- update_cell_types(rna, "mast_cells", args).
	cluster.ids <- update_cluster.ids(rna, cluster.ids, cell_type)
    } # for
  } # if


  # cell.type
  df.tmp <- as.data.frame(cluster.ids)
  levels(Idents(rna)) <- df.tmp$cluster.ids

  #rna$cell.type <- Idents(rna)
  #rna$cell.type <- paste0(rna$RNA_snn_res.0.7, "-", rna$cell.type)

  cell.type <- paste0(fac_cluster, "-", Idents(rna))

  cell.type


} # get_most_common_cell_type_for_each_cluster











# expand_epi_cells
#
# input:
#   idx_epi:
#   vec_cell_types:
#
# output:
#   list_out$idx_epi: more epithelial cells
#   list_out$vec_cell_types: vec_cell_types[idx_epi] <- "Epithelial cells"
#
# usage:
# list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
# idx_epi <- list_out$idx_epi
# vec_cell_types <- list_out$vec_cell_types
expand_epi_cells <- function(rna, idx_epi, vec_cell_types, args) {

  col_cell_types <- get_column_name_for_cell_types(rna, args)

  # expand epithelial cells
  switch(tolower(args$cancer_subtype),
	"mbc"={
		cat(sprintf("\t\targs$cancer_subtype: %s\n", args$cancer_subtype))
		cat(sprintf("\t\tusing SingleR.HTAPP_toolbox as well as %s\n", col_cell_types))
		if (!"SingleR.HTAPP_toolbox" %in% colnames(rna@meta.data)) {
			rna <- identify_cell_types(rna, "singler_htapp_toolbox", args)
		}
		idx_epi_htapp_mbc <- grep(pattern_epi, rna$SingleR.HTAPP_toolbox)
		cat(sprintf("\t\tn_epi_cells_htapp_mbc: %d\n", length(idx_epi_htapp_mbc)))
		idx_epi <- union(idx_epi, idx_epi_htapp_mbc)

	},
	"test"={
		idx_epi_e <- grep("^epi", rna$epi_krt_epcam)
		cat(sprintf("\t\tn_epi_cells_epi_krt_epcam: %d\n", length(idx_epi_e)))
		idx_epi <- union(idx_epi, idx_epi_e)

		if ("celltype_garnett.major" %in% colnames(rna@meta.data)) {
			idx_epi_g <- grep(pattern_epi, rna$celltype_garnett.major)
			cat(sprintf("\t\tn_epi_celltype_garnett.major: %d\n", length(idx_epi_g)))
			idx_epi <- union(idx_epi, idx_epi_g)
		}
		if ("celltype_td.major" %in% colnames(rna@meta.data)) {
			idx_epi_t <- grep(pattern_epi, rna$celltype_td.major)
			cat(sprintf("\t\tn_epi_celltype_td.major: %d\n", length(idx_epi_t)))
			idx_epi <- union(idx_epi, idx_epi_t)
		}
		if ("celltype_markercount" %in% colnames(rna@meta.data)) {
			idx_epi_m <- grep(pattern_epi, rna$celltype_markercount)
			cat(sprintf("\t\tn_epi_celltype_markercount: %d\n", length(idx_epi_m)))
			idx_epi <- union(idx_epi, idx_epi_m)
		}
	},
	{}
  ) # switch

  method_to_update_epi <- "all"
  switch(method_to_update_epi,

	"all"={
		cat(sprintf("\t\tn_epi_cells: %d\n", length(idx_epi)))
		vec_cell_types[idx_epi] <- "Epithelial cells"
	},

	"singler_miss"={
		# try to fix potential misclassifications
		switch(args$method_to_identify_cell_types,
			"singler_htapp_toolbox"={ 
			},
			"singler_hpca_celldex"={
			},
			"singler_blueprint_encode"={
				pattern_cell_types_potentially_misclassified <- "(^|-)Adipocytes|Chondrocytes|(^|-)Fibroblasts|(^|-)Pericytes"
				idx_potentially_misclassified <- grep(pattern_cell_types_potentially_misclassified, vec_cell_types)
				idx.epi.miss <- intersect(idx_potentially_misclassified, idx_epi)
				cat(sprintf("\t\tn_epi_new: %d\n", length(idx.epi.miss)))
				vec_cell_types[idx.epi.miss] <- "Epithelial cells"
			},
			{}
		) # switch
	},

	"within_obs_clusters"={
		if ((args$method_to_identify_reference_clusters == "none") || (args$type_infercnv_argset == "none")) {
			# do nothing
		} else {
			reference.clusters <- identify_reference_clusters(rna, args)
			clusters <- levels(Idents(rna))
			for (cluster in clusters) {
				# consider only obsrvation clusters
				if (cluster %in% reference.clusters) next
				cat(sprintf("\t\tobservation cluster: %s\n", cluster))

				idx_cell_cluster <- which(Idents(rna) == cluster)
				tb_cluster <- sort(table(vec_cell_types[idx_cell_cluster]), decreasing=T)
				cat(sprintf("\t\t\tn_cells: %d\n", length(idx_cell_cluster)))
				if (grepl(pattern_epi, names(tb_cluster)[1])) {
					# epithelial cell is the dominant cell type in this cluster where there may be misclassification.
					f.epi.cluster <- grepl(pattern_epi, vec_cell_types[idx_cell_cluster])
					idx.new.epi <- intersect(idx_cell_cluster[!f.epi.cluster], idx_epi)
					
					cat(sprintf("\t\t\tn_epi_cells: %d, n_new_epi_cells: %d\n", length(which(f.epi.cluster)), length(idx.new.epi)))
					vec_cell_types[idx.new.epi] <- "Epithelial cells"
				} # if
			} # for
		} # if
	},
	{}
  ) # switch


  list_out <- list()
  list_out$idx_epi <- idx_epi
  list_out$vec_cell_types <- vec_cell_types

  list_out

} # expand_epi_cells















# update_epi_cell_types
#
# input:
#   method: {epi_normal_tumor, normal_epi_type, tumor_epi_type}
#
# output:
#   when method="epi_normal_tumor"
#     rna$epi_normal_tumor: {nonepi, epi, normal_epi}
#   when method="normal_epi_type"
#     rna$normal_epi_type: {non-Ep, Ep, LEp_prog, BEp, LEp}
#     rna$normal_epi_type: {non-Ep, Ep, LEp_prog, BEp, BEp_myo, LEp_secretory, LEp_hormone}
#   when method="tumor_epi_type"
#     rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB}
#
update_epi_cell_types <- function(rna, method, args) {

  cat(sprintf("\tupdate_epi_cell_types for %s\n", method))
  cat(sprintf("\tupdate_epi_cell_types for %s\n", method), file=stderr())

  #list_cancer_type_specific_info <- get_cancer_type_specific_info(args)
  list_cancer_type_specific_info <- args$list_cancer_type_specific_info
  vec_cell_types <- get_vec_cell_types(rna, args)

  switch(method,

	"epi_normal_tumor"={

		# rna$epi_normal_tumor={"nonepi", "epi", "normal_epi"}
		# set rna$epi_normal_tumor <- "normal_epi" when normal epi. was called with high confidence with gene expression without copy number info.
		# currently, there is no valid method to call normal epi. without copy number info.
		rna$epi_normal_tumor <- "nonepi"
		th_score_epi <- 1.0
		th_score_normal <- 1.0

		# epithelialcell_signature.4
		#score_epi <- rna$epithelialcell_signature.4
		#f_epi <- (score_epi > th_score_epi)

		# using SingleR
		f_epi <- grepl(pattern_epi, vec_cell_types)

		idx_epi <- which(f_epi)

		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		# rna$epi_normal_tumor
		rna$epi_normal_tumor[idx_epi] <- "epi"

		type_calling_normal_epi <- "unknown"
		switch(type_calling_normal_epi,
			"test"={
				# this is not a valid method, but a test.
				# calling normal epithelial cells with normal_luminalprog.1, normal_basal.2, and normal_matureluminal.3 > 1.0 tends to miss many normal epithelial cells since some normal epithelial cells > 0.5. but, the objective here is to call normal epithelail cells with higher confidence so that they can be used for infercnv reference cells.
				# normal_luminalprog.1
				score_normal <- rna$normal_luminalprog.1
				idx_normal <- which(f_epi & (score_normal > th_score_normal))
				rna$epi_normal_tumor[idx_normal] <- "normal_epi"

				# normal_basal.2
				score_normal <- rna$normal_basal.2
				idx_normal <- which(f_epi & (score_normal > th_score_normal))
				rna$epi_normal_tumor[idx_normal] <- "normal_epi"

				# normal_matureluminal.3
				score_normal <- rna$normal_matureluminal.3
				idx_normal <- which(f_epi & (score_normal > th_score_normal))
				rna$epi_normal_tumor[idx_normal] <- "normal_epi"
			},
			{}
		) # switch

	},

	"normal_epi_type"={

		rna$normal_epi_type <- "non-Ep"

		# using epithelialcell_signature.4 and SingleR
		#score_epi <- rna$epithelialcell_signature.4
		#f_epi <- (score_epi > 1.0) | (grepl(pattern_epi, vec_cell_types))
		#f_epi <- grepl("Epithelial", vec_cell_types)

		# using SingleR
		# Keratinocytes are epithelial cells that form the superficial layer of the skin  https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/keratinocytes
		f_epi <- grepl(pattern_epi, vec_cell_types)

		idx_epi <- which(f_epi)
		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		# rna$normal_epi_type
		rna$normal_epi_type[idx_epi] <- "Ep"


		switch(list_cancer_type_specific_info$type_normal_epi_type,
			"peroulab"={
				types <- c("LEp_prog", "BEp", "LEp")
				mtx <- cbind(rna$normal_luminalprog.1, rna$normal_basal.2, rna$normal_matureluminal.3)
				mtx <- mtx[idx_epi,]
				# ties.method = c("random", "first", "last")
				idx_col <- max.col(mtx, ties.method="last")

				for (i in 1:length(types)) {
		  			idx_normal <- idx_epi[which(idx_col == i)]
		  			rna$normal_epi_type[idx_normal] <- types[i]
				}
			},
			"swarbricklab"={
				rna$normal_epi_type[idx_epi] <- rna$breast_normal_epi_type[idx_epi]
			},
			"kessenbrocklab"={
				rna$normal_epi_type[idx_epi] <- rna$breast_normal_epi_type[idx_epi]
			},
			{}
		) # switch

	},

	"tumor_epi_type"={

		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}

		rna$tumor_epi_type <- "non-Ep"
		# update tumor_epi_type for all epithelaial cells
		# normal epithelial cells will be determind later by update_cell_types_after_infercnv()
		f_epi <- grepl(pattern_epi, vec_cell_types)

		idx_epi <- which(f_epi)
		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		# rna$tumor_epi_type
		rna$tumor_epi_type[idx_epi] <- "Ep"


		switch(args$method_to_identify_subtypes,
			"scbcsubtype"={
				# swarbricklab_breast_cancer_epi
				rna$tumor_epi_type[idx_epi] <- rna$breast_cancer_epi_type[idx_epi]
				# assign "Ep" when subtype was not available.
				f.na <- is.na(rna$tumor_epi_type[idx_epi])
				rna$tumor_epi_type[idx_epi[f.na]] <- "Ep"
			},
			"brca_pam50_tcga_centroids"={
				rna$tumor_epi_type[idx_epi] <- rna$brca_pam50[idx_epi]
			},
			{
				cols <- list_cancer_type_specific_info$cols
				subtypes <- list_cancer_type_specific_info$subtypes
				mtx <- NULL
				cols <- paste0(cols, ".", 1:length(cols))
				if (!all(cols %in% colnames(rna@meta.data))) {
					return(rna)
				}
				for (col in cols) {
				  mtx <- cbind(mtx, rna@meta.data[,col])
				}

				mtx <- mtx[idx_epi,]
				# ties.method = c("random", "first", "last")
				idx_col <- max.col(mtx, ties.method="last")

				for (i in 1:length(subtypes)) {
				  idx_epi_cancer <- idx_epi[which(idx_col == i)]
		  		  rna$tumor_epi_type[idx_epi_cancer] <- subtypes[i]
				} # for
			}
		) # switch

	},

	{}
  ) # switch

  rna

} # update_epi_cell_types























# update_cell_types
#
# input:
#   args$method_to_identify_cell_types
#   args$method_to_update_cell_types: 
#	normal_epithelial_cell_types: for normal samples
#	cancer_epithelial_cell_types: for cancer samples
#
# output:
#   when args$method_to_identify_cell_types == "singler_blueprint_encode"
#       update rna$SingleR.BED
#
# usage:
# rna <- update_cell_types(rna, args$method_to_update_cell_types, args)
update_cell_types <- function(rna, method, args) {

  cat(sprintf("\tupdate_cell_types for %s\n", method))
  cat(sprintf("\tupdate_cell_types for %s\n", method), file=stderr())

  #list_cancer_type_specific_info <- get_cancer_type_specific_info(args)
  list_cancer_type_specific_info <- args$list_cancer_type_specific_info
  col_cell_types <- get_column_name_for_cell_types(rna, args)
  vec_cell_types <- get_vec_cell_types(rna, args)

  # update vec_cell_types
  switch(method,

	"mast_cells"={

		if ("mast_cells.1" %in% colnames(rna@meta.data)) {
  			# try to annotate mast cells with 3 mast gene signature, but do not touch epithelial cells
			score_mast <- rna$mast_cells.1
			f_mast <- (score_mast > 1.0) & !grepl(pattern_epi, vec_cell_types)
			idx_mast <- which(f_mast)
			cat(sprintf("\t\t# of mast cells: %d\n", length(idx_mast)))
			vec_cell_types[idx_mast] <- "Mast cells"
		}

	},

	"epithelial_cell_types"={

		# epithelial cells with multiple cell type callers

		# using epithelialcell_signature.4 and SingleR
		#score_epi <- rna$epithelialcell_signature.4
		#f_epi <- (score_epi > 1.0) | (grepl(pattern_epi, vec_cell_types))
		#f_epi <- grepl("Epithelial", vec_cell_types)

		# using SingleR
		# Keratinocytes are epithelial cells that form the superficial layer of the skin  https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/keratinocytes
		f_epi <- grepl(pattern_epi, vec_cell_types)
		idx_epi <- which(f_epi)
		vec_cell_types[idx_epi] <- "Epithelial cells"
		cat(sprintf("\t\tn_epi_%s: %d\n", col_cell_types, length(idx_epi)))


	},
		
	"normal_epithelial_cell_types"={

		f_epi <- grepl(pattern_epi, vec_cell_types)
		idx_epi <- which(f_epi)

		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		vec_cell_types <- list_out$vec_cell_types

		switch(list_cancer_type_specific_info$type_normal_epi_type,
			"peroulab"={
				types <- c("LEp_prog", "BEp", "LEp")
				mtx <- cbind(rna$normal_luminalprog.1, rna$normal_basal.2, rna$normal_matureluminal.3)
				mtx <- mtx[idx_epi,]
				# ties.method = c("random", "first", "last")
				idx_col <- max.col(mtx, ties.method="last")

				for (i in 1:length(types)) {
		  			idx_epi_normal <- idx_epi[which(idx_col == i)]
		  			vec_cell_types[idx_epi_normal] <- types[i]
				}
			},
			"swarbricklab"={
				vec_cell_types[idx_epi] <- rna$breast_normal_epi_type[idx_epi]
			},
			"kessenbrocklab"={
				vec_cell_types[idx_epi] <- rna$breast_normal_epi_type[idx_epi]
			},
			{}
		) # switch

	},

	"cancer_epithelial_cell_types"={

		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}

		# update tumor_epi_type for all epithelaial cells
		# normal epithelial cells will be determind later by update_cell_types_after_infercnv()
		f_epi <- grepl(pattern_epi, vec_cell_types)
		idx_epi <- which(f_epi)

		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		vec_cell_types <- list_out$vec_cell_types

		switch(args$method_to_identify_subtypes,

			"scbcsubtype"={
				# swarbricklab_breast_cancer_epi
				if ("breast_cancer_epi_type" %in% colnames(rna@meta.data)) {
					vec_cell_types[idx_epi] <- rna$breast_cancer_epi_type[idx_epi]
					f.na <- is.na(vec_cell_types[idx_epi])
					vec_cell_types[idx_epi[f.na]] <- "Epithelial cells"
				}
			},
			"brca_pam50_tcga_centroids"={
				vec_cell_types[idx_epi] <- rna$brca_pam50[idx_epi]
			},
			{

				cols <- list_cancer_type_specific_info$cols
				subtypes <- list_cancer_type_specific_info$subtypes
				mtx <- NULL
				cols <- paste0(cols, ".", 1:length(cols))
				if (!all(cols %in% colnames(rna@meta.data))) {
					return(rna)
				}
				for (col in cols) {
				  mtx <- cbind(mtx, rna@meta.data[,col])
				} # for

				mtx <- mtx[idx_epi,]
				# ties.method = c("random", "first", "last")
				idx_col <- max.col(mtx, ties.method="last")

				for (i in 1:length(subtypes)) {
				  idx_epi_cancer <- idx_epi[which(idx_col == i)]
		 		  vec_cell_types[idx_epi_cancer] <- subtypes[i]
				} # for
			}
		) # switch

	},

	"cancer_normal_epithelial_cell_types"={

		# update vec_cell_types when tumor/normal cells can be predicted with gene expression without copy number info
		# currently, there is no valid method to call tumor/normal cells without copy number info.
		rna$epi_normal_tumor <- "nonepi"
		df_sig <- list_cancer_type_specific_info$df_sig
		if (is.null(df_sig)) {
			return(rna)
		}


		# update tumor_epi_type for all epithelaial cells
		# normal epithelial cells will be determind later by update_cell_types_after_infercnv()
		f_epi <- grepl(pattern_epi, vec_cell_types)
		idx_epi <- which(f_epi)

		# expand epithelial cells
		list_out <- expand_epi_cells(rna, idx_epi, vec_cell_types, args)
		idx_epi <- list_out$idx_epi
		vec_cell_types <- list_out$vec_cell_types

		cols <- list_cancer_type_specific_info$cols
		subtypes <- list_cancer_type_specific_info$subtypes
		mtx <- NULL
		cols <- paste0(cols, ".", 1:length(cols))
		if (!all(cols %in% colnames(rna@meta.data))) {
			return(rna)
		}
		for (col in cols) {
		  mtx <- cbind(mtx, rna@meta.data[,col])
		} # for

		type_calling_tumor_normal_cells <- "unknown"
		switch(list_cancer_type_specific_info$type_normal_epi_type,
			"test"={
				# this is not a valid method, but a test.
				subtypes <- c(subtypes, "LEp_prog", "BEp", "LEp")
				mtx <- cbind(mtx, rna$normal_luminalprog.1, rna$normal_basal.2, rna$normal_matureluminal.3)
				mtx <- mtx[idx_epi,]
				# ties.method = c("random", "first", "last")
				idx_col <- max.col(mtx, ties.method="last")

				for (i in 1:length(subtypes)) {
		  			idx_epi_type <- idx_epi[which(idx_col == i)]
		  			vec_cell_types[idx_epi_type] <- subtypes[i]
				} # for
			},
			{}
		) # switch
	},
	{}

  ) # swtich

  # update rna${SingleR.HTAPP_toolbox, SingleR.HPCA, SingleR.BED} depending on args$method_to_identify_cell_types
  #rna <- set_vec_cell_types(rna, vec_cell_types, args)

  # update rna$cell.type
  rna$cell.type <- vec_cell_types

  rna

} # update_cell_types
























# identify_cell_types_with_multiple_methods
#
# input:
#   rna: seurat obj
#   args:
#   f_update_cancer_subtype_after_infercnv:
#
# output:
#   rna$epi_krt_epcam
#   rna${normal_luminalprog.1, normal_basal.2, normal_matureluminal.3, epithelialcell_signature.4}
#   rna${Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4}
#   rna$epi_normal_tumor
#   rna$normal_epi_type
#   rna$tumor_epi_type
#   update rna$SingleR.BED when args$method_to_identify_cell_types == "singler_blueprint_encode"
#
# usage:
# find_clusters_after_doublet_removal.R:  rna <- identify_cell_types_with_multiple_methods(rna, args)
#
identify_cell_types_with_multiple_methods <- function(rna, args, f_update_cancer_subtype_after_infercnv=FALSE) {

  cat(sprintf("\tidentify_cell_types_with_multiple_methods(rna, args, f_update_cancer_subtype_after_infercnv=%s)\n", f_update_cancer_subtype_after_infercnv))

  if (!f_update_cancer_subtype_after_infercnv) {

	# rna <- identify_cell_types(rna, "panglaodb", args)

	# args$method_to_identify_cell_type: {"singler_htapp_toolbox", "singler_hpca_cellidx", ["singler_blueprint_encode"]}
	rna <- identify_cell_types(rna, args$method_to_identify_cell_types, args)

	# rna$epi_krt_epcam
	rna <- identify_cell_types(rna, "epi_krt_epcam", args)

	# rna${normal_luminalprog.1, normal_basal.2, normal_matureluminal.3, epithelialcell_signature.4} when list_cancer_type_specific_info$type_normal_epi_type="peroulab"
	# rna$breast_normal_epi_type when list_cancer_type_specific_info$type_normal_epi_type="kessenbrocklab"
	rna <- identify_cell_types(rna, "normal_epithelial_cells", args)

	switch(tolower(args$cancer_subtype),
		"test"={
			# rna${celltype_garnett.major, celltype_garnett.major.ext, celltype_garnett.mincor, celltype_garnett.minor.ext, celltype_garnett.subset, celltype_garnett.subset.ext}
			rna <- identify_cell_types(rna, "celltype_garnett", args)

			# rna${celltype_td.major, celltype_td.major.score.max, celltype_td.minor, celltype_td.minor.score.max, celltype_td.subset, celltype_td.subset.score.max}
			#rna <- identify_cell_types(rna, "celltype_td", args)

			# rna$celltype_markercount
			rna <- identify_cell_types(rna, "celltype_markercount", args)
		},
		{}
	) # switch

  } # if





  # subtype
  switch(args$cancer_type_standard,
	"bc"={

		switch(args$method_to_identify_subtypes,

			"none"={
			},

			"brca_cell_atlas_scsubtyper"={
				# rna${Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4}
				if (f_update_cancer_subtype_after_infercnv) {
					# performa subtyping after infercnv
					rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
				} else {
					# performa subtyping after infercnv
					args$f_update_cancer_subtype_after_infercnv <<- TRUE
				}
			},

			"brca_cell_atlas_scsubtyper_with_signatures"={
				# rna${Basal_SC.1, Her2E_SC.2, LumA_SC.3, LumB_SC.4}
				if (f_update_cancer_subtype_after_infercnv) {
					# performa subtyping after infercnv
					rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
				} else {
					# performa subtyping after infercnv
					args$f_update_cancer_subtype_after_infercnv <<- TRUE
				}
			},

			"scbcsubtype"={
				# swarbricklab_breast_cancer_epi
				# rna$breast_cancer_epi_type
				if (f_update_cancer_subtype_after_infercnv) {
					# performa subtyping after infercnv
					rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
				} else {
					# performa subtyping after infercnv
					args$f_update_cancer_subtype_after_infercnv <<- TRUE
				}
			},

			"brca_pam50_tcga_centroids"={
				# rna$brca_pam50
				rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
			},

			"brca_pam50_cor"={
				# brca_pam50_cor
				# rna${Basal.1, Her2E.2, LumA.3, LumB.4, Normal-like.5}
				rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
			},

			"brca_pam50_genefu"={
				# brca_pam50_genefu
				# rna${Basal.1, Her2E.2, LumA.3, LumB.4, Normal-like.5}
				rna <- identify_cell_types(rna, args$method_to_identify_subtypes, args)
			},
			{}
		) # switch
	},
	{}
  ) # switch



  ### update epithelial cell types

  if (!f_update_cancer_subtype_after_infercnv) {
	# rna$epi_normal_tumor: {nonepi, epi, normal_epi}
	rna <- update_epi_cell_types(rna, "epi_normal_tumor", args)
	# rna$normal_epi_type: {non-Ep, Ep, LEp_prog, BEp, LEp}
	# rna$normal_epi_type: {non-Ep, Ep, LEp_prog, BEp, BEp_myo, LEp_secretory, LEp_hormone}
	rna <- update_epi_cell_types(rna, "normal_epi_type", args)
  } # if


  switch(args$cancer_type_standard,
	"bc"={
		switch(args$method_to_identify_subtypes,
			"none"={
			},
			"brca_cell_atlas_scsubtyper"={
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB}
				if (f_update_cancer_subtype_after_infercnv) {
					rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
				}
			},
			"brca_cell_atlas_scsubtyper_with_signatures"={
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB}
				if (f_update_cancer_subtype_after_infercnv) {
					rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
				}
			},
			"scbcsubtype"={
				# swarbricklab_breast_cancer_epi
				# rna$breast_cancer_epi_type
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB}
				if (f_update_cancer_subtype_after_infercnv) {
					rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
				}
			},
			"brca_pam50_tcga_centroids"={
				# rna$brca_pam50
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, CLow, Her2E, LumA, LumB, Normal-like, Normal}
				rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
			},
			"brca_pam50_cor"={
				# brca_pam50_cor
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB, Normal-like}
				rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
			},
			"brca_pam50_genefu"={
				# brca_pam50_genefu
				# rna$tumor_epi_type: {non-Ep, Ep, Basal, Her2E, LumA, LumB, Normal-like}
				rna <- update_epi_cell_types(rna, "tumor_epi_type", args)
			},
			{}
		) # switch
	},
	{}
  ) # switch



  ### update cell types
  # update rna$SingleR.XXX
  # update rna$SingleR.BED when args$method_to_identify_cell_types == "singler_blueprint_encode"


  # Mast cells were not defined by HTAPP, HPCA, or BED.
  # try to annotate mast cells with 3 mast gene signature, but do not touch epithelial cells
  rna <- update_cell_types(rna, "mast_cells", args)


  # using SingleR
  # Keratinocytes are epithelial cells that form the superficial layer of the skin  https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/keratinocytes
  rna <- update_cell_types(rna, "epithelial_cell_types", args)

  if (args$method_to_update_cell_types != "epithelial_cell_types") {
	# method_to_update_cell_types: {[epithelial_cell_types], normal_epithelial_cell_types, cancer_epithelial_cell_types, cancer_normal_epithelial_cell_types}
	# normal_epithelial_cell_types for normal samples
	# cancer_epithelial_cell_types for tumor samples
	rna <- update_cell_types(rna, args$method_to_update_cell_types, args)
  } # if


  rna


} # identify_cell_types_with_multiple_methods














# update_cell_types_after_infercnv
#
# input:
#   args$f_infercnv_pos_notpos: classify tumor and non-tumor with infercnv, otherwise classify tumor, non-tumor, and unassigned (default)
#
# output:
#   update rna$SingleR.BED when args$method_to_identify_cell_types == "singler_blueprint_encode"
#   rna$epi_type
#     args$method_to_update_cell_types=="normal_epithelial_cell_types": {non-Ep, Ep, LEp_prog, BEp, LEp}
#     args$method_to_update_cell_types=="cancer_epithelial_cell_types": {non-Ep, Ep, Basal, Her2E, LumA, LumB}
#     args$method_to_update_cell_types=="cancer_normal_epithelial_cell_types": {non-Ep, Ep, LEp_prog, BEp, LEp, Basal, Her2E, LumA, LumB}
#
#
# usage:
# rna <- update_cell_types_after_infercnv(rna, args)
update_cell_types_after_infercnv <- function(rna, args) {

  cat(sprintf("\tupdate_cell_types_after_infercnv()\n\n"))

  rna$epi_type <- "non-Ep"

  # col_cell_types: {"SingleR.HTAPP_toolbox", "SingleR.HPCA", "SingleR.BED", ...}
  #col_cell_types <- get_column_name_for_cell_types(rna, args)
  col_cell_types <- "cell.type"

  n <- length(rna$CNV.Pos)
  f_cnv_pos <- rep(NA, n)
  f_cnv_not_pos <- rep(NA, n)
  f_cnv_neg <- rep(NA, n)
  f_cnv_unassigned <- rep(NA, n)

  if ("CNV.Pos" %in% colnames(rna@meta.data)) {

 	cat(sprintf("\tusing CNV positive cells\n"))
	f_cnv_pos <- grepl("11", rna$CNV.Pos)
	f_cnv_not_pos <- grepl("01|10|00", rna$CNV.Pos)
	f_cnv_neg <- grepl("00", rna$CNV.Pos)
	f_cnv_unassigned <- grepl("01|10", rna$CNV.Pos)

	if (!is.null(args$f_update_cancer_subtype_after_infercnv)) {
		# update cell types with CNV.Pos
		rna <- identify_cell_types_with_multiple_methods(rna, args, f_update_cancer_subtype_after_infercnv=TRUE)
	}
	vec_cell_types <- rna@meta.data[,col_cell_types]

  } else {

	vec_cell_types <- rna@meta.data[,col_cell_types]
	f_epi <- grepl(pattern_epi, vec_cell_types)
	f_cnv_pos[f_epi] <- TRUE
        f_cnv_not_pos[f_epi] <- FALSE
        f_cnv_neg[f_epi] <- FALSE
        f_cnv_unassigned[f_epi] <- FALSE

  } # if

  if (args$f_infercnv_pos_notpos) {
	cat(sprintf("\t# of cells of cnv_pos: %d\n", length(which(f_cnv_pos))))
	cat(sprintf("\t# of cells of cnv_not_pos: %d\n", length(which(f_cnv_not_pos))))
  } else {
	cat(sprintf("\t# of cells of cnv_pos: %d\n", length(which(f_cnv_pos))))
	cat(sprintf("\t# of cells of cnv_neg: %d\n", length(which(f_cnv_neg))))
	cat(sprintf("\t# of cells of cnv_unassigned: %d\n", length(which(f_cnv_unassigned))))
  }

  f_epi <- grepl(pattern_epi, vec_cell_types)
  rna$epi_type[f_epi] <- "Ep"

  f_tumor_epi <- grepl(pattern_tumor_epi, vec_cell_types)
  f_normal_epi <- grepl(pattern_normal_epi, vec_cell_types)

  cat(sprintf("\t# of cells of tumor_epi: %d\n", length(which(f_tumor_epi))))
  cat(sprintf("\t# of cells of normal_epi: %d\n", length(which(f_normal_epi))))

  # tumor
  idx_tumor_epi_cnv_pos <- which(f_tumor_epi & f_cnv_pos)
  idx_tumor_epi_cnv_not_pos <- which(f_tumor_epi & f_cnv_not_pos)
  idx_tumor_epi_cnv_neg <- which(f_tumor_epi & f_cnv_neg)
  idx_tumor_epi_cnv_unassigned <- which(f_tumor_epi & f_cnv_unassigned)
  if (args$f_infercnv_pos_notpos) {
	cat(sprintf("\t# of cells of tumor_epi & cnv_pos: %d\n", length(idx_tumor_epi_cnv_pos)))
	cat(sprintf("\t# of cells of tumor_epi & cnv_not_pos: %d\n", length(idx_tumor_epi_cnv_not_pos)))
  } else {
	cat(sprintf("\t# of cells of tumor_epi & cnv_pos: %d\n", length(idx_tumor_epi_cnv_pos)))
	cat(sprintf("\t# of cells of tumor_epi & cnv_neg: %d\n", length(idx_tumor_epi_cnv_neg)))
	cat(sprintf("\t# of cells of tumor_epi & cnv_unassigned: %d\n", length(idx_tumor_epi_cnv_unassigned)))
  }


  # normal
  idx_normal_epi_cnv_pos <- which(f_normal_epi & f_cnv_pos)
  idx_normal_epi_cnv_not_pos <- which(f_normal_epi & f_cnv_not_pos)
  idx_normal_epi_cnv_neg <- which(f_normal_epi & f_cnv_neg)
  idx_normal_epi_cnv_unassigned <- which(f_normal_epi & f_cnv_unassigned)
  if (args$f_infercnv_pos_notpos) {
	cat(sprintf("\t# of cells of normal_epi & cnv_pos: %d\n", length(idx_normal_epi_cnv_pos)))
	cat(sprintf("\t# of cells of normal_epi & cnv_not_pos: %d\n", length(idx_normal_epi_cnv_not_pos)))
  } else {
	cat(sprintf("\t# of cells of normal_epi & cnv_pos: %d\n", length(idx_normal_epi_cnv_pos)))
	cat(sprintf("\t# of cells of normal_epi & cnv_neg: %d\n", length(idx_normal_epi_cnv_neg)))
	cat(sprintf("\t# of cells of normal_epi & cnv_unassigned: %d\n", length(idx_normal_epi_cnv_unassigned)))
  }

  # rna${cycling_score_epi, cycling_score_nonepi, celltype.cycling}
  rna <- identify_cell_types(rna, "celltype_cycling", args)

  switch(args$method_to_update_cell_types,

	"normal_epithelial_cell_types"={
		# all normal samples
		if (args$f_infercnv_pos_notpos) {
  			rna$epi_type[idx_normal_epi_cnv_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_pos]
  			rna$epi_type[idx_normal_epi_cnv_not_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_not_pos]
		} else {
  			rna$epi_type[idx_normal_epi_cnv_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_pos]
  			rna$epi_type[idx_normal_epi_cnv_neg] <- rna$normal_epi_type[idx_normal_epi_cnv_neg]
  			rna$epi_type[idx_normal_epi_cnv_unassigned] <- rna$normal_epi_type[idx_normal_epi_cnv_unassigned]
		}
	},

	"epithelial_cell_types"={
		# set tumor/normal epithelial cells with CNV.Pos
		if (args$f_infercnv_pos_notpos) {
			# Ep/non-Ep --> Epi. Tumor/Epi. Non-tumor
			rna$epi_type[idx_normal_epi_cnv_pos] <- "Epi. Tumor"
			rna$epi_type[idx_normal_epi_cnv_not_pos] <- "Epi. Non-tumor"
		} else {
			# Ep/non-Ep --> Epi. Tumor/Epi. Non-tumor/Epi. Unassigned
			rna$epi_type[idx_normal_epi_cnv_pos] <- "Epi. Tumor"
			rna$epi_type[idx_normal_epi_cnv_neg] <- "Epi. Non-tumor"
			rna$epi_type[idx_normal_epi_cnv_unassigned] <- "Epi. Unassigned"
		}

		if (args$f_infercnv_pos_notpos) {
			# Epithelial cells --> Epi. Tumor/Epi. Non-tumor
			# epithelial cells outside epi. clusters are annotated as non-tumor.
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_normal_epi_cnv_pos, col_cell_types] <- "Epi. Tumor"
			rna@meta.data[idx_normal_epi_cnv_not_pos, col_cell_types] <- "Epi. Non-tumor"
		} else {
			# Epithelial cells --> Epi. Tumor/Epi. Non-tumor/Epi. Unassigned
			# epithelial cells outside epi. clusters are annotated as unassigned.
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Unassigned"
			rna@meta.data[idx_normal_epi_cnv_pos, col_cell_types] <- "Epi. Tumor"
			rna@meta.data[idx_normal_epi_cnv_neg, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_normal_epi_cnv_unassigned, col_cell_types] <- "Epi. Unassigned"
		}
	},

	"cancer_epithelial_cell_types"={
		# all tumor samples
		if (args$f_infercnv_pos_notpos) {
			# Ep/non-Ep --> LumA/LumB/Her2E/Basal/Epi. Tumor/Epi. Non-tumor
			rna$epi_type[idx_tumor_epi_cnv_pos] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna$epi_type[idx_tumor_epi_cnv_not_pos] <- "Epi. Non-tumor" # LumA/LumB/Her2E/Basal/Epi. Tumor
			rna$epi_type[idx_normal_epi_cnv_not_pos] <- "Epi. Non-tumor" # Epi. Non-tumor
		} else {
			# Ep/non-Ep --> LumA/LumB/Her2/Basal/Epi. Tumor/Epi. Non-tumor/Epi. Unassigned
			rna$epi_type[idx_tumor_epi_cnv_pos] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna$epi_type[idx_tumor_epi_cnv_neg] <- "Epi. Non-tumor"
			rna$epi_type[idx_normal_epi_cnv_neg] <- "Epi. Non-tumor"
			rna$epi_type[idx_tumor_epi_cnv_unassigned] <- "Epi. Unassigned"
		}

		if (args$f_infercnv_pos_notpos) {
			# LumA/LumB/Her2E/Basal/Epi. Tumor --> LumA/LumB/Her2/Basal/Epi. Tumor/Epi. Non-tumor
			# epithelial cells outside epi. clusters are annotated as non-tumor
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_tumor_epi_cnv_pos, col_cell_types] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna@meta.data[idx_tumor_epi_cnv_not_pos, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_normal_epi_cnv_not_pos, col_cell_types] <- "Epi. Non-tumor"
		} else {
			# LumA/LumB/Her2E/Basal/Epi. Tumor --> LumA/LumB/Her2/Basal/Epi. Tumor/Epi. Non-tumor/Epi. Unassigned
			# epithelial cells outside epi. clusters are annotated as unassigned.
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Unassigned"
			rna@meta.data[idx_tumor_epi_cnv_pos, col_cell_types] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna@meta.data[idx_tumor_epi_cnv_neg, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_normal_epi_cnv_neg, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_tumor_epi_cnv_unassigned, col_cell_types] <- "Epi. Unassigned"
			#rna@meta.data[idx_tumor_epi_cnv_unassigned, col_cell_types] <- paste(rna$tumor_epi_type[idx_tumor_epi_cnv_unassigned], "Unassigned")
		}
	},

	"cancer_normal_epithelial_cell_types"={
		# normal + tumor samples

		# update rna$epi_type with rna$normal_epi_type
		if (args$f_infercnv_pos_notpos) {
  			rna$epi_type[idx_normal_epi_cnv_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_pos]
  			rna$epi_type[idx_normal_epi_cnv_not_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_not_pos]
		} else {
  			rna$epi_type[idx_normal_epi_cnv_pos] <- rna$normal_epi_type[idx_normal_epi_cnv_pos]
  			rna$epi_type[idx_normal_epi_cnv_neg] <- rna$normal_epi_type[idx_normal_epi_cnv_neg]
  			rna$epi_type[idx_normal_epi_cnv_unassigned] <- rna$normal_epi_type[idx_normal_epi_cnv_unassigned]
		}

		# update rna$epi_type with idx_tumor_epi
		if (args$f_infercnv_pos_notpos) {
			rna$epi_type[idx_tumor_epi_cnv_pos] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna$epi_type[idx_tumor_epi_cnv_not_pos] <- rna$normal_epi_type[idx_tumor_epi_cnv_not_pos] # use normal epi type for cnv-
		} else {
			rna$epi_type[idx_tumor_epi_cnv_pos] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			rna$epi_type[idx_tumor_epi_cnv_neg] <- rna$normal_epi_type[idx_tumor_epi_cnv_neg] # use normal epi type for cnv-
			rna$epi_type[idx_tumor_epi_cnv_unassigned] <- "Epi. Unassigned"
		}

		if (args$f_infercnv_pos_notpos) {
			# LumA/LumB/Her2E/Basal/Epi. Tumor --> LumA/LumB/Her2/Basal/Epi. Tumor/Epi. Non-tumor
			# epithelial cells outside epi. clusters are annotated as non-tumor.
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Non-tumor"
			rna@meta.data[idx_tumor_epi_cnv_pos, col_cell_types] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			# use normal epi type for cnv-
			rna@meta.data[idx_tumor_epi_cnv_not_pos, col_cell_types] <- rna$normal_epi_type[idx_tumor_epi_cnv_not_pos]
			rna@meta.data[idx_normal_epi_cnv_not_pos, col_cell_types] <- rna$normal_epi_type[idx_normal_epi_cnv_not_pos]
		} else {
			# LumA/LumB/Her2E/Basal/Epi. Tumor --> LumA/LumB/Her2/Basal/Epi. Tumor/Epi. Non-tumor/Epi. Unassigned
			# epithelial cells outside epi. clusters are annotated as unassigned.
			rna@meta.data[f_epi, col_cell_types] <- "Epi. Unassigned"
			rna@meta.data[idx_tumor_epi_cnv_pos, col_cell_types] <- rna$tumor_epi_type[idx_tumor_epi_cnv_pos]
			# use normal epi type for cnv-
			rna@meta.data[idx_tumor_epi_cnv_neg, col_cell_types] <- rna$normal_epi_type[idx_tumor_epi_cnv_neg]
			rna@meta.data[idx_normal_epi_cnv_neg, col_cell_types] <- rna$normal_epi_type[idx_normal_epi_cnv_neg]
			# unassigned
			rna@meta.data[idx_tumor_epi_cnv_unassigned, col_cell_types] <- "Epi. Unassigned"
			#rna@meta.data[idx_tumor_epi_cnv_unassigned, col_cell_types] <- paste(rna$tumor_epi_type[idx_tumor_epi_cnv_unassigned], "Unassigned")
		}
	},
	{}
  ) # switch

  rna

} # update_cell_types_after_infercnv


















### utilities


# get_genes_from_signatures
get_genes_from_signatures <- function(args) {

  switch(args$type_signatures,
	"panglaodb"={
		genes <- unlist(panglaodb)
	},
	"panglaodb_plus_others"={
		genes <- unlist(panglaodb)
		genes <- union(genes, mast_cells_11genes)
		genes <- union(genes, unlist(list_markers_peroulab))
		genes <- union(genes, unlist(list_markers_bartlettlab))
		genes <- union(genes, unlist(list_markers_kessenbrocklab_paper))
		genes <- union(genes, unlist(list_markers_kessenbrocklab))
	},
	"peroulab_normal_epi"={
		genes <- unlist(list_markers_peroulab)
	},
	"swarbricklab_normal_epi"={
		#genes <- unlist(list_markers_kessenbrocklab)
		filename_rds_centroids <- "reference/centroids_for_singler/swarbricklab_breast_normal_epi_3centroids.rds"
		cat(sprintf("\t\treadRDS('%s')\n", filename_rds_centroids))
		se_centroid <- readRDS(filename_rds_centroids)
		genes <- rownames(se_centroid)
	},
	"kessenbrocklab_normal_epi"={
		#genes <- unlist(list_markers_kessenbrocklab)
		filename_rds_centroids <- "reference/centroids_for_singler/kessenbrocklab_breast_normal_epi_5centroids.rds"
		cat(sprintf("\t\treadRDS('%s')\n", filename_rds_centroids))
		se_centroid <- readRDS(filename_rds_centroids)
		genes <- rownames(se_centroid)
	},
	{
		genes <- NULL
	}
  ) # switch

  genes <- unname(genes)
  genes <- unique(genes)
  genes

} # get_genes_from_signatures




# select_genes_with_seurat_obj 
select_genes_with_seurat_obj <- function(rna, genes, f_intersect=TRUE) {

  genes <- genes[nchar(genes) > 0]
  if (f_intersect) {
    genes <- intersect(genes, rownames(rna))
  }

  genes

} # select_genes_with_seurat_obj







# make_se_centroid
#
# input:
#   rna: seurat obj
#   col_celltype <- "celltype_minor"
#   list_markers: list
#   filename_prefix_dge:
#   filename_rds_centroids:
#   nv_cell_type_conversion_table:
#   th_ncount: "med" or numerical value
#   th_nfeature: "med" or numerical value
#   genes_addition: biomarkers to add
#   list_rules: 
#	list_rules <- list("Her2E"=list("min_expr"=c("ERBB2"=1.0)))
#   list_biomarkers:
#	list_biomarkers <- list("CLow"=list("expressed"=c("ALDH1A1", "ZEB2"), "not expressed"=c("KRT18", "KRT19", "ESR1", "GATA3", "ERBB2", "CLDN3", "CLDN4", "CLDN7", "CDH1", "OCLN", "CD24", "EPCAM", "MUC1"), "subtype of"="Basal"))
#
# usage:
# filename_prefix_dge <- "tsv/swarbricklab_breast_normal_epi_3centroids"
# filename_rds_centroids <- "reference/centroids_for_singler/swarbricklab_breast_normal_epi_3centroids.rds"
# nv_cell_type_conversion_table <- c( "Mature Luminal"="LEp", "Luminal Progenitors"="LEp_prog", "Myoepithelial"="BEp")
# se_centroid <- make_se_centroid(rna, "celltype_minor", list_markers, filename_prefix_dge=filename_prefix_dge, filename_rds_centroids=filename_rds_centroids, nv_cell_type_conversion_table=nv_cell_type_conversion_table)
make_se_centroid <- function(rna, col_celltype, list_markers=NULL, filename_prefix_dge=NULL, filename_rds_centroids=NULL, nv_cell_type_conversion_table=NULL, th_percent.mt=NULL, th_ncount=NULL, th_nfeature=NULL, df_centroid_to_check=NULL, genes_addition=NULL, list_rules=NULL, list_biomarkers=NULL) {

  suppressPackageStartupMessages(library(SummarizedExperiment))

  cat(sprintf("\n\tmake_se_centroid()\n\n"))
  groups <- names(list_markers)

  # filter cells
  f_cells <- rep(TRUE, ncol(rna))

  f_celltype <- rna@meta.data[,col_celltype] %in% groups
  f_cells <- f_cells & f_celltype

  if (!is.null(th_percent.mt)) {
  	if (th_percent.mt == "med") {
		th_percent.mt <- median(rna$percent.mt, na.rm = TRUE) - 2*stats::mad(rna$percent.mt, na.rm = TRUE)
 	}
  	cat(sprintf("\t\tth_percent.mt: %g\n", th_percent.mt))
	f_cells <- f_cells & (rna$percent.mt < th_percent.mt)
  }

  if (!is.null(th_ncount)) {
  	if (th_ncount == "med") {
		th_ncount <- median(rna$nCount_RNA, na.rm = TRUE)
 	 }
  	cat(sprintf("\t\tth_ncount: %g\n", th_ncount))
	f_cells <- f_cells & (rna$nCount_RNA > th_ncount)
  }

  if (!is.null(th_nfeature)) {
	if (th_nfeature == "med") {
		th_nfeature <- median(rna$nFeature_RNA, na.rm = TRUE)
  	}
	cat(sprintf("\t\tth_nfeature: %g\n", th_nfeature))
	f_cells <- f_cells & (rna$nFeature_RNA > th_nfeature)
  }

  if (!all(f_cells)) {
	rna <- rna[,f_cells]
  }

  # mtx_logcounts
  mtx_logcounts <- GetAssayData(object = rna, assay="RNA", slot = "data")

  # genes_marker
  genes_marker <- unique(unname(unlist(list_markers)))
  if (length(genes_marker) == 0) {
	#genes_marker <- rownames(rna)

	log_obj(table(rna@meta.data[,col_celltype]), tab=2)

	list_df <- list()
	genes_marker <- c()
	for (group1 in groups) {

		# find_markers()
		cat(sprintf("\t\t%s\n", group1))
		list_out <- find_markers(rna, cell_type1=group1, cell_type_ref=NULL, group_name1=group1, group_name_ref="others", col_cluster_types=NULL, col_cell_types=col_celltype, th_log2fc=0.25, th_padj=0.01, min.pct=0.25, min.diff.pct=-Inf, max.cells.per.ident=Inf, method_dge="seurat_findmarkers_enricher", str_condition=NULL, n_log=1)
		df <- list_out$markers
		df <- df[order(-df$avg_log2FC),]
		# df_out
		df_out <- df
		cols <- colnames(df_out); idx <- grep('.', cols)
		df_out[,cols[idx]] <- sapply(cols[idx], function(x) {
			formatC(df_out[,x], digits=2)
		})
		df_out$comparison <- sprintf("%s vs. others", group1)
		filename_dge <- sprintf("%s_%s_vs_others.tsv", filename_prefix_dge, tolower(group1))
		filename_dge <- gsub(" ", "_", filename_dge)
 		write.table(df_out, file=filename_dge, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
		df <- list_out$markers
		#         p_val   avg_log2FC  pct.1  pct.2  p_val_adj
		f_select <- ((df$avg_log2FC > 1) & (df$pct.1 > 0.5) & (df$pct.2 < 0.2) & (df$p_val_adj < 1e-12))
		# KRT18   2.4e-61 0.88      1       1     6.6e-57 0.96    2.7     4.9     Cancer LumA SC vs. others
		# KRT8    3.9e-14 0.66    0.99    0.99    1.1e-09   0       2     4.4     Cancer LumA SC vs. others
		# KRT19   1.4e-252        1.2       1     0.97    3.8e-248        1.3     3.8     4.7     Cancer LumB SC vs. others
		# KRT18   2.4e-86 0.3       1       1     6.5e-82 1.1     2.6     3.6     Cancer LumB SC vs. others
		# ERBB2   4.2e-172   3.97      0.99  0.632  1.17e-167
		f_select <- f_select | ((df$avg_log2FC > 3) & (df$pct.1 > 0.5) & (df$pct.1 - df$pct.2 > 0.3) & (df$p_val_adj < 1e-12))
		# check markers with df_centroid_to_check
		if (!is.null(df_centroid_to_check)) {
			group1_ <- group1
  			if (!is.null(nv_cell_type_conversion_table)) {
				group1_ <- nv_cell_type_conversion_table[group1]
			}
			syms_include <- rownames(df_centroid_to_check)[(df_centroid_to_check$subtype == group1_)]
			syms_include <- syms_include[syms_include %in% rownames(df_out)]
			log_txt("\t\tdf_centroid_to_check\n")
			log_obj(df_out[syms_include,], tab=3)
			#f_select <- f_select | rownames(df_out) %in% syms_include
		}
		# exclude too lowly expressed markers
		f_select <- f_select & (df$expr_max > 2.0)
		# df_select
		df_select <- df[f_select,]
		df_select <- df_select[order(-df_select$expr_mean),]
		list_df[[group1]] <- df_select
		genes_marker <- c(genes_marker, rownames(df_select))
	} # for
	tb <- table(genes_marker)
	genes_dup <- names(tb[tb > 1])
	cat(sprintf("\t\tgenes_dup: %s\n", paste(genes_dup, collapse=", ")))
	
	# select
	genes_marker <- c()
	df_out_select <- data.frame()
	for (group1 in groups) {
		cat(sprintf("\t\t%s\n", group1))
		df <- list_df[[group1]]
		# df_select
		f_select <- !rownames(df) %in% genes_dup
		df_select <- df[f_select,]
		df_select <- head(df[f_select,], 10)
		# df_out
		df_out <- df_select
		cols <- colnames(df_out); idx <- grep('.', cols)
		df_out[,cols[idx]] <- sapply(cols[idx], function(x) {
			 formatC(df_out[,x], digits=2)
		})
		#filename_dge <- sprintf("%s_%s_vs_others_selected.tsv", filename_prefix_dge, tolower(group1))
		#filename_dge <- gsub(" ", "_", filename_dge)
 		#write.table(df_out, file=filename_dge, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
		df_out_select <- rbind(df_out_select, df_out)
		log_obj(df_out, tab=2)
		genes_marker <- union(genes_marker, rownames(df_select))
	} # for

	filename_dge <- sprintf("%s_dge.tsv", filename_prefix_dge)
	filename_dge <- gsub(" ", "_", filename_dge)
 	write.table(df_out_select, file=filename_dge, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

  } # if


  genes_marker_dge <- genes_marker
  if (!is.null(genes_addition)) {
	cat(sprintf("\t\tgenes_addition: %s\n", paste(genes_addition, collapse=", ")))
	genes_marker <- union(genes_marker, genes_addition)
  }

  genes_marker_x <- genes_marker[!genes_marker %in% rownames(mtx_logcounts)]
  if (length(genes_marker_x) > 0) {
	cat(sprintf("\t\tgenes_marker_x: %s\n", paste(genes_marker_x, collapse=", ")))
  }


  genes_marker <- genes_marker[genes_marker %in% rownames(mtx_logcounts)]

  # mtx_centroid
  mtx_centroid <- matrix(0, length(genes_marker), length(list_markers),
                    dimnames=list(genes_marker, groups))

  for (group1 in groups) {
    idx.cells <- which(rna@meta.data[,col_celltype] == group1)
    mtx_centroid[,group1] <- rowMeans(mtx_logcounts[genes_marker, idx.cells])
  } # for

  rownames(mtx_centroid) <- genes_marker

  if (!is.null(nv_cell_type_conversion_table)) {
	cols <- colnames(mtx_centroid)
	idx <- match(cols, names(nv_cell_type_conversion_table))
	f <- !is.na(idx); idx <- idx[f]
	if (any(f)) {
		cols[f] <- nv_cell_type_conversion_table[idx]	
		colnames(mtx_centroid) <- cols
	}
  } # if

  cat(sprintf("\t\tcentroid rows %d: %s, ...\n", dim(mtx_centroid)[1], paste(head(rownames(mtx_centroid)), collapse=", ")))
  cat(sprintf("\t\tcentroid columns: %s\n", paste(colnames(mtx_centroid), collapse=", ")))

  # df_rowdata
  subtypes <- colnames(mtx_centroid)
  subtypes_syms_centroid <- subtypes[max.col(mtx_centroid, ties.method="first")]
  df_rowdata <- data.frame(sym=genes_marker, subtype=subtypes_syms_centroid, f_dge=FALSE, min_expr=NA)
  rownames(df_rowdata) <- genes_marker
  if (!is.null(genes_addition)) {
	df_rowdata[genes_marker_dge, "f_dge"] <- TRUE
	df_rowdata[!df_rowdata$f_dge, "subtype"] <- NA
  }

  # df_coldata
  df_coldata <- data.frame(label.main=colnames(mtx_centroid))
  rownames(df_coldata) <- colnames(mtx_centroid)

  # metadata
  metadata <- list()
  if (!is.null(list_rules)) {
	metadata[["list_rules"]] <- list_rules
	for (subtype in names(list_rules)) {
		list_info <- list_rules[[subtype]]
		if ("min_expr" %in% list_info) {
			nv_min_expr <- list_info[["min_expr"]]	
			for (gene in names(nv_min_expr)) {
				min_v <- nv_min_expr[gene]
				df_rowdata[gene, "min_expr"] <- min_v
			} # for gene
		}
	} # for subtype
  }

  if (!is.null(list_biomarkers)) {
	metadata[["list_biomarkers"]] <- list_biomarkers
  }

  se_centroid <- SummarizedExperiment(assays=list(logcounts=mtx_centroid), rowData = df_rowdata, colData = df_coldata, metadata=metadata)
  if (!is.null(filename_rds_centroids)) {
  	cat(sprintf("\t\tsave %s\n", filename_rds_centroids))
	saveRDS(se_centroid, filename_rds_centroids)
  }

  cat(sprintf("\n"))

  filename_dge <- sprintf("%s.tsv", filename_prefix_dge)
  filename_dge <- gsub(" ", "_", filename_dge)
  cat(sprintf("\t\twrite.table('%s')\n", filename_dge))
  write.table(mtx_centroid, file=filename_dge, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

  se_centroid

} # make_se_centroid






