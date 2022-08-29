#
# motif_enrichment_analysis.R
# author: H. Kim
# date created: 2021, Dec.
# date last modified: 2021, Dec.
#
#
#
# contents:
# add_archr_anno_enrichemnt()
# add_archr_anno_enrichemnt_encodetfbs()
#
#
#



suppressPackageStartupMessages(library(latex2exp))















# add_archr_anno_enrichemnt
# 
# input:
#   proj.archr:
#   p2g.df.sub.plot:
#   markerPeaks.all:
#   peakName_selected:
#   markerPeaks.name:
#   peakName_selected.name:
#   str_motif_set:
#   df_output:
#   args:
#
# output:
#   list_out$df_output: updated df_output
#
# usage:
# list_out <- add_archr_anno_enrichemnt(proj.archr, markerPeaks.all, peakName.non_cancer_specific, "normal-vs-cancer", "non-cancer-specific", str_motif_set="cisbp", df_output, args)
add_archr_anno_enrichemnt <- function(proj.archr, p2g.df.sub.plot, markerPeaks.all, peakName_selected, markerPeaks.name, peakName_selected.name, str_motif_set="cisbp", df_output, args) {

  # filename_appendix
  # e.g. normal-vs-cancer_cancer-specific_enhancer
  filename_appendix <- sprintf("%s_%s_%s", markerPeaks.name, peakName_selected.name, peaktype)


  cat(sprintf("\n\t------------------------------------\n"))
  cat(sprintf("\t%s\n\n", filename_appendix))

  # peakset
  peakset <- proj.archr@peakSet
  peakset.peakName <- sprintf("%s:%d-%d", seqnames(peakset), start(peakset), end(peakset))

  # peakName.markerpeaks.all
  df_rowdata <- rowData(markerPeaks.all)
  peakName.markerpeaks.all <- sprintf("%s:%d-%d", df_rowdata$seqnames, df_rowdata$start, df_rowdata$end)

  # peakName.markerpeaks.for_comparison
  idx <- match(peakName.markerpeaks.all, peakName_selected)
  f.peaks_for_comparison <- !is.na(idx)
  markerPeaks.for_comparison <- markerPeaks.all[f.peaks_for_comparison,]
  peakName.markerpeaks.for_comparison <- peakName.markerpeaks.all[f.peaks_for_comparison]

  # p2g.coords.sub
  f <- p2g.df.sub.plot$peakName %in% peakName_selected
  p2g.coords.sub <- p2g.df.sub.plot[f,]
  

  # update df_output
  n_markerpeaks.for_comparison <- nrow(markerPeaks.for_comparison)
  var <- sprintf("n_markerpeaks.for_comparison.%s", filename_appendix)
  cat(sprintf("\t%s: %d\n", var, n_markerpeaks.for_comparison))
  df_output[var, "info"] <- n_markerpeaks.for_comparison

  # idx_match
  idx_match <- match(peakName.markerpeaks.for_comparison, peakset.peakName)
  f.na <- is.na(idx_match)
  if (any(f.na)) {
        stop(sprintf("there are unmappable peaks in peakName.markerpeaks.for_comparison"))
  }

  # https://www.archrproject.com/reference/addPeakSet.html
  cat(sprintf("\taddPeakSet\n"))
  proj.archr.for_comparison <- addPeakSet(
                ArchRProj = proj.archr,
                peakSet = peakset[idx_match,],
                genomeAnnotation = getGenomeAnnotation(proj.archr),
                force = TRUE
        ) # addPeakSet

  # https://www.archrproject.com/reference/addMotifAnnotations.html
  cat(sprintf("\taddMotifAnnotations\n"))
  proj.archr.for_comparison <- addMotifAnnotations(
                ArchRProj = proj.archr.for_comparison,
                motifSet = str_motif_set, # The motif set to be used for annotation. Options include: (i) "JASPAR2016", "JASPAR2018", "JASPAR2020" which gives the 2016, 2018 or 2020 version of JASPAR motifs or (ii) one of "cisbp", "encode", or "homer" which gives the corresponding motif sets from the chromVAR package.
                name = "Motif",
                force = TRUE,
                logFile = createLogFile("addMotifAnnotations", logDir=dir_log)
        )

  # save
  saveRDS(markerPeaks.for_comparison, sprintf("%s/markerpeaks.for_comparison.%s.rds", dir_rds, filename_appendix))
  saveRDS(proj.archr.for_comparison, sprintf("%s/proj.archr.for_comparison.%s.rds", dir_rds, filename_appendix))






  ### up

  # markerPeaks.for_comparison <- readRDS("output_p2g_male-bc/rds/markerpeaks.for_comparison.enhancer.rds"); proj.archr.for_comparison <- readRDS("output_p2g_male-bc/rds/proj.archr.for_comparison.enhancer.rds"); peaktype <- "enhancer"; df_output <- data.frame(); dir_log <- "log" # for debug
  cat(sprintf("\tpeakAnnoEnrichment for up-regulated motifs in %s\n", filename_appendix))
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr.for_comparison,
    peakAnnotation = "Motif",
    #cutOff = "FDR < 0.1 & Log2FC > 0.5",
    cutOff = "FDR < 0.01 & Log2FC > 1.0",
    logFile = createLogFile("peakAnnoEnrichmentUp", logDir=dir_log)
  )


  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  filename_tsv <- sprintf("%s/peakannoenrichment_%s_motif_up.%s.tsv", dir_tsv, str_motif_set, filename_appendix)
  cat(sprintf("\t\twrite.table(df, '%s')\n", filename_tsv))
  write.table(df, filename_tsv, row.names = FALSE, col.names = TRUE, sep = "\t", quote = F)

  n_motifs.up <- nrow(df)
  var <- sprintf("n_motifs.up.%s", filename_appendix)
  cat(sprintf("\t\t%s: %d\n", var, n_motifs.up))
  df_output[var, "info"] <- n_motifs.up


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
    scale_color_gradientn(
        name = TeX(r'(-$log_{10} (P_{adj})$)'),
        colors = paletteContinuous(set = "comet"))

  ggsave(sprintf("%s/scatterplot_%s_motif_up.%s.pdf", dir_pdf, str_motif_set, filename_appendix), width = 5, height = 5.5, plot=ggUp)
  #ggsave(sprintf("%s/scatterplot_%s_motif_up.%s.png", dir_png, str_motif_set, filename_appendix), width = 5, height = 5.5, plot=ggUp)

  cutOff_mlog10Padj <- 5
  if (length(which(df$mlog10Padj > cutOff_mlog10Padj)) > 0) {
    heatmapEM <- plotEnrichHeatmap(motifsUp, n = 20,
                cutOff = cutOff_mlog10Padj, # default=20
                transpose = TRUE,
                returnMatrix = FALSE,
                logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log)
    )

    # ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    plotPDF(heatmapEM, name = sprintf("heatmap_%s_motif_up.%s", str_motif_set, filename_appendix), width = 8, height = 5, ArchRProj = proj.archr.for_comparison, addDOC = FALSE)
  } # if







  ### down


  # https://www.archrproject.com/reference/peakAnnoEnrichment.html
  # This function will perform hypergeometric enrichment of a given peak annotation within the defined marker peaks.
  cat(sprintf("\tpeakAnnoEnrichment for down-regulated motifs in %s\n", filename_appendix))
  motifsDo <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr.for_comparison,
    peakAnnotation = "Motif",
    #cutOff = "FDR < 0.1 & Log2FC < -0.5",
    cutOff = "FDR < 0.01 & Log2FC < -1.0",
    logFile = createLogFile("peakAnnoEnrichment_Do", logDir=dir_log)
    )


  df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  filename_tsv <- sprintf("%s/peakannoenrichment_%s_motif_down.%s.tsv", dir_tsv, str_motif_set, filename_appendix)
  cat(sprintf("\t\twrite.table(df, '%s')\n", filename_tsv))
  write.table(df, filename_tsv, row.names = FALSE, col.names = TRUE, sep = "\t", quote = F)

  n_motifs.down <- nrow(df)
  var <- sprintf("n_motifs.down.%s", filename_appendix)
  cat(sprintf("\t\t%s: %d\n", var, n_motifs.down))
  df_output[var, "info"] <- n_motifs.down

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
    scale_color_gradientn(
        name = TeX(r'(-$log_{10} (P_{adj})$)'),
        colors = paletteContinuous(set = "comet"))

  ggsave(sprintf("%s/scatterplot_%s_motif_down.%s.pdf", dir_pdf, str_motif_set, filename_appendix), width = 5, height = 5.5, plot=ggDo)
  #ggsave(sprintf("%s/scatterplot_%s_motif_down.%s.png", dir_png, str_motif_set, filename_appendix), width = 5, height = 5.5, plot=ggDo)

  if (length(which(df$mlog10Padj > cutOff_mlog10Padj)) > 0) {
    heatmapEM <- plotEnrichHeatmap(motifsDo, n = 20,
                cutOff = cutOff_mlog10Padj, # default=20
                transpose = TRUE,
                returnMatrix = FALSE,
                logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log)
    )

    # ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    plotPDF(heatmapEM, name = sprintf("heatmap_%s_motif_down.%s", str_motif_set, filename_appendix), width = 8, height = 5, ArchRProj = proj.archr.for_comparison, addDOC = FALSE)
  } # if






  ### up_down

  plotPDF(ggUp, ggDo, name = sprintf("scatterplots_%s_motif_updown.%s", str_motif_set, filename_appendix), width = 5, height = 5.5, ArchRProj = proj.archr.for_comparison, addDOC = FALSE)




  if (args$f_make_tfsee_input) {

        source("./r/make_tfsee_input.R")
	dir_bam_subset <- sprintf("%s/bam_subset_archr.%s", dir_output, filename_appendix)
        out <- make_tfsee_input(proj.archr.for_comparison, dir_bam_subset, p2g.coords.sub=p2g.coords.sub, markerpeaks=NULL)

  } # if





  list_out <- list()
  #list_out$proj.archr.for_comparison <- proj.archr.for_comparison
  list_out$df_output <- df_output

  list_out





} # add_archr_anno_enrichemnt


















# add_archr_anno_enrichemnt_encodetfbs
# 
# input:
#   proj.archr.for_comparison:
#   df_output:
#
# output:
#   list_out$df_output: updated df_output
#
#
add_archr_anno_enrichemnt_encodetfbs <- function(proj.archr.for_comparison, df_output) {


  # 12.3 ArchR enrichment
  str_collection <- "EncodeTFBS"
  #str_collection <- "ATAC"
  #str_collection <- "Codex"


  # https://www.archrproject.com/reference/addArchRAnnotations.html
  cat(sprintf("\taddArchRAnnotations for %s\n", str_collection))
  proj.archr.for_comparison <- addArchRAnnotations(
        ArchRProj = proj.archr.for_comparison,
        collection = str_collection,
        force = TRUE,
        logFile = createLogFile("addArchRAnnotations", logDir=dir_log)
    )




  ### up

  cat(sprintf("\tpeakAnnoEnrichment for up-regulated motifs with %s\n", str_collection))
  motifsUp <- peakAnnoEnrichment(
    seMarker = markerPeaks.for_comparison,
    ArchRProj = proj.archr.for_comparison,
    peakAnnotation = str_collection,
    #cutOff = "FDR < 0.1 & Log2FC > 0.5",
    cutOff = "FDR < 0.01 & Log2FC > 1.0",
    logFile = createLogFile(sprintf("peakAnnoEnrichment_%s_Up", str_collection), logDir=dir_log)
  )

  df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
  df <- df[order(df$mlog10Padj, decreasing = TRUE),]
  df$rank <- seq_len(nrow(df))

  write.table(df, sprintf("%s/peakannoenrichment_%s_motif_up.cancer_clusters.%s.tsv", dir_tsv, str_collection, peaktype), row.names = TRUE, col.names = NA, sep = "\t", quote = F)

  n_motifs.up.EncodeTFBS <- nrow(df)
  var <- sprintf("n_motifs.up.%s.EncodeTFBS", peaktype)
  cat(sprintf("\t\t%s: %d\n", var, n_motifs.up.EncodeTFBS))
  df_output[var, "info"] <- n_motifs.up.EncodeTFBS

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
    scale_color_gradientn(
	name = TeX(r'(-$log_{10} (P_{adj})$)'),
	colors = paletteContinuous(set = "comet")) 

  ggsave(sprintf("%s/scatterplot_%s_motif_up.cancer_clusters.%s.pdf", dir_pdf, str_collection, peaktype), width = 5, height = 5.5, plot=ggUp)
  #ggsave(sprintf("%s/scatterplot_%s_motif_up.cancer_clusters.%s.png", dir_png, str_collection, peaktype), width = 5, height = 5.5, plot=ggUp)


  cutOff_mlog10Padj <- 5
  if (length(which(df$mlog10Padj > cutOff_mlog10Padj)) > 0) {
    heatmapEM <- plotEnrichHeatmap(motifsUp, n = 20,
                cutOff = cutOff_mlog10Padj, # default=20
                transpose = TRUE,
                returnMatrix = FALSE,
                logFile = createLogFile("plotEnrichHeatmap", logDir=dir_log)
    )

    # ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
    plotPDF(heatmapEM, name = sprintf("heatmap_%s_motif_up.cancer_clusters.%s", str_collection, peaktype), width = 8, height = 6, ArchRProj = proj.archr.for_comparison, addDOC = FALSE)

  } # if






  ### down


  # not implemented yet

  list_out <- list()
  #list_out$proj.archr.for_comparison <- proj.archr.for_comparison
  list_out$df_output <- df_output

  list_out

} # add_archr_anno_enrichemnt_encodetfbs




