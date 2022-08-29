#
# plotPeak2GeneHeatmap.distal.R
# modified by H. Kim
# last modified at 2021, Oct.
#
# input:
#   k: (default=25) An integer describing the number of k-means clusters to group peak-to-gene links prior to plotting heatmaps.
#   returnMatrices: (default=FALSE) A boolean value that determines whether the matrices should be returned with kmeans id versus plotting.
#
# output:
#   returnMatrices=FALSE, heatmap
#   returnMatrices=TRUE
#     out <- SimpleList(
#       ATAC = SimpleList(
#         matrix = mATAC[kDF[,2],colOrder],
#         kmeansId = kDF[,3],
#         colData = cD[colOrder,,drop=FALSE]
#       ),
#       RNA = SimpleList(
#         matrix = mRNA[kDF[,2],colOrder],
#         kmeansId = kDF[,3],
#         colData = cD[colOrder,,drop=FALSE]
#       ),
#       Peak2GeneLinks = p2g[kDF[,2],]
#     )
#
#
#
# reference:
# This function plots side by side heatmaps of linked ATAC and Gene regions from addPeak2GeneLinks.  https://www.archrproject.com/reference/plotPeak2GeneHeatmap.html
# 





plotPeak2GeneHeatmap.distal <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  FDRCutOff = 0.0001,
  peaks,
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  k = 25,
  nPlot = 25000,
  limitsATAC = c(-2, 2),
  limitsRNA = c(-2, 2),
  groupBy = "Clusters",
  palGroup = NULL,
  palATAC = paletteContinuous("solarExtra"),
  palRNA = paletteContinuous("blueYellow"),
  verbose = TRUE,
  returnMatrices = FALSE,
  seed = 1,
  # begin of addition by H. Kim
  legend_ncol = 4,
  # end of addition
  # begin of modification by H. Kim
  #logFile = createLogFile("plotPeak2GeneHeatmap")
  logFile = createLogFile("plotPeak2GeneHeatmap", logDir=dir_log, useLogs=TRUE)
  # end of modification
) {
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  .validInput(input = FDRCutOff, name = "FDRCutOff", valid = c("numeric"))
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = limitsATAC, name = "limitsATAC", valid = c("numeric"))
  .validInput(input = limitsRNA, name = "limitsRNA", valid = c("numeric"))
  .validInput(input = groupBy, name = "groupBy", valid = c("character"))
  .validInput(input = palGroup, name = "palGroup", valid = c("palette", "null"))
  .validInput(input = palATAC, name = "palATAC", valid = c("palette", "null"))
  .validInput(input = palRNA, name = "palRNA", valid = c("palette", "null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = returnMatrices, name = "returnMatrices", valid = c("boolean"))
  .validInput(input = seed, name = "seed", valid = c("integer"))
  .validInput(input = logFile, name = "logFile", valid = c("character"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "peak2GeneHeatmap Input-Parameters", logFile = logFile)
  
  if(is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    stop("No Peak2GeneLinks Found! Try addPeak2GeneLinks!")
  }
  
  #########################################
  # Get Inputs
  #########################################
  ccd <- getCellColData(ArchRProj, select = groupBy)
  p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks
  p2g$idx <- paste0(p2g$idxATAC,"-",p2g$idxRNA)
  p2g <- p2g[p2g$idx %in% peaks$idx,]
  # if(!is.null(varCutOffATAC)){
  #   p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
  # }
  # 
  # if(!is.null(varCutOffRNA)){
  #   p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
  # }
  
  if(nrow(p2g) == 0){
    stop("No peak2genelinks found with your cutoffs!")
  }
  
  if(!file.exists(metadata(p2g)$seATAC)){
    stop("seATAC does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  if(!file.exists(metadata(p2g)$seRNA)){
    stop("seRNA does not exist! Did you change paths? If this does not work, please try re-running addPeak2GeneLinks!")
  }
  mATAC <- readRDS(metadata(p2g)$seATAC)[p2g$idxATAC, ]
  mRNA <- readRDS(metadata(p2g)$seRNA)[p2g$idxRNA, ]
  p2g$peak <- paste0(rowRanges(mATAC))
  p2g$gene <- rowData(mRNA)$name
  gc()
  
  mATAC <- assay(mATAC)
  mRNA <- assay(mRNA)
  
  #########################################
  # Determine Groups from KNN
  #########################################
  .logDiffTime(main="Determining KNN Groups!", t1=tstart, verbose=verbose, logFile=logFile)
  KNNList <- as(metadata(readRDS(metadata(p2g)$seRNA))$KNNList, "list")
  KNNGroups <- lapply(seq_along(KNNList), function(x){
    KNNx <- KNNList[[x]]
    names(sort(table(ccd[KNNx, 1, drop = TRUE]), decreasing = TRUE))[1]
  }) %>% unlist
  cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups)
  pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
  if(!is.null(palGroup)){
    pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
  }
  colorMap <- list(groupBy = pal)
  attr(colorMap[[1]], "discrete") <- TRUE
  
  #########################################
  # Organize Matrices
  #########################################
  mATAC <- .rowZscores(mATAC)
  mRNA <- .rowZscores(mRNA)
  rownames(mATAC) <- NULL
  rownames(mRNA) <- NULL
  colnames(mATAC) <- paste0("K_", seq_len(ncol(mATAC)))
  colnames(mRNA) <- paste0("K_", seq_len(ncol(mRNA)))
  rownames(mATAC) <- paste0("P2G_", seq_len(nrow(mATAC)))
  rownames(mRNA) <- paste0("P2G_", seq_len(nrow(mRNA)))
  rownames(p2g) <- paste0("P2G_", seq_len(nrow(p2g)))
  
  .logDiffTime(main="Ordering Peak2Gene Links!", t1=tstart, verbose=verbose, logFile=logFile)


  if(!is.null(seed)){
    set.seed(seed)
  }

  k1 <- kmeans(mATAC, k)

  if(nrow(mATAC) > nPlot){
    nPK <- nPlot * table(k1$cluster) / length(k1$cluster) 
    splitK <- split(seq_len(nrow(mATAC)), k1$cluster)
    kDF <- lapply(seq_along(splitK), function(x){
      idx <- sample(splitK[[x]], floor(nPK[x]))
      k <- rep(x, length(idx))
      DataFrame(k = k, idx = idx)
    }) %>% Reduce("rbind", .)
  }else{
    kDF <- DataFrame(k = k1$cluster, idx = seq_len(nrow(mATAC)))
  }


  bS <- .binarySort(t(.groupMeans(t(mATAC[kDF[,2],]), kDF[,1])),  clusterCols = TRUE, cutOff = 1)
  rowOrder <- rownames(bS[[1]])
  colOrder <- colnames(bS[[1]])
  kDF[,3] <- as.integer(mapLabels(paste0(kDF[,1]), newLabels = paste0(seq_along(rowOrder)), oldLabels = rowOrder))
  
  if (returnMatrices){
    
    out <- SimpleList(
      ATAC = SimpleList(
        matrix = mATAC[kDF[,2],colOrder],
        kmeansId = kDF[,3],
        colData = cD[colOrder,,drop=FALSE]
      ),
      RNA = SimpleList(
        matrix = mRNA[kDF[,2],colOrder],
        kmeansId = kDF[,3],
        colData = cD[colOrder,,drop=FALSE]
      ),
      Peak2GeneLinks = p2g[kDF[,2],]
    )
    
    return(out)
    
  }
  
  #Log Info
  .logThis(colorMap, "colorMap", logFile = logFile)
  .logThis(colOrder, "colOrder", logFile = logFile)
  .logThis(kDF, "kDF", logFile = logFile)
  .logThis(mATAC, "mATAC", logFile = logFile)
  .logThis(mRNA, "mRNA", logFile = logFile)
  .logThis(cD[colOrder,,drop=FALSE], "cD", logFile = logFile)
  .logThis(mATAC[kDF[,2],colOrder], "mATAC2", logFile = logFile)
  .logThis(mRNA[kDF[,2],colOrder], "mRNA2", logFile = logFile)
  
  #########################################
  # Plot Heatmaps
  #########################################
  .logDiffTime(main="Constructing ATAC Heatmap!", t1=tstart, verbose=verbose, logFile=logFile)
  htATAC <- .ArchRHeatmap(
    mat = mATAC[kDF[,2],colOrder],
    scale = FALSE,
    limits = limitsATAC,
    color = palATAC, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    # begin of addition by H. Kim
    useRaster = TRUE,
    rasterQuality = 1, # default=5
    rasterDevice = "png",
    show_annotation_name = FALSE,
    annotation_legend_param = list(ncol = legend_ncol),
    # end of addition
    draw = FALSE,
    name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks")
  )
  
  .logDiffTime(main = "Constructing RNA Heatmap!", t1 = tstart, verbose = verbose, logFile = logFile)
  htRNA <- .ArchRHeatmap(
    mat = mRNA[kDF[,2],colOrder], 
    scale = FALSE,
    limits = limitsRNA,
    color = palRNA, 
    colData = cD[colOrder,,drop=FALSE],
    colorMap = colorMap,
    clusterCols = FALSE,
    clusterRows = FALSE,
    split = kDF[,3],
    labelRows = FALSE,
    labelCols = FALSE,
    # begin of addition by H. Kim
    useRaster = TRUE,
    rasterQuality = 1, # default=5
    rasterDevice = "png",
    show_annotation_name = FALSE,
    annotation_legend_param = list(ncol = legend_ncol),
    # end of addition
    draw = FALSE,
    name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
  )
  
  .endLogging(logFile = logFile)
  
  htATAC + htRNA
  
} # plotPeak2GeneHeatmap.distal











#######################
# ArchRHeatmap.R




# .ArchRHeatmap
.ArchRHeatmap <- function(
  mat = NULL, 
  scale = FALSE,
  limits = c(min(mat), max(mat)),
  colData = NULL, 
  color = paletteContinuous(set = "solarExtra", n = 100),
  clusterCols = TRUE,
  clusterRows = FALSE,
  labelCols = FALSE,
  labelRows = FALSE,
  colorMap = NULL,
  useRaster = TRUE,
  rasterQuality = 5,
  split = NULL,
  fontSizeRows = 10,
  fontSizeCols = 10,
  fontSizeLabels = 8,
  colAnnoPerRow = 4,
  showRowDendrogram = FALSE,
  showColDendrogram = FALSE,
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.75,
  rasterDevice = "png",
  padding = 45,
  borderColor = NA,
  draw = TRUE,
  name = "Heatmap",
  # begin of addition by H. Kim
  show_annotation_name = TRUE,
  annotation_legend_param = NULL,
  heatmap_legend_param = NULL
  # end of addition
  ) {
  
  #Packages
  .requirePackage("ComplexHeatmap", source = "bioc")
  .requirePackage("circlize", source = "cran")
  
  #Z-score
  if (scale) {
    message("Scaling Matrix..")
    mat <- .rowZscores(mat, limit = FALSE)
    name <- paste0(name," Z-Scores")
  }
  
  #Get A Color map if null
  if (is.null(colorMap)) {
    colorMap <- .colorMapAnno(colData)
  }
  
  #Prepare ColorMap format for Complex Heatmap
  if (!is.null(colData)){
    # begin of modification by H. Kim
    #colData = data.frame(colData)
    colData = data.frame(colData, check.names=FALSE)
    # end of modifiation
    colorMap <- .colorMapForCH(colorMap, colData) #change
    showLegend <- .checkShowLegend(colorMap[match(names(colorMap), colnames(colData))]) #change
  }else {
    colorMap <- NULL
    showLegend <- NULL
  }
  
  #Prepare Limits if needed
  breaks <- NULL
  if (!is.null(limits)) {
    mat[mat > max(limits)] <- max(limits)
    mat[mat < min(limits)] <- min(limits)
  }else{
    limits <- c(round(min(mat),2), round(max(mat),2))
  }

  #Scale Values 0 - 1
  mat <- (mat - min(limits)) / (max(limits) - min(limits))
  breaks <- seq(0, 1, length.out = length(color))
  color <- circlize::colorRamp2(breaks, color)

  if(exists('anno_mark', where='package:ComplexHeatmap', mode='function')){
    anno_check_version_rows <- ComplexHeatmap::anno_mark
    anno_check_version_cols <- ComplexHeatmap::anno_mark
  }else{
    anno_check_version_rows <- ComplexHeatmap::row_anno_link
    anno_check_version_cols <- ComplexHeatmap::column_anno_link
  }

  #Annotation Heatmap
  # begin of addition by H. Kim
  if (is.null(annotation_legend_param)) {
        annotation_legend_param <- list( nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1)) )
  } else {
	# custom parameters
  }
  # end of addition
  if(!is.null(colData) & !is.null(customColLabel)){
    message("Adding Annotations..")
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      #show_annotation_name = TRUE,
      show_annotation_name = show_annotation_name,
      gp = gpar(col = "NA"),
      # begin of modification by H. Kim
      #annotation_legend_param =
      #  list(
      #    nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
      #  ),
      annotation_legend_param = annotation_legend_param,
      # end of modification
      foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels))
    )

  } else if(!is.null(colData)){
    message("Adding Annotations..")
    ht1Anno <- HeatmapAnnotation(
      df = colData,
      col = colorMap, 
      show_legend = showLegend,
      #show_annotation_name = TRUE,
      show_annotation_name = show_annotation_name,
      gp = gpar(col = "NA"),
      # begin of modification by H. Kim
      #annotation_legend_param =
      #  list(
      #    nrow = min(colAnnoPerRow, max(round(nrow(colData)/colAnnoPerRow), 1))
      #  )
      annotation_legend_param = annotation_legend_param
      # end of modification
    )
  } else if(is.null(colData) & !is.null(customColLabel)){
    if(is.null(customColLabelIDs)){
      customColLabelIDs <- colnames(mat)[customColLabel]
    }
    message("Adding Annotations..")
    #print(customColLabel)
    #print(customColLabelIDs)
    #ht1Anno <- columnAnnotation(foo = anno_check_version_cols(
    #   at = customColLabel, labels = customColLabelIDs),
    #   width = unit(customLabelWidth, "cm") + max_text_width(customColLabelIDs))
    #ht1Anno <- HeatmapAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 1097:1100), labels = month.name[1:10]))
    ht1Anno <- HeatmapAnnotation(foo = anno_check_version_cols(at = customColLabel, labels = customColLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)))
  }else{
    ht1Anno <- NULL
  }

  # begin of addition by H. Kim
  heatmap_legend_param_org <- list(
           at = c(0, 1),
           labels = c(round(min(limits),2), round(max(limits),2)),
           color_bar = "continuous", 
           legend_direction = "horizontal",
           legend_width = unit(3, "cm")
      )
  idx <- match(names(heatmap_legend_param_org), names(heatmap_legend_param))
  f <- !is.na(idx); idx <- idx[f]
  if (length(idx) > 0) {
	heatmap_legend_param_org[f] <- heatmap_legend_param[idx]	
	heatmap_legend_param[idx] <- NULL
  }
  heatmap_legend_param <- c(heatmap_legend_param_org, heatmap_legend_param)
  # end of addition
  message("Preparing Main Heatmap..")
  ht1 <- Heatmap(
    
    #Main Stuff
    matrix = as.matrix(mat),
    name = name,
    col = color, 
    
    #Heatmap Legend
    # begin of modification by H. Kim
    #heatmap_legend_param = 
    #  list(
    #       at = c(0, 1),
    #       labels = c(round(min(limits),2), round(max(limits),2)),
    #       color_bar = "continuous", 
    #       legend_direction = "horizontal",
    #       legend_width = unit(3, "cm"),
    #  ), 
    heatmap_legend_param = heatmap_legend_param,
    # end of modification

    rect_gp = gpar(col = borderColor), 
    
    #Column Options
    show_column_names = labelCols,
    cluster_columns = clusterCols, 
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = fontSizeCols), 
    column_names_max_height = unit(100, "mm"),
    
    #Row Options
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows), 
    cluster_rows = clusterRows, 
    show_row_dend = showRowDendrogram, 
    clustering_method_rows = "ward.D2",
    split = split, 
    
    #Annotation
    top_annotation = ht1Anno, 

    #Raster Info
    use_raster = useRaster, 
    raster_device = rasterDevice, 
    raster_quality = rasterQuality
  )

  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    ht1 <- ht1 + rowAnnotation(link = 
        anno_check_version_rows(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSizeLabels)),
        width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs))
  }

  if(draw){
    draw(ht1, 
      padding = unit(c(padding, padding, padding, padding), "mm"), 
      heatmap_legend_side = "bot", 
      annotation_legend_side = "bot")
  }else{
    ht1
  }

} # .ArchRHeatmap






