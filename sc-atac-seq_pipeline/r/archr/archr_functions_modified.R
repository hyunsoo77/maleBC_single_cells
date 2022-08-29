#
# archr_functions_modified.R
# modified by H. Kim
# last modified on 2022, Mar.
#
# content:
# IntegrativeAnalysis.R
#	 .constructGR()
# 	getPeak2GeneLinks()
#	plotPeak2GeneHeatmap()
# ArchRBrowser.R
#	plotBrowserTrack()
#	.bulkTracks()
#	.geneTracks
#	.featureTracks()
#	.loopTracks()
# ArchRHeatmap.R
#	.ArchRHeatmap()
#	.binarySort()
# ArrowRead.R
# ReproduciblePeakSet.R
#	.plotPeakCallSummary()


source("./r/archr/ArchR-1.0.1/R/ArrowUtils.R")
source("./r/archr/ArchR-1.0.1/R/HiddenUtils.R")
source("./r/archr/ArchR-1.0.1/R/LoggerUtils.R")
source("./r/archr/ArchR-1.0.1/R/ValidationUtils.R")

source("./r/archr/ArchR-1.0.1/R/ArchRHeatmap.R")




#######################
# IntegrativeAnalysis.R






# .constructGR
.constructGR <- function(
  seqnames = NULL,
  start = NULL,
  end = NULL,
  ignoreStrand = TRUE
){
  .validInput(input = seqnames, name = "seqnames", valid = c("character", "rleCharacter"))
  #.validInput(input = start, name = "start", valid = c("integer"))
  #.validInput(input = end, name = "end", valid = c("integer"))
  .validInput(input = ignoreStrand, name = "ignoreStrand", valid = c("boolean"))
  df <- data.frame(seqnames, start, end)
  idx <- which(df[,2] > df[,3])
  df[idx,2:3] <-  df[idx,3:2]
  if(!ignoreStrand){
    strand <- rep("+",nrow(df))
    strand[idx] <- "-"
  }else{
    strand <- rep("*",nrow(df))
  }
  gr <- GRanges(df[,1], IRanges(df[,2],df[,3]), strand = strand)
  return(gr)


} # .constructGR







# getPeak2GeneLinks
# Make modified getP2G function:
getPeak2GeneLinks <- function(
  ArchRProj = NULL, 
  corCutOff = 0.45, 
  # begin of modification
  #FDRCutOff = 1e-4,
  FDRCutOff = NULL, # default=1e-4
  PvalCutOff = NULL, # default=1e-4
  RawPValCutOff = 1e-12,
  peakName = NULL,
  f_include_negative_cor = FALSE,
  # end of modification
  varCutOffATAC = 0.25,
  varCutOffRNA = 0.25,
  resolution = 1, 
  returnLoops = TRUE
){
  
  .validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProject")
  .validInput(input = corCutOff, name = "corCutOff", valid = "numeric")
  #.validInput(input = PvalCutOff, name = "PvalCutOff", valid = "numeric")
  #.validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = "numeric")
  .validInput(input = varCutOffATAC, name = "varCutOffATAC", valid = c("numeric", "null"))
  #.validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = "numeric")
  .validInput(input = varCutOffRNA, name = "varCutOffRNA", valid = c("numeric", "null"))
  .validInput(input = resolution, name = "resolution", valid = c("integer", "null"))
  .validInput(input = returnLoops, name = "returnLoops", valid = "boolean")
  
  if (is.null(ArchRProj@peakSet)){

    return(NULL)

  }
  


  if (is.null(metadata(ArchRProj@peakSet)$Peak2GeneLinks)){
    
    return(NULL)
    
  } else {
    
    p2g <- metadata(ArchRProj@peakSet)$Peak2GeneLinks

    # begin of modification by H. Kim
    p2g$geneName <- mcols(metadata(p2g)$geneSet)$name[p2g$idxRNA]
    p2g$peakName <- (metadata(p2g)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2g$idxATAC]

    #p2g <- p2g[which(p2g$Correlation >= corCutOff & p2g$FDR <= FDRCutOff), ,drop=FALSE]


    if (f_include_negative_cor) {

	if (!is.null(FDRCutOff)) {
		p2g <- p2g[which((abs(p2g$Correlation) >= corCutOff) & (p2g$FDR <= FDRCutOff)), ,drop=FALSE]
	} else if (!is.null(PvalCutOff)) {
		p2g <- p2g[which((abs(p2g$Correlation) >= corCutOff) & (p2g$Pval <= PvalCutOff)), ,drop=FALSE]
	} else {
		p2g <- p2g[which((abs(p2g$Correlation) >= corCutOff) & (p2g$RawPVal <= RawPValCutOff)), ,drop=FALSE]
	} # if

    } else {

	if (!is.null(FDRCutOff)) {
		p2g <- p2g[which((p2g$Correlation >= corCutOff) & (p2g$FDR <= FDRCutOff)), ,drop=FALSE]
	} else if (!is.null(PvalCutOff)) {
		p2g <- p2g[which((p2g$Correlation >= corCutOff) & (p2g$Pval <= PvalCutOff)), ,drop=FALSE]
	} else {
		p2g <- p2g[which((p2g$Correlation >= corCutOff) & (p2g$RawPVal <= RawPValCutOff)), ,drop=FALSE]
	} # if

    } # if

    if (!is.null(peakName)) {
	idx <- match(p2g$peakName, peakName)
	f <- !is.na(idx)
	p2g <- p2g[which(f),,drop=FALSE]
    }

    # end of modification
    
    if(!is.null(varCutOffATAC)){
      p2g <- p2g[which(p2g$VarQATAC > varCutOffATAC),]
    }
    
    if(!is.null(varCutOffRNA)){
      p2g <- p2g[which(p2g$VarQRNA > varCutOffRNA),]
    }

    if (returnLoops) {
      
      peakSummits <- resize(metadata(p2g)$peakSet, 1, "center")
      geneStarts <- resize(metadata(p2g)$geneSet, 1, "start")
      
      if(!is.null(resolution)){
        summitTiles <- floor(start(peakSummits) / resolution) * resolution + floor(resolution / 2)
        geneTiles <- floor(start(geneStarts) / resolution) * resolution + floor(resolution / 2)
      }else{
        summitTiles <- start(peakSummits)
        geneTiles <- start(geneTiles)
      }
      
      loops <- .constructGR(
        seqnames = seqnames(peakSummits[p2g$idxATAC]),
        start = summitTiles[p2g$idxATAC],
        end = geneTiles[p2g$idxRNA]
      )
      mcols(loops)$value <- p2g$Correlation
      mcols(loops)$FDR <- p2g$FDR
      # begin of modification by H. Kim
      mcols(loops)$geneName <- p2g$geneName
      #mcols(loops)$peakName <- p2g$peakName
      # end of modification
      
      loops <- loops[order(mcols(loops)$value, decreasing=TRUE)]
      loops <- unique(loops)
      loops <- loops[width(loops) > 0]
      loops <- sort(sortSeqlevels(loops))
      
      loops <- SimpleList(Peak2GeneLinks = loops)
      
      return(loops)
      
    }else{
      
      return(p2g)
      
    }
    
  }
  
} # getPeak2GeneLinks











# plotPeak2GeneHeatmap
# This function plots side by side heatmaps of linked ATAC and Gene regions from addPeak2GeneLinks.  https://www.archrproject.com/reference/plotPeak2GeneHeatmap.html
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
  # begin of modification by H. Kim
  pal = NULL,
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


  # begin of modification by H. Kim
  #cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups)
  #pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
  #if(!is.null(palGroup)){
  #  pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
  #}

  if (is.null(pal)) {
	cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = KNNGroups)
	pal <- paletteDiscrete(values=gtools::mixedsort(unique(ccd[,1])))
	if(!is.null(palGroup)){
		pal[names(palGroup)[names(palGroup) %in% names(pal)]] <- palGroup[names(palGroup) %in% names(pal)]
	}
  } else {
	cD <- DataFrame(row.names=paste0("K_", seq_len(ncol(mATAC))), groupBy = factor(KNNGroups, levels=names(pal)))
  } # if
  # end of modification

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
  
  if(returnMatrices){
    
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
    # end of addition
    draw = FALSE,
    name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
  )
  
  .endLogging(logFile = logFile)
  
  htATAC + htRNA
  
} # plotPeak2GeneHeatmap.distal















#######################
# ArchRBrowser.R



####################################################################
# Signal Track Plotting Methods
####################################################################

#' Launch ArchR Genome Browser
#' 
#' This function will open an interactive shiny session in style of a browser track. It allows for normalization of the signal which
#' enables direct comparison across samples.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument
#' can be used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param browserTheme A `shinytheme` from shinythemes for viewing the ArchR Browser. If not installed this will be NULL.
#' To install try devtools::install_github("rstudio/shinythemes").
#' @param threads The number of threads to use for parallel execution.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
ArchRBrowser <- function(
  ArchRProj = NULL,
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  minCells = 25,
  baseSize = 10,
  borderWidth = 0.5,
  tickWidth = 0.5,
  facetbaseSize = 12,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  browserTheme = "cosmo",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("ArchRBrowser")
  ){

  #.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  #.validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  #.validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  #.validInput(input = minCells, name = "minCells", valid = c("integer"))
  #.validInput(input = baseSize, name = "baseSize", valid = c("integer"))
  #.validInput(input = borderWidth, name = "borderWidth", valid = c("numeric"))
  #.validInput(input = tickWidth, name = "tickWidth", valid = c("numeric"))
  #.validInput(input = facetbaseSize, name = "facetbaseSize", valid = c("numeric"))
  #geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  #.validInput(input = browserTheme, name = "browserTheme", valid = c("character"))
  #.validInput(input = threads, name = "threads", valid = c("integer"))
  #.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  #.validInput(input = logFile, name = "logFile", valid = c("character"))

  #.startLogging(logFile=logFile)
  #.logThis(mget(names(formals()),sys.frame(sys.nframe())), "ArchRBrowser Input-Parameters", logFile = logFile)

  #.requirePackage("shiny", installInfo = 'install.packages("shiny")')
  #.requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')

  #Determine Grouping Methods
  ccd <- getCellColData(ArchRProj)
  discreteCols <- lapply(seq_len(ncol(ccd)), function(x){
    .isDiscrete(ccd[, x])
  }) %>% unlist %>% {colnames(ccd)[.]}
  if("Clusters" %in% discreteCols){
    selectCols <- "Clusters"
  }else{
    selectCols <- "Sample"
  }

  #Extend where upstream can be negative for browser
  extendGR2 <-  function(gr = NULL, upstream = NULL, downstream = NULL){
    #.validInput(input = gr, name = "gr", valid = c("GRanges"))
    #.validInput(input = upstream, name = "upstream", valid = c("integer"))
    #.validInput(input = downstream, name = "downstream", valid = c("integer"))
    #Get Info From gr
    st <- start(gr)
    ed <- end(gr)
    #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
    isMinus <- BiocGenerics::which(strand(gr) == "-")
    isOther <- BiocGenerics::which(strand(gr) != "-")
    #Forward
    st[isOther] <- st[isOther] - upstream
    ed[isOther] <- ed[isOther] + downstream
    #Reverse
    ed[isMinus] <- ed[isMinus] + upstream
    st[isMinus] <- st[isMinus] - downstream
    #If Any extensions now need to be flipped.
    end(gr) <- pmax(st, ed)
    start(gr) <- pmin(st, ed)
    return(gr)
  }


  #####################
  #Shiny App UI
  #####################
  if(!requireNamespace("shinythemes", quietly = TRUE)){
    #.logMessage("shinythemes not found! To see a nice theme use :\n\tinstall.packages('shinythemes')\nContinuing wihtout shinythemes!", verbose = verbose, logFile = logFile)
    theme <- NULL
  }else{
    theme <- shinythemes::shinytheme(browserTheme)
  }

  ui <- fluidPage(
    theme = theme,
    titlePanel(
        h1(div(HTML(paste0("<b>ArchR Browser v1 : nCells = ", formatC(nCells(ArchRProj), format="f", big.mark = ",", digits=0), "</b>"))), align = "left")
    ),
    sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style(HTML("hr {border-top: 1px solid #000000;}"))
        ),
        tags$head(
          tags$style(HTML('#exitButton{background-color:#D60000}'))
        ), 
        tags$head(
          tags$style(HTML('#restartButton{background-color:#02A302}'))
        ),
        tags$head(
          tags$style(HTML('#plot_height{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#plot_width{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#ymax{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#tile_size{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_min{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_max{height: 35px}'))
        ),
        actionButton(inputId = "exitButton", label = "Exit Browser", icon = icon("times-circle")),
        br(),
        br(),
        actionButton(inputId = "restartButton", label = "Plot Track!", icon = icon("play-circle")),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("exitButton", label = "Exit Browser", icon = icon("paper-plane"))),
        #br(),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("restartButton", label = "Plot Track!", icon = icon("paper-plane"))),
        br(),
        br(),
        selectizeInput("name",
                       label = "Gene Symbol",
                       choices = as.vector(geneAnnotation$genes$symbol)[!is.na(as.vector(geneAnnotation$genes$symbol))],
                       multiple = FALSE,
                       options = list(placeholder = 'Select a Center'),
                       selected = "CD4"
                      ),
        selectizeInput("grouping",
           label = "groupBy",
           choices = discreteCols,
           multiple = FALSE,
           options = list(placeholder = 'Select Grouping'),
           selected = selectCols
        ),
        sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
          numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
        ),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
          numericInput("ymax", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
        ),
        hr(),
        downloadButton(outputId = "down", label = "Download the Track!"),
        br(),
        br(),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
          numericInput("plot_height", "Height", min = 0, max = 250, value = 12)
        ),
        width = 2, height = "750px", position = "left"
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Plot",
            plotOutput(outputId = "ATAC", width= "800px", height = "725px")
          ),
          tabPanel("Additional Params",
            br(),
            br(),
            selectizeInput("normATAC",
               label = "normMethod",
               choices = c("ReadsInTSS", "ReadsInPromoter", "nFrags", "None"),
               multiple = FALSE,
               options = list(placeholder = 'Select NormMethod'),
               selected = "ReadsInTSS"
            ),
            br(),
            div(HTML("<b>Group Metadata</b>")),
            rHandsontableOutput("Metadata", width= "1085px",height = "800px")
          )
        )
      )
    )
  )

  #####################
  #Shiny App Server
  #####################
  server <- function(
        input = input, 
        output = output, 
        session = session
    ){

    output$Metadata <- renderRHandsontable({
        groups <- gtools::mixedsort(unique(ccd[,input$grouping]))
        mdata <- data.frame(
          groupBy = input$grouping,
          include = rep(TRUE,length(groups)), 
          group = groups, 
          color = paletteDiscrete(values = groups)[groups], 
          nCells = as.vector(table(ccd[,input$grouping])[groups]),
          medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
          medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
          stringsAsFactors = FALSE
        )
        rownames(mdata) <- NULL
        rhandsontable(mdata)
      })

    #Update Sliders
    observeEvent(input$range_min, {
      updateSliderInput(session, "range",
                        value = c(input$range_min,max(input$range)))
    })
    
    observeEvent(input$range_max, {
      updateSliderInput(session, "range",
                        value = c(input$range_min,input$range_max))
    })

    observeEvent(input$range , {

      updateNumericInput(session, "range_min", value = min(input$range))
      updateNumericInput(session, "range_max", value = max(input$range))

    }, priority = 200)

    output$checkbox <- renderUI({
      choice <- gtools::mixedsort(unique(ccd[,input$grouping,drop=TRUE]))
      checkboxGroupInput("checkbox","Select Groups", choices = choice, selected = choice)
    })    

    #################################
    # Inputs that cause re-plotting
    #################################
    #toListen <- reactive({
    #  list(input$restartButton, input$name)
    #})

    restartFN <- observeEvent(input$restartButton, {

      if (input$name == ""){
          
          output$ATAC <- renderPlot({
            p <- ggplot() +
                xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
            print(p)
          })

      }else{
        output$ATAC <- renderPlot({

          withProgress(message = 'Plotting', style = "notification", value = 0, {

            #Get Region if Gene Symbol
            region <- geneAnnotation$genes

            if(tolower(input$name) %ni% tolower(mcols(region)$symbol)){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
              return(print(p))
            }

            region <- region[which(tolower(mcols(region)$symbol) %in% tolower(input$name))]
            region <- region[order(match(tolower(mcols(region)$symbol), tolower(input$name)))]
            region1 <- resize(region, 1, "start")
            strand(region1) <- "*"

            #Extend Region
            #region <- extendGR(region, upstream = -min(input$range) * 1000, downstream = max(input$range) * 1000)
            #Pre-Load full window for even faster plotting
            region <- extendGR2(region1, upstream = 250000, downstream = 250000)
            tmpArchRRegion <<- extendGR2(region1, 
              upstream = -min(isolate(input$range)) * 1000, 
              downstream = max(isolate(input$range)) * 1000
            )
            region <- tmpArchRRegion

            setProgress(0.1)

            #User Inputs
            groupBy <- isolate(input$grouping)

            groupDF <- tryCatch({
              isolate(hot_to_r(input$Metadata))
            },error=function(x){
              groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
              mdata <- data.frame(
                groupBy = input$grouping,
                include = rep(TRUE,length(groups)), 
                group = groups, 
                color = paletteDiscrete(values = groups)[groups], 
                nCells = as.vector(table(ccd[,input$grouping])[groups]),
                medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
                medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
                stringsAsFactors = FALSE
              )
              rownames(mdata) <- NULL
              mdata
            })

            if(groupDF$groupBy[1] != groupBy){
              groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
              groupDF <- data.frame(
                groupBy = groupBy,
                include = rep(TRUE,length(groups)), 
                group = groups, 
                color = paletteDiscrete(values = groups)[groups], 
                stringsAsFactors = FALSE
              )
              rownames(groupDF) <- NULL
            }

            useGroups <- groupDF[groupDF[,"include"],"group"]


            if(!all(.isColor(groupDF[groupDF[,"include"], "color"]))){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
              return(print(p))               
            }

            if(any(useGroups %ni% ccd[, groupBy])){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
              return(print(p)) 
            }

            pal <- groupDF[groupDF[,"include"], "color"]
            names(pal) <- groupDF[groupDF[,"include"], "group"]

            ylim <- c(0, isolate(input$ymax))
            normMethod <- isolate(input$normATAC)
            tileSize <- isolate(input$tile_size)

            p <- .bulkTracks(
                ArchRProj = ArchRProj, 
                region = region, 
                tileSize = tileSize, 
                useGroups = useGroups,
                groupBy = groupBy,
                threads = threads, 
                minCells = minCells,
                ylim = ylim,
                baseSize = baseSize,
                borderWidth = borderWidth,
                tickWidth = tickWidth,
                facetbaseSize = facetbaseSize,
                normMethod = normMethod,
                geneAnnotation = geneAnnotation,
                title = "",
                pal = pal, 
                tstart = NULL,
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))

            #p <- p + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            tmpArchRP <<- p

            setProgress(0.5)

            if(!is.null(features)){

              f <- .featureTracks(
                  features = features, 
                  region = tmpArchRRegion,
                  facetbaseSize = facetbaseSize, 
                  hideX = TRUE, 
                  title = "PeakSets",
                  logFile = logFile
                ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

              #f <- f + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            }

            if(!is.null(loops)){

              l <- .loopTracks(
                loops = loops, 
                region = tmpArchRRegion, 
                facetbaseSize = facetbaseSize,
                hideX = TRUE, 
                hideY = TRUE,
                title = "Loops",
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
            }

            setProgress(0.6)

            g <- .geneTracks(
              geneAnnotation = geneAnnotation, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              labelSize = 3,
              title = "Genes",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

            #g <- .suppressAll(g + scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            setProgress(0.8)

            if(!is.null(loops)){
              if(!is.null(features)){
                suppressWarnings(print(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = TRUE)))
              }else{
                suppressWarnings(print(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = TRUE)))
              }
            }else{
              if(!is.null(features)){
                suppressWarnings(print(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = TRUE)))
              }else{
                suppressWarnings(print(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = TRUE)))
              }
            }

            setProgress(1)

          })

        })
      }    

    })

    #######################################
    # When Download Is Initiated
    #######################################

    # downloadHandler contains 2 arguments as functions, namely filename, content
    output$down <- downloadHandler(

      filename <- function(){
        paste0("ArchRBrowser-",input$name,"-",seqnames(tmpArchRRegion)[1],":",start(tmpArchRRegion)[1],"-",end(tmpArchRRegion)[1],".pdf")
      },

      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        
        withProgress(message = 'Creating PDF', style = "notification", value = 0, {
         
          if(!exists("tmpArchRP")){

            #User Inputs
            groupBy <- isolate(input$grouping)
            groupDF <- isolate(hot_to_r(input$Metadata))
            useGroups <- groupDF[groupDF[,"include"],"group"]

            .isColor <- function(x = NULL) {
                unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), 
                  error = function(e) FALSE)))
            }

            if(!all(.isColor(groupDF[groupDF[,"include"], "color"]))){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
              return(print(p))               
            }

            if(any(useGroups %ni% ccd[, groupBy])){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
              return(print(p)) 
            }

            pal <- groupDF[groupDF[,"include"], "color"]
            names(pal) <- groupDF[groupDF[,"include"], "group"]

            ylim <- c(0, isolate(input$ymax))
            normMethod <- isolate(input$normATAC)
            tileSize <- isolate(input$tile_size)

            p <- .bulkTracks(
                ArchRProj = ArchRProj, 
                region = tmpArchRRegion, 
                tileSize = tileSize, 
                useGroups = useGroups,
                groupBy = groupBy,
                threads = threads, 
                minCells = minCells,
                ylim = ylim,
                baseSize = baseSize,
                borderWidth = borderWidth,
                tickWidth = tickWidth,
                facetbaseSize = facetbaseSize,
                normMethod = normMethod,
                geneAnnotation = geneAnnotation,
                title = "",
                pal = pal, 
                tstart = NULL,
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))

          }else{

            print("Using previous ggplot")

            p <- tmpArchRP

          }

          setProgress(0.5)

          if(!is.null(features)){

            f <- .featureTracks(
                features = features, 
                region = tmpArchRRegion,
                facetbaseSize = facetbaseSize, 
                hideX = TRUE, 
                title = "PeakSets",
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

          }

          setProgress(0.6)

          if(!is.null(loops)){

            l <- .loopTracks(
              loops = loops, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              hideX = TRUE, 
              hideY = TRUE,
              title = "Loops",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          }

          setProgress(0.7)

          g <- .geneTracks(
            geneAnnotation = geneAnnotation, 
            region = tmpArchRRegion, 
            facetbaseSize = facetbaseSize,
            labelSize = 3,
            title = "Genes",
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

          setProgress(0.8)

          pdf(file = file, width = input$plot_width, height = input$plot_height)
          
          if(!is.null(loops)){
            if(!is.null(features)){
              a <- suppressWarnings(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = FALSE))
            }else{
              a <- suppressWarnings(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = FALSE))
            }
          }else{
            if(!is.null(features)){
              a <- suppressWarnings(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = FALSE))
            }else{
              a <- suppressWarnings(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = FALSE))
            }
          }
          
          suppressWarnings(grid::grid.draw(a))
          dev.off()

          setProgress(1)

        })

    })


    exitFN <- observeEvent(input$exitButton, {
      if(exists("tmpArchRRegion")){
        .suppressAll(rm(tmpArchRRegion))
      }
      if(exists("tmpArchRP")){
        .suppressAll(rm(tmpArchRP))
      }
      shiny::stopApp()
    })

  }

  shiny::runGadget(ui, server)

}

#' @export
ArchRBrowserTrack <- function(...){
    .Deprecated("plotBrowserTrack")
    plotBrowserTrack(...)
}


















#' Plot an ArchR Region Track
#' 
#' This function will plot the coverage at an input region in the style of a browser track. It allows for normalization of the signal
#' which enables direct comparison across samples.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param region A `GRanges` region that indicates the region to be plotted. If more than one region exists in the `GRanges` object,
#' all will be plotted. If no region is supplied, then the `geneSymbol` argument can be used to center the plot window at the
#' transcription start site of the supplied gene.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in
#' `cellColData`. This limits the groups to be plotted.
#' @param plotSummary A character vector containing the features to be potted. Possible values include "bulkTrack" (the ATAC-seq signal),
#' "scTrack" (scATAC-seq signal), "featureTrack" (i.e. the peak regions), "geneTrack" (line diagrams of genes with introns and exons shown. 
#' Blue-colored genes are on the minus strand and red-colored genes are on the plus strand), and "loopTrack" (links between a peak and a gene).
#' @param sizes A numeric vector containing up to 3 values that indicate the sizes of the individual components passed to `plotSummary`.
#' The order must be the same as `plotSummary`.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param geneSymbol If `region` is not supplied, plotting can be centered at the transcription start site corresponding to the gene symbol(s) passed here.
#' @param useMatrix If supplied geneSymbol, one can plot the corresponding GeneScores/GeneExpression within this matrix. I.E. "GeneScoreMatrix"
#' @param log2Norm If supplied geneSymbol, Log2 normalize the corresponding GeneScores/GeneExpression matrix before plotting.
#' @param upstream The number of basepairs upstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param downstream The number of basepairs downstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument can be
#' used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality.
#' @param threads The number of threads to use for parallel execution.
#' @param ylim The numeric quantile y-axis limit to be used for for "bulkTrack" plotting. If not provided, the y-axis limit will be c(0, 0.999).
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override coloring for groups.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param scTileSize The width of the tiles in scTracks. Larger numbers may make cells overlap more. Default is 0.5 for about 100 cells.
#' @param scCellsMax The maximum number of cells for scTracks.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param title The title to add at the top of the plot next to the plot's genomic coordinates.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
plotBrowserTrack <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  # begin of addition by H. Kim
  pal_strip = NULL,
  bulktracks_scale=NULL,
  bulktracks_hjust=NULL,
  bulktracks_vjust=0,
  bulktracks_angle=0,
  # end of addition
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
  ) {
  
  #.validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  #.validInput(input = region, name = "region", valid = c("granges","null"))
  #.validInput(input = groupBy, name = "groupBy", valid = "character")
  #.validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  #.validInput(input = plotSummary, name = "plotSummary", valid = "character")
  #.validInput(input = sizes, name = "sizes", valid = "numeric")
  #.validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  #.validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  #.validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  #.validInput(input = useMatrix, name = "useMatrix", valid = c("character", "null"))
  #.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  #.validInput(input = upstream, name = "upstream", valid = c("integer"))
  #.validInput(input = downstream, name = "downstream", valid = c("integer"))
  #.validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  #.validInput(input = minCells, name = "minCells", valid = c("integer"))
  #.validInput(input = normMethod, name = "normMethod", valid = c("character"))
  #.validInput(input = threads, name = "threads", valid = c("integer"))
  #.validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  #.validInput(input = pal, name = "pal", valid = c("palette", "null"))
  #.validInput(input = baseSize, name = "baseSize", valid = "numeric")
  #.validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
  #.validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
  #.validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  #.validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  #.validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  #geneAnnotation <- .validGeneAnnotation(geneAnnotation)
  #.validInput(input = title, name = "title", valid = "character")

  tstart <- Sys.time()
  #.startLogging(logFile=logFile)
  #.logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotBrowserTrack Input-Parameters", logFile = logFile)

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  #.logDiffTime("Validating Region", t1=tstart, verbose=verbose, logFile=logFile)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      #print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  #region <- .validGRanges(region)
  #.logThis(region, "region", logFile = logFile)

 
  # begin of addition by H. Kim
  f_loop_inside_region <- FALSE
  for (loop in names(loops)) {
	gr_loop <- loops[[loop]]
	if (length(gr_loop) > 0) {
		if (any(gr_loop %within% region)) {
			f_loop_inside_region <- TRUE
		}
	}
  } # for
  # end of addition

  



  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }

  if(!is.null(useMatrix)){
    featureMat <- .getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }



  # for multiple region
  ggList <- lapply(seq_along(region), function(x){

    plotList <- list()


    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){

      #.logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = minCells,
        pal = pal,
	# begin of addition by H. Kim
	pal_strip = pal_strip,
	f_loop_inside_region = f_loop_inside_region,
	bulktracks_scale = bulktracks_scale,
	bulktracks_hjust = bulktracks_hjust,
	bulktracks_vjust = bulktracks_vjust,
	bulktracks_angle = bulktracks_angle,
	# end of addition
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }

    
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){

      #.logDiffTime(sprintf("Adding SC Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$sctrack <- .scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){

      if(!is.null(features)){
        #.logDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$featuretrack <- .featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "PeakSets",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){

      if(!is.null(loops)){
        #.logDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$looptrack <- .loopTracks(
            loops = loops, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){

      #.logDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }

    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]

    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]

    if(!is.null(useMatrix)){

      suppressWarnings(.combinedFeaturePlot(
        plotList = plotList,
        log2Norm = log2Norm,
        featureMat = featureMat,
        feature = region[x]$symbol[[1]],
        useMatrix = useMatrix,
        pal = pal,
        sizes = sizes,
        baseSize = baseSize,
        facetbaseSize = facetbaseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth
      ))

    }else{

      #.logThis(names(plotList), sprintf("(%s of %s) names(plotList)",x,length(region)), logFile=logFile)
      #.logThis(sizes, sprintf("(%s of %s) sizes",x,length(region)), logFile=logFile)
      #.logThis(nullSummary, sprintf("(%s of %s) nullSummary",x,length(region)), logFile=logFile)
      #.logDiffTime("Plotting", t1=tstart, verbose=verbose, logFile=logFile)
      
      tryCatch({
        suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))
      }, error = function(e){
        #.logMessage("Error with plotting, diagnosing each element", verbose = TRUE, logFile = logFile)
        for(i in seq_along(plotList)){
          tryCatch({
            print(plotList[[i]])
          }, error = function(f){
            #.logError(f, fn = names(plotList)[i], info = "", errorList = NULL, logFile = logFile)
          })
        }
        #.logError(e, fn = "ggAlignPlots", info = "", errorList = NULL, logFile = logFile)
      })

    }

  })

  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

  #.endLogging(logFile=logFile)

  ggList

} # plotBrowserTrack



    















#######################################################
# Bulk Aggregated ATAC Track Methods
#######################################################
.bulkTracks <- function(
  ArchRProj = NULL, 
  region = NULL, 
  tileSize = 100, 
  minCells = 25,
  groupBy = "Clusters",
  useGroups = NULL,
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  # begin of addition by H. Kim
  pal_strip = NULL,
  f_loop_inside_region = TRUE,
  bulktracks_scale=NULL,
  bulktracks_hjust=NULL,
  bulktracks_vjust=0,
  bulktracks_angle=0,
  # end of addition
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ) {

  #.requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  df <- .groupRegionSumArrows(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    normMethod = normMethod,
    useGroups = useGroups,
    minCells = minCells,
    region = region, 
    tileSize = tileSize, 
    threads = threads,
    verbose = verbose,
    logFile = logFile
  )

  #.logThis(split(df, df[,3]), ".bulkTracks df", logFile = logFile)

  ######################################################
  # Plot Track
  ######################################################
  if(!is.null(ylim)){
    ylim <- quantile(df$y, ylim)
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }else{
    ylim <- c(0,quantile(df$y, probs=c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }

  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
  }
  df$group <- factor(df$group, levels = uniqueGroups)
  # begin of modification by H. Kim
  #title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  title <- paste0(as.character(seqnames(region)),":", format(start(region)-1, big.mark=",", scientific=F), "-", format(end(region), big.mark=",", scientific=F), " ", title)
  # end of modification

  allGroups <- gtools::mixedsort(unique(getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)))

  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = allGroups))
  }
  
  # begin of addition by H. Kim
  levels(df$group) <- names(pal)
  f <- (df$x >= start(region) & df$x <= end(region))
  df <- df[f,]
  # end of addition

  # Plot Track
  p <- ggplot(df, aes_string("x", "y", color = "group", fill = "group")) + 
    geom_area(stat = "identity") + 
    facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
    ylab(sprintf("Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    # begin of modification by H. Kim
    #scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
    scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE), limits = c(start(region), end(region)), expand = c(0,0)) +
    # end of modification
    scale_y_continuous(limits = ylim, expand = c(0,0)) +
    theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) +
    #guides(fill = FALSE, colour = FALSE) + ggtitle(title)
    guides(fill = "none", colour = "none") + ggtitle(title)


  # begin of modification by H. Kim
  # p
  g <- ggplot_gtable( ggplot_build(p) )
  stripr <- which( grepl( "strip-r", g$layout$name ))
  k <- 1
  for ( i in stripr ) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- pal_strip[k]
    k <- k + 1
  }

  if (f_loop_inside_region) {
  	if (is.null(bulktracks_scale)) {
		# 3 peak sets and 2 loop tracks
		bulktracks_scale=1.21
		bulktracks_hjust=0.082
		bulktracks_vjust=0
		bulktracks_angle=0
	}
  } else {
		# 1 peak set
		#bulktracks_scale=1.167
		#bulktracks_hjust=0.065
		# 3 peak sets
		bulktracks_scale=1.181
		bulktracks_hjust=0.071
		bulktracks_vjust=0
		bulktracks_angle=0
  } # if

  ggplotify::as.ggplot(g, scale=bulktracks_scale,
		 hjust=bulktracks_hjust,
		 vjust=bulktracks_vjust,
		 angle=bulktracks_angle)

  # end of modification




} # .bulkTracks


















################################################################
# Create Average Tracks from Arrows
################################################################
.groupRegionSumArrows <- function(
  ArchRProj = NULL,
  useGroups = NULL,
  groupBy = NULL,
  region = NULL,
  tileSize = NULL,
  normMethod = NULL,
  verbose = FALSE,
  minCells = 25,
  maxCells = 500,
  threads = NULL,
  logFile = NULL
  ){

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    #.logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      #print(errorList)
      #.logError(e, fn = ".groupRegionSumArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("+" , .)

  #Plot DF
  df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
  df$y <- groupMat[cbind(df[,1], df[,2])]

  #Minus 1 Tile Size
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0

  #Plus 1 Size
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0

  #Create plot DF
  df <- rbind(df, dfm1, dfp1)
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0,]
  df$x <- regionTiles[df$row]
  df$group <- uniqueGroups[df$col]

  #Add In Ends
  dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
  )

  dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
  )
  
  #Final output
  plotDF <- rbind(df,dfs,dfe)
  plotDF <- df[order(df$group,df$x),]
  plotDF <- df[,c("x", "y", "group")]
  
  #Normalization 
  g <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "readsinpromoter"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else if(tolower(normMethod) == "none"){
      groupNormFactors <- rep(10^4, length(g))
      names(groupNormFactors) <- g
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }

  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors
  matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
  plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])

  return(plotDF)

}

.regionSumArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)

  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = c(matchID[ids], matchID[ide]),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  #Create Group Matrix
  groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
  colnames(groupMat) <- uniqueGroups
  uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
  for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
  }

  return(groupMat)

}












#######################################################
# Gene Tracks
#######################################################
.geneTracks <- function(
  geneAnnotation = NULL, 
  region = NULL, 
  baseSize = 9, 
  borderWidth = 0.4, 
  title = "Genes",
  geneWidth = 2, 
  exonWidth = 4, 
  labelSize = 2,
  facetbaseSize,
  colorMinus = "dodgerblue2",
  colorPlus = "red",
  logFile = NULL
  ){

  #.requirePackage("ggplot2", source = "cran")
  #.requirePackage("ggrepel", source = "cran")

  #only take first region
  #region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  genes <- sort(sortSeqlevels(geneAnnotation$genes), ignore.strand = TRUE)
  exons <- sort(sortSeqlevels(geneAnnotation$exons), ignore.strand = TRUE)
  genesO <- data.frame(subsetByOverlaps(genes, region, ignore.strand = TRUE))

  if(nrow(genesO) > 0){

    #Identify Info for Exons and Genes
    exonsO <- data.frame(subsetByOverlaps(exons, region, ignore.strand = TRUE))
    exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol),]
    genesO$facet = title
    genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, start(region)))
    genesO$end <- matrixStats::rowMins(cbind(genesO$end, end(region)))

    #Collapse Iteratively
    #backwards iteration so that the last value chosen is the lowest cluster possible to fit in.
    genesO$cluster <- 0
    for(i in seq_len(nrow(genesO))){
      if(i==1){
        genesO$cluster[i] <- 1
      }else{
        for(j in seq_len(max(genesO$cluster))){
          jEnd <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], rev(genesO$cluster))]
          if(genesO$start[i] > jEnd + median(genesO$width)){
            genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j],rev(genesO$cluster))]
          }
        }
        if(genesO$cluster[i]==0){
          genesO$cluster[i] <- genesO$cluster[i-1] + 1
        }
      }
    }
    exonsO$cluster <- genesO$cluster[match(exonsO$symbol, genesO$symbol)]
    pal <- c("-"=colorMinus,"+"=colorPlus,"*"=colorPlus)
    
    p <- ggplot(data = genesO, aes(color = strand, fill = strand)) +
      facet_grid(facet~.) +
      #################################################
      #Limits
      #################################################
      ylim(c(0.5, max(genesO$cluster) + 0.5)) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
      #################################################
      #Segment for Not Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)!="-"),], 
        aes(x = start, xend = end, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segment for Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)=="-"),], 
        aes(x = end, xend = start, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segement for Exons
      #################################################
      geom_segment(data = exonsO, aes(x = start, xend = end, y = cluster, 
        yend = cluster, color = strand),size=exonWidth) +
      #################################################
      #Colors
      #################################################
      scale_color_manual(values = pal, guide = FALSE) + 
      scale_fill_manual(values = pal) +
      #################################################
      #Theme
      #################################################
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(size = facetbaseSize, angle = 0)) +
      #guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = FALSE) + 
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = "none") + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand!="-"),], 
	# begin of modification by H. Kim
        #aes(x = start, y = cluster, label = symbol, color = strand, fill = NA), 
        aes(x = start, y = cluster, label = symbol, color = strand), 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          #size = labelSize, direction = "x")
          size = labelSize, direction = "x", inherit.aes=FALSE)
	# end of modification
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand=="-"),], 
	# begin of modification by H. Kim
        #aes(x = end, y = cluster, label = symbol, color = strand, fill = NA), 
        aes(x = end, y = cluster, label = symbol, color = strand), 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          #size = labelSize, direction = "x")
          size = labelSize, direction = "x", inherit.aes=FALSE)
	# end of modification
    }

    p <- p + theme(legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA), legend.position="none")

  }else{

    #create empty plot
    df <- data.frame(facet = "GeneTrack", start = 0, end = 0, strand = "*", symbol = "none")
    pal <- c("*"=colorPlus)
    p <- ggplot(data = df, aes(start, end, fill = strand)) + geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(!is.ggplot(p)){
    #.logError("geneTrack is not a ggplot!", fn = ".geneTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

} # .geneTracks















#######################################################
# Feature Tracks
#######################################################
.featureTracks <- function(
  features = NULL, 
  region = NULL, 
  title = "FeatureTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = NULL,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  #.requirePackage("ggplot2", source = "cran")

  #only take first region
  #region <- .validGRanges(region)

  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!.isGRList(features)){
      #features <- .validGRanges(features)
      featureList <- SimpleList(FeatureTrack = features)
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }
    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
      mcols(featurex) <- NULL
      sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
      if(length(sub) > 0){
        data.frame(sub, name = namex)
      }else{
        empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
        data.frame(empty, name = namex)
      }

    })

    featureO <- Reduce("rbind", featureO)
    
    #.logThis(featureO, "featureO", logFile = logFile)

    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }
    
    featureO$name <- factor(paste0(featureO$name), levels=names(featureList))

    # begin of addition by H. Kim
    baseSize_grlist <- baseSize
    if (.isGRList(features)) {
	baseSize_grlist <- baseSize_grlist - 2
    }
    # end of addition

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme(legend.text = element_text(size = baseSize)) + 
      #theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme_ArchR(baseSize = baseSize_grlist, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      #guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())
      guides(color = "none", fill = "none") + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

  }else{

    #create empty plot
    df <- data.frame(facet = "FeatureTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    #.logError("featureTrack is not a ggplot!", fn = ".featureTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

} # .featureTracks



























#######################################################
# Loop Tracks
#######################################################
.loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){





  getArchDF <- function(lp, r = 100){

    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx

    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }

    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)

    return(df)

  } # getArchDF






  if(!is.null(loops)){

    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    } else if(.isGRList(loops)){

    } else {
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }

    # begin of modification by H. Kim
    #valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    #valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))

    valueMin <- min(unlist(lapply(loops, function(x) {
		if (length(x) == 0) return(0.45)
		min(x$value)
	})))
    valueMax <- max(unlist(lapply(loops, function(x) {
		if (length(x) == 0) return(0.45)
		max(x$value)
	})))
    valueMin_local <- min(unlist(lapply(loops, function(x) {
		if (length(x) == 0) return(0.45)
		f <- x %within% region
		if (!any(f)) return(0.45)
		min(x[x %within% region]$value)
	 })))
    valueMax_local <- max(unlist(lapply(loops, function(x) {
		if (length(x) == 0) return(0.45)
		f <- x %within% region
		if (!any(f)) return(0.45)
		max(x[x %within% region]$value)
	 })))
    if (valueMin_local > 0) {
		valueMin <- 0.45
    } else if (valueMin_local < 0) {
		abs_max <- max(abs(valueMin_local), abs(valueMax_local))
		valueMin <- -abs_max
		valueMax <- abs_max
      		if (is.null(pal)) {
			cs <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")
        		pal <- colorRampPalette(c(cs[2], cs[3], cs[4], "#E6E7E8","#3A97FF","#8816A7","black"))(100)
      		}
    } # if
    # end of modification


    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 

       if(length(subLoops)>0){
         dfx <- getArchDF(subLoops)
         dfx$name <- Rle(paste0(names(loops)[x]))
	 # begin of modification by H. Kim
         #dfx
         as.data.frame(dfx)
	 # end of modification
       } else {
         NULL
       }
    }) %>% Reduce("rbind",.)
    #.logThis(loopO, "loopO", logFile = logFile)

 
    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if (testDim) {

      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
      }

      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line() +
        facet_grid(name ~ .) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        theme(legend.text = element_text(size = baseSize)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3))

    } else {

      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        # begin of addition by H. Kim
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        # end of addition
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

    }

  } else {

    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      # begin of addition by H. Kim
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
      # end of addition

      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  } # if



  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    #.logError("loopTracks is not a ggplot!", fn = ".loopTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

} # .loopTracks



















# .subsetSeqnamesGR
.subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  #.validInput(input = gr, name = "gr", valid = c("GRanges"))
  #.validInput(input = names, name = "names", valid = c("character"))
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))

  return(gr)

} # .subsetSeqnamesGR











#######################################################
# scATAC Track Methods
#######################################################

.scTracks <- function(
  ArchRProj = NULL,
  region = NULL,
  tileSize = 100,
  minCells = 5,
  maxCells = 100,
  groupBy = "Clusters",
  useGroups = NULL,
  threads = 1,
  baseSize = 7,
  scTileSize = 0.5,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  #.requirePackage("ggplot2", source = "cran")

  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    #.logMessage(sprintf("Getting Region From Arrow Files %s of %s", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSCArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
      #.logError(e, fn = ".groupRegionSCArrows", info = .sampleName(ArrowFiles[i]), errorList = errorList, logFile = logFile)
    })
  }, threads = threads) %>% Reduce("cbind" , .)

  groupDF <- data.frame(Matrix::summary(groupMat))
  groupDF$group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[colnames(groupMat)[groupDF$j], 1]
  groupDF <- lapply(split(groupDF, groupDF$group), function(z){
    nz <- tabGroups[z$group[1]]
    nc <- length(unique(z$j))
    idx <- sort(sample(seq_len(nz), nc))
    idx[1] <- 1
    idx[length(idx)] <- nz
    z$y <- idx[match(z$j, unique(z$j))]
    z
  }) %>% Reduce("rbind", .)
  groupDF$bp <- regionTiles[groupDF$i]
  
  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = names(tabGroups)))
  }

  nn <- paste0(names(tabGroups), ":", tabGroups)
  names(nn) <- names(tabGroups)
  groupDF$group2 <- nn[groupDF$group]
  names(pal) <- nn[names(pal)]

  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  
  #Re-Order
  groupDF$group2 <- factor(
    paste0(groupDF$group2), 
    levels = gtools::mixedsort(unique(paste0(groupDF$group2)))
  )

  p <- ggplot(groupDF, aes(x=bp, y=y, width = tileSize, fill = group2, color = group2)) + 
      geom_tile(size = scTileSize) + 
      facet_grid(group2 ~ ., scales="free_y") + 
      theme_ArchR() + 
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      ylab("Binarized SC Coverage") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                  baseRectSize = borderWidth,
                  baseLineSize = tickWidth,
                  legendPosition = "right",
                  axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text = element_text(
              size = facetbaseSize, 
              color = "black", 
              margin = margin(0,0.35,0,0.35, "cm")),
              strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color="black")) +
      #guides(fill = FALSE, colour = FALSE) + ggtitle(title)
      guides(fill = "none", colour = "none") + ggtitle(title)

    p

}

.regionSCArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = as.vector(c(matchID[ids], matchID[ide])),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  return(mat)

}



####################################
# Combined Feature Plot
####################################

.combinedFeaturePlot <- function(
  plotList = NULL,
  useMatrix = NULL,
  featureMat = NULL,
  log2Norm = FALSE,
  feature = NULL,
  pal = NULL,
  sizes = NULL,
  baseSize = NULL,
  facetbaseSize = NULL,
  borderWidth = NULL,
  tickWidth = NULL
  ){

  #.requirePackage("patchwork", installInfo = "devtools::install_github('thomasp85/patchwork')")

  if(is.null(pal)){
    pal <- paletteDiscrete(values=featureMat$Group, set = "stallion")
  }

  if(log2Norm){
    title <- paste0("Log2 ", useMatrix, " : ", feature)
  }else{
    title <- paste0("Raw ", useMatrix, " : ", feature) 
  }

  featurePlot <- ggGroup(
      x = featureMat$Group,
      y = featureMat[,feature],
      groupOrder = gtools::mixedsort(paste0(unique(featureMat$Group))),
      pal = pal
    ) + 
    facet_wrap(x~., ncol=1,scales="free_y",strip.position="right") +
    #guides(fill = FALSE, colour = FALSE) +
    guides(fill = "none", colour = "none") +
    theme_ArchR(baseSize = baseSize,
              baseRectSize = borderWidth,
              baseLineSize = tickWidth,
              legendPosition = "right",
              axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(
          size = facetbaseSize, 
          color = "black", 
          margin = margin(0,0.35,0,0.35, "cm")),
          strip.text.y = element_text(angle = 0),
        strip.background = element_rect(color="black")) +
    theme(plot.margin = unit(c(0.35, 0.15, 0.35, 0.15), "cm")) +
    ggtitle(title)

  if(any(tolower(names(plotList)) %in% "bulktrack")){

    idx <- which(tolower(names(plotList)) == "bulktrack")
    
    p <- plotList[[idx]] + featurePlot + plot_spacer()
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }else{


    idx <- which(tolower(names(plotList)) == "sctrack")
    
    p <- plotList[[idx]] + featurePlot + plot_spacer()
    
    plotList[idx] <- NULL
    
    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }
    
    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2), 
      heights = sizes
    )

  }

  p

} # .combinedFeaturePlot












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








# .checkShowLegend
# comment: 
# no modification, just for debugging
.checkShowLegend <- function(colorMap = NULL, max_discrete = 30){

  show <- lapply(seq_along(colorMap), function(x){
      if(attr(colorMap[[x]],"discrete") && length(unique(colorMap[[x]])) > max_discrete){
        sl <- FALSE
      }else{
        sl <- TRUE
      }
      return(sl)
    }) %>% unlist

  names(show) <- names(colorMap)

  return(show)

} # .checkShowLegend












#######################
# ArrowRead.R


















####################################################################
# Reading fragments from Arrow Files
####################################################################

#' Get the fragments from an ArchRProject 
#' 
#' This function retrieves the fragments from a given ArchRProject as a GRangesList object.
#'
#' @param ArchRProject An `ArchRProject` object to get fragments from.
#' @param subsetBy A Genomic Ranges object to subset fragments by.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromProject <- function(
  ArchRProj = NULL,
  subsetBy = NULL,
  cellNames = NULL,
  verbose = FALSE,
  logFile = createLogFile("getFragmentsFromProject")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = subsetBy, name = "subsetBy", valid = c("GRanges", "null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))

  ArrowFiles <- getArrowFiles(ArchRProj)

  if(!is.null(subsetBy)){
    chr <- paste0(unique(seqnames(subsetBy)))
  }else{
    chr <- NULL
  }

  ArchR:::.startLogging(logFile = logFile)

  FragmentsList <- lapply(seq_along(ArrowFiles), function(x){
    message(sprintf("Reading ArrowFile %s of %s", x, length(ArrowFiles)))
    fragx <- getFragmentsFromArrow(
      ArrowFile = ArrowFiles[x], 
      chr = chr, 
      cellNames = cellNames, 
      verbose = verbose,
      logFile = logFile
    )
    if(!is.null(subsetBy)){
      fragx <- subsetByOverlaps(fragx, subsetBy, ignore.strand = TRUE)
    }
    fragx
  }) %>% SimpleList

  names(FragmentsList) <- names(ArrowFiles)

  FragmentsList

}

#' Get the fragments from an ArrowFile 
#' 
#' This function retrieves the fragments from a given ArrowFile as a GRanges object.
#'
#' @param ArrowFile The path to the ArrowFile from which fragments should be obtained.
#' @param chr A name of a chromosome to be used to subset the fragments `GRanges` object to a specific chromsome if desired.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells
#' from the provided ArrowFile using `getCellNames()`.
#' @param verbose A boolean value indicating whether to use verbose output during execution of this function. Can be set to `FALSE` for a cleaner output.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getFragmentsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  cellNames = NULL, 
  verbose = TRUE,
  logFile = createLogFile("getFragmentsFromArrow")
  ){

  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = chr, name = "chr", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getFragmentsFromArrow Input-Parameters", logFile = logFile)

  ArrowFile <- .validArrow(ArrowFile)

  if(is.null(chr)){
    chr <- .availableSeqnames(ArrowFile, subGroup = "Fragments")
  }

  if(any(chr %ni% .availableSeqnames(ArrowFile, subGroup = "Fragments"))){
    stop("Error Chromosome not in ArrowFile!")
  }
  
  out <- lapply(seq_along(chr), function(x){
    .logDiffTime(sprintf("Reading Chr %s of %s", x, length(chr)), t1 = tstart, verbose = verbose, logFile = logFile)
    .getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = chr[x], 
      out = "GRanges", 
      cellNames = cellNames, 
      method = "fast"
    )
  })


  .logDiffTime("Merging", tstart, t1 = tstart, verbose = verbose, logFile = logFile)

  out <- tryCatch({

    o <- .suppressAll(unlist(GRangesList(out, compress = FALSE)))

    if(.isGRList(o)){
      stop("Still a GRangesList")
    }

    o

  }, error = function(x){
    
    o <- c()
    
    for(i in seq_along(out)){
      if(!is.null(out[[i]])){
        if(i == 1){
          o <- out[[i]]
        }else{
          o <- c(o, out[[i]])
        }
      }
    }
    
    o

  })

  out

}

.getFragsFromArrow <- function(
  ArrowFile = NULL, 
  chr = NULL, 
  out = "GRanges", 
  cellNames = NULL, 
  method = "fast"
  ){

  if(is.null(chr)){
    stop("Need to provide chromosome to read!")
  }

  o <- h5closeAll()
  #ArrowFile <- .validArrow(ArrowFile)
  
  avSeq <- .availableSeqnames(ArrowFile)
  if(chr %ni% avSeq){
    stop(paste0("Chromosome ", chr ," not in ArrowFile! Available Chromosomes are : ", paste0(avSeq, collapse=",")))
  }

  #Get Sample Name
  sampleName <- .h5read(ArrowFile, paste0("Metadata/Sample"), method = method)

  o <- h5closeAll()
  nFrags <- sum(.h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method))

  if(nFrags==0){
    if(tolower(out)=="granges"){
      output <- GRanges(seqnames = chr, IRanges(start = 1, end = 1), RG = "tmp")
      output <- output[-1,]
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- output[-1,]
    }
    return(output)
  }

  if(is.null(cellNames) | tolower(method) == "fast"){
    
    output <- .h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), method = method) %>% 
      {IRanges(start = .[,1], width = .[,2])}
    mcols(output)$RG <- Rle(
      values = paste0(sampleName, "#", .h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues"), method = method)), 
      lengths = .h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths"), method = method)
    )
    if(!is.null(cellNames)){
      output <- output[BiocGenerics::which(mcols(output)$RG %bcin% cellNames)]
    }

  }else{
    
    if(!any(cellNames %in% .availableCells(ArrowFile))){

      stop("None of input cellNames are in ArrowFile availableCells!")

    }else{

      barRle <- Rle(h5read(ArrowFile, paste0("Fragments/",chr,"/RGValues")), h5read(ArrowFile, paste0("Fragments/",chr,"/RGLengths")))
      barRle@values <- paste0(sampleName, "#", barRle@values)
      idx <- BiocGenerics::which(barRle %bcin% cellNames)
      if(length(idx) > 0){
        output <- h5read(ArrowFile, paste0("Fragments/",chr,"/Ranges"), index = list(idx, 1:2)) %>% 
          {IRanges(start = .[,1], width = .[,2])}
        mcols(output)$RG <- barRle[idx]
      }else{
        output <- IRanges(start = 1, end = 1)
        mcols(output)$RG <- c("tmp")
        output <- output[-1,]
      }
    }

  }
  
  o <- h5closeAll()

  if(tolower(out)=="granges"){
    if(length(output) > 0){
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)    
    }else{
      output <- IRanges(start = 1, end = 1)
      mcols(output)$RG <- c("tmp")
      output <- GRanges(seqnames = chr, ranges(output), RG = mcols(output)$RG)
      output <- output[-1,]
    }
  }

  return(output)
}

####################################################################
# Reading Matrices/Arrays from Arrow Files
####################################################################

#' Get a data matrix stored in an ArchRProject
#' 
#' This function gets a given data matrix from an `ArchRProject`.
#'
#' @param ArchRProj An `ArchRProject` object to get data matrix from.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromProject <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
  ){

  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  .validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))
  .validInput(input = threads, name = "threads", valid = c("integer"))
  
  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromProject Input-Parameters", logFile = logFile)

  ArrowFiles <- getArrowFiles(ArchRProj)

  cellNames <- ArchRProj$cellNames

  avMat <- getAvailableMatrices(ArchRProj)
  if(useMatrix %ni% avMat){
    stop("useMatrix is not in Available Matrices see getAvailableMatrices")
  }

  seL <- .safelapply(seq_along(ArrowFiles), function(x){

    .logDiffTime(paste0("Reading ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
      t1 = tstart, verbose = FALSE, logFile = logFile)

    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]

    if(length(allCells) != 0){

      o <- getMatrixFromArrow(
        ArrowFile = ArrowFiles[x],
        useMatrix = useMatrix,
        useSeqnames = useSeqnames,
        cellNames = allCells, 
        ArchRProj = ArchRProj,
        verbose = FALSE,
        binarize = binarize,
        logFile = logFile
      )

      .logDiffTime(paste0("Completed ", useMatrix," : ", names(ArrowFiles)[x], "(",x," of ",length(ArrowFiles),")"), 
        t1 = tstart, verbose = FALSE, logFile = logFile)

      o

    }else{

      NULL
      
    }

  }, threads = threads) 

  #ColData
  .logDiffTime("Organizing colData", t1 = tstart, verbose = verbose, logFile = logFile)
  cD <- lapply(seq_along(seL), function(x){
    colData(seL[[x]])
  }) %>% Reduce("rbind", .)
  
  #RowData
  .logDiffTime("Organizing rowData", t1 = tstart, verbose = verbose, logFile = logFile)
  rD1 <- rowData(seL[[1]])
  rD <- lapply(seq_along(seL), function(x){
    identical(rowData(seL[[x]]), rD1)
  }) %>% unlist %>% all
  if(!rD){
    stop("Error with rowData being equal for every sample!")
  }

  #RowRanges
  .logDiffTime("Organizing rowRanges", t1 = tstart, verbose = verbose, logFile = logFile)
  rR1 <- rowRanges(seL[[1]])
  rR <- lapply(seq_along(seL), function(x){
    identical(rowRanges(seL[[x]]), rR1)
  }) %>% unlist %>% all
  if(!rR){
    stop("Error with rowRanges being equal for every sample!")
  }

  #Assays
  nAssays <- names(assays(seL[[1]]))
  asy <- lapply(seq_along(nAssays), function(i){
    .logDiffTime(sprintf("Organizing Assays (%s of %s)", i, length(nAssays)), t1 = tstart, verbose = verbose, logFile = logFile)
    m <- lapply(seq_along(seL), function(j){
      assays(seL[[j]])[[nAssays[i]]]
    }) %>% Reduce("cbind", .)
    m
  }) %>% SimpleList()
  names(asy) <- nAssays
  
  .logDiffTime("Constructing SummarizedExperiment", t1 = tstart, verbose = verbose, logFile = logFile)
  if(!is.null(rR1)){
    se <- SummarizedExperiment(assays = asy, colData = cD, rowRanges = rR1)
  }else{
    se <- SummarizedExperiment(assays = asy, colData = cD, rowData = rD1)
  }
  rm(seL)
  gc()

  .logDiffTime("Finished Matrix Creation", t1 = tstart, verbose = verbose, logFile = logFile)

  se
  
}

#' Get a data matrix stored in an ArrowFile
#' 
#' This function gets a given data matrix from an individual ArrowFile.
#'
#' @param ArrowFile The path to an ArrowFile from which the selected data matrix should be obtained.
#' @param useMatrix The name of the data matrix to retrieve from the given ArrowFile. Options include "TileMatrix", "GeneScoreMatrix", etc.
#' @param useSeqnames A character vector of chromosome names to be used to subset the data matrix being obtained.
#' @param cellNames A character vector indicating the cell names of a subset of cells from which fragments whould be extracted.
#' This allows for extraction of fragments from only a subset of selected cells. By default, this function will extract all cells from
#' the provided ArrowFile using `getCellNames()`.
#' @param ArchRProj An `ArchRProject` object to be used for getting additional information for cells in `cellColData`.
#' In some cases, data exists within the `ArchRProject` object that does not exist within the ArrowFiles. To access this data, you can
#' provide the `ArchRProject` object here.
#' @param verbose A boolean value indicating whether to use verbose output during execution of  this function. Can be set to FALSE for a cleaner output.
#' @param binarize A boolean value indicating whether the matrix should be binarized before return. This is often desired when working with insertion counts.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export
getMatrixFromArrow <- function(
  ArrowFile = NULL, 
  useMatrix = "GeneScoreMatrix",
  useSeqnames = NULL,
  cellNames = NULL, 
  ArchRProj = NULL,
  verbose = TRUE,
  binarize = FALSE,
  logFile = createLogFile("getMatrixFromArrow")
  ){

  .validInput(input = ArrowFile, name = "ArrowFile", valid = "character")
  .validInput(input = useMatrix, name = "useMatrix", valid = "character")
  .validInput(input = useSeqnames, name = "useSeqnames", valid = c("character","null"))
  .validInput(input = cellNames, name = "cellNames", valid = c("character","null"))
  .validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj","null"))
  .validInput(input = verbose, name = "verbose", valid = c("boolean"))
  .validInput(input = binarize, name = "binarize", valid = c("boolean"))

  tstart <- Sys.time()
  .startLogging(logFile = logFile)
  .logThis(mget(names(formals()),sys.frame(sys.nframe())), "getMatrixFromArrow Input-Parameters", logFile = logFile)

  ArrowFile <- .validArrow(ArrowFile)
  sampleName <- .sampleName(ArrowFile)

  seqnames <- .availableSeqnames(ArrowFile, subGroup = useMatrix)
  featureDF <- .getFeatureDF(ArrowFile, subGroup = useMatrix)
  .logThis(featureDF, paste0("featureDF ", sampleName), logFile = logFile)

  if(!is.null(useSeqnames)){
    seqnames <- seqnames[seqnames %in% useSeqnames]
  }

  if(length(seqnames) == 0){
    stop("No seqnames available!")
  }

  featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames), ]

  .logDiffTime(paste0("Getting ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
    t1 = tstart, verbose = verbose, logFile = logFile)

  if(!is.null(cellNames)){
    allCells <- .availableCells(ArrowFile = ArrowFile, subGroup = useMatrix)
    if(!all(cellNames %in% allCells)){
      stop("cellNames must all be within the ArrowFile!!!!")
    }
  }

  mat <- .getMatFromArrow(
    ArrowFile = ArrowFile, 
    featureDF = featureDF, 
    cellNames = cellNames, 
    useMatrix = useMatrix,
    binarize = binarize,
    useIndex = FALSE
  )
  .logThis(mat, paste0("mat ", sampleName), logFile = logFile)

  .logDiffTime(paste0("Organizing SE ",useMatrix," from ArrowFile : ", basename(ArrowFile)), 
    t1 = tstart, verbose = verbose, logFile = logFile)
  matrixClass <- h5read(ArrowFile, paste0(useMatrix, "/Info/Class"))

  if(matrixClass == "Sparse.Assays.Matrix"){
    rownames(mat) <- paste0(featureDF$name)
    splitIdx <- split(seq_len(nrow(mat)), featureDF$seqnames)
    mat <- lapply(seq_along(splitIdx), function(x){
      mat[splitIdx[[x]], , drop = FALSE]
    }) %>% SimpleList
    names(mat) <- names(splitIdx)
    featureDF <- featureDF[!duplicated(paste0(featureDF$name)), ,drop = FALSE]
    featureDF <- featureDF[,which(colnames(featureDF) %ni% "seqnames"), drop=FALSE]
    rownames(featureDF) <- paste0(featureDF$name)
  }else{
    mat <- SimpleList(mat)
    names(mat) <- useMatrix    
  }

  colData <- .getMetadata(ArrowFile)
  colData <- colData[colnames(mat[[1]]),,drop=FALSE]

  if(!is.null(ArchRProj)){
    projColData <- getCellColData(ArchRProj)[rownames(colData), ]
    colData <- cbind(colData, projColData[ ,colnames(projColData) %ni% colnames(colData)])
  }

  rowData <- tryCatch({
    makeGRangesFromDataFrame(featureDF, keep.extra.columns = TRUE)
  }, error = function(x){
    featureDF
  })

  se <- SummarizedExperiment(
    assays = mat,
    rowData = rowData,
    colData = colData
  )
  .logThis(se, paste0("se ", sampleName), logFile = logFile)

  se

}

.getMatFromArrow <- function(
  ArrowFile = NULL, 
  featureDF = NULL, 
  binarize = NULL, 
  cellNames = NULL,
  useMatrix = "TileMatrix", 
  useIndex = FALSE,
  threads = 1
  ){

  if(is.null(featureDF)){
    featureDF <- .getFeatureDF(ArrowFile, useMatrix)
  }

  if(any(c("seqnames","idx") %ni% colnames(featureDF))){
    stop("Need to provide featureDF with columns seqnames and idx!")
  }

  #Add RowNames for Check at the end
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))

  o <- h5closeAll()

  matClass <- h5read(ArrowFile, paste0(useMatrix,"/Info/Class"))
  if(matClass %ni% c("Sparse.Binary.Matrix", "Sparse.Integer.Matrix", "Sparse.Double.Matrix", "Sparse.Assays.Matrix")){
    stop("Arrow Mat is not a valid Sparse Matrix!")
  }
  if(is.null(binarize)){
    if(matClass == "Sparse.Binary.Matrix"){
      binarize <- TRUE
    }else{
      binarize <- FALSE
    }
  }
  if(matClass == "Sparse.Binary.Matrix"){
    if(!binarize){
      stop("Sparse Matrix in Arrow is Binarized! Set binarize = TRUE to use matrix!")
    }
  }

  matColNames <- paste0(.sampleName(ArrowFile), "#", h5read(ArrowFile, paste0(useMatrix,"/Info/CellNames")))
  if(!is.null(cellNames)){
    idxCols <- which(matColNames %in% cellNames)
  }else{
    idxCols <- seq_along(matColNames)
  }

  seqnames <- unique(featureDF$seqnames)

  mat <- .safelapply(seq_along(seqnames), function(x){

    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex),]
    idxRows <- featureDFx$idx

    j <- Rle(
      values = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jValues")), 
      lengths = h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/jLengths"))
      )

    #Match J
    matchJ <- S4Vectors::match(j, idxCols, nomatch = 0)
    idxJ <- BiocGenerics::which(matchJ > 0)
    if(useIndex){
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"), index = list(idxJ, 1))
    }else{
      i <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/i"))[idxJ]
    }
    j <- matchJ[idxJ]

    #Match I
    matchI <- match(i, idxRows, nomatch = 0)
    idxI <- which(matchI > 0)
    i <- i[idxI]
    j <- j[idxI]
    i <- matchI[idxI]
    
    if(!binarize){
      x <- h5read(ArrowFile, paste0(useMatrix,"/",seqnamex,"/x"))[idxJ][idxI]
    }else{
      x <- rep(1, length(j))
    }

    mat <- Matrix::sparseMatrix(
      i=as.vector(i),
      j=j,
      x=x,
      dims = c(length(idxRows), length(idxCols))
    )
    rownames(mat) <- rownames(featureDFx)

    rm(matchI, idxI, matchJ, idxJ, featureDFx, idxRows)

    return(mat)

  }, threads = threads) %>% Reduce("rbind", .)

  o <- h5closeAll()

  colnames(mat) <- matColNames[idxCols]

  #Double Check Order!
  mat <- mat[rownames(featureDF), , drop = FALSE]
  rownames(mat) <- NULL

  if(!is.null(cellNames)){
    mat <- mat[,cellNames,drop=FALSE]
  }

  return(mat)

}

####################################################################
# Helper read functioning
####################################################################
.getGroupMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  groupList = NULL,
  threads = 1, 
  useIndex = FALSE, 
  verbose = TRUE, 
  useMatrix = "TileMatrix",
  asSparse = FALSE,
  tstart = NULL
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################
  seqnames <- unique(featureDF$seqnames)
  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  cellNames <- unlist(groupList, use.names = FALSE) ### UNIQUE here? doublet check QQQ

  allCellsList <- lapply(seq_along(ArrowFiles), function(x){
    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]
    if(length(allCells) != 0){
      allCells
    }else{
      NULL
    }
  })

  mat <- .safelapply(seq_along(seqnames), function(x){

    .logDiffTime(sprintf("Constructing Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)

    #Construct Matrix
    seqnamex <- seqnames[x]
    featureDFx <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnamex), ]

    matChr <- matrix(0, nrow = nrow(featureDFx), ncol = length(groupList))
    colnames(matChr) <- names(groupList)
    rownames(matChr) <- rownames(featureDFx)

    for(y in seq_along(ArrowFiles)){

      allCells <- allCellsList[[y]]
      
      if(!is.null(allCells)){

        maty <- .getMatFromArrow(
          ArrowFile = ArrowFiles[y], 
          useMatrix = useMatrix,
          featureDF = featureDFx, 
          cellNames = allCells, 
          useIndex = useIndex
        )

        for(z in seq_along(groupList)){

          #Check Cells In Group
          cellsGroupz <- groupList[[z]]
          idx <- BiocGenerics::which(colnames(maty) %in% cellsGroupz)

          #If In Group RowSums
          if(length(idx) > 0){
            matChr[,z] <- matChr[,z] + Matrix::rowSums(maty[,idx,drop=FALSE])
          }

        }

        rm(maty)

      }
     

      if(y %% 20 == 0 | y %% length(ArrowFiles) == 0){
        gc()
      } 

    }

    if(asSparse){
      matChr <- as(matChr, "dgCMatrix")
    }

    .logDiffTime(sprintf("Finished Group Matrix %s of %s", x, length(seqnames)), tstart, verbose = verbose)
    
    matChr

  }, threads = threads) %>% Reduce("rbind", .)

  mat <- mat[rownames(featureDF), , drop = FALSE]
  
  .logDiffTime("Successfully Created Group Matrix", tstart, verbose = verbose)

  gc()

  return(mat)
  
}

.getPartialMatrix <- function(
  ArrowFiles = NULL, 
  featureDF = NULL, 
  cellNames = NULL, 
  progress = TRUE, 
  threads = 1, 
  useMatrix = "TileMatrix",
  doSampleCells = FALSE, 
  sampledCellNames = NULL, 
  tmpPath = .tempfile(pattern = paste0("tmp-partial-mat")), 
  useIndex = FALSE,
  tstart = NULL,
  verbose = TRUE
  ){

  #########################################
  # Time Info
  #########################################
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #########################################
  # Construct Matrix
  #########################################

  mat <- .safelapply(seq_along(ArrowFiles), function(x){
    
    .logDiffTime(sprintf("Getting Partial Matrix %s of %s", x, length(ArrowFiles)), tstart, verbose = verbose)

    allCells <- .availableCells(ArrowFile = ArrowFiles[x], subGroup = useMatrix)
    allCells <- allCells[allCells %in% cellNames]

    if(length(allCells) == 0){
      if(doSampleCells){
        return(list(mat = NULL, out = NULL))
      }else{
        return(NULL)
      }
    }

    o <- h5closeAll()
    matx <- .getMatFromArrow(
      ArrowFile = ArrowFiles[x], 
      featureDF = featureDF, 
      cellNames = allCells,
      useMatrix = useMatrix, 
      useIndex = useIndex
    )
   
    if(doSampleCells){

      #Save Temporary Matrix
      outx <- paste0(tmpPath, "-", .sampleName(ArrowFiles[x]), ".rds")
      .safeSaveRDS(matx, outx, compress = FALSE)     

      #Sample Matrix 
      matx <- matx[, which(colnames(matx) %in% sampledCellNames),drop = FALSE]
      
      return(list(mat = matx, out = outx))

    }else{
      
      return(matx)

    }

  }, threads = threads)

  gc()
  

  if(doSampleCells){

    matFiles <- lapply(mat, function(x) x[[2]]) %>% Reduce("c", .)
    mat <- lapply(mat, function(x) x[[1]]) %>% Reduce("cbind", .)
    if(!all(cellNames %in% colnames(mat))){
      .logThis(sampledCellNames, "cellNames supplied", logFile = logFile)
      .logThis(colnames(mat), "cellNames from matrix", logFile = logFile)
      stop("Error not all cellNames found in partialMatrix")
    }
    mat <- mat[,sampledCellNames, drop = FALSE]
    mat <- .checkSparseMatrix(mat, length(sampledCellNames))

    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(list(mat = mat, matFiles = matFiles))

  }else{

    mat <- Reduce("cbind", mat)
    if(!all(cellNames %in% colnames(mat))){
      .logThis(cellNames, "cellNames supplied", logFile = logFile)
      .logThis(colnames(mat), "cellNames from matrix", logFile = logFile)
      stop("Error not all cellNames found in partialMatrix")
    }
    mat <- mat[,cellNames, drop = FALSE]
    mat <- .checkSparseMatrix(mat, length(cellNames))
    
    .logDiffTime("Successfully Created Partial Matrix", tstart, verbose = verbose)

    return(mat)

  }


}

.checkSparseMatrix <- function(x, ncol = NULL){
  isSM <- is(x, 'sparseMatrix')
  if(!isSM){
    if(is.null(ncol)){
      stop("ncol must not be NULL if x is not a matrix!")
    }
    cnames <- tryCatch({
      names(x)
    }, error = function(e){
      colnames(x)
    })
    if(length(cnames) != ncol){
      stop("cnames != ncol!")
    }
    x <- Matrix::Matrix(matrix(x, ncol = ncol), sparse=TRUE)
    colnames(x) <- cnames
  }
  x
}

########################################################################
# Compute Summary Statistics!
########################################################################

.getRowSums <- function(
  ArrowFiles = NULL,
  useMatrix = NULL,
  seqnames = NULL,
  verbose = TRUE,
  tstart = NULL,
  filter0 = FALSE,
  threads = 1,
  addInfo = FALSE
  ){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
    
  if(is.null(seqnames)){
    seqnames <- .availableSeqnames(ArrowFiles, useMatrix)
  }

  #Compute RowSums
  summaryDF <- .safelapply(seq_along(seqnames), function(x){
    o <- h5closeAll()
    for(y in seq_along(ArrowFiles)){
      if(y == 1){
        sumy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
      }else{
        sumy1 <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/rowSums"))
        if(length(sumy1) != length(sumy)){
          stop("rowSums lengths do not match in ArrowFiles for a seqname!")
        }else{
          sumy <- sumy + sumy1
        }
      }
    }
    #Return Setup In Feature DF Format (seqnames, idx columns)
    DataFrame(seqnames = Rle(seqnames[x], lengths = length(sumy)), idx = seq_along(sumy), rowSums = as.vector(sumy))
  }, threads = threads) %>% Reduce("rbind", .)
  
  if(addInfo){
    featureDF <- .getFeatureDF(ArrowFiles, useMatrix)
    rownames(featureDF) <- paste0(featureDF$seqnames, "_", featureDF$idx)
    rownames(summaryDF) <- paste0(summaryDF$seqnames, "_", summaryDF$idx)
    featureDF <- featureDF[rownames(summaryDF), , drop = FALSE]
    featureDF$rowSums <- summaryDF[rownames(featureDF), "rowSums"]
    summaryDF <- featureDF
    rownames(summaryDF) <- NULL
    remove(featureDF)
  }

  if(filter0){
    summaryDF <- summaryDF[which(summaryDF$rowSums > 0), ,drop = FALSE]
  }

  return(summaryDF)

}

.getRowVars <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  useLog2 = FALSE,
  threads = 1
  ){
  
  .combineVariances <- function(dfMeans = NULL, dfVars = NULL, ns = NULL){

    #https://rdrr.io/cran/fishmethods/src/R/combinevar.R

    if(ncol(dfMeans) != ncol(dfVars) | ncol(dfMeans) != length(ns)){
      stop("Means Variances and Ns lengths not identical")
    }

    #Check if samples have NAs due to N = 1 sample or some other weird thing.
    #Set it to min non NA variance
    dfVars <- lapply(seq_len(nrow(dfVars)), function(x){
      vx <- dfVars[x, ]
      if(any(is.na(vx))){
        vx[is.na(vx)] <- min(vx[!is.na(vx)])
      }
      vx
    }) %>% Reduce("rbind", .)

    combinedMeans <- rowSums(t(t(dfMeans) * ns)) / sum(ns)
    summedVars <- rowSums(t(t(dfVars) * (ns - 1)) + t(t(dfMeans^2) * ns))
    combinedVars <- (summedVars - sum(ns)*combinedMeans^2)/(sum(ns)-1)

    data.frame(combinedVars = combinedVars, combinedMeans = combinedMeans)

  }

  featureDF <- .getFeatureDF(ArrowFiles, useMatrix)

  if(!is.null(seqnames)){
    featureDF <- featureDF[BiocGenerics::which(featureDF$seqnames %bcin% seqnames),]
  }

  rownames(featureDF) <- paste0("f", seq_len(nrow(featureDF)))
  fnames <- rownames(featureDF)

  featureDF <- S4Vectors::split(featureDF, as.character(featureDF$seqnames))

  ns <- lapply(seq_along(ArrowFiles), function(y){
    length(.availableCells(ArrowFiles[y], useMatrix))
  }) %>% unlist

  #Compute RowVars
  summaryDF <- .safelapply(seq_along(featureDF), function(x){
    
    o <- h5closeAll()
    seqx <- names(featureDF)[x]
    meanx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))
    varx <- matrix(NA, ncol = length(ArrowFiles), nrow = nrow(featureDF[[x]]))

    for(y in seq_along(ArrowFiles)){

      if(useLog2){
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeansLog2"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVarsLog2")) 
      }else{
        meanx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowMeans"))
        varx[, y] <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqx, "/rowVars"))
      }

    }

    cbind(featureDF[[x]], DataFrame(.combineVariances(meanx, varx, ns)))

  }, threads = threads) %>% Reduce("rbind", .)

  summaryDF <- summaryDF[fnames, , drop = FALSE]
  
  return(summaryDF)

}

.getColSums <- function(
  ArrowFiles = NULL,
  seqnames = NULL,
  useMatrix = NULL,
  verbose = TRUE,
  tstart = NULL,
  threads = 1
  ){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Compute ColSums
  cS <- .safelapply(seq_along(seqnames), function(x){
    
    lapply(seq_along(ArrowFiles), function(y){
      
      o <- h5closeAll()
      cSy <- h5read(ArrowFiles[y], paste0(useMatrix, "/", seqnames[x], "/colSums"))
      rownames(cSy) <- .availableCells(ArrowFiles[y], useMatrix)
      cSy
      
    }) %>% Reduce("rbind", .)
    
  }, threads = threads) %>% Reduce("cbind", .) %>% rowSums

  .logDiffTime("Successfully Computed colSums", tstart, verbose = verbose)

  return(cS)

}

# h5read implementation for optimal reading
.h5read <- function(
  file = NULL,
  name = NULL,
  method = "fast",
  index = NULL,
  start = NULL,
  block = NULL,
  count = NULL
  ){

  if(tolower(method) == "fast" & is.null(index) & is.null(start) & is.null(block) & is.null(count)){
    fid <- H5Fopen(file)
    dapl <- H5Pcreate("H5P_DATASET_ACCESS")
    did <- .Call("_H5Dopen", fid@ID, name, dapl@ID, PACKAGE='rhdf5')
    res <- .Call("_H5Dread", did, NULL, NULL, NULL, TRUE, 0L, FALSE, fid@native, PACKAGE='rhdf5')
    invisible(.Call("_H5Dclose", did, PACKAGE='rhdf5'))   
  }else{
    res <- h5read(file = file, name = name, index = index, start = start, block = block, count = count)
  }
  o <- h5closeAll()
  return(res)
}

.getMatrixClass <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
  ){

  threads <- min(length(ArrowFiles), threads)

  matrixClass <- .safelapply(seq_along(ArrowFiles), function(i){
    h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Class"))
  }, threads = threads) %>% unlist %>% unique

  if(length(matrixClass) != 1){
    stop("Not all matrix classes are the same!")
  }

  matrixClass

}

.getMatrixUnits <- function(
  ArrowFiles = NULL, 
  useMatrix = NULL,
  threads = getArchRThreads()
  ){

  threads <- min(length(ArrowFiles), threads)

  matrixUnits <- .safelapply(seq_along(ArrowFiles), function(i){
    tryCatch({ #This handles backwards compatibility!
      h5read(ArrowFiles[i], paste0(useMatrix, "/Info/Units"))
    }, error = function(x){
      "None"
    })
  }, threads = threads) %>% unlist %>% unique

  if(length(matrixUnits) != 1){
    stop("Not all matrix units are the same!")
  }

  matrixUnits

}







#######################
# ReproduciblePeakSet.R


# plotPeakCallSummary
#
# input:
#    type: {["freq"], "percent"}
#
# usage: plotPDF(.plotPeakCallSummary(ArchRProj), name = "Peak-Call-Summary", width = 8, height = 5, ArchRProj = ArchRProj, addDOC = FALSE)
plotPeakCallSummary <- function(
        ArchRProj = NULL,
        pal = NULL,
	type = "freq",
	cell.type.order=NULL
        ){

  peakDF <- metadata(ArchRProj@peakSet)$PeakCallSummary

  if(is.null(peakDF)){
    stop("Error no Peak Call Summary available are you sure these peaks were called with CreateReproduciblePeakSet?")
  }

  # begin of addition by H. Kim
  switch(type,
	"percent"={
		peakDF <- peakDF %>% group_by(Group) %>% mutate(percent = 100*(Freq / sum(Freq)))
		var_y <- "percent"
		ylab <- "Percent"
	},
	{
		var_y <- "Freq"
		ylab <- "Number of Peaks (x10^3)"
	}
  ) # switch

  if (!is.null(cell.type.order)) {
    group_u <- unique(peakDF$Group)
    cell.type.order <- c(cell.type.order, "UnionPeaks")
    idx_order <- c()
    for (p in cell.type.order) {
	p <- gsub("\\+", "\\\\\\+", p)
 	p <- sprintf("^%s", p)
	idx <- grep(p, group_u)
	idx_order <- c(idx_order, idx)
    }
    peakDF$Group <- factor(peakDF$Group, levels=group_u[idx_order])
  }
  # end of addition

  if(is.null(pal)){
    pal <- paletteDiscrete(values=peakDF$Var1)
  }
  pal <- c("Distal" = "#60BA64", "Exonic" = "#73C6FF", "Intronic" = "#620FA3", "Promoter" = "#FFC554", pal)
  pal <- pal[!duplicated(names(pal))]

  #lengthMax <- split(peakDF$Freq, peakDF$Group) %>% lapply(sum) %>% unlist %>% max
  lengthMax <- split(peakDF[,var_y], peakDF$Group) %>% lapply(sum) %>% unlist %>% max

  colnames(peakDF)[colnames(peakDF)=="Var1"] <- "PeakAnno"

 # p <- ggplot(peakDF, aes(x=Group, y=Freq, fill=PeakAnno)) +
 p <- ggplot(peakDF, aes_string(x="Group", y=var_y, fill="PeakAnno")) +
    geom_bar(stat = "identity") +
    theme_ArchR(xText90 = TRUE) +
    #ylab("Number of Peaks (x10^3)") +
    ylab(ylab) +
    xlab("") +
    theme(legend.position = "bottom",
      legend.key = element_rect(size = 2),
      legend.box.background = element_rect(color = NA),
      # begin of addition by H. Kim
      legend.title=element_text(size=12), 
      legend.text=element_text(size=12),
      plot.margin = margin(0.5, 0.5, 0.5, 2, "cm"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
      # end of addition
    ) +
    scale_fill_manual(values=pal)

  # begin of addition by H. Kim
  switch(type,
	"percent"={
  	},
	{
    		p <- p +scale_y_continuous(
  		    breaks = seq(0, lengthMax * 2,50),
 		     limits = c(0, lengthMax * 1.1),
    		  expand = c(0,0)
     		 )
	}
  ) # switch
  # end of addition

  attr(p, "ratioYX") <- 0.5

  return(p)

} # plotPeakCallSummary












