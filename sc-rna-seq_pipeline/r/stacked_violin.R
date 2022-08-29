#
# stacked_violin.R
# modified by H. Kim
# last modified on 2021, Oct.
#
# motivation:
# modified for R >= 4.0.2
#
# reference:
# https://rpubs.com/DarrenVan/628853
# https://satijalab.org/seurat/reference/vlnplot
#

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))








# modify_vlnplot
# remove the x-axis text and tick
# plot.margin to adjust the white space between each plot.
# ... pass any arguments to VlnPlot in Seurat
modify_vlnplot <- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {

  # begin of addition by H. Kim
  # B_cells.8 --> B_cells
  ylab_feature <- feature
  ylab_feature <- gsub("\\..*$", "", ylab_feature)
  ylab_feature <- gsub("_", "\n", ylab_feature)
  # end of addition



  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    theme_nothing() +
    xlab("") + ylab(ylab_feature) + ggtitle("") + 
    theme(rect = element_blank(),
	  axis.text.x = element_text(angle = 90), 
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_text(size = rel(0.9), angle = 90), 
          axis.text.y = element_text(size = rel(1)), 
	  axis.line = element_line(size = 0.9, colour = "black"),
	  axis.ticks.length = unit(0.2, "cm"),
	  legend.position = "none",
	  legend.margin = margin(0, 0, 0, 0),
	  legend.box.margin =  margin(0, 0, 0, 0),
          plot.margin = plot.margin
   )



  # begin of addition by H. Kim
  x_labs <- ggplot_build(p)$layout$panel_params[[1]]$x$get_labels()
  y_labs <- ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()

  cluster_num <- gsub("-.*", "", x_labs)
  breaks = seq(max(min(y_labs, na.rm=T), 0), max(y_labs, na.rm=T), length.out = 3)

  if (grepl("Smooth", feature)) {

    p <- p + scale_x_discrete(labels=cluster_num)
    suppressMessages( p <- p + scale_y_continuous(breaks=breaks) )

  } else {

    p <- p + scale_x_discrete(labels=NULL)
    #p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
    suppressMessages( p <- p + scale_y_continuous(breaks=breaks) )

  } # if

  # end of addition



  return(p)

} # modify_vlnplot







# extract the max value of the y axis
extract_max <- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
} # extract_max














# main function
StackedVlnPlot <- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  


  # begin of modification by H. Kim 
  #plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x,  ...))

  idx <- which((features %in% colnames(obj@meta.data)) | (features %in% rownames(obj)))
  if (length(idx) == 0) {
	return(NULL)
  }
  features <- features[idx]

  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, pt.size=pt.size, plot.margin=plot.margin,  ...))
  # end of modification
  



  # add back x-axis title to bottom plot. patchwork is going to support this?
  #plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] + theme(axis.text.x=element_text(), axis.ticks.x = element_line(), axis.ticks.y = element_line())

  for (i in 1:length(plot_list)) {
  	plot_list[[i]]<- plot_list[[i]] + theme(axis.text.x=element_text(), axis.ticks.x = element_line(), axis.ticks.y = element_line())
  } # for
  



  # change the y-axis tick to only max value 
  #ymaxs<- purrr::map_dbl(plot_list, extract_max)
  # plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
  #        scale_y_continuous(breaks = c(y)) + 
  #        expand_limits(y = y))
  # 


  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)

  return(p)

} # StackedVlnPlot


