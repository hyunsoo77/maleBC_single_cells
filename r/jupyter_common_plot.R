#
# jupyter_common_plot.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Mar.
#
# content:
# print_ggpubr()
# get_volcano_plot()
# theme_geometry()
# 





### ggpubr



# print_ggpubr
#
# input:
#
#    obj: data.frame or seurat obj
#    method_comparison: {"t.test", ["wilcox.test"]}
#    list_comparision:
#    symnum.args: list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*",  "ns"))
#    test: anova, kruskal.test
#    ylimits: limits=c(-1, 6),
#    ybreaks: pretty( c(0, 6) , n = 3),
#    yexpand: c(0, 0)
#
#
print_ggpubr <- function(obj, x, y, color, fill="white", pattern_cell_type=NULL, col_cell_types=NULL, order_x=NULL, facet.by=NULL, facet.ncol=NULL, plot_type="boxplot", method_comparison="wilcox.test", list_comparison=NULL, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*",  "ns")), step.increase = -0.03, label.y = 5.0, vjust = 0, tip.length=c(0.01, 0.01), type_printing_mean=NULL, test=NULL, test_label.x=1, test_label.y=1, test_vjust=0.5, color_manual=NULL, fill_manual=NULL, ylimits=NULL, ybreaks = pretty(c(0, 6), n = 3), yexpand = c(0, 0), yintercept=NULL, xlim=NULL, ylim=NULL, title=NULL, xlab=NULL, ylab=NULL, angle=0, font_size=font_size_sc_clusters, width=8, height=8, str_condition_plot="", n_log=0) {

  suppressPackageStartupMessages(library(ggpubr))

  switch(class(obj),
	"data.frame"={
		df <- obj
	},
	"ArchRProject"={
		df <- NULL
		meta.data <- as.data.frame(getCellColData(obj))
	},
	"Seurat"={
		df <- NULL
		meta.data <- obj@meta.data
	},
	{}
  ) # switch

  if (is.null(df) && (length(y) > 1)) {

	df <- meta.data[,y[y %in% colnames(meta.data)]]

	# add gene expression
	mtx_logcounts <- GetAssayData(object = rna, assay="RNA", slot = "data")
	genes <- y[y %in% rownames(mtx_logcounts)]
	if (length(genes) == 1) {
		df[, genes] <- mtx_logcounts[genes,]
	} else if (length(genes) > 1) {
		df[, genes] <- t(mtx_logcounts[genes,])
	}

	suppressMessages(suppressWarnings(
		df <- reshape2::melt(df)
	))

	df$variable <- as.character(df$variable)
	if (!is.null(names(order_x))) {
		for (name in names(order_x)) {
			idx <- which(df$variable == name)
			df$variable[idx] <- order_x[name]
		}
	}
	x <- "variable"
	y <- "value"
	color <- "variable"
	order_x <- unname(order_x)

  } else if (is.null(df) && (!is.null(y) && (y == "cnt"))) {

	df <- meta.data %>% group_by_at(x) %>% dplyr::summarise(cnt = n())
	if (!is.null(order_x)) {
		if ((length(order_x) == 1) && (order_x %in% colnames(df))) {
			order_x <- df[[x]][order(df[[order_x]])]
		}
		f <- df[[x]] %in% order_x
		#f <- f & (df[[x]] > 0.25)
		df <- df[f,]
		df[[x]] <- factor(df[[x]], levels=order_x)
	}

  } else if (is.null(df)) {
	vec <- c(x, y, facet.by, col_cell_types)
	f <- vec %in% colnames(meta.data)
	df <- data.frame()
	if (any(f)) {
		df <- meta.data[, vec[f], drop=F]
		cols_info <- setdiff(vec[f], x)
		genes <- NULL
		if (any(!f)) {
			genes <- vec[!f]
			genes <- genes[genes %in% rownames(rna)]
			if (length(genes) > 0) {
				mtx_logcounts <- GetAssayData(object = rna, assay="RNA", slot = "data")
				if (length(genes) == 1) {
					df[, genes] <- mtx_logcounts[genes,]
				} else {
					df[, genes] <- t(mtx_logcounts[genes,])
				}
			}
		}

		if (!is.null(pattern_cell_type)) {
			f <- grepl(pattern_cell_type, df[[col_cell_types]])
			df <- df[f,,drop=F]
			log_obj(table(df[[col_cell_types]]), tab=2)
		}	

		if (n_log > 0) {
			for (col in cols_info) {
				if (!is.numeric(df[[col]])) next
				df %>% group_by_at(x) %>% dplyr::summarize_at(col, mean, na.rm=TRUE) %>% arrange(desc(!!sym(col))) -> df1
				df %>% group_by_at(x) %>% dplyr::summarise(nzero = sum(!!sym(col) != 0), cnt=n()) %>% mutate(perc = formattable::percent(nzero / cnt)) %>% arrange(desc(perc)) -> df2
				df_summary <- left_join(df1, df2, by=x)
				switch(plot_type,
					"scatterplot"={
						log_obj(head(df_summary), tab=2)
					},
					{
						log_obj(df_summary, tab=2)
					}
				) # switch
				if (!is.null(order_x) && (length(order_x) == 1) && (order_x %in% colnames(df_summary))) {
					order_x <- df_summary[[x]][order(df_summary[[order_x]])]
				}
			}
		}

		if (!is.null(order_x)) {
			df[[x]] <- as.character(df[[x]])
			if (!is.null(names(order_x))) {
				for (name in names(order_x)) {
					idx <- which(df[[x]] == name)
					df[[x]][idx] <- order_x[name]
				}
			}

			f <- df[[x]] %in% order_x
			#f <- f & (df[[x]] > 0.25)
			df <- df[f,,drop=F]
			df[[x]] <- factor(df[[x]], levels=order_x)
		}

		if (n_log > 0) {
			for (gene in genes) {
				if (!is.numeric(df[[gene]])) next
				df1 <- df %>% group_by_at(x) %>% dplyr::summarize_at(gene, mean, na.rm=TRUE) %>% arrange(desc(!!sym(gene)))
				switch(plot_type,
					"scatterplot"={
						log_obj(head(df1), tab=2)
					},
					{
						log_obj(df1, tab=2)
					}
				) # switch

			}
		}

	
		if (is.null(y)) {
			df <- reshape2::melt(df, id.vars=x)
			y <- "value"
			facet.by <- "variable"
		}
	} # if
  } # if

  #df <- df[complete.cases(df),]

  if (n_log > 0) {
	log_obj(head(df), tab=2)
  }

  switch(plot_type,
	"barplot"={
		suppressMessages(suppressWarnings(
		gg <- ggbarplot(df, x = x, y = y,
			#ylim = c(0, 45)
			#add = "dotplot",
			#add = "jitter",
			add = "mean",
			#add = "mean_se",
			#add.params = list(color='black'),
			#add = c("dotplot", "mean_se"),
			#add = c("jitter", "mean_se"),
			#palette = "jco"
			color = color,
			fill = fill,
			facet.by = facet.by,
			ncol = facet.ncol,
			panel.labs.font = list(face = NULL, color = NULL, size = font_size, angle = NULL),
			ggtheme = theme_bw() # theme_pubr()
		)
		))
	},
	"boxplot"={
		suppressMessages(suppressWarnings(
		gg <- ggboxplot(df, x = x, y = y,
			notch = TRUE,
			#add = "jitter",
			add = "mean",
			#palette = "jco",
			color = color,
			fill = fill,
			facet.by = facet.by,
			ncol = facet.ncol,
			panel.labs.font = list(face = NULL, color = NULL, size = font_size, angle = NULL),
			ggtheme = theme_bw() # theme_pubr()
		)
		))



	},
	"scatterplot"={
		suppressMessages(suppressWarnings(
		# https://rpkgs.datanovia.com/ggpubr/reference/ggscatter.html
		gg <- ggscatter(df, x = x, y = y,
			#add = "jitter",
			#palette = "jco",
			color = color,
			fill = fill,
			facet.by = facet.by,
			ncol = facet.ncol,
			panel.labs.font = list(face = NULL, color = NULL, size = font_size, angle = NULL),
			ggtheme = theme_bw() # theme_pubr()
		)
		))
	},
	"violinplot"={
		suppressMessages(suppressWarnings(
		gg <- ggviolin(df, x = x, y = y,
			color = color,
			fill = fill,
			alpha = 0.8,
			add = "boxplot",
			#add = "dotplot",
			#add = c("boxplot", "dotplot"),
			add.params = list(fill = "white", size = 0.5, alpha = 0.7),
			facet.by = facet.by,
			ncol = facet.ncol,
			panel.labs.font = list(face = NULL, color = NULL, size = font_size, angle = NULL),
			ggtheme = theme_bw() # theme_pubr()
		)
		))
	},
	{}
  ) # switch

  if (!is.null(facet.by) & FALSE) {
	# https://rpkgs.datanovia.com/ggpubr/reference/facet.html
	gg <- facet(gg + theme_bw(),
		facet.by = facet.by,
		ncol = facet.ncol,
		#scales = "free", # fixed, free, free_x, free_y
		#short.panel.labs = FALSE,   # Allow long labels in panels
		#panel.labs.background = list(fill = "steelblue", color = "steelblue")
		panel.labs.font = list(face = NULL, color = NULL, size = font_size, angle = NULL)
	)
	#gg <- gg + scale_x_continuous(breaks = get_breaks(n=3))
  }


  if (!is.null(list_comparison)) {
	gg <- gg + stat_compare_means(method = method_comparison, # t.test, wilcox.test
		comparisons = list_comparison,
		symnum.args = symnum.args,
		hide.ns = FALSE,
		step.increase = step.increase,
		label.y = label.y,
		vjust = vjust,
		tip.length = tip.length
	)

	if (!is.null(type_printing_mean)) {
		switch(type_printing_mean,
			"mean1"={
				# https://stackoverflow.com/questions/50202895/add-mean-value-on-boxplot-with-ggpubr
				gg <- gg + stat_summary(fun.data = function(x) data.frame(y=label.y*0.95, label = paste("Mean=",round(mean(x),3))), geom="text") 
			},
			{}
		) # switch
	} # if

  }


  if (!is.null(test)) {
	gg <- gg + stat_compare_means(method = test, # anova, kruskal.test
		label.x = test_label.x,
		label.y = test_label.y,
		#label = "p.signif"
		vjust = test_vjust,
	)
  }

  if (!is.null(color_manual)) {
	gg <- gg + scale_color_manual(values = color_manual)
  }

  if (!is.null(fill_manual)) {
	gg <- gg + scale_fill_manual(values = fill_manual)
  }

  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  if (!is.null(ylim)) gg <- gg + ylim(ylim)

  if (!is.null(ylimits)) {
	gg <- gg + scale_y_continuous(limits=ylimits,
		     breaks = ybreaks,
		     expand = yexpand)
  }

  if (!is.null(yintercept)) {
	gg <- gg + geom_hline(yintercept=yintercept, size=0.8, linetype=2, color="black", alpha=0.5)
  }

  #gg <- gg + annotate("text", x=1.5, y=0, label=c("","","Private","Private"))

  if (!is.null(title)) gg <- gg + ggtitle(title)
  if (!is.null(xlab)) gg <- gg + xlab(xlab)
  if (!is.null(ylab)) gg <- gg + ylab(ylab)

  if (angle == 90) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle, vjust = 0.5, hjust=1))
  } else if (angle > 0) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle, vjust = 1, hjust=1))
  }

  gg <- gg + theme(plot.title=element_text(size=font_size, face = "bold"),
	axis.text.x=element_text(size=font_size),
	axis.text.y=element_text(size=font_size, lineheight=0.9),
	axis.title.x=element_text(size=font_size),
	axis.title.y=element_text(size=font_size),
	#axis.title.y=element_blank(),
	legend.position = "none",
	plot.margin = margin(t=1, r=0, b=1, l=0, "pt"),
	panel.spacing = unit(1, "lines"))


  if (height <= 0) {
	if (!is.null(facet.by)) {
		height <- (length(unique(df[[facet.by]]))/facet.ncol) * 4.0
	}
  }

  str_condition_plot <- tolower(gsub(" ", "_", str_condition_plot))
  filename_pdf <- sprintf("%s_%s", plot_type, str_condition_plot)
  print_figure(gg, width=width, height=height, file=filename_pdf)

  gg

} # print_ggpubr







### volacno plot


# get_volcano_plot
get_volcano_plot <- function(df, th_log2fc=1.0, th_padj=0.05, str_condition="", label.select=NULL, col_log2fc="log2FoldChange", col_padj="padj", col_symbol="rownames_") {

  df$sig <- "NS"
  idx_up <- which(df[,col_log2fc] > th_log2fc & df[,col_padj] < th_padj)
  df[idx_up, "sig"] <- "up"
  idx_dn <- which(df[,col_log2fc] < -th_log2fc & df[,col_padj] < th_padj)
  df[idx_dn, "sig"] <- "dn"

  vec_log2fc <- abs(df[,col_log2fc])
  xmax <- ceiling_dec(max(vec_log2fc, na.rm=T), 0)

  df[,col_padj] <- pmax(df[,col_padj], 1e-99)
  vec_mlog10padj <- -log10(df[,col_padj])
  ymax <- max(vec_mlog10padj, na.rm=T)
  if (ymax < 5) {
    ymax <- ymax_ + 1 
  } else {
    #ymax <- ceiling_dec(ymax, -1) 
    ymax <- plyr::round_any(ymax, 5, f = ceiling) 
    if (ymax == 100) ymax <- 105
  }

  str_title <- sprintf("%s  (up: %d, down: %d)",
		     str_condition, length(idx_up), length(idx_dn))
  ylab_tex <- r'($-log_{10}$ (adjusted pvalue))'
  if (grepl("FDR", col_padj)) {
  	ylab_tex <- r'($-log_{10}$ FDR)'
  }

  gg <- ggplot(data=df,
    aes_string(x=col_log2fc, y=sprintf("-log10(%s)", col_padj), colour="sig")) +
    geom_point(size=2.0, alpha=0.5) +
    ggtitle(str_title) +
    xlim(c(-xmax,xmax)) + ylim(c(0,ymax)) +
    xlab(TeX(r'($log_2$ (fold change))')) +
    ylab(TeX(ylab_tex)) +
    scale_colour_manual(values = c("NS"="gray30", "up"="#ee5500", "dn"="#0055ee"))

  #Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
  #abline(v=0, col="black", lty=3, lwd=1.0)
  gg <- gg + geom_vline(xintercept=-th_log2fc, col="black", lty="22", lwd=1.0)+
	geom_vline(xintercept=th_log2fc, col="black", lty="22", lwd=1.0) +
	#geom_hline(yintercept=-log10(max(df[df[,col_padj]<0.05,col_padj], na.rm=TRUE)), col="black", lty="22", lwd=1.0)
	geom_hline(yintercept=-log10(th_padj), col="black", lty="22", lwd=1.0)

  gg <- gg +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=14),
	axis.title=element_text(size=14,face="bold"),
	legend.title=element_text(size=14), 
	legend.text=element_text(size=14))


  df$nudge_x <- 0.5
  df$nudge_x[idx_dn] <- -0.5 
  df$nudge_y <- -log10(df[,col_padj])/30 + 1		 

  # change name
  if (col_symbol == "rownames_") {
  	df$symbol <- rownames(df)
	col_symbol <- "symbol"
  }

  df[,col_symbol] <- mgsub::mgsub(df[,col_symbol],
	c('replication-dependent',', acting on a heme group of donors'),
	c('rep.-dep.',''))
  df[,col_symbol] <- str_wrap(df[,col_symbol], width=15)

  if (!is.null(label.select)) {
    f <- df[,col_symbol] %in% label.select
    #f <- df[,col_symbol] %in% c("PALM3", "ANKRD22", "BMF")
    #f <- f | df[col_symbol,] %in% c("NIPAL2", "SLC40A1", "TFAP2C", "RET")
    idx <- which(f)
    df <- df[idx,,drop=F]
  }

  if (nrow(df) > 0) {
    gg <- gg +   
      geom_text_repel(data=subset(df, sig=='dn'),
	   aes_string(label=col_symbol),
	   size=2.5,
	   colour='black',
	   nudge_x=subset(df, sig=='dn')$nudge_x,
	   nudge_y=subset(df, sig=='dn')$nudge_y,
	   force=1, point.padding=0.5, min.segment.length=0,
	   segment.color='#555555', segment.alpha=0.5, hjust=1) +    
      geom_text_repel(data=subset(df, sig=='up'),
	   aes_string(label=col_symbol),
	   size=2.5,
	   colour='black',
	   nudge_x=subset(df, sig=='up')$nudge_x,
	   nudge_y=subset(df, sig=='up')$nudge_y, 
	   force=1, point.padding=0.5, min.segment.length=0,
	   segment.color='#555555', segment.alpha=0.5, hjust=0)
  }

  gg

} # get_volcano_plot





# theme_geometry
# original function from https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2
# modified by H. Kim
# input:
#   xvals: values of x that will be plotted
#   yvals: values of y that will be plotted
#   xgeo: x intercept value for y axis
#   ygeo: y intercept value for x axis
#   color: default color for axis
#   linesize: line size for axis
#   xlab: label for x axis
#   ylab: label for y axis
#   labsize: label size for axis
#   labgap: gap betwen axis and axis label
#   ticks: number of ticks to add to plot in each axis
#   textsize: size of text for ticks
#   xlimit: limit value for x axis
#   ylimit: limit value for y axis
#   epsilon: parameter for small space
#   gap_tick_label: gap betwen axis and tick label
#
# usage:
# see fig2.ipynb or r_ggplot_scatterplot.txt
#
theme_geometry <- function(xvals, yvals, xgeo = 0, ygeo = 0, 
	color = "black", linesize = 1, 
	xlab = "x", ylab = "y", labsize=3.5, labgap=2,
	ticks = 10,
	textsize = 3, 
	xlimit = max(abs(xvals),abs(yvals)),
	ylimit = max(abs(yvals),abs(xvals)),
	epsilon = max(xlimit,ylimit)/50, gap_tick_label=3){

  #Create axis 
  xaxis <- data.frame(x_ax = c(-xlimit, xlimit), y_ax = rep(ygeo,2))
  yaxis <- data.frame(x_ax = rep(xgeo, 2), y_ax = c(-ylimit, ylimit))

  #Add axis
  theme.list <- 
  list(
    theme_void(), #Empty the current theme
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = linesize, data = xaxis),
    geom_line(aes(x = x_ax, y = y_ax), color = color, size = linesize, data = yaxis),
    # begin of modification by H. Kim
    #annotate("text", x = xlimit + 2*epsilon, y = ygeo, label = xlab, size = 2*textsize),
    #annotate("text", x = xgeo, y = ylimit + 4*epsilon, label = ylab, size = 2*textsize),
    annotate("text", x = 0, y = -ylimit-labgap, label = xlab, size = labsize),
    annotate("text", x = -xlimit-labgap, y = 0, label = ylab, size = labsize, angle=90),
    # end of modification
    xlim(-xlimit - 7*epsilon, xlimit + 7*epsilon), #Add limits to make it square
    ylim(-ylimit - 7*epsilon, ylimit + 7*epsilon)  #Add limits to make it square
  )

  #Add ticks programatically
  ticks_x <- round(seq(-xlimit, xlimit, length.out = ticks),2)
  ticks_y <- round(seq(-ylimit, ylimit, length.out = ticks),2)

  #Add ticks of x axis
  nlist <- length(theme.list)
  for (k in 1:ticks){

    #Create data frame for ticks in x axis
    xtick <- data.frame(xt = rep(ticks_x[k], 2), 
			yt = c(xgeo + epsilon, xgeo - epsilon))

    #Create data frame for ticks in y axis
    ytick <- data.frame(xt = c(ygeo + epsilon, ygeo - epsilon), 
			yt = rep(ticks_y[k], 2))

    #Add ticks to geom line for x axis
    theme.list[[nlist + 4*k-3]] <- geom_line(aes(x = xt, y = yt), 
					 data = xtick, size = linesize, 
					 color = color)

    #Add labels to the x-ticks
    theme.list[[nlist + 4*k-2]] <- annotate("text", 
					    x = ticks_x[k], 
					    y = ygeo - gap_tick_label*epsilon,
					    size = textsize,
					    label = paste(ticks_x[k]))


    #Add ticks to geom line for y axis
    theme.list[[nlist + 4*k-1]] <- geom_line(aes(x = xt, y = yt), 
					     data = ytick, size = linesize, 
					     color = color)

    #Add labels to the y-ticks
    theme.list[[nlist + 4*k]] <- annotate("text", 
					    x = xgeo - gap_tick_label*epsilon, 
					    y = ticks_y[k],
					    size = textsize,
					    label = paste(ticks_y[k]))
  }

  #Add theme
  #theme.list[[3]] <- 
  return(theme.list)

} # theme_geometry







