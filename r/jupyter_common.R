#
# jupyter_common.R
# author: H. Kim
# created: 2019, May
# last modified: 2020, Aug.
#
# contents:
# print_figure()
# 
#
#

# silencing the warning message of dplyr.summarise
options(dplyr.summarise.inform=F) 


# memory
# use pryr::object_size(), pryr:mem_used() rather than loading pryr.
#suppressPackageStartupMessages(library(pryr)))

# data structure 

# string
suppressPackageStartupMessages(library(stringr))
#suppressPackageStartupMessages(library(stringi)))
# use qdap::mgsub() instead of loading qdap package
#suppressPackageStartupMessages(library(qdap)))

# factor
#suppressPackageStartupMessages(library(forcats)))

# matrix
# caution: the order of loading matrixStats after Rfast was intended since rowMaxs() needs to be overlayed.
# solution: try to avoid Rfast, but to use matrixStats, e.g. matrixStats::rowMins(), matrixStats::rowMaxs()
#suppressPackageStartupMessages(library(Rfast))) # for colMinsMaxs()
#suppressPackageStartupMessages(library(matrixStats))) # for rowMins(), rowMaxs()
#suppressPackageStartupMessages(library(scales)))

# data.frame
# do not use dplyr since pathview/eg2id function did not work with dplyr
# consider using magrittr for pipe %>% and dplyr::function() instead of loading dplyr
suppressPackageStartupMessages(library(dplyr)) # filter, mutate, rename
suppressPackageStartupMessages(library(magrittr)) # for pipe %>%
suppressPackageStartupMessages(library(tidyr)) # gather
suppressPackageStartupMessages(library(reshape2))


# tibble
#suppressPackageStartupMessages(library(tibble))

# data.table
suppressPackageStartupMessages(library(data.table))

# id conversion
#suppressPackageStartupMessages(library(org.Hs.eg.db)))
#suppressPackageStartupMessages(library(AnnotationDbi)))

# distance
#suppressPackageStartupMessages(library(amap)))

# basic graph
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(ggrepel))
#suppressPackageStartupMessages(library(ggvenn))
suppressPackageStartupMessages(library(ggpubr))

# color
suppressPackageStartupMessages(library(RColorBrewer))

# visualization
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize)) # chordDiagram(), colorRamp2()
#suppressPackageStartupMessages(library(dittoSeq)))

# enrichment analysis
suppressPackageStartupMessages(library(clusterProfiler))
#suppressPackageStartupMessages(library(enrichplot))

# sequence
#suppressPackageStartupMessages(library(Biostrings)))

# display
suppressPackageStartupMessages(library(IRdisplay))



### load subroutines
source("r/jupyter_message.R")
source("r/jupyter_common_plot.R")
source("r/jupyter_common_heatmap.R")


### common parameters

figure_format <- "png"






### color
c1<-brewer.pal(n=9, name = "Set1")
c1[6]<-'#fee287' # light yellow --> yellow
c2<-brewer.pal(n=8, name = "Set2")
c3<-brewer.pal(n=12, name = "Set3")
p1<-brewer.pal(n=9, name = "Pastel1")
p2<-brewer.pal(n=8, name = "Pastel2")
spectral <- brewer.pal(n=11, name = "Spectral")
npg_nrc <- ggsci::pal_npg("nrc", alpha=0.7)(9)
#scales::show_col(npg_nrcl)









nv_row_annot_color <- c(

   # unc project
   "Female BC"=npg_nrc[1],
   "Male BC"=npg_nrc[2],
   "tamoxifen"=npg_nrc[1],
   "control"=npg_nrc[2],

   '49758L'="#A6CEE3",
   '49758L_12hrSUS_E2'="#FDBF6F",
   '49758L_12hrSUS_E2_TAM'="#B2DF8A",

   #'unt48'=c1[2], 'tgfb48'=c1[1],
   #'tgfbWO'=c1[3], 'tgfbCX5461100nm'=c1[4],
   'unt48'="#00BA38",
   'tgfb48'="#F8766D", 
   'tgfbCX5461100nm'="#619CFF",

   # stjude.org project
   'unt'="#00BA38",
   'tgfb'="#F8766D", 
   'tgfbCX'="#619CFF",

   'GNPs'=c1[4],
   'SHH model'="#ff0000",
    
   'Retro-Myc model'="#999900", # dd yellow
   'Retro-Myc model 19568-Cas9-Tumor'="#cccc00", # d yellow
   'CRISPR-Myc model'="#eeee00", # ~yellow
    
   'Utx-E3-Cre+p53DN+Myc'=c1[3], # green
   'Utx-E3-Cre+p53DN+Mycn'="#000000", # black
    
   # dark orange to light orange
   'Retro-Myc+Gfi1 model'=c1[5], # ff7f00,
   'WT-NS-Myc+Gfi1'="#ffaf30",
   'WT-NS-Myc+Gfi1 (secondary from 1232)'="#ffd700", # gold
    
   'WT-Prom1-Myc+Gfi1'=c1[7], # a65628
   'WT-Prom1-Myc+Gfi1 (secondary from 1266)'="#c07050",
   'WT-Prom1-Myc+Gfi1 (secondary from 57163)'="#d09070",
   'WT-NS or Prom1-Myc+Gfi1B'="#cccccc"

) # nv_row_annot_color






# map_row_annot_color
map_row_annot_color <- function(group_u, type_group="default", verbose=F) {
  group_u <- unique(as.character(group_u))

  group_u <- gsub('HALLMARK_', '', group_u)
  col_group <- nv_row_annot_color[group_u]
  names(col_group) <- group_u;
  idx_na <- which(is.na(col_group));
  n_na <- length(idx_na)
  if (n_na > 0) {
    # when any color was not mapped, use palette
    switch(type_group,
	'up'={col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Reds"))(n_na+3)[(n_na+1):2]},
	'dn'={col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Blues"))(n_na+3)[2:(n_na+1)]},
	{

		mtx <- stringr::str_split(group_u[idx_na], pattern=" ", simplify=T)
		rownames(mtx) <- group_u[idx_na]	
		for (j in 1:ncol(mtx)) {
			idx <- which(nchar(mtx[,j]) > 0)
			col_group[idx_na[idx]] <- nv_row_annot_color[mtx[idx,j]]
			idx_na <- which(is.na(col_group));
			if (length(idx_na) == 0) break
		}
		idx_na <- which(is.na(col_group));
		if (length(idx_na) > 0) {
			col_group[idx_na] <- colorRampPalette(brewer.pal(9,"Blues"))(n_na+3)[2:(n_na+1)]
		}
	}

    )
  } # if

  if (verbose) {
    display(col_group)
  } # if

  return(col_group)

} # map_row_annot_color








# add_ortholog_info
# add columns of ortholog information (e.g., HomoloGene.ID, mouse.sym, mouse.eid, human.sym, human.eid)
add_ortholog_info <- function(df, from="mouse", to="human") {

  df_from <- orthodf[orthodf$specie==from,]
  idx <- match(rownames(df), df_from[,"Symbol"])
  df$HomoloGene.ID <- df_from[idx,"HomoloGene.ID"]
  df[, sprintf("%s.sym",from)] <- df_frome[idx,"Symbol"]
  df[, sprintf("%s.eid",from)] <- df_from[idx,"EntrezGene.ID"]

  df_to <- orthodf[orthodf$specie==to,]
  idx <- match(df$HomoloGene.ID , df_to$HomoloGene.ID)
  df[, sprintf("%s.sym",to)] <- df_to[idx,"Symbol"]
  df[, sprintf("%s.eid",to)] <- df_to[idx,"EntrezGene.ID"]

  return(df)

} # add_ortholog_info















# print_figure
#
# input:
#   padding: unit of padding values correspond to the bottom, left, top and right e.g. unit(c(2, 5, 2, 2), "mm")
#
# usage:
# print_figure(list_out$list_ht, width=8, height=8, padding = unit(c(2, 20, 2, 2), "mm"), file=sprintf("heatmap_%s_%s",str_condition, list_out$type_col_short))
print_figure <- function(obj, width, height, padding=NULL, file=NULL, annotation_legend_side="bottom", heatmap_legend_side="right", figure_format_=figure_format, resolution=300, f_display2screen=T) {

  if (is.null(obj) || (width <= 0) || (height <= 0)) return()

  type_obj=class(obj)[1]

  if (!is.null(file) && (length(file) > 0)) {

    # save figure to a file
    switch(figure_format_,
	"pdf"={
    		filename <- sprintf('pdf/%s.pdf', file);
    		cairo_pdf(filename, width=width, height=height)      
	},
	"png"={ 
    		filename <- sprintf('png/%s.png', file);
    		png(filename, width=width, height=height, units='in', res=resolution)
	},
	"tiff"={ 
    		filename <- sprintf('tiff/%s.tiff', file);
    		tiff(filename, width=width, height=height, units='in', res=resolution, compression='lzw', type='cairo')
	},
	{}
    ) # switch

    switch(type_obj,

      "Heatmap"={
	if (is.null(padding)) {
		obj <- draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side)
	} else {
		obj <- draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side, padding=padding)
	}
      },

      "HeatmapList"={
	if (is.null(padding)) {
		obj <- draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side)
	} else {
		obj <- draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side, padding=padding)
	}
      },

      "list"={
         # print gg
         obj <- obj$gg
      },
      { }
    ) # switch

    suppressMessages(suppressWarnings(print(obj)))
    dev.off()  
  } 

  if (f_display2screen) {
    options(repr.plot.width=width, repr.plot.height=height)
    switch(type_obj,
      "Heatmap"={
         draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side)
      },
      "HeatmapList"={
         draw(obj, annotation_legend_side=annotation_legend_side, heatmap_legend_side=heatmap_legend_side)
      },
      "list"={
         grid.draw(obj$gt)
      },
      { 
	suppressMessages(suppressWarnings(print(obj)))
      }
    )
  } # f_display2screen

  #IRdisplay::display_html('')  
  cat(sprintf(""))
}



##### dge

# usage:
# list_venn <- list()
# list_venn[['ours']] <- sym_up; list_venn[['elife2018']] <- sym_up_nmumg
# vec_color <- c(c1[2],c1[1],c1[4])
# gg <- get_euler_diagram_from_list(list_venn, shape='circle', vec_color=vec_color, fills=list(fill=vec_color, alpha=0.6), fill_alpha=0.6, edges=list(), labels=list(cex=0.8), quantities=list(cex=0.7), strips=list(), legend=FALSE, label_pos_x=NULL, label_pos_y=NULL, main='up genes in in uint48 vs. tgfb48 (transcription)', main.gp=list(cex=0.75))
# print_figure(gg, width=3.1, height=2.8, file='euler_diagram_up_genes')
#
get_euler_diagram_from_list <- function(list_venn, shape='circle', vec_color=NULL, fills=NULL, fill_alpha=0.6, edges=list(), labels=list(cex=0.8), quantities=list(cex=0.7), strips=list(), legend=FALSE, label_pos_x=NULL, label_pos_y=NULL, quantities_pos_x=NULL, quantities_pos_y=NULL, main=NULL, main.gp=list(cex=0.75), ...) {

  mtx <- conv_list2mtx(list_venn)
  if (is.null(fills)) {
    vec_color <- get_colors_from_list(list_venn, vec_color)
    fills <- list(fill=vec_color, alpha=fill_alpha)
  }

  require(eulerr)
  fit <- euler(mtx, shape=shape)
  gg <- plot(fit, fills = fills,
           edges = edges,
           labels = labels,
           quantities = quantities,
           strips = strips, legend = legend, main = main, ...)

  if (!is.null(main)) {
    t <- getGrob(gg,'main.grob')
    gg <- editGrob(gg, 'main.grob', gp=main.gp)
  }

  if (!is.null(label_pos_x)) {
    #grid.ls(gg)
    t <- getGrob(gg,'labels.grob')
    x <- t$x; y <- t$y
    for (i in 1:length(label_pos_x)) {
      x[i] <- t$x[i] + unit(label_pos_x[i], "mm")
      y[i] <- t$y[i] + unit(label_pos_y[i], "mm")
    }
    gg <- editGrob(gg,"labels.grob", x=x, y=y)
  }

  if (!is.null(quantities_pos_x)) {
    #grid.ls(gg)
    t <- getGrob(gg,'quantities.grob')
    x <- t$x; y <- t$y
    for (i in 1:length(quantities_pos_x)) {
      x[i] <- t$x[i] + unit(quantities_pos_x[i], "mm")
      y[i] <- t$y[i] + unit(quantities_pos_y[i], "mm")
    }
    gg <- editGrob(gg,"quantities.grob", x=x, y=y)
  }


  return(gg)
}



#https://github.com/cran/VennDiagram/blob/master/R/hypergeometric.test.R
# This function performs the hypergeometric test on the two categories. Taken from package BoutrosLab.statistics.general
# reference:
# http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
#
# usage:
# (1) under-representation: calculate.overlap.and.pvalue(sym_mrna_up, sym_mrna_up_elife, total.size=n.total.coding.genes, lower.tail = TRUE, adjust=FALSE) # under-representation, use lower.tal=TRUE (default) and adjust=FALSE, since probabilities are P[X ≤ x].
# (2) over-representation: calculate.overlap.and.pvalue(sym_mrna_up, sym_mrna_up_elife, total.size=n.total.coding.genes, lower.tail = FALSE, adjust=TRUE) # over-representation, use lower.tail=FALSE, adjust=TRUE, subtract x by 1, when P[X ≥ x] is needed.
calculate.overlap.and.pvalue = function(list1, list2, total.size, lower.tail = TRUE, adjust = FALSE) {

        # calculate actual overlap
        actual.overlap <- length(intersect(list1, list2));

        # calculate expected overlap
        # need to cast to avoid integer overflow when length(list1) * length(list2) is extremely large
        expected.overlap <- as.numeric(length(list1)) * length(list2) / total.size;

        adjust.value <- 0;

        # adjust actual.overlap to reflect P[X >= x]
        if (adjust & !lower.tail) {
                adjust.value <- 1;
                #warning('Calculating P[X >= x]');
                }

        # calculate significance of the overlap
        overlap.pvalue <- phyper(
                q = actual.overlap - adjust.value,
                m = length(list1),
                n = total.size - length(list1),
                k = length(list2),
                lower.tail = lower.tail
                );

        # return values
        return( c(actual.overlap, expected.overlap, overlap.pvalue) );

}


# https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
fisher.test.overlap.and.pvalue = function(list1, list2, total.size) {

    #contingency table:
    mtx <- matrix(c(total.size - length(union(list1,list2)),
                    length(setdiff(list1,list2)),
                    length(setdiff(list2,list1)),
                    length(intersect(list1,list2))), nrow=2)

    return( fisher.test(mtx, alternative="greater"))

}







##### utility


# https://stackoverflow.com/questions/35807523/r-decimal-ceiling
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)

# https://stackoverflow.com/questions/35807523/r-decimal-ceiling
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)







# row.match
# input:
#   nomatch: the value to be returned in the case when no match is found.
# output:
#   see help(match)
# reference:
# https://github.com/tagteam/prodlim/blob/master/R/row.match.R
#
row.match <- function(x, table, nomatch=NA){

  if (class(table)=="matrix") table <- as.data.frame(table)
  if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
  cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
  ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
  match(cx, ct, nomatch=nomatch)

}

conv_list2mtx <- function(list_strings) {

  sym <- Reduce(union, list_strings)
  n_elements <- length(list_strings)
  mtx <- matrix(FALSE, length(sym), n_elements,
           dimnames=list(sym,names(list_strings)))
  for (i in 1:n_elements) {
     f <- sym %in% list_strings[[i]]
     mtx[,i] <- f
  }

  return(mtx)
}


##' This function retrieves the indices of non-zero elements in sparse matrices
##' of class dgCMatrix from package Matrix. This function is largely inspired from 
##' the package \code{Ringo}
##' 
##' @title Retrieve the indices of non-zero elements in sparse matrices
##' @param x A sparse matrix of class dgCMatrix
##' @return A two-column matrix
##' @author Samuel Wieczorek
##' @examples
##' library(Matrix)
##' mat <- Matrix(c(0,0,0,0,0,1,0,0,1,1,0,0,0,0,1),nrow=5, byrow=TRUE, sparse=TRUE)
##' res <- nonzero(mat)
nonzero <- function(x){
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix
    
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
                      dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
}

# matchlast
# reference:
# https://stackoverflow.com/questions/37404084/return-last-match-from-vector
matchlast <- function (needles, haystack) {
   length(haystack) + 1L - match(needles, rev(haystack))
}



# split_string_with_length
# reference: 
# https://stackoverflow.com/questions/11619616/how-to-split-a-string-into-substrings-of-a-given-length
split_string_with_length <- function(text, n) {
    substring(text, seq(1, nchar(text)-(n-1), n), seq(n, nchar(text), n))
}



# conv_codon2aa
conv_codon2aa <- function(vec_codon, f_aaa=FALSE, f_append_aa=FALSE, f_append_codon=FALSE) {

  mycode <- GENETIC_CODE
  #Selenocysteine  Sec/U   UGA     https://en.wikipedia.org/wiki/Selenocysteine
  mycode[["TGA"]] <- "U"
  #Pyrrolysine     Pyl/O   UAG     https://en.wikipedia.org/wiki/Pyrrolysine
  mycode[["TAG"]] <- "O"

  suppressWarnings( vec_aa  <-  Biostrings::translate( DNAStringSet(vec_codon), genetic.code=mycode, no.init.codon=FALSE, if.fuzzy.codon="error" ) )

  if (f_aaa) {
    vec_aa <- conv_aa2aaa(vec_aa, f_append_aa=f_append_aa)
  } 
  if (f_append_codon) {
    paste(vec_aa, vec_codon, sep="_")
  } else {
    vec_aa
  }

} # conv_codon2aa


# conv_aa2aaa
conv_aa2aaa <- function(vec_aa, f_append_aa=FALSE) {

  vec_aa <- as.character(vec_aa)
  # convet amino-acid one-letter code into the three-letter one.
  suppressWarnings( vec_aaa <- aaa(vec_aa) )

  # 21st and 22nd amino acids, alternative stop codons
  idx <- match("U", vec_aa)
  if (length(idx) > 0) vec_aaa[idx] <- "Sec"
  idx <- match("O", vec_aa)
  if (length(idx) > 0) vec_aaa[idx] <- "Pyl"

  if (f_append_aa) {
    paste(vec_aaa, vec_aa, sep="")
  } else {
    vec_aaa
  }

} # conv_aa2aaa



# order_by_pattern
order_by_pattern <- function(vec, pattern) {

  idx <- c()
  items <- strsplit(pattern, "\\|")[[1]]
  for (item in items) {
	idx <- c(idx, grep(item, vec))
  }

  idx

} # order_by_pattern


# verb
verb <- function(...) cat(sprintf(...), sep='', file=stdout())






