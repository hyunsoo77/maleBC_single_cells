#
# enrichment_analysis.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2022, Jan.
#
# contents:
# open_xlsx():
# add_worksheet_gsea_results()
# add_worksheet_enricher_results()
# write_xlsx()
#
# execute_enrichment_analysis()
#
# usage:
# source("enrichment_analysis.R")
#
# reference:
#


suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(IRdisplay)) # qplots::space(), S4Vectors::space()





"H	hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
C1	positional gene sets  for each human chromosome and cytogenetic band.
C2	curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.
C3	regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.
C4	computational gene sets  defined by mining large collections of cancer-oriented microarray data.
C5	ontology gene sets  consist of genes annotated by the same ontology term.
C6	oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.
C7	immunologic signature gene sets  represent cell states and perturbations within the immune system.
C8	cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue."



dir_gmt <- "./reference/gmt"

collections <- c("h", "c2", "c5", "c6")


list_gmt <- list()
for (collection in collections) {
    fname_gmt <- sprintf("%s/%s.all.v6.1.symbols.gmt", dir_gmt, collection)
    list_gmt[[collection]] <- read.gmt(fname_gmt)

    #do.call(log_cmd, list(head(list_gmt[[collection]])) )
}



n_sheet <- 0

# open_xlsx
open_xlsx <- function() {

  wb <- createWorkbook()
  n_sheet <<- 0

  wb

} # open_xlsx







seed_gsea <- 51
pattern_enrich_result_exclude <- "^ID$|Count"



# add_worksheet_gsea_results
#
# input:
#   res: data.frame, cols={"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
#   collections: names(list_gmt), list_gmt[[term1]] <- read.gmt(filename_gmt)
#   log_cmd: {"display", "print"}
#   log_cmd_txt: {"display_html", "print"}
#
# output:
#   list_gsea: list
#     df_up:
#     df_dn:
#     collections:
#
# usage:
# wb <- open_xlsx()
# list_gsea <- add_worksheet_gsea_results(wb, markers, col_log2fc="avg_log2FC", col_pvalue="p_val", col_padj="p_val_adj", th_log2fc=0.25, th_padj=0.05)
# write_xlsx(wb, sprintf("xlsx/%s.xlsx", str_condition))
add_worksheet_gsea_results <- function(wb, res, collections=NULL, col_log2fc="log2FoldChange", col_pvalue="pvalue", col_padj="padj", th_log2fc=log2(2), th_padj=0.05, th_gsea_padj=0.05, log_cmd="display", log_cmd_txt="display_html") {

  list_out <- list()

  idx_up <- which((res[,col_log2fc] > th_log2fc) & (res[,col_padj] < th_padj))
  df_up <- as.data.frame(res[idx_up,])
  df_up <- df_up[order(df_up[,col_padj]),]

  do.call(log_cmd_txt, list(sprintf("df_up")) )
  do.call(log_cmd, list(head(df_up)) )
  do.call(log_cmd, list(dim(df_up)) )

  idx_dn <- which((res[,col_log2fc] < -th_log2fc) & (res[,col_padj] < th_padj))
  df_dn <- as.data.frame(res[idx_dn,])
  df_dn <- df_dn[order(df_dn[,col_padj]),]

  do.call(log_cmd_txt, list(sprintf("df_dn")) )
  do.call(log_cmd, list(head(df_dn)) )
  do.call(log_cmd, list(dim(df_dn)) )


  list_out[["df_up"]] <- df_up
  list_out[["df_dn"]] <- df_dn


  # take care of pvalue=0
  pvalues <- res[,col_pvalue]
  min_pvalues <- min(pvalues[pvalues > 0], na.rm=T)
  pvalues[pvalues == 0] <- min_pvalues*0.1

  nv <- sign(res[,col_log2fc])*(-log10(pvalues))
  names(nv) <- rownames(res);
  nv <- nv[!is.na(nv)]

  #head(nv)
  #tail(nv)
  #length(nv)


  sheet_name <- "up"
  if (!sheet_name %in% names(wb)) {
    n_sheet <<- n_sheet + 1
    addWorksheet(wb, sheetName = sheet_name, gridLines = TRUE)
    writeDataTable(wb, sheet = n_sheet, x = df_up, colNames = TRUE, rowNames = TRUE)
    #setColWidths(wb, sheet = 1, cols = 1:(n + 1), widths = c(15, 50, 25, rep(15, 7), rep(12, 3), rep(15, 7), rep(12, 3)))
    # style1 <- createStyle(valign = 'center', wrapText = TRUE) addStyle(wb, sheet =
    # 1, style1, rows = 1:(m+1), cols = n+1, gridExpand = TRUE) save xlsx
  }

  sheet_name <- "down"
  if (!sheet_name %in% names(wb)) {
    n_sheet <<- n_sheet + 1
    addWorksheet(wb, sheetName = sheet_name, gridLines = TRUE)
    writeDataTable(wb, sheet = n_sheet, x = df_dn, colNames = TRUE, rowNames = TRUE)
  }

  if (log_cmd_txt == "display_html") {
	display_html("<hr style=\"height:1px\">")
  } 

  if (is.null(collections)) {
	collections <- names(list_gmt)
  }

  for (collection in collections) {
    
    set.seed(seed_gsea)
    
    sheet_name_collection <- sprintf("gsea %s", collection)
    do.call(log_cmd_txt, list(sheet_name_collection) )

    
    suppressMessages(suppressWarnings(
	# https://rdrr.io/bioc/clusterProfiler/man/GSEA.html
	gsea_result <- GSEA(sort(nv[!is.na(nv)], decreasing=T),
		  exponent =1,
                  minGSSize = 10, maxGSSize = 500,   
                  pvalueCutoff = th_gsea_padj, # adjusted pvalue cutoff
		  pAdjustMethod = "BH", 
                  TERM2GENE=list_gmt[[collection]],
                  TERM2NAME=NA,
                  seed=FALSE, by="fgsea")
	))

    if (!is.null(gsea_result) && (nrow(gsea_result) > 0)) {

      do.call(log_cmd, list(dim(gsea_result)) )
      f_up <- (gsea_result$NES > 0)
      f_dn <- (gsea_result$NES < 0)

      log_txt("\tup")
      do.call(log_cmd, list(head(gsea_result[f_up,])) )
      do.call(log_cmd, list(dim(gsea_result[f_up,])) )

      log_txt("\tdown")
      do.call(log_cmd, list(head(gsea_result[f_dn,])) )
      do.call(log_cmd, list(dim(gsea_result[f_dn,])) )

    } # if

    list_out[[collection]] <- gsea_result

    if (!sheet_name_collection %in% names(wb)) {
	# gsea h
	addWorksheet(wb, sheetName = sheet_name_collection, gridLines = TRUE)
	n_sheet <<- n_sheet + 1

	if (!is.null(gsea_result)) {
		df_gsea <- as.data.frame(gsea_result)
		f_include <- !grepl("ID", colnames(df_gsea))
		df_gsea <- df_gsea[,f_include]
		writeDataTable(wb, sheet = n_sheet, x = df_gsea, colNames = TRUE, rowNames = FALSE)
		setColWidths(wb, sheet = n_sheet, cols = 1, widths = 50)
	}

    } # if 
    
  } # for


  list_out$wb <- wb

  list_out


} # add_worksheet_gsea_results







# add_worksheet_enricher_results
#
# input:
#   res: data.frame, cols={"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
#   collections: names(list_gmt), list_gmt[[term1]] <- read.gmt(filename_gmt)
#   th_log2fc: abs(log2fc) > th_log2fc
#   th_padj: pajd < th_padj
#   log_cmd: {"display", "print"}
#   log_cmd_txt: {"display_html", "print"}
#
# output:
#   list_out: list
#     df_up:
#     df_dn:
#     collections:
#
# usage:
# wb <- open_xlsx()
# list_enrich_result <- add_worksheet_enricher_results(wb, markers, col_log2fc="avg_log2FC", col_pvalue="p_val", col_padj="p_val_adj", th_log2fc=0.25, th_padj=0.05)
# write_xlsx(wb, sprintf("xlsx/%s.xlsx", str_condition))
add_worksheet_enricher_results <- function(wb, res, collections=NULL, col_log2fc="log2FoldChange", col_pvalue="pvalue", col_padj="padj", th_log2fc=log2(2), th_padj=0.05, th_enricher_padj=0.05, th_enricher_qval=0.05, log_cmd="display", log_cmd_txt="display_html") {


  list_out <- list()

  idx_up <- which((res[,col_log2fc] > th_log2fc) & (res[,col_padj] < th_padj))
  df_up <- as.data.frame(res[idx_up,])
  df_up <- df_up[order(df_up[,col_padj]),]

  do.call(log_cmd_txt, list(sprintf("df_up")) )
  do.call(log_cmd, list(head(df_up)) )
  do.call(log_cmd, list(dim(df_up)) )

  idx_dn <- which((res[,col_log2fc] < -th_log2fc) & (res[,col_padj] < th_padj))
  df_dn <- as.data.frame(res[idx_dn,])
  df_dn <- df_dn[order(df_dn[,col_padj]),]

  do.call(log_cmd_txt, list(sprintf("df_dn")) )
  do.call(log_cmd, list(head(df_dn)) )
  do.call(log_cmd, list(dim(df_dn)) )

  list_out[["df_up"]] <- df_up
  list_out[["df_dn"]] <- df_dn


  # up
  sheet_name <- "up"
  gene_symbols <- rownames(df_up)

  if (log_cmd_txt == "display_html") {
	display_html("<hr style=\"height:1px\">")
  } 

  if (!sheet_name %in% names(wb)) {
    # up
    addWorksheet(wb, sheetName = sheet_name, gridLines = TRUE)
    n_sheet <<- n_sheet + 1
    #writeDataTable(wb, sheet = n_sheet, x = data.frame(symbols=gene_symbols), colNames = TRUE, rowNames = FALSE)
     df_up <- tibble::rownames_to_column(df_up, "gene")
    writeDataTable(wb, sheet = n_sheet, x = df_up, colNames = TRUE, rowNames = FALSE)
  }

  if (is.null(collections)) {
	collections <- names(list_gmt)
  }

  for (collection in collections) {
    
    set.seed(seed_gsea)
    
    sheet_name_collection <- sprintf("%s %s", sheet_name, collection)
    do.call(log_cmd_txt, list(sheet_name_collection) )

    if (length(gene_symbols) > 5) {

      # https://rdrr.io/bioc/clusterProfiler/man/enricher.html
      enrich_result <- enricher(gene_symbols,
                   pvalueCutoff = th_enricher_padj, # adjusted pvalue cutoff on enrichment tests to report
                   pAdjustMethod = "BH",
                   universe = NULL,
                   minGSSize = 10,
                   maxGSSize = 500,
                   qvalueCutoff = th_enricher_qval, # qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
                   TERM2GENE = list_gmt[[collection]],
                   TERM2NAME = NA)

      if (is.null(enrich_result)) next

      if (!is.null(enrich_result) && (nrow(enrich_result) > 0)) {
        do.call(log_cmd, list(head(enrich_result)) )
        do.call(log_cmd, list(dim(enrich_result)) )
      }

      list_out[[sheet_name_collection]] <- enrich_result

      if (!sheet_name_collection %in% names(wb)) {
	# up h
        addWorksheet(wb, sheetName = sheet_name_collection, gridLines = TRUE)
        n_sheet <<- n_sheet + 1

        if (!is.null(enrich_result)) {
          df_enrich_result <- as.data.frame(enrich_result)
          f_include <- !grepl(pattern_enrich_result_exclude, colnames(df_enrich_result))    
          df_enrich_result <- df_enrich_result[,f_include]
          writeDataTable(wb, sheet = n_sheet, x = df_enrich_result, colNames = TRUE, rowNames = FALSE)
          setColWidths(wb, sheet = n_sheet, cols = 1, widths = 50)  
        }

      } # if

    } # if
    
  } # for




  # down
  sheet_name <- "down"
  gene_symbols <- rownames(df_dn)

  if (log_cmd_txt == "display_html") {
	display_html("<hr style=\"height:1px\">")
  } 

  if (!sheet_name %in% names(wb)) {
    # down
    addWorksheet(wb, sheetName = sheet_name, gridLines = TRUE)
    n_sheet <<- n_sheet + 1
    #writeDataTable(wb, sheet = n_sheet, x = data.frame(symbols=gene_symbols), colNames = TRUE, rowNames = FALSE)
     df_dn <- tibble::rownames_to_column(df_dn, "gene")
    writeDataTable(wb, sheet = n_sheet, x = df_dn, colNames = TRUE, rowNames = FALSE)
  }

  for (collection in collections) {
    
    set.seed(seed_gsea)
    
    sheet_name_collection <- sprintf("%s %s", sheet_name, collection)
    do.call(log_cmd_txt, list(sheet_name_collection) )

    if (length(gene_symbols) > 5) {

      enrich_result <- enricher(gene_symbols,
                   pvalueCutoff = th_enricher_padj, # adjusted pvalue cutoff on enrichment tests to report
                   pAdjustMethod = "BH",
                   universe = NULL,
                   minGSSize = 10,
                   maxGSSize = 500,
                   qvalueCutoff = th_enricher_qval, # qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
                   TERM2GENE = list_gmt[[collection]],
                   TERM2NAME = NA)

      if (is.null(enrich_result)) next

      if (!is.null(enrich_result) && (nrow(enrich_result) > 0)) {
        do.call(log_cmd, list(head(enrich_result)) )
        do.call(log_cmd, list(dim(enrich_result)) )
      }

      list_out[[sheet_name_collection]] <- enrich_result

      if (!sheet_name_collection %in% names(wb)) {
	# down h
        addWorksheet(wb, sheetName = sheet_name_collection, gridLines = TRUE)
        n_sheet <<- n_sheet + 1

        if (!is.null(enrich_result)) {
          df_enrich_result <- as.data.frame(enrich_result)
          f_include <- !grepl(pattern_enrich_result_exclude, colnames(df_enrich_result))    
          df_enrich_result <- df_enrich_result[,f_include]
          writeDataTable(wb, sheet = n_sheet, x = df_enrich_result, colNames = TRUE, rowNames = FALSE)
          setColWidths(wb, sheet = n_sheet, cols = 1, widths = 50)  
        }

      } # if

    } # if
    
  } # for


  list_out$wb <- wb

  list_out

} # add_worksheet_enricher_results







# write_xlsx
write_xlsx <- function(wb, filename_xlsx) {

  log_txt(sprintf("filename_xlsx: %s\n", filename_xlsx))
  saveWorkbook(wb, filename_xlsx, overwrite = TRUE)

} # write_xlsx















### wrappers


# execute_enrichment_analysis
#
# input:
#   th_log2fc: abs(log2fc) > th_log2fc
#   th_padj: pajd < th_padj
#
# output:
#   list_ea$df_up:
#   list_ea$df_dn:
#   gsea:
#      h, c2, c5, c6
#   enricher:
#      up h, up c2, up c5, up c6
#      down h, down c2, down c5, down c6
# usage:
# list_ea <- execute_enrichment_analysis(markers, str_condition_tmp)
#
execute_enrichment_analysis <- function(markers, str_condition=NULL, col_log2fc="avg_log2FC", col_pvalue="p_val", col_padj="p_val_adj", th_log2fc=th_log2fc_global, th_padj=th_padj_global, th_gsea_padj=0.05, th_enricher_padj=0.05, th_enricher_qval=0.05, method_dge="seurat_findmarkers_enricher", dir_xlsx="xlsx", n_log=0, log_cmd="display", log_cmd_txt="display_html") {

  if (is.null(markers)) {
	return(NULL)
  }

  if (nrow(markers) == 0) {
	return(NULL)
  }

  cols <- colnames(markers)
  if (!col_log2fc %in% cols) {
	idx <- grep("log[0-9]*FC", cols)
	col_log2fc <- cols[idx]
  }

  if (!col_pvalue %in% cols) {
	idx <- grep("[Pp][-]*val(ue)*", cols)
	col_pvalue <- cols[idx]
  }

  if (!col_padj %in% cols) {
	idx <- grep("adj", cols)
	col_padj <- cols[idx]
  }

  if (n_log > 0) {
    log_txt(sprintf("%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n\t\t%s\n",
		method_dge,
    		sprintf("\t\tcol_log2fc=%s", col_log2fc),
    		sprintf("\t\tcol_pvalue=%s", col_pvalue),
    		sprintf("\t\tcol_padj=%s", col_padj),
    		sprintf("\t\tth_log2fc=%g", th_log2fc),
    		sprintf("\t\tth_padj=%g", th_padj)), log_cmd_txt)
  } # if

  wb <- open_xlsx()

  list_ea <- list()
  if (grepl("gsea", method_dge)) {
    list_gsea_result <- add_worksheet_gsea_results(wb, markers,
                    col_log2fc=col_log2fc,
                    col_pvalue=col_pvalue,
                    col_padj=col_padj,
                    th_log2fc=th_log2fc,
		    th_padj=th_padj,
		    th_gsea_padj=th_gsea_padj)
    list_ea <- c(list_ea, list_gsea_result)
  }

  if (grepl("enricher", method_dge)) {
    list_enrich_result <- add_worksheet_enricher_results(wb, markers,
                    col_log2fc=col_log2fc,
                    col_pvalue=col_pvalue,
                    col_padj=col_padj,
                    th_log2fc=th_log2fc,
		    th_padj=th_padj,
		    th_enricher_padj=th_enricher_padj,
		    th_enricher_qval=th_enricher_qval)
    list_ea <- c(list_ea, list_enrich_result)
  } 

  if (!is.null(str_condition)) {
    dir.create(dir_xlsx, showWarnings = FALSE, recursive = TRUE)
    str_condition <- gsub(" ", "_", str_condition)
    write_xlsx(wb, sprintf("%s/%s.xlsx", dir_xlsx, str_condition))
  }

  list_ea

} # execute_enrichment_analysis










### other utilities




# GSA.read.gmt_modified
# modified by H. Kim
# date last modified: 2021, Dec.
#
# modification log:
# (1) comment.char
# (2) slient mode not to log messages in scan()
# (3) names(genesets) <- geneset.names
#
# reference:
# https://rdrr.io/cran/GSA/src/R/GSA.read.gmt.R
GSA.read.gmt_modified <- function(filename, comment.char = "") {

  # Read in and parse a gmt file (gene set file) from the  Broad institute
  # this is tricky, because each lines (geneset) has a variable length
  #  I read the file twice, first to pick up the geneset name and description
  # in the   first two  columns, then I read it all in as a long string

  # The beginning and end of each gene set in the string
  # is determined by matching
  # BOTH on  geneset name and description (since geneset names sometimes
  # occur as genenames elsewhere in the file)

  # begin of modification by H. Kim
  # https://rdrr.io/r/base/scan.html
  #a=scan(filename,what=list("",""),sep="\t", quote=NULL, fill=T, flush=T,multi.line=F)
  a <- scan(filename,what=list("",""), sep="\t", quote=NULL, fill=T, flush=T, quiet=TRUE, multi.line=F, comment.char = comment.char)

  geneset.names <- a[1][[1]]
  geneset.descriptions <- a[2][[1]]

  #dd=scan(filename, what="", sep="\t", quote=NULL)
  dd <- scan(filename, what="", sep="\t", quote=NULL, quiet=TRUE, comment.char = comment.char)
  # end of modification

  nn <- length(geneset.names)
  n <- length(dd)
  ox <- rep(NA, nn)

  ii <- 1
  for(i in 1:nn) {
    #cat(i) # removed by H. Kim
    while((dd[ii] != geneset.names[i]) | (dd[ii+1] != geneset.descriptions[i]) ){
	ii <- ii + 1
    } # while
    ox[i] <- ii
    ii <- ii + 1
  } # for

  genesets <- vector("list",nn)

  for (i in 1:(nn-1)) {
	# cat(i,fill=T)
	i1 <- ox[i]+2
	i2 <- ox[i+1]-1
	geneset.descriptions[i] <- dd[ox[i]+1]
	genesets[[i]] <- dd[i1:i2]
  } # for

  geneset.descriptions[nn] <- dd[ox[nn]+1]
  genesets[[nn]] <- dd[(ox[nn]+2):n]

  # begin of addition by H. Kim
  names(genesets) <- geneset.names
  # end of addition

  out <- list(genesets=genesets, geneset.names=geneset.names, geneset.descriptions=geneset.descriptions)
  class(out) <- "GSA.genesets"

  return(out)

} # GSA.read.gmt_modified





