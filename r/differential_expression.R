#
# differential_expression.R
# modified by H. Kim
# date last modified: 2021, Oct
#
# contents:
#
# usage:
# source("differential_expression.R")
#
# refernce:
#

suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(sva)))



# https://github.com/vd4mmind/RNASeq/blob/master/svaBatchCor.R
svaBatchCor <- function(dat, mmi, mm0,n.sv=NULL){
  dat <- as.matrix(dat)
  Y <- t(dat)
  if(is.null(n.sv))   n.sv <- num.sv(dat,mmi,method="leek")
  o <- svaseq(dat,mmi,mm0,n.sv=n.sv)
  W <- o$sv
  alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
  o$corrected <- t(Y - W %*% alpha)
  return(o)
} # svaBatchCor






#
# select_most_variable_genes
# input:
#    cutoff_perc: select 25% of genes with highest relative standard deviation
#    ntop_genes: the number of the most variable genes
#
#
select_most_variable_genes <- function(mtx, method, str_condition, cutoff_perc=0.25, ntop_genes=3000) {

  list_out <- list()

  switch(method,
	"rowvars"={
		# variance estimates for each row in a matrix
		rv <- rowVars(mtx)
		o <- order(rv, decreasing=TRUE)[1:ntop_genes]
		mtx <- mtx[o, ]
		if (!grepl("top", str_condition)) {
			str_condition <- sprintf("%s.top_%d_most_variable_genes", str_condition, ntop_genes)
		}
	},
	"rsvd"={
		vec_rstd <- apply(mtx, 1, function(x) sd(x) / mean(x))
		#summary(vec_rstd)
		#summary(vec_rstd[rank(vec_rstd) / length(vec_rstd) > 1 - cutoff])

		f_select <- (rank(vec_rstd) / length(vec_rstd)) > (1 - cutoff)
		#min(vec_rstd[f_select])
		mtx <- mtx[f_select,]
		if (!grepl("top", str_condition)) {
			str_condition <- sprintf("%s.top_%d%%%%_most_variable_genes", str_condition, cutoff*100)
		}
	},
	{
	}
  ) # switch

  list_out$mtx <- mtx
  list_out$str_condition <- str_condition

  list_out

}







