#
# convert_symbol.R
# author: H. Kim
# date created: 2021, Oct.
# date last modified: 2021, Dec.
#
#
# called by make_sc-rna-seq_seurat_obj.R
#
#

# file format: ENSG chromosome (1, 2, 3,..., 21, MT, KI270713.1) start end
#GRCH38.annotations <- "reference/genome_annotation/Homo_sapiens.GRCh38.86.txt"
#gtf <- read.delim(GRCH38.annotations, header = F)
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))

# convert_symbol
# usage:
# df_gtf <- read.delim(GRCH38.annotations, header = F)
# df_gtf <- convert.symbol(df_gtf)
# write.table(df_gtf, file_name_gene_order, sep = "\t", row.names = FALSE, col.names = FALSE)
convert.symbol = function(df_gtf) {

	ensembls <- df_gtf$V1
	ensembls <- gsub("\\.[0-9]*$", "", ensembls)
	geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembls, keytype = "GENEID", columns = "SYMBOL")
	df_gtf <- cbind.fill(df_gtf, geneIDs1, fill = NA)
	df_gtf <- na.omit(df_gtf)
	df_gtf$feature <- df_gtf$SYMBOL
	df_gtf.new <- data.frame(df_gtf$SYMBOL,df_gtf$V2,df_gtf$V3,df_gtf$V4)

	# chr1, ..., chr21, chrMT, ..., chrKI270713.1
	df_gtf.new$df_gtf.V2 <- paste("chr",df_gtf.new$df_gtf.V2,sep = "")
	df_gtf.new$df_gtf.SYMBOL <- make.unique(df_gtf.new$df_gtf.SYMBOL)

	# begin of addition by H. Kim
	df_gtf <- df_gtf.new
	colnames(df_gtf) <- c("symbol", "chrom", "start", "end")
	chrOrder <- paste0("chr", c(1:22, "X", "Y", "MT"))
	df_gtf <- df_gtf[df_gtf$chrom %in% chrOrder,]
	df_gtf$chrom <- factor(df_gtf$chrom, chrOrder, ordered=TRUE)
	df_gtf.new <- df_gtf[do.call(order, df_gtf[, c("chrom", "start")]), ]
	# end of addition

      return(df_gtf.new)

} # convert.symbol





