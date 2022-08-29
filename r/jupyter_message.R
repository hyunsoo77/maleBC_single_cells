#
# jupyter_message.R
# author: H. Kim
# created: 2019, May
# last modified: 2020, Aug.
#
# contents:
# .checkJupyter()
# log_obj()
# log_txt()
# 
#




# .checkJupyter
# reference:
#   ArchR-1.0.1/R/GlobalDefaults.R
.checkJupyter <- function(){
  tryCatch({
    sysID <- Sys.getenv("JPY_PARENT_PID")
    if(!is.character(sysID)){
      return(FALSE)
    }
    if(sysID == ""){
      FALSE
    }else{
      TRUE
    }
  },error= function(e){
    FALSE
  })
} # .checkJupyter






# log_obj
# 
# input:
#   obj:
#   log_cmd: {["display"], "print"}
#
# usage:
# log_obj(head(df))
# log_obj(head(df), tab=2)
# log_obj(head(df), indent=2*8)
log_obj <- function(obj, log_cmd="display", indent=0, tab=0, width=50, ...) {

   if (!.checkJupyter() && (log_cmd=="display")) {
	log_cmd <- "print"
   }

   if (log_cmd == "display") {

	display(obj)

   } else if (log_cmd == "print") {

	width_org <- getOption("width")
	options(width = width)

	str_indent <- ""
	if (indent > 0) {
		str_indent <- paste(rep(" ", indent), sep="", collapse="")
	}

	str_tab <- ""
	if (tab > 0) {
		str_tab <- paste(rep("\t", tab), sep="", collapse="")
	}

	if ((nchar(str_indent) > 0) || (nchar(str_tab) > 0)) {
		tcdesc <- ""	
		tc_obj <- textConnection("tcdesc", "w", local=TRUE)
		sink(tc_obj); try(print(obj, ...)); sink(); close(tc_obj)
		if (length(tcdesc) > 1) {
			if (nchar(tcdesc[1]) == 0) {
				# e.g. log_obj(table(df[,1]), tab=2)
				tcdesc <- tcdesc[-1]
			}
		}

		cat(paste0(str_tab, str_indent, tcdesc), sep="\n")
	} else {
		do.call(log_cmd, list(obj) )
	}

	options(width = width_org)

   } else {

	width_org <- getOption("width")
	options(width = width)

	do.call(log_cmd, list(obj) )

	options(width = width_org)

   }

} # log_obj





# log_txt
#
# input:
#   msg:
#   log_cmd_txt: {["display_html"], "cat"]
#
# usage:
# log_txt(sprintf("%s\n\t\t%s\n", "\ttitle", "\t\tline1"), "display_html")
log_txt <- function(msg, log_cmd_txt="display_html", width=50) {

   if (length(msg) == 0) return(NULL)
   if (!.checkJupyter() && (log_cmd_txt == "display_html")) {
	log_cmd_txt <- "cat"
   }

   if (log_cmd_txt == "display_html") {

	if (grepl("\\t", msg)) {
		msg <- gsub("\\t", "&nbsp;&nbsp;", msg)
		msg <- gsub("\\n", "<br>", msg)
	}
	display_html(msg)

   } else {

	width_org <- getOption("width")
	options(width = width)

	do.call(log_cmd_txt, list(msg))

	options(width = width_org)

   } # if

} # log_txt






