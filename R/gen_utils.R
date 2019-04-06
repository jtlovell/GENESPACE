#' @title GENESPACE utility functions
#' @description
#' \code{gen_utils} Some general utility functions
#' @name gen_utils
#'
#' @param expr expression to be evaluated
#' @param error warnings / error to return
#' @param gff.file path to gff file
#' @param str2drop string to gsub out
#' @param str2parse single character to split on
#' @param whichAttr numeric of length 1, indicating which element to use
#'
#' @param verbose logical, should updates be printed?
#' @param ... not currently in use
#'
#' @note \code{gen_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{gen_utils} returns its own arguments.
#'
#' @title stop_withMessage
#' @description
#' \code{stop_withMessage} stop_withMessage
#' @rdname gen_utils
#' @export
stop_withMessage <- function(expr, error){
  if (!all(expr))
    stop(paste0(error,"/n"), call. = FALSE)
}

#' @title parse_gff
#' @description
#' \code{parse_gff} parse_gff
#' @rdname gen_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @importFrom Biostrings width readAAStringSet
#' @export
parse_gff <- function(gff.file,
                      str2drop = "Name=",
                      str2parse = ";",
                      whichAttr = 2){
  g <- suppressWarnings(
    fread(gff.file,
          showProgress = F,
          verbose = F))
  g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
  g$V9 <- sapply(g$V9, function(x)
    gsub(str2drop, "",
         strsplit(x, str2parse)[[1]][whichAttr]))
  setnames(g, c("id", "chr", "start", "end", "strand"))
  return(g)
}
