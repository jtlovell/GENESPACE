#' @title Import gffs
#'
#' @description
#' \code{import_gff} Import gff annotation files as a data.table, culling
#' to gene entries and parsing the attribute column.
#'
#' @param gff.dir path to directory storing the gff annotations
#' @param genomeIDs character vector of genome identifiers
#' @param str2drop character, string in attribute column of gff file to be dropped
#' @param str2parse character, string in attribute column of gff file to use as the separator
#' @param whichAttr numeric, which attribute should be returned in the
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return culled gff data.table, concatenated across genomeIDs
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
import_gff <- function(gff.dir,
                       genomeIDs,
                       verbose = T,
                       str2drop = "Name=",
                       str2parse = ";",
                       whichAttr = 2){
  #######################################################
  # -- Error checking
  stop_withMessage(dir.exists(gff.dir) | missing(gff.dir),
                   paste(gff.dir, "does not exist"))
  stop_withMessage(c(is.character(genomeIDs),
                     length(genomeIDs) > 1,
                     missing(genomeIDs)),
                   "genomeIDs must be a character vector of length > 1")
  stop_withMessage(is.numeric(whichAttr),
                   "whichAttr must be an integer")
  stop_withMessage(c(nchar(str2parse) == 1,
                     is.character(str2parse)),
                   "str2parse must be a single character")
  stop_withMessage(c(is.character(str2drop),
                     length(str2drop) == 1),
                   "str2drop must be a character vector of length 1")
  stop_withMessage(c(is.character(str2drop),
                     length(str2drop) == 1),
                   "str2drop must be a character vector of length 1")
  stop_withMessage(is.logical(verbose),
                   "verbose must be logical")
  whichAttr <- as.integer(whichAttr)
  #######################################################
  gff.files <- file.path(gff.dir,
                         paste0(genomeIDs, ".gff3"))

  stop_withMessage(file.exists(gff.files),
                   paste("some gff files do not exist in", gff.fir))

  names(gff.files) <- genomeIDs
  #######################################################
  gff <- rbindlist(lapply(names(gff.files), function(i){
    if (verbose)
      cat("\tReading",i,"... ")
    tmp <- parse_gff(gff.file = gff.files[[i]],
                     str2drop = str2drop,
                     str2parse = str2parse,
                     whichAttr = whichAttr)
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    if (verbose)
      cat("n. genes =", nrow(tmp),"\n")
    return(tmp)
  }))
  #######################################################

  #######################################################
  setkey(gff, "genome", "id")
  return(gff)
}

