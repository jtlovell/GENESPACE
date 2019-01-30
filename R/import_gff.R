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
#' @importFrom compiler cmpfun
#' @export
import_gff <- function(gff.dir,
                       genomeIDs,
                       verbose = T,
                       str2drop = "Name=",
                       str2parse = ";",
                       whichAttr = 2){
  #######################################################
  #######################################################
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
  #######################################################
  #######################################################
  parse_gff <- cmpfun(parse_gff)
  ########################################################
  ########################################################

  #######################################################
  if (verbose)
    cat("Importing gff3 annotation files\n")
  gff.files <- file.path(gff.dir,
                         paste0(genomeIDs, ".gff3"))
  names(gff.files) <- genomeIDs
  #######################################################

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
  if (verbose)
    cat("\tDone!\n")
  return(gff)
}

