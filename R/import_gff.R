#' @title Import a gff annotation
#'
#' @description
#' \code{import_gff} Parse and simplify a gff
#'
#' @param gff.dir file.path, to the subdirectory containing
#' the parsed gff annotation files
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param use character, the string to match in the 3rd gff column
#' @param str2drop character, string in attributes column to drop
#' @param str2parse character, string in attributes column to split on
#' @param whichAttr numeric, specifying the index of attribute column to pull
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#'
#' @details Imports and parses the gff annotation
#'
#' @return a data.table containing relevant information in the gff.
#'
#' @examples
#' \dontrun{
#' gff <- import_gff(
#' gff.dir = dir.locs$gff,
#' genomeIDs = genomeIDs)
#' }
#' @import data.table
#' @export
import_gff <- function(gff.dir,
                       genomeIDs,
                       use = "gene",
                       verbose = T,
                       str2drop = "Name=",
                       str2parse = ";",
                       whichAttr = 2){

  # -- check the params
  ################################################
  ################################################
  ################################################
  if (!is.character(genomeIDs) | length(genomeIDs) == 1)
    stop("genomeIDs must be a character vector of length > 1\n")

  if (!is.character(use) | length(use) != 1)
    stop("use must be a single string to search for in the 3rd gff column\n")

  if (!is.character(str2drop) | length(str2drop) != 1)
    stop("str2drop must be a single string to remove in the 9th gff column\n")

  if (!is.character(str2parse) | length(str2parse) != 1)
    stop("str2parse must be a single string
         to use as the demiliter in the 9th gff column\n")

  if (!is.numeric(whichAttr) | length(str2parse) != 1)
    stop("whichAttr must be a single number,
         specifying which parsed field in the 9th gff column contains the gene identifier\n")
  if (!dir.exists(gff.dir))
    stop("can't find gff.dir ... \n")
  ################################################
  ################################################
  ################################################

  ################################################
  parse_gff <- function(gff.file,
                        str2drop = "Name=",
                        str2parse = ";",
                        use = "gene",
                        whichAttr = 2){
    g <- suppressWarnings(
      fread(gff.file,
            showProgress = F,
            verbose = F, skip = use))
    g <- g[g$V3 == use, c(9, 1, 4, 5, 7)]
    g$V9 <- sapply(g$V9, function(x)
      gsub(str2drop, "",
           strsplit(x, str2parse)[[1]][whichAttr]))
    setnames(g, c("id", "chr", "start", "end", "strand"))
    return(g)
  }
  #######################################################
  #######################################################

  gff.files <- file.path(gff.dir,
                         paste0(genomeIDs, ".gff3"))
  names(gff.files) <- genomeIDs

  gff <- rbindlist(lapply(names(gff.files), function(i){
    if (verbose)
      cat("\tReading",i,"... ")
    tmp <- parse_gff(gff.file = gff.files[[i]],
                     str2drop = str2drop,
                     str2parse = str2parse,
                     use = use,
                     whichAttr = whichAttr)

    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    if (verbose)
      cat("n. genes =", nrow(tmp),"\n")
    return(tmp)
  }))

  setkey(gff, "genome", "id")
  return(gff)
}
