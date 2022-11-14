#' @title Tabulate copy number variation
#'
#' @description
#' \code{tabulate_cnv} Count ortholog copy number by gene and genome.
#'
#' @param gsParam GENESPACE parameter list.
#' @details Txxx
#'
#' @return xxx
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#
#' @import data.table
#' @export
tabulate_cnv <- function(gsParam){

  nGeneGlobOG <- globHOG <- genome <- nGeneOg <- og <- NULL
  # read in the bed file
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  # 1. Add count of genes within genome for each orthogroup
  if(gsParam$params$useHOGs || is.na(gsParam$params$useHOGs)){
    bed[,nGeneGlobOG := .N, by = c("genome", "globHOG")]
    ggw <- dcast(bed, globHOG ~ genome, value.var = "nGeneGlobOG", fun.aggregate = length)
  }else{
    bed[,nGeneGlobOG := .N, by = c("genome", "globHOG")]
    ggw <- dcast(bed, globOG ~ genome, value.var = "nGeneGlobOG", fun.aggregate = length)
  }

  bed[,nGeneOg := .N, by = c("genome", "og")]
  ogw <- dcast(bed, og ~ genome, value.var = "nGeneOg", fun.aggregate = length)

  return(list(global = ggw, syntenic = ogw, bed = bed))
}



