#' @title Simple orthogroup-constrained tandem array inference.
#'
#' @description
#' \code{find_arrayClusters} Simple orthogroup-constrained tandem array inference
#' @param map the blast dataset to screen for syntenic hits
#' @param min.hits1 The minimum amount of hits in either genome
#' to be counted as an array
#' @param min.hits2 The minimum amount of hits in the over-represented
#' genome to be counted as an array
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
find_tandemArrays <- function(map,
                              min.hits1 = 1,
                              min.hits2 = 3,
                              verbose = T,
                              ...){
  #######################################################
  if (verbose)
    cat("Clustering blast hits within orthogroups and blocks ... \n")
  map$array.id <- with(map, paste0(block.id, "|", og1))

  map$unique.genome <- with(map, paste(genome1, genome2))
  spl <- split(map, "unique.genome")

  out <- rbindlist(lapply(spl, function(x){
    if (verbose)
      cat(paste0("\t",x$genome1[1]),
          "<-->", x$genome2[1], "... ")

    x <- data.table(x)
    x[ ,n.hits1 := length(unique(id1)),
      by = list(array.id)]

    x[ ,n.hits2 := length(unique(id2)),
      by = list(array.id)]

    x$array.id[with(x,
                    (n.hits1 <= min.hits1 &
                       n.hits2 <= min.hits2) |
                      (n.hits1 <= min.hits2 &
                         n.hits2 <= min.hits1))] <- NA
    if (verbose)
      cat("found",
          length(unique(x$og1[!is.na(x$array.id)])),
          "tandem arrays across",
          sum(!is.na(x$array.id)),
          "blast hits\n")
    return(x)
  }))

  if (verbose)
    cat("Done! Found",length(unique(out$array.idss)),
        "tandem arrays across",
        sum(!is.na(out$array.id)),"blast hits\n")
  #######################################################
  return(out)
}
