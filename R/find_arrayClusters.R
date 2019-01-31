#' @title Simple orthogroup-constrained tandem array inference.
#'
#' @description
#' \code{find_arrayClusters} Simple orthogroup-constrained tandem array inference
#' @param map the blast dataset to screen for syntenic hits
#' @param rerank Logical, should ranks be remade before each step?
#' @param clean.radius Passed on to clean_blocks
#' @param clean.mappings Passed on to clean_blocks
#' @param merge.buffer Passed on to merge_blocks
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
cluster_tandemArrays <- function(map,
                                 rerank = T,
                                 clean.radius = 250,
                                 clean.mappings = 5,
                                 merge.buffer = 10,
                                 verbose = T,
                                 min.hits1 = 2,
                                 min.hits2 = 3,
                                 ...){
  if (verbose)
    cat("Dropping intra-genomic hits\n")
  map <- map[with(map, genome1 != genome2),]
  if (rerank){
    map[,rank1 := frank(start1,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method =  "dense"),
        by = list(genome1, genome2, chr2)]
  }

  if (rerank) {
    if (verbose)
      cat("Using the gff gene ranks as positions\n")
    rr <- with(map,
               rerank_fromIDs(id1 = id1,
                              id2 = id2,
                              gff = gff))
    map$rank1 <- NULL
    map$rank2 <- NULL
    map <- merge(rr[,c("id1", "id2",
                       "rank1", "rank2")],
                 map,
                 by = c("id1", "id2"))

  }
  cols2keep <- c(colnames(map),
                 "n.hits1",
                 "n.hits2",
                 "tandemarray.id")
  map$unique <- with(map,
                     paste(genome1, genome2,
                           block.id, og1))
  map[,n.hits1 := length(unique(id1)),
      by = list(unique)]
  map[,n.hits2 := length(unique(id2)),
      by = list(unique)]
  mo <- map
  mo$tandemarray.id <- with(mo,
                            paste(og1,
                                  unique.genome,
                                  block.id))
  mo <- map[with(map,
                 (n.hits1 >= min.hits1 &
                    n.hits2 >= min.hits2) |
                   (n.hits1 >= min.hits2 &
                      n.hits2 >= min.hits1)),
            c("unique", "id1", "id2","tandemarray.id")]
  n.genes <- nrow(mo)
  mo$id1 <- NULL
  mo$id2 <- NULL
  mo <- mo[!duplicated(mo),]


  out <- merge(map,
               mo[, c("unique", "tandemarray.id")],
               by = "unique",
               all = T)

  if (verbose)
    cat("Found",
        nrow(mo),
        "arrays that contain",
        n.genes,
        "blast hits\nDone!\n")

  return(out[, cols2keep, with = F])
}
