#' @title Simple orthogroup-constrained tandem array inference.
#'
#' @description
#' \code{find_arrayClusters} Simple orthogroup-constrained tandem array inference
#' @param blast the blast dataset to screen for syntenic hits
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
find_arrayClusters = function(blast,
                              clean.radius = 250,
                              clean.mappings = 5,
                              merge.buffer = 10,
                              verbose = T,
                              min.hits1 = 2,
                              min.hits2 = 3,
                              ...){
  if(verbose)
    cat("Dropping intra-genomic hits\n")
  blast <- blast[with(blast, genome1 != genome2),]
  if(verbose)
    cat("Forming blocks from blast\n")
  syn.clean = clean_blocks(blast,
                           radius = clean.radius,
                           n.mappings = clean.mappings,
                           verbose = F)
  if(verbose)
    cat("Merging adjacent blocks\n")
  syn.merge = with(syn.clean,
                   merge_blocks(map = map,
                                blk = block,
                                buffer = merge.buffer,
                                verbose = F))

  if(verbose)
    cat("Counting hits and assigning array ids\n")
  map = syn.merge$map

  bl.og <- blast[,c("id1","id2","og1")]
  setnames(bl.og, "og1","og")
  m <- merge(bl.og, map, by = c("id1","id2"))

  m$unique = with(m, paste(genome1, genome2, block.id, og))
  m[,n.hits1 := length(unique(id1)),
    by = list(unique)]
  m[,n.hits2 := length(unique(id2)),
    by = list(unique)]
  mo <- m[(m$n.hits1 >= min.hits1 & m$n.hits2 >= min.hits2) |
            (m$n.hits1 >= min.hits2 & m$n.hits2 >= min.hits1) ,
          c("unique","id1","id2")]
  n.genes = length(unique(c(mo$id1, mo$id2)))
  mo$id1<-NULL
  mo$id2<-NULL
  mo <- mo[!duplicated(mo),]
  mo$tandemarray.id <- paste(mo$unique, 1:nrow(mo))

  out <- merge(m, mo, by = "unique", all = T)
  out$og1 <- out$og
  out$og2 <- out$og
  if(verbose)
    cat("Found", nrow(mo), "arrays that contain", n.genes, "unique genes\nDone!\n")
  return(out)
}
