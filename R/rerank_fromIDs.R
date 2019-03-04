#' @title Cull blast to syntenic regions
#'
#' @description
#' \code{cull_syntenicBlast} An internal function, designed to be called
#' by extend_blocks and find_syntenicOrthogs.
#'
#' @param id1 Character vector of ids
#' @param id2 Character vector of ids
#' @param gff gff data.table
#' @param ties.method passed on to data.table::frank
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
rerank_fromIDs <- function(map,
                           gff,
                           ...){

  m.gff1 <- map[,c("genome1","id1", "chr1", "start1","end1")]
  m.gff2 <- map[,c("genome2","id2", "chr2", "start2","end2")]
  gff <- gff[,c("genome","id","chr","start","end")]
  setnames(m.gff1, colnames(gff))
  setnames(m.gff2, colnames(gff))
  g <- rbind(gff, m.gff1, m.gff2)
  g <- g[!duplicated(g),]

  g$unique.chr <- with(g, paste(genome, chr))
  spl <- split(g, "unique.chr")
  g <- rbindlist(lapply(spl, function(x){
    x$rank <- frank(x, start, end, id, ties.method = "dense")
    return(x)
  }))
  g$unique.chr <- NULL

  g1 <- data.table(g)
  g2 <- data.table(g)

  setnames(g1, paste0(colnames(g1), "1"))
  setnames(g2, paste0(colnames(g2), "2"))
  setkeyv(g1, cols = colnames(g1)[-6])
  setkeyv(g2, cols = colnames(g2)[-6])

  setkeyv(map, cols = colnames(g2)[-6])
  if("rank1" %in% colnames(map))
    map$rank1 <- NULL
  if("rank2" %in% colnames(map))
    map$rank2 <- NULL
  if("rank" %in% colnames(map))
    map$rank <- NULL
  if("n" %in% colnames(map))
    map$n <- NULL

  m1 <- merge(g2, map)
  setkeyv(m1, cols = colnames(g1)[-6])
  out <- merge(g1, m1)

  return(out)
}
