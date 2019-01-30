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
rerank_fromIDs <- function(id1,
                           id2,
                           gff,
                           ties.method = "dense",
                           ...){
  u.id <- unique(c(id1, id2))
  if ("rank" %in% colnames(gff))
    gff$rank <- NULL

  gff[, rank := frank(start, ties.method = ties.method),
      by = list(genome, chr)]
  gff <- gff[gff$id %in% u.id, ]

  gff1 <- data.table(gff)
  gff2 <- data.table(gff)

  setnames(gff1, paste0(colnames(gff1), "1"))
  setnames(gff2, paste0(colnames(gff2), "2"))
  setkey(gff1, id1)
  setkey(gff2, id2)

  id.dt <- data.table(id1 = id1,
                      id2 = id2,
                      stringsAsFactors = F)

  setkey(id.dt, id2)
  m1 <- merge(gff2, id.dt)
  setkey(m1, id1)
  out <- merge(gff1, m1)

  out[, unique.genome := paste(genome1, genome2)]
  out[, unique.chr := paste(genome1, genome2,
                            chr1, chr2)]

  return(out)
}
