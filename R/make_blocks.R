#' @title make blocks from mappings
#'
#' @description
#' \code{make_blocks} Internal function to build blocks
#'
#' @param map The map object (data.frame or data.table)
#' @param rerank Should the ranks of gene order be re-calculated?
#' @param rename.blocks Logical, should the block be re-named according to
#' order of the genome, chromosome and block IDs?
#' @param drop.NAs Logical, if TRUE, only retain complete cases of the map.
#' @param add.metadata Logical, should metadata about blocks be added?
#' @param clean.columns Logical, should un-necessary columns be dropped?
#' @param ties.method passed on to data.table::frank ties.method.
#' @param ... not currently in use
#' @details Nothing yet
#' @return List with blocks and mappings
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
make_blocks <- function(map,
                        rerank = T,
                        drop.NAs = F,
                        rename.blocks = F,
                        add.metadata = F,
                        clean.columns = T,
                        ties.method = "dense",
                        ...){


  if(clean.columns){
    cols2keep <- c("block.id", "orthogroup",
                   "genome1", "genome2",
                   "id1", "id2",
                   "og1", "og2",
                   "chr1", "start1", "end1", "strand1", "order1",
                   "chr2", "start2", "end2", "strand2", "order2",
                   "rank1", "rank2")
    cols2keep <- cols2keep[cols2keep %in% colnames(map)]
    map <- map[,cols2keep, with = F]
  }

  map <- data.table(map)
  if (rename.blocks) {
    map$block.id <- as.numeric(as.factor(with(map, paste(genome1, genome2, block.id))))
  }

  setkey(map, chr1, chr2, start1, start2)

  if (rerank) {
    map[,rank1 := frank(start1,
                        ties.method = ties.method),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = ties.method),
        by = list(genome1, genome2, chr2)]
  }

  if (drop.NAs) {
    map <- map[complete.cases(map),]
  }

  if (!add.metadata) {
    out.blk <- map[,list(chr1 = chr1[1],
                         chr2 = chr2[1],
                         start1 = min(start1),
                         start2 = min(start2),
                         end1 = max(end1),
                         end2 = max(end2),
                         rankstart1 = min(rank1),
                         rankstart2 = min(rank2),
                         rankend1 = max(rank1),
                         rankend2 = max(rank2),
                         n.mapping = length(start1),
                         orient = ifelse(length(start1) <= 1, "+",
                                         ifelse(cor(jitter(start1),
                                                    jitter(start2)) > 0,"+", "-"))),
                   by = list(block.id, genome1, genome2)]
  }else{
    out.blk <- map[,list(chr1 = chr1[1],
                         chr2 = chr2[1],
                         start1 = min(start1),
                         start2 = min(start2),
                         end1 = max(end1),
                         end2 = max(end2),
                         rankstart1 = min(rank1),
                         rankstart2 = min(rank2),
                         rankend1 = max(rank1),
                         rankend2 = max(rank2),
                         n.mapping = length(start1),
                         n.unique.map1 = length(unique(id1)),
                         n.unique.map2 = length(unique(id2)),
                         orient = ifelse(length(start1) <= 1, "+",
                                         ifelse(cor(jitter(start1),
                                                    jitter(start2)) > 0,"+", "-"))),
                   by = list(block.id, genome1, genome2)]
    out.blk$width1 <- with(out.blk, rankend1 - rankstart1) + 1
    out.blk$width2 <- with(out.blk, rankend2 - rankstart2) + 1
    out.blk$prop.map1 <- with(out.blk, n.unique.map1 / width1)
    out.blk$prop.map2 <- with(out.blk, n.unique.map2 / width2)
    out.blk$prop.ratio1 <- with(out.blk, prop.map1 / prop.map2)
    out.blk$prop.ratio2 <- with(out.blk, prop.map2 / prop.map1)
  }


  map <- data.table(map,
                    stringsAsFactors = F)
  blk <- data.table(out.blk,
                    stringsAsFactors = F)

  return(list(block = blk,
              map = map))
}
