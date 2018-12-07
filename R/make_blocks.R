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
                        rename.blocks = T){
  map <- data.table(map)
  if(rename.blocks){
    map$block.id <- as.numeric(as.factor(with(map, paste(genome1, genome2, chr1, chr2, block.id))))
  }

  setkey(map, chr1, chr2, start1, start2)
  if(rerank){
    map[,rank1 := frank(start1,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr2)]
  }
  if(drop.NAs){
    map <- map[complete.cases(map),]
  }

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
                       n.mapping = length(score),
                       orient = ifelse(cor(start1, start2) > 0,"+","-")),
                 by = list(block.id, genome1, genome2)]

  map <- data.table(map,
                    stringsAsFactors = F)
  blk <- data.table(out.blk,
                    stringsAsFactors = F)

  return(list(block = blk,
              map = map))
}
