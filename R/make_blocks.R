#' @title make blocks from mappings
#'
#' @description
#' \code{make_blocks} Internal function to build blocks
#'
#' @param map The map object (data.frame or data.table)
#' @details Nothing yet
#' @return List with blocks and mappings
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
make_blocks<-function(map){
  map = data.table(map)

  map$rank1 = frank(map, chr1, start1, ties.method = "dense")
  map$rank2 = frank(map, chr2, start2, ties.method = "dense")
  setkey(map, chr1, chr2, start1, start2)
  out.blk <- map[,list(genome1 = genome1[1],
                       genome2 = genome2[1],
                       chr1 = chr1[1],
                       chr2 = chr2[1],
                       start1 = min(start1),
                       start2 = min(start2),
                       end1 = max(end1),
                       end2 = max(end2),
                       rankstart1 = min(rank1),
                       rankstart2 = min(rank2),
                       rankend1 = max(rank1),
                       rankend2 = max(rank2),
                       medianscore = median(score),
                       n.mapping = length(score)),
                 by=list(block.id)]

  out.blk$density = with(out.blk,
                         n.mapping/((abs(rankend1-rankstart1)+abs(rankend2-rankstart2))/2))

  orient = sapply(split(map, map$block.id), function(x) cor(x$start1, x$start2))
  orient = data.table(block.id = names(orient), orient = ifelse(orient>0,"+","-"))
  orient$block.id<-as.numeric(orient$block.id)
  setkey(out.blk, block.id)
  setkey(orient, block.id)
  out.blk = merge(out.blk, orient)
  out.blk$first1 = with(out.blk, ifelse(orient == "+", start1, end1))
  out.blk$first2 = with(out.blk, ifelse(orient == "+", start2, end2))
  out.blk$last1 = with(out.blk, ifelse(orient == "+", end1, start1))
  out.blk$last2 = with(out.blk, ifelse(orient == "+", end2, start2))
  out.blk$firstrank1 = with(out.blk, ifelse(orient == "+", rankstart1, rankend1))
  out.blk$firstrank2 = with(out.blk, ifelse(orient == "+", rankstart2, rankend2))
  out.blk$lastrank1 = with(out.blk, ifelse(orient == "+", rankend1, rankstart1))
  out.blk$lastrank2 = with(out.blk, ifelse(orient == "+", rankend2, rankstart2))
  map = data.frame(map[order(map$genome1, map$chr1, map$start1),], stringsAsFactors = F)
  blk = data.frame(out.blk, stringsAsFactors = F)
  return(list(block = blk, map = map))
}
