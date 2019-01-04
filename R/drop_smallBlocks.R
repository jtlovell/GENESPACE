#' @title Drop small blocks
#' @description
#'  \code{drop_smallBlocks} Function to limit blocks to a minimum size
#' @param blk a data.table or data.frame containing the block information
#' @param map a data.table or data.frame containing the map information
#' @param min.block.size numeric, what is the smallest block to retain?
#' @param ... Not currently in use
#' @details Nothing yet
#' @return new list of block and map data.tables
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
drop_smallBlocks <- function(blk,
                             map,
                             min.block.size = 5){
  blk <- data.table(blk)
  map <- data.table(map)
  wh2keep <- blk$block.id[blk$n.mapping >= min.block.size]
  map <- map[map$block.id %in% wh2keep, ]
  out <- make_blocks(map,
                     rerank = T,
                     rename.blocks = T)
  return(out)
}
