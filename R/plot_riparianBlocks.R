#' @title plot_riparianBlocks
#'
#' @description
#' \code{plot_riparianBlocks} plot_riparianBlocks
#'
#' @param map data.table, format similar to blast
#' @param rip.out list, output from plot_riparian
#' @param genes2plot character of gene IDs
#' @param blocks2plot character of block IDs
#' @param blk.col character, to be coerced to a color of the polygons that
#' represent syntenic blocks in riparian plot
#' @param blk.border logical (or character coerced to color), of the block
#' borders in riparian plot
#' @param ... not currently in use
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
plot_riparianBlocks <- function(map,
                                rip.out,
                                genes2plot = NULL,
                                blocks2plot = NULL,
                                blk.border,
                                radius = 5,
                                n.mapping = 5,
                                blk.col){

  ####################################################################
  # 1. get all unique genes in a block
  m <- data.table(map)
  if(is.null(genes2plot) == is.null(blocks2plot))
    stop("must specify a region, defined by EITHER genes OR blocks\n")
  if(is.null(blocks2plot)){
    b <- subset(m, id1 %in% genes2plot | id2 %in% genes2plot)
  }else{
    b <- subset(m, block.id %in% blocks2plot)
  }

  genes.in.blk <- unique(unlist(subset(
    m,
    id1 %in% unique(b$id1, b$id2) |
      id2 %in% unique(b$id1, b$id2))[,c("id1","id2")]))
  mb <- subset(m, id1 %in% genes.in.blk & id2 %in% genes.in.blk)

  ####################################################################
  # 2. re-make blocks
  mb[,rank1 := frank(start1, ties.method = "dense"),
     by = c("genome1","genome2","chr1","chr2")]
  mb[,rank2 := frank(start2, ties.method = "dense"),
     by = c("genome1","genome2","chr1","chr2")]
  mb[,block.id := dbscan(frNN(cbind(rank1, rank2),
                              eps = radius),
                         minPts = n.mapping)$cluster,
     by = c("genome1", "genome2","chr1","chr2")]

  ####################################################################
  # 3. re-project the blocks
  mb2 <- data.table(rip.out$map)
  mb2[, block.id := NULL]
  mb2 <- merge(mb2, mb[,c("id1","id2","block.id")], by = c("id1", "id2"))
  blk <- mb2[,list(start1 = pos1[1],
                   start2 = pos2[1],
                   end1 = pos1[length(pos1)],
                   end2 = pos2[length(pos2)]),
             by = list(block.id,
                       y1, y2,
                       genome1, genome2,
                       chr1, chr2)]

  ####################################################################
  # 4. draw the polygon
  draw_blkPolygon(
    blk = blk,
    chr.buffer = rip.out$rip.param$chr.buff,
    blk.border = blk.border,
    blk.col = blk.col,
    simplify.poly = rip.out$rip.param$simplify.poly,
    points.per.curve = rip.out$rip.param$points.per.curve)

  return(list(blk = blk, map = mb))
}
