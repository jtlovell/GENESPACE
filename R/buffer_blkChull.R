#' @title Create buffer around blocks
#'
#' @description
#' \code{buffer_blkChull} Generate a convex hull and buffer it
#' by gene rank order. Used primarily for culling blast hits.
#'
#' @param rank1 Numeric, ranks (genome1) of hits in blocks
#' @param rank2 Numeric, ranks (genome2) of hits in blocks
#' @param block.id Character, specifying block identity for
#' gene rank hits
#' @param rank.buffer Numeric, length 1, specifying the number
#' of ranks to extend block convex hull.
#' @param plotit Logical, should a plot be made?
#' @param plot.title Character, length 1, specifying plot title.
#' @param ... Additional arguments passed to plot, controlling
#' plotting of rank1 v. rank2.
#' @details None yet

#' @return Output from rgeos::gBuffer, a SpatialPolygons object
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#'
#' @import data.table
#' @importFrom sp CRS SpatialPointsDataFrame over Polygon Polygons SpatialPolygons
#' @importFrom rgeos gBuffer
#' @importFrom grDevices chull
#' @export
buffer_blkChull <- function(rank1,
                            rank2,
                            block.id,
                            rank.buffer = 500,
                            plotit = F,
                            plot.title = NULL){

  m <- data.frame(block.id = as.character(block.id),
                  rank1 = as.numeric(rank1),
                  rank2 = as.numeric(rank2),
                  stringsAsFactors = F)
  xy <- data.matrix(data.frame(rank1 = as.numeric(rank1),
                   rank2 = as.numeric(rank2),
                   stringsAsFactors = F))

  crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

  chulls <- lapply(split(m, m$block.id), function(x){
    block.id <- x$block.id[1]
    ch <- chull(jitter(data.matrix(x[,c("rank1","rank2")])))
    ch <- c(ch,ch[1])
    hull.df <- x[ch, ]
    hull.xy <- data.frame(hull.df[,c("rank1","rank2")])
    poly <- Polygon(hull.xy, hole = F)
    return(Polygons(list(poly), block.id))
  })

  spoly <- SpatialPolygons(chulls,
                           proj4string = crs)
  gpoly <- gBuffer(spoly,
                   width = rank.buffer)

  if(plotit){
    plot(rank1, rank2,
         xlab = "rank, genome1",
         ylab = "rank, genome2",
         main = plot.title)
    plot(spoly,
         add = T,
         col = rgb(1,0,0,.5))
    plot(gpoly,
         add = T,
         col = rgb(0,0,1,.2),
         lty = 2)
  }
  return(gpoly)
}
