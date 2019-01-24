#' @title Cull blast hits by a SpatialPolygons object
#'
#' @description
#' \code{cull_blastInPolygon} Cull blast hits by a SpatialPolygons object
#'
#' @param map The map data.frame or data.table
#' @param blast Some blast data that needs to be culled.
#' Has the same format as the map data, but does not need to have
#' a 'block.id' column.
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param rank.buffer The size of the buffer (in gene order ranks)
#' to extend the chull around each block.
#' @param plotit Should buffer_blkChull plots be made?
#' @param ... Plotting arguments passed on to buffer_blkChull and
#' then to plot.
#' @details None yet

#' @return A blast object, culled by the gBuffer of the map.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom sp CRS SpatialPointsDataFrame over Polygon Polygons SpatialPolygons
#' @importFrom rgeos gBuffer
#' @importFrom grDevices chull
#' @export
cull_blastInPolygon <- function(map,
                                gff,
                                blast,
                                rank.buffer,
                                plotit = F){

  #######################################################
  #######################################################

  #######################################################
  #######################################################

  tg1 <- gff[gff$id %in% blast$id1, c("id","rank")]
  tg2 <- gff[gff$id %in% blast$id2, c("id","rank")]

  setnames(tg1, c("id1","rank1"))
  setnames(tg2, c("id2","rank2"))
  setkey(tg1, id1)
  setkey(tg2, id2)
  tm0 <- map[, c("id1", "id2", "block.id")]

  setkey(tm0, id2)
  m <- merge(tg2, tm0)
  setkey(m, id1)
  m <- merge(tg1, m)

  out <- with(m,
              buffer_blkChull(rank1 = rank1,
                            rank2 = rank2,
                            block.id = block.id,
                            plotit = plotit,
                            rank.buffer = rank.buffer,
                            plot.title = paste(map$genome1[1], map$genome2[1],
                                               map$chr1[1], map$chr2[1])))

  tg1 <- gff[gff$id %in% blast$id1, c("id", "rank")]
  tg2 <- gff[gff$id %in% blast$id2, c("id", "rank")]
  setnames(tg1, c("id1", "rank1"))
  setnames(tg2, c("id2", "rank2"))
  setkey(tg1, id1)
  setkey(tg2, id2)
  tm0 <- blast[,c("id1","id2")]
  setkey(tm0, id2)
  m <- merge(tg2, tm0)
  setkey(m, id1)
  m <- merge(tg1, m)
  m <- data.frame(m, stringsAsFactors = F)

  xy <- data.matrix(data.frame(m[, c("rank1", "rank2")],
                   stringsAsFactors = F))

  crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  spdf <- SpatialPointsDataFrame(coords = xy,
                                 data = m,
                                 proj4string = crs)

  test <- over(spdf,
               out,
               fn = NULL)
  bl.cull <- m[!is.na(test),]
  if(plotit)
    points(bl.cull$rank1,
           bl.cull$rank2,
           col = "green3",
           cex = .2)
  return(bl.cull)
}

