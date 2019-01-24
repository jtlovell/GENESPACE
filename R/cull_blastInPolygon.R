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

  gff1.idrank <- gff[gff$id %in% map$id1, c("id","rank")]
  gff2.idrank <- gff[gff$id %in% map$id2, c("id","rank")]

  setnames(gff1.idrank, c("id1","rank1"))
  setnames(gff2.idrank, c("id2","rank2"))
  setkey(gff1.idrank, id1)
  setkey(gff2.idrank, id2)
  map.id.blk <- map[, c("id1", "id2", "block.id")]

  setkey(map.id.blk, id2)
  merge.gff2.map.id.blk <- merge(gff2.idrank, map.id.blk)
  setkey(merge.gff2.map.id.blk, id1)
  merge.gff.map.id.blk <- merge(gff1.idrank,
                                merge.gff2.map.id.blk)
  totest<<-merge.gff.map.id.blk
  chulls <- with(merge.gff.map.id.blk,
              buffer_blkChull(rank1 = rank1,
                            rank2 = rank2,
                            block.id = block.id,
                            plotit = plotit,
                            rank.buffer = rank.buffer,
                            plot.title = paste(map$genome1[1], map$genome2[1],
                                               map$chr1[1], map$chr2[1])))

  blast.ids <- blast[,c("id1","id2")]
  setkey(blast.ids, id2)
  merge.blast.gff <- merge(gff2.idrank, blast.ids)
  setkey(merge.blast.gff, id1)
  merge.blast.gff <- merge(gff1.idrank, merge.blast.gff)
  merge.blast.gff <- data.frame(merge.blast.gff, stringsAsFactors = F)

  xy.blast <- data.matrix(data.frame(merge.blast.gff[, c("rank1", "rank2")],
                   stringsAsFactors = F))

  crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  spdf <- SpatialPointsDataFrame(coords = xy.blast,
                                 data = merge.blast.gff,
                                 proj4string = crs)

  test <- over(spdf,
               chulls,
               fn = NULL)
  merge.blast.gff.cull <- merge.blast.gff[!is.na(test),]
  if(plotit)
    points(merge.blast.gff.cull$rank1,
           merge.blast.gff.cull$rank2,
           col = "green3",
           cex = .2)
  return(merge.blast.gff.cull)
}

