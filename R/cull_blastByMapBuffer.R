#' @title cull blast by map hits
#'
#' @description
#' \code{cull_blastByMapBuffer} reduce blast hits to those in proximity to
#' hits in map object.
#'
#' @param map The map object (data.frame or data.table)
#' @param blast The blast object, which is like a map object, but
#' does not necessarily contain a 'block.id' column.
#' @param gff The gff data frame, containing concatenated annotations
#' across all genomes.
#' @param rerank.gff Should the ranks of gene order be re-calculated?
#' @param buffer Distance in gene rank order from hits in blocks to retain
#' blast hits
#' @param n.cores Number of parallel processes to run
#' @param plotit Logical, should plots be generated?
#' @details Nothing yet
#' @return A blast object culled to the regions.
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
cull_blastByMapBuffer = function(gff,
                                 blast,
                                 map,
                                 rerank.gff = T,
                                 buffer = 1000,
                                 n.cores = 1,
                                 plotit = T,
                                 verbose = T){

  #######################################################
  #######################################################
  cull_blastbyMap = function(map,
                             gff,
                             blast,
                             buffer,
                             plotit = F){

    x = map
    tg1 = gff[gff$id %in% x$id1,c("id","rank")]
    tg2 = gff[gff$id %in% x$id2,c("id","rank")]
    setnames(tg1, c("id1","rank1"))
    setnames(tg2, c("id2","rank2"))
    setkey(tg1,id1)
    setkey(tg2, id2)
    tm0 = x[,c("id1","id2","block.id")]
    setkey(tm0, id2)
    m = merge(tg2, tm0)
    setkey(m, id1)
    m = merge(tg1, m)

    out = buffer_chulls(x = m$rank1,
                        y = m$rank2,
                        id = m$block.id,
                        plotit = plotit,
                        buffer = buffer,
                        plot.title = paste0(x$genome1[1]," (",x$chr1[1],") vs. ",
                                            x$genome2[1]," (",x$chr2[1],")"))

    x = blast
    tg1 = gff[gff$id %in% x$id1,c("id","rank")]
    tg2 = gff[gff$id %in% x$id2,c("id","rank")]
    setnames(tg1, c("id1","rank1"))
    setnames(tg2, c("id2","rank2"))
    setkey(tg1,id1)
    setkey(tg2, id2)
    tm0 = x[,c("id1","id2")]
    setkey(tm0, id2)
    m = merge(tg2, tm0)
    setkey(m, id1)
    m = merge(tg1, m)
    m = data.frame(m, stringsAsFactors = F)
    xy = data.frame(m[,c("rank1","rank2")], stringsAsFactors = F)
    crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    spdf <- SpatialPointsDataFrame(coords = xy,
                                   data = m,
                                   proj4string = crs)

    test = over(spdf, out, fn = NULL)
    bl.cull = m[!is.na(test),]
    if(plotit)
      points(bl.cull$rank1,
             bl.cull$rank2,
             col = "green3", cex = .2)
    return(bl.cull)
  }
  #######################################################
  #######################################################
  buffer_chulls <- function(x,y,id, buffer = 500, plotit = F, plot.title = NULL){
    m = data.frame(id, x, y, stringsAsFactors = F)
    xy = data.frame(x,y, stringsAsFactors = F)

    crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    spdf <- SpatialPointsDataFrame(coords = xy,
                                   data = m,
                                   proj4string = crs)

    chulls = lapply(split(m, m$id), function(x){
      id = x$id[1]
      ch = chull(x[,c("x","y")])
      ch = c(ch,ch[1])
      hull.df = x[ch, ]
      xy = data.frame(hull.df[,c("x","y")])
      poly = Polygon(xy, hole=as.logical(NA))
      return(Polygons(list(poly), id))
    })
    crs <- CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    spoly = SpatialPolygons(chulls, proj4string = crs)
    gpoly = gBuffer(spoly, width = buffer)

    if(plotit){
      plot(x,y, xlab = "rank, genome1", ylab = "rank, genome2",
           main = plot.title)
      plot(spoly, add = T, col = rgb(1,0,0,.5))
      plot(gpoly, add = T, col = rgb(0,0,1,.2), lty = 2)
    }
    return(gpoly)
  }
  #######################################################
  #######################################################

  #######################################################
  #######################################################

  #######################################################
  if(verbose)
    cat("Parsing gff, blast and map objects\n")
  if(rerank.gff){
    gff[,rank := frank(start, ties.method = "dense"),
        by = list(genome, chr)]
  }

  map$unique.genome = with(map,
                           paste(genome1, genome2))
  blast$unique.genome = with(blast,
                             paste(genome1, genome2))
  map$unique = with(map,
                    paste(genome1, genome2,
                          chr1, chr2))
  blast$unique = with(blast,
                      paste(genome1, genome2,
                            chr1, chr2))
  maps = split(map, "unique.genome")
  blasts = split(blast, "unique.genome")

  if(verbose)
    cat("Culling blasts by pairwise genome comparisons ...\n\t")

  out.all = rbindlist(lapply(names(maps), function(j){
    map.x = maps[[j]]
    blast.x = blasts[[j]]
    if(verbose)
      cat(map.x$genome1[1],"-->", map.x$genome2[1],
          paste0("(", nrow(map.x),"/",nrow(blast.x),") ... "))

    spl.map = split(map.x, "unique")
    spl.blast = split(blast.x, "unique")
    #######################################################

    #######################################################
    out <- rbindlist(mclapply(names(spl.map), mc.cores = n.cores, function(i){
      blast.x = spl.blast[[i]]
      map.x = spl.map[[i]]
      cull.blast = cull_blastbyMap(map = map.x,
                                   blast = blast.x,
                                   gff = gff,
                                   buffer = buffer,
                                   plotit = plotit)
      return(cull.blast)
    }))
    #######################################################

    #######################################################
    setkey(blast, id1, id2)
    out.id = out[,c("id1","id2")]
    setkey(out.id, id1, id2)
    out2 = merge(blast,out.id)
    #######################################################

    #######################################################
    out2$orthogroup = out2$og1
    if(verbose)
      cat("culled to", nrow(out2),"\n\t")
    return(out2)
  }))

  if(verbose)
    cat("Done!\n")
  return(out.all)
}


