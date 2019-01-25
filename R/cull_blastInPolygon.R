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
#' @importFrom sp CRS SpatialPointsDataFrame over Polygon Polygons SpatialPolygons coordinates proj4string
#' @importFrom raster extent buffer
#' @export
buffer_blast <- function(map,
                         gff,
                         blast,
                         rank.buffer = 200,
                         round2nearest = rank.buffer/3,
                         n.cores = 1,
                         verbose = T,
                         plotit = F){
  #######################################################
  #######################################################
  mround <- function(x,base){
    base * round(x / base)
  }
  #######################################################
  #######################################################
  cull_blastBuffer <- function(map.chr,
                               blast.chr,
                               rank.buffer,
                               round2nearest,
                               plotit = F,
                               crs = CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")){
    # -- make a culled map dataset by nearby points
    map.pts <- with(map.chr,
                    data.frame(
                      rank1 = mround(rank1, round2nearest),
                      rank2 = mround(rank2, round2nearest),
                      stringsAsFactors = F))
    map.pts <- map.pts[!duplicated(map.pts),]

    coordinates(map.pts) <- ~ rank1 + rank2
    proj4string(map.pts) <- crs

    # -- make the buffer polygon
    buf <- buffer(map.pts, rank.buffer)

    # -- reformat blast
    blast.pts <- with(blast.chr,
                      data.frame(
                        rank1 = rank1,
                        rank2 = rank2,
                        stringsAsFactors = F))
    coordinates(blast.pts) <- ~ rank1 + rank2
    proj4string(blast.pts) <- crs

    # -- check which hits are in the buffer
    is.in.buf <- !is.na(over(blast.pts,
                             buf,
                             fn = NULL))

    # -- keep only those inside the buffer
    blast.chr <- blast.chr[is.in.buf,]

    # -- if desired, make the plot
    if (plotit) {
      ylims <- extent(buf)[3:4]
      xlims <- extent(buf)[1:2]
      plot(jitter(blast.pts$rank1),
           jitter(blast.pts$rank2),
           pch = ".",
           xlab = paste0("gene rank (",
                         blast.chr$genome1[1], ", ",
                         blast.chr$chr1[1], ")"),
           ylab = paste0("gene rank (",
                         blast.chr$genome2[1], ", ",
                         blast.chr$chr2[1], ")"),
           ylim = ylims,
           xlim = xlims)

      plot(buf, add = T,
           col = rgb(0,0,1,.2),
           lty = 2)

      points(map.chr$rank1,
             map.chr$rank2,
             col = "darkorange",
             cex = .2)
    }

    return(blast.chr)
  }
  #######################################################
  #######################################################

  # -- Re-calculate gff ranks
  gff[, rank := as.numeric(frank(start, ties.method = "dense")),
      by = list(genome, chr)]
  tmp.gff <- data.table(gff[gff$id %in% unique(c(blast$id1, blast$id2)),])
  tmp.map <- data.table(map[,c("genome1", "genome2",
                             "chr1", "chr2",
                             "id1", "id2", "block.id")])

  # -- Double gff to be merged with map / blast
  gff1 <- data.table(tmp.gff[, c("id","rank")])
  gff2 <- data.table(tmp.gff[, c("id","rank")])
  setnames(gff1, c("id1","rank1"))
  setnames(gff2, c("id2","rank2"))
  setkey(gff1, id1)
  setkey(gff2, id2)

  # -- merge with map

  setkey(tmp.map, id2)
  m.map <- merge(gff2, tmp.map)
  setkey(m.map, id1)
  m.map <- merge(gff1, m.map)
  m.map$unique.genome <- with(m.map, paste(genome1, genome2))
  m.map$unique.chr <- with(m.map, paste(chr1, chr2))

  # -- merge with blast
  tmp.blast <- data.table(blast)
  setkey(tmp.blast, id2)
  m.blast <- merge(gff2, tmp.blast)
  setkey(m.blast, id1)
  m.blast <- merge(gff1, m.blast)
  m.blast <- data.frame(m.blast, stringsAsFactors = F)
  m.blast$unique.genome <- with(m.blast, paste(genome1, genome2))
  m.blast$unique.chr <- with(m.blast, paste(chr1, chr2))

  # -- cull to genomes and chromosomes shared by both map and blast
  ugm <- with(m.map,
             paste(genome1, genome2,
                   chr1, chr2))
  ugb <- with(m.blast,
             paste(genome1, genome2,
                          chr1, chr2))

  gug <- intersect(unique(ugb), unique(ugm))
  wh.map <- which(ugm %in% gug)
  wh.blast <- which(ugb %in% gug)

  m.map <- data.table(m.map[wh.map,])
  m.blast <- data.table(m.blast[wh.blast,])

  # -- split the blast and map by chr, within genomes
  spl.g.map <- split(m.map, "unique.genome")
  spl.g.blast <- split(m.blast, "unique.genome")

  # -- do the work
  out.chr.blast <- lapply(names(spl.g.map), function(i){

    x.map <- spl.g.map[[i]]
    x.blast <- spl.g.blast[[i]]

    if(verbose)
      cat(x.map$genome1[1], "-->",
          x.map$genome2[1],"... ")

    spl.x.map <- split(x.map, "unique.chr")
    spl.x.blast <- split(x.blast, "unique.chr")
    if(verbose)
      cat("n.map/blast =",
          paste0("(",nrow(x.map),"/",
                 nrow(x.blast),")"),
          "... ")

    ns <- names(spl.x.map)
    blast.cull <- mclapply(ns,  mc.cores = n.cores, function(j){

      spl.y.map <- spl.x.map[[j]]
      spl.y.blast <- spl.x.blast[[j]]

      blast.out <- cull_blastBuffer(
        map.chr = spl.y.map,
        blast.chr = spl.y.blast,
        rank.buffer = rank.buffer,
        plotit = plotit & n.cores == 1,
        round2nearest =round2nearest)
      return(blast.out)
    })

    blast.cull <- rbindlist(blast.cull)

    if(verbose)
      cat("culled to", nrow(blast.cull), "\n")
    return(blast.cull)
  })

  return(rbindlist(out.chr.blast))
}

