#' @title Cull blast to syntenic regions
#'
#' @description
#' \code{cull_syntenicBlast} An internal function, designed to be called
#' by extend_blocks and find_syntenicOrthogs.
#'
#' @param map the map data.table or data.frame. This is used to infer
#' syntenic regions.
#' @param blast the blast dataset to screen for syntenic hits
#' @param gff gff data.table
#' @param rank.buffer The buffer, in rank order of gff annotation
#' data.table.
#' @param plotit Logical, should plots be made? If so, and n.cores
#' is 1, then plotting is done.
#' @param n.cores The number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @import sp
#' @importFrom raster buffer extent
#' @import rgeos
#' @export
cull_syntenicBlast <- function(map,
                               blast,
                               gff,
                               rank.buffer = 250,
                               verbose = T,
                               plotit = F,
                               n.cores = 1){
  #######################################################
  #######################################################
  rerank_fromIDs <- function(id1,
                             id2,
                             gff){
    u.id <- unique(c(id1, id2))
    if ("rank" %in% colnames(gff))
      gff$rank <- NULL

    gff[, rank := as.numeric(frank(start,
                                   ties.method = "dense")),
        by = list(genome, chr)]
    gff <- gff[gff$id %in% u.id, ]

    gff1 <- data.table(gff)
    gff2 <- data.table(gff)
    setnames(gff1, paste0(colnames(gff1), "1"))
    setnames(gff2, paste0(colnames(gff2), "2"))
    setkey(gff1, id1)
    setkey(gff2, id2)

    id.dt <- data.table(id1 = id1,
                        id2 = id2,
                        stringsAsFactors = F)

    setkey(id.dt, id2)
    m1 <- merge(gff2, id.dt)
    setkey(m1, id1)
    out <- merge(gff1, m1)

    out$unique.genome <- with(out, paste(genome1, genome2))
    out$unique.chr <- with(out, paste(genome1, genome2,
                                      chr1, chr2))

    return(out)
  }
  #######################################################
  #######################################################
  mround <- function(x,base){
    base * round(x / base)
  }
  #######################################################
  #######################################################
  buffer_xy <- function(x,
                        y,
                        buffer,
                        roundto = buffer / 3,
                        crs = CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")){
    xy <- data.frame(
      x = mround(x, roundto),
      y = mround(y, roundto),
      stringsAsFactors = F)

    xy <- xy[!duplicated(xy), ]

    coordinates(xy) <- ~ x + y
    proj4string(xy) <- crs

    buf <- buffer(xy, buffer)
    return(buf)
  }
  #######################################################
  #######################################################
  find_xy2keep <- function(buffer.sp,
                           x,
                           y,
                           plotit = T,
                           xlab = NULL,
                           ylab = NULL){

    xy.pts <- data.frame(x = x,
                         y = y,
                         stringsAsFactors = F)

    coordinates(xy.pts) <- ~ x + y
    proj4string(xy.pts) <- proj4string(buffer.sp)

    # -- check which hits are in the buffer
    is.in.buf <- !is.na(over(xy.pts,
                             buffer.sp,
                             fn = NULL))

    # -- if desired, make the plot
    if (plotit) {
      xy2keep <- xy.pts[is.in.buf,]
      ylims <- extent(buffer.sp)[3:4]
      xlims <- extent(buffer.sp)[1:2]

      plot(jitter(xy.pts$x),
           jitter(xy.pts$y),
           pch = ".",
           ylim = ylims,
           xlim = xlims,
           xlab = xlab,
           ylab = ylab)

      plot(buffer.sp,
           add = T,
           col = rgb(0, 0, 1, .2),
           lty = 2)
    }
    return(is.in.buf)
  }
  #######################################################
  #######################################################
  buffer_withinGenome <- function(map,
                                  blast,
                                  buffer,
                                  n.cores = 1,
                                  verbose = TRUE,
                                  plotit = TRUE){
    if (verbose)
      cat(map$genome1[1], "-->",
          map$genome2[1], "... ")

    s.m <- split(map, "unique.chr")
    s.b <- split(blast, "unique.chr")

    if (verbose)
      cat("n.map/blast =",
          nrow(map),
          "/", nrow(blast), "... ")

    ns <- names(s.m)
    blast.cull <- mclapply(ns,  mc.cores = n.cores, function(i){

      tmp.map <- s.m[[i]]
      tmp.blast <- s.b[[i]]

      i.xlab <- paste0(tmp.map$genome1[1],
                       " (", tmp.map$chr1[1], ")")
      i.ylab <- paste0(tmp.map$genome2[1],
                       " (", tmp.map$chr2[1], ")")
      buf <- with(tmp.map,
                  buffer_xy(x = rank1,
                            y = rank2,
                            buffer = buffer))
      wh2keep <- with(tmp.blast,
                      find_xy2keep(x = rank1,
                                   y = rank2,
                                   buffer.sp = buf,
                                   plotit = plotit,
                                   xlab = i.xlab,
                                   ylab = i.ylab))

      return(list(blast = tmp.blast[wh2keep,],
                  buf = buf))
    })

    names(blast.cull) <- ns
    buf <- sapply(ns, USE.NAMES = T, simplify = F, function(i)
      blast.cull[[i]]$buf)
    blast.cull <- rbindlist(lapply(blast.cull, function(x)
      x$blast))

    if (verbose)
      cat("culled to",
          nrow(blast.cull), "\n\t")
    return(list(blast = blast.cull,
                buffer.sp = buf))
  }
  #######################################################
  #######################################################
  buffer_all <- function(map,
                         blast,
                         buffer,
                         plotit = FALSE,
                         n.cores = 1,
                         verbose = TRUE){
    spl.map <- split(map, "unique.genome")
    spl.blast <- split(blast, "unique.genome")
    ns <- unique(names(spl.blast))
    ns <- ns[ns %in% unique(names(spl.map))]
    out.all <- lapply(ns, function(i){
      i.map <- spl.map[[i]]
      i.blast <- spl.blast[[i]]
      out <- buffer_withinGenome(map = i.map,
                                 blast = i.blast,
                                 buffer = buffer,
                                 plotit = plotit,
                                 n.cores = n.cores)
      return(out)
    })

    blast.out <- rbindlist(lapply(out.all, function(x) x$blast))
    buf.out <- lapply(out.all, function(x) x$buf)
    buf.out <- unlist(buf.out, recursive = F)
    return(list(blast = blast.out,
                buffer.sp = buf.out))
  }
  #######################################################
  if (verbose)
    cat("Making new blast and map with gff-based ranks ... ")

  r.map <- with(map,
                rerank_fromIDs(id1 = id1,
                               id2 = id2,
                               gff = gff))
  r.blast <- with(blast,
                  rerank_fromIDs(id1 = id1,
                                 id2 = id2,
                                 gff = gff))
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Splitting map and blast by chr / genome combinations ... ")
  ugm <- with(r.map,
              paste(genome1, genome2,
                    chr1, chr2))
  ugb <- with(r.blast,
              paste(genome1, genome2,
                    chr1, chr2))
  gug <- intersect(unique(ugb), unique(ugm))
  wh.map <- which(ugm %in% gug)
  wh.blast <- which(ugb %in% gug)
  m.map <- data.table(r.map[wh.map,])
  m.blast <- data.table(r.blast[wh.blast,])


  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Culling by pairwise genome comparison ...\n\t")
  out <- buffer_all(map = m.map,
                    blast = m.blast,
                    plotit = plotit,
                    verbose = verbose,
                    n.cores = n.cores,
                    buffer = rank.buffer)

  if (verbose)
    cat("Done!\n")
  out.blast <- out$blast
  re.out <- with(out.blast,
                 rerank_fromIDs(id1 = id1,
                                id2 = id2,
                                gff = gff))
  #######################################################
  return(re.out)
}
