#' @title plotting utility functions
#' @description
#' \code{plot_utils} Functions that allow selection stat calculation in GENESPACE
#' @name genespace_utils
#'
#' @param alpha numeric, length
#' @param bg.col character, length
#' @param blk block-formatted data.table
#' @param blk.border numeric, length
#' @param blk.ids character, length
#' @param border.color character, length
#' @param chr character, length
#' @param chr.abbrev.fun function to
#' @param chr.bg.cex numeric, length
#' @param chr.bg.col character, length
#' @param chr.bg.pch numeric, length
#' @param chr.buff numeric, length
#' @param chr.buffer numeric, length
#' @param chr.id.cex numeric, length
#' @param chr.id.col character, length
#' @param chr.lab.buff numeric, length
#' @param chr.list list of chromosomes to plot
#' @param chr.segm.col character, length
#' @param chr.segm.lwd numeric, length
#' @param col character, length
#' @param cols character, length
#' @param do.cumulative logical, should
#' @param dodge.geneIDs numeric, length
#' @param dodge.x numeric, length
#' @param e1 numeric, length
#' @param e2 numeric, length
#' @param end numeric, length
#' @param fais fai-like data.table
#' @param fill.color character, length
#' @param forCircos logical, should
#' @param gap.prop numeric, length
#' @param gene.colors character, length
#' @param geneID.abbrev.fun function to
#' @param geneid.cex numeric, length
#' @param geneid.offset numeric, length
#' @param genes2plot character, length
#' @param genomeIDs character, length
#' @param genomes character, length
#' @param gff gff-like data.table
#' @param lab.chr logical, should
#' @param lab.chr.1only logical, should
#' @param map map-formatted data.table
#' @param n.out numeric, length
#' @param n.reps numeric, length
#' @param n.sample numeric, length
#' @param num numeric, length
#' @param ortho.col character, length
#' @param points.per.curve numeric, length
#' @param return.start logical, should
#' @param s1 numeric, length
#' @param s2 numeric, length
#' @param scale.it logical, should
#' @param scale2dodge numeric, length
#' @param simplify.poly logical, should
#' @param start numeric, length
#' @param use.rank logical, should
#' @param verbose logical, should
#' @param x numeric or character, length
#' @param y numeric, length
#' @param y.end numeric, length
#' @param y.start numeric, length
#'
#' @note \code{plot_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{plotutils} returns its own arguments.
#'
#'


#' @title format_gffChrlist
#' @description
#' \code{format_gffChrlist} rename and simplify gff
#' @rdname genespace_utils
#' @export
format_gffChrlist <- function(gff,
                              genomes,
                              chr.list,
                              use.rank,
                              gap.prop,
                              do.cumulative = T){
  gff <- rbindlist(lapply(1:length(genomes), function(i){
    tmp <- subset(gff, genome == genomes[i] &
                    chr %in% chr.list[[i]])
    tmp$genome <- names(chr.list)[i]
    return(tmp)
  }))

  if (use.rank) {
    gff[,start := frank(start, ties.method = "random"),
        by = list(genome, chr)]
    gff[,end := frank(end, ties.method = "random"),
        by = list(genome, chr)]
  }
  gff$chr <- as.character(gff$chr)
  if (do.cumulative) {
    gff <- convert_gff2coords(gff = gff,
                              chr.list,
                              gap.prop = gap.prop)
  }else{
    gff$start.p <- gff$start
    gff$end.p <- gff$end
  }

  return(gff)
}

#' @title format_mapChrlist
#' @description
#' \code{format_mapChrlist} rename and simplify map
#' @rdname genespace_utils
#' @export
format_mapChrlist <- function(map,
                              gff,
                              chr.list,
                              genomes,
                              forCircos = F,
                              do.cumulative = T){
  if(length(genomes) == 1){
    map <- subset(map, genome1 == genomes &
                    chr1 %in% chr.list[[1]] &
                    genome2 == genomes &
                    chr2 %in% chr.list[[1]])
    map$genome1 <- names(chr.list)
    map$genome2 <- names(chr.list)
  }else{
    map <- rbindlist(lapply(1:(length(genomes) - 1), function(i){
      tmp <- subset(map, genome1 == genomes[i] &
                      chr1 %in% chr.list[[i]] &
                      genome2 == genomes[i + 1] &
                      chr2 %in% chr.list[[i + 1]])
      tmp$genome1 <- names(chr.list)[i]
      tmp$genome2 <- names(chr.list)[i + 1]
      return(tmp)
    }))
  }

  gi <- genomes
  genomes <- as.character(names(chr.list))

  gf1 <- with(gff,
              data.table(id1 = id,
                         chr1 = chr,
                         pos1 = start.p))
  gf2 <- with(gff,
              data.table(id2 = id,
                         chr2 = chr,
                         pos2 = start.p))

  map <- map[!duplicated(map[,c("block.id",
                                "genome1","genome2",
                                "id1","id2")]),]
  gf1 <- gf1[!duplicated(gf1), ]
  gf2 <- gf2[!duplicated(gf2), ]

  map <- merge(gf2,
               map[,c("block.id",
                      "genome1","genome2",
                      "id1","id2")],
               by = c("id2"))
  map <- merge(gf1,
               map, by = c("id1"))

  if (forCircos) {
    if(length(genomes) > 1){
      combns <- data.table(t(combn(genomes,2)))
      setnames(combns, c("genome1", "genome2"))
      map <- merge(combns,
                   map,
                   by = c("genome1", "genome2"))
    }
  }else{
    combs <- data.table(y1 = 0:(length(genomes) - 2),
                        y2 = 1:(length(genomes) - 1),
                        genome1 = genomes[1:(length(genomes) - 1)],
                        genome2 = genomes[2:length(genomes)],
                        stringsAsFactors = F)
    map <- rbind(map,
                 with(map,
                      data.table(block.id = block.id,
                                 genome1 = genome2,
                                 genome2 = genome1,
                                 id1 = id2,
                                 id2 = id1,
                                 chr1 = chr2,
                                 chr2 = chr1,
                                 pos1 = pos2,
                                 pos2 = pos1,
                                 stringsAsFactors = F)))
    map <- map[!duplicated(map),]
    map <- merge(combs,
                 map,
                 by = c("genome1", "genome2"))
  }

  return(map)
}

#' @title convert_gff2coords
#' @description
#' \code{convert_gff2coords} makes blocks and merges by block.id
#' @rdname genespace_utils
#' @export
convert_gff2coords <- function(gff,
                               chr.list,
                               gap.prop){

  spl <- split(gff, by = "genome")
  out <- rbindlist(lapply(names(spl), function(i){
    x <- spl[[i]]
    y <- chr.list[[i]]
    x$chr <- factor(x$chr, levels = y)
    setkey(x, chr)
    x[, start.l := calc_linearCoord(chr = chr,
                                    start = start,
                                    end = end,
                                    gap.prop = gap.prop,
                                    scale.it = F,
                                    return.start = T)]
    x[, end.l := calc_linearCoord(chr = chr,
                                  start = start,
                                  end = end,
                                  gap.prop = gap.prop,
                                  scale.it = F,
                                  return.start = F)]
    x[, start.p := calc_linearCoord(chr = chr,
                                    start = start,
                                    end = end,
                                    gap.prop = gap.prop,
                                    scale.it = T,
                                    return.start = T)]
    x[, end.p := calc_linearCoord(chr = chr,
                                  start = start,
                                  end = end,
                                  gap.prop = gap.prop,
                                  scale.it = T,
                                  return.start = F)]
    x$chr <- as.character(x$chr)
    return(x)
  }))
  return(out)
}


#' @title calc_linearCoord
#' @description
#' \code{calc_linearCoord} convert start / end positions into
#' linear coordinates..
#' @rdname genespace_utils
#' @export
calc_linearCoord <- function(chr,
                             start,
                             end,
                             gap.prop,
                             scale.it,
                             return.start){
  ml <- data.table(chr = chr,
                   start = start,
                   end = end,
                   stringsAsFactors = F)
  fais <- ml[, list(start = 0,
                    end = max(end)),
             by = list(chr)]
  setkey(fais, chr)
  fais <- add_gap(fais = fais,
                  gap.prop = gap.prop)
  maxf <- max(fais$endl)
  if (!scale.it)
    maxf <- 1

  bt <- merge(fais[, c("chr","startl")], ml, by = "chr")

  if (return.start) {
    return((bt$start + bt$startl) / maxf)
  }else{
    return((bt$end + bt$startl) / maxf)
  }
}

#' @title add_gap
#' @description
#' \code{add_gap} add gap to linear chromosome x axis.
#' @rdname genespace_utils
#' @export
add_gap <- function(fais,
                    gap.prop){
  gap <- sum(fais$end) * gap.prop
  fais$endl <- (cumsum(fais$end + gap)) - gap
  fais$startl <- fais$endl - fais$end
  return(fais)
}

#' @title make_fais
#' @description
#' \code{make_fais} make data.table akin to fasta index.
#' @rdname genespace_utils
#' @export
make_fais <- function(genomes,
                      gff){
  ypos <- data.table(y = 0:(length(genomes) - 1),
                     genome = genomes)

  fais <- gff[,list(start = min(start.p),
                    end = max(start.p)),
              by = list(genome, chr)]
  fais <- merge(fais,
                ypos,
                by = "genome")
  fais$genome <- factor(fais$genome, levels = genomes)
  setkey(fais, genome)
  fais$genome <- as.character(fais$genome)
  return(fais)
}

#' @title draw_blkPolygon
#' @description
#' \code{draw_blkPolygon} draw blk Polygon
#' @rdname plot_utils
#' @export
draw_blkPolygon <- function(blk,
                            chr.buffer,
                            blk.border,
                            blk.col,
                            simplify.poly,
                            points.per.curve){
  blk[,unique := paste(genome1, genome2)]
  spl <- split(blk, by = "unique")
  for (j in names(spl)) {
    tmp <- spl[[j]]
    ys <- chr.buffer + tmp$y1[1]
    ye <- (1 - chr.buffer) + tmp$y1[1]
    y <- cos_y(y.start = ys,
               y.end = ye,
               n.out = points.per.curve)

    for (i in 1:nrow(tmp)) {
      with(tmp[i,], make_polygon(
        s1 = start1, e1 = end1,
        s2 = start2, e2 = end2,  y = y,
        fill.color = blk.col,
        simplify.poly = simplify.poly,
        border.color = blk.border))
    }
  }
}

#' @title draw_genePolygon
#' @description
#' \code{draw_genePolygon} draw gene Polygon
#' @rdname genespace_utils
#' @export
draw_genePolygon <- function(map,
                             genes2plot,
                             chr.buff,
                             gene.col,
                             simplify.poly,
                             points.per.curve){

  map[,unique := paste(genome1, genome2)]
  map <- subset(map, id1 %in% genes2plot | id2 %in% genes2plot)
  spl <- split(map, by = "unique")
  for (j in names(spl)) {
    tmp <- spl[[j]]
    ys <- chr.buff + tmp$y1[1]
    ye <- (1 - chr.buff) + tmp$y1[1]
    y <- cos_y(y.start = ys,
               y.end = ye,
               n.out = points.per.curve)

    for (i in 1:nrow(tmp)) {
      with(tmp[i,], make_polygon(
        s1 = pos1,
        e1 = pos1,
        s2 = pos2,
        e2 = pos2,
        y = y,
        fill.color = gene.col,
        simplify.poly = simplify.poly,
        border.color = gene.col))
    }
  }
}

#' @title annotate_riparian
#' @description
#' \code{annotate_riparian} annotate_riparian
#' @rdname genespace_utils
#' @export
annotate_riparian <- function(map,
                              genomes,
                              chr.buff,
                              dodge.geneIDs,
                              scale2dodge,
                              geneid.cex,
                              geneid.offset,
                              gene.colors,
                              geneID.abbrev.fun,
                              genes2plot){
  mapg1 <-  with(subset(map,
                        genome1 == genomes[1] &
                          !duplicated(id1) &
                          id1 %in% genes2plot),
                 data.table(genome = genome1,
                            id = id1,
                            y = y1,
                            x = pos1))
  mapg2 <-  with(subset(map,
                        genome2 == genomes[length(genomes)] &
                          !duplicated(id2) &
                          id2 %in% genes2plot),
                 data.table(genome = genome2,
                            id = id2,
                            y = y2,
                            x = pos2))
  mapg1$x1 <- get_xClus(x = mapg1$x,
                        n.reps = 5,
                        dodge.x = dodge.geneIDs,
                        scale2dodge = scale2dodge)
  mapg2$x1 <- get_xClus(x = mapg2$x,
                        n.reps = 5,
                        dodge.x = dodge.geneIDs,
                        scale2dodge = scale2dodge)

  mapg1$y1 <- mapg1$y - geneid.offset
  mapg1$y <- mapg1$y - chr.buff

  mapg2$y1 <- mapg2$y + geneid.offset
  mapg2$y <- mapg2$y + chr.buff

  with(rbind(mapg1, mapg2),
       segments(y0 = y,
                y1 = y1,
                x0 = x,
                x1 = x1))


  mapg1$lab1 <- geneID.abbrev.fun(mapg1$id)
  # for (i in 1:nrow(mapg1)) {
  #   if (mapg1$id[i] %in% genes2plot) {
  #     if (is.null(names(genes2plot))) {
  #       mapg1$lab1[i] <- mapg1$id[i]
  #     }else{
  #       mapg1$lab1[i] <- names(genes2plot)[genes2plot == mapg1$id[i]]
  #     }
  #   }
  # }
  mapg2$lab2 <- geneID.abbrev.fun(mapg2$id)
  # for (i in 1:nrow(mapg2)) {
  #   if (mapg2$id[i] %in% genes2plot) {
  #     if (is.null(names(genes2plot))) {
  #       mapg2$lab2[i] <- mapg2$id[i]
  #     }else{
  #       mapg2$lab2[i] <- names(genes2plot)[genes2plot == mapg2$id[i]]
  #     }
  #   }
  # }

  with(mapg1,
       text(lab1,
            y = y1,  x = x1, srt = 90,
            cex = geneid.cex, adj = c(1.05,.5),
            col = gene.colors[id]))

  with(mapg2,
       text(lab2,
            y = y1,  x = x1, srt = 90,
            cex = geneid.cex, adj = c(-0.05,.5),
            col = gene.colors[id]))

}

#' @title label_riparian
#' @description
#' \code{label_riparian} label_riparian
#' @rdname genespace_utils
#' @export
label_riparian <- function(fais,
                           chr.segm.lwd,
                           chr.segm.col,
                           genomes,
                           chr.lab.buff,
                           chr.bg.col,
                           chr.bg.cex,
                           chr.bg.pch,
                           chr.id.col,
                           chr.id.cex,
                           lab.chr,
                           lab.chr.1only,
                           chr.abbrev.fun){

  if (length(chr.segm.col) == length(genomes)) {
    names(chr.segm.col) <- genomes
    chr.segm.col <- chr.segm.col[match(fais$genome, names(chr.segm.col))]
  }else{
    chr.segm.col <- chr.segm.col[1]
  }

  if (length(chr.segm.lwd) == length(genomes)) {
    names(chr.segm.lwd) <- genomes
    chr.segm.lwd <- chr.segm.lwd[match(fais$genome, names(chr.segm.lwd))]
  }else{
    chr.segm.lwd <- chr.segm.lwd[1]
  }

  # add in chromosomes and IDs
  with(fais, segments(x0 = start,
                      x1 = end,
                      y0 = y,
                      y1 = y,
                      lwd = chr.segm.lwd,
                      col = chr.segm.col))

  fais$class <- with(fais,
                     ifelse(genome == genomes[1],
                            "bottom",
                            ifelse(genome == genomes[length(genomes)],
                                   "top","middle")))
  fais$lab.y <- with(fais,
                     ifelse(class == "middle", y,
                            ifelse(class == "top", y - chr.lab.buff,
                                   y + chr.lab.buff)))

  if (!is.null(chr.bg.col) & lab.chr) {
    points(x = rowMeans(fais[,c("start","end")]),
           y = fais$lab.y,
           col = chr.bg.col,
           pch = chr.bg.pch,
           cex = chr.bg.cex)
  }
  if (lab.chr) {
    if (lab.chr.1only) {
      fais2 <- fais[fais$genome == genomes[1],]
    }else{
      fais2 <- fais
    }
    text(chr.abbrev.fun(fais2$chr),
         x = rowMeans(fais2[,c("start","end")]),
         y = fais2$lab.y,
         col = chr.id.col,
         cex = chr.id.cex)
  }
}

#' @title add_alpha
#' @description
#' \code{add_alpha} add transparency to a color.
#' @rdname genespace_utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col,
                      alpha = 1){
  if (missing(col))
    stop("Please provide a vector of colours.")

  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1],
              x[2],
              x[3],
              alpha = alpha))
}


#' @title make_polygon
#' @description
#' \code{make_polygon} given a y string and start/end coords, draw a polygon
#' @rdname genespace_utils
#' @importFrom graphics polygon
#' @export
make_polygon <- function(s1,
                         s2,
                         e1,
                         e2,
                         y,
                         fill.color,
                         border.color,
                         simplify.poly){
  y <- y[!is.na(y)]
  ye <- max(y)
  ys <- min(y)

  x1 <- seq(from = s1,
            to = s2,
            length.out = length(y))
  x3 <- seq(from = e2,
            to = e1,
            length.out = length(y))

  x2 <- c(s2,e2)
  x4 <- c(e1,s1)

  y1 <- y
  y2 <- c(ye,ye)
  y3 <- rev(y)
  y4 <- c(ys,ys)
  if (!is.null(simplify.poly)) {
    qs <- qnorm(seq(.0001,
                    .9999,
                    length.out = length(y)))
    wh <- which(!duplicated(round_toAny(qs, simplify.poly)))
    y1 <- y1[wh]
    y3 <- y3[wh]
    x1 <- x1[wh]
    x3 <- x3[wh]
  }

  polygon(c(x1,x2,x3,x4),
          c(y1,y2,y3,y4),
          col = fill.color,
          border = border.color)
}

#' @title round_toAny
#' @description
#' \code{round_toAny} round_toAny
#' @rdname genespace_utils
#' @importFrom dbscan frNN dbscan
#' @export
round_toAny <- function(x,
                        num){
  num * round(x / num)
}

#' @title get_xClus
#' @description
#' \code{get_xClus} for a set of x,
#' @rdname genespace_utils
#' @importFrom dbscan frNN dbscan
#' @export
get_xClus <- function(x,
                      n.reps,
                      dodge.x,
                      scale2dodge = 2){
  for (i in n.reps) {
    clus <- run_dbs(map = data.table(rank1 = as.numeric(x),
                                   rank2 = as.numeric(x)),
                    radius = dodge.x*scale2dodge,
                    n.mappings = 1)$cluster
    out <- data.table(x = x,
                      clus = clus)
    out[,round := mean(x), by = list(clus)]
    out[, rank := as.numeric(frank(x, ties.method = "random")),
        by = list(clus)]
    out[,n.inclus := .N,
        by = list(clus)]
    out[,n.inclus := .N,
        by = list(clus)]
    out[, rank := ifelse(n.inclus == 1, 1,(rank - mean(rank)) + 1),
        by = list(clus)]
    out[, x1 := (rank * dodge.x + round) - dodge.x,
        by = list(clus)]
    x <- out$x1
  }
  return(x)
}

#' @title cos_y
#' @description
#' \code{cos_y} given a range of x, return a cosine curve.
#' @rdname genespace_utils
#' @export
cos_y <- function(y.start,
                  y.end,
                  n.out,
                  n.sample = n.out){
  rs <- seq(from = 0,
            to = pi,
            length.out = n.sample)
  xt <- 1 - cos(rs)
  xo <- seq(from = min(xt),
            to = max(xt),
            length.out = n.out)
  yvals <- sapply(xo, function(x)
    which.min(abs(x - xt)))
  y <- y.start + (yvals - min(yvals)) *
    ((y.end - y.start) / (max(yvals) - min(yvals)))
  return(y)
}
