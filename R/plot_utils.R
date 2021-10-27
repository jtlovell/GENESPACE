#' @title Functions to make genespace plots
#' @description
#' \code{plot_utils} Internal ad hoc functions for plotting routines
#'
#' @name plot_utils
#'
#' @param gsParam a list containing all parameters for a GENESPACE run. See
#' init_genespace
#' @param genomeIDs an optional vector of genomeIDs to consider. If not
#' specified (default) taken from gsParam$genomeIDs$genomeIDs
#' @param hits data.table containing annotated blast-format pairwise hits
#' @param blkSize integer of length 1 specifying the minimum size for a syntenic
#' block and the -s 'size' MCScanX parameter
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param useBlks logical, for cross compatibility with plot_hits
#' @param refGenome character string matching one of the genomeIDs in gsParam
#' @param y y position of the scale bar
#' @param interpTails logical, should unbounded tails be interpolated?
#' @param xspan amount of span on the x axis
#' @param yspan amount of span on the y axis
#' @param label scale bar label
#' @param cex scale bar label character expansion
#' @param lwd line thickness for scale bar
#' @param xleft numeric, specifying the coordinate of left x position
#' @param ybottom numeric, specifying the coordinate of lower y position
#' @param xright numeric, specifying the coordinate of right x position
#' @param ytop numeric, specifying the coordinate of upper y position
#' @param start1 numeric, specifying the coordinate of blk start in genome1
#' @param end1 numeric, specifying the coordinate of  blk end in genome1
#' @param start2 numeric, specifying the coordinate of  blk start in genome2
#' @param end2 numeric, specifying the coordinate of blk end in genome2
#' @param y1 numeric, specifying the coordinate of y position in genome1
#' @param y2 numeric, specifying the coordinate of  y position in genome1
#' @param hitsRef a hits data table against the reference
#' @param hitsRip a hits data table along the riparian daisy chain
#' @param chrd data.table of chromosomes for re-ordering
#' @param chrl list of chromosomes for re-ordering
#' @param minrl minimum run length for block breaking
#' @param minprp minimum proportion of overlapping hits
#' @param minGenesOnChr integer specifying the minimum number of hits on a
#' chromosome for it to be displayed
#' @param plotRegions logical, should the regions be plot
#' @param refHits a hits data table against the reference
#' @param gapProp number (0-1) specifying the proportion of the the total
#' plot that should be gaps.
#' @param minHits2plot integer specifying the minimum number of hits on a
#' chromosome for it to be displayed
#' @param lightChrFill color for chromosome backgrounds with the fewest hits
#' @param darkChrFill color for chromosome backgrounds with the most hits
#' @param emptyChrFill color for chromosome backgrounds without any hits
#' @param x vector
#' @param col integer or character coercible to an R color
#' @param alpha numeric (0-1) specifying the transparency

#' @note \code{plot_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{plot_utils} returns its own arguments.


#' @title split_blksByRef
#' @description
#' \code{split_blksByRef} split_blksByRef
#' @rdname plot_utils
#' @export
split_blksByRef <- function(hitsRef,
                            hitsRip,
                            gff,
                            chrd,
                            chrl,
                            refGenome,
                            minrl = 5,
                            minprp = .9){
  refChr <- NULL
  gen1 <- gen2 <- chr1 <- chr2 <- chr <- genome <- start <- end <- ord2 <- NULL
  ord1 <- n <- refChro <- nu <- ofID1 <- ofID2 <- rl <- blkID <- NULL

  chro <- hitsRef[,list(genome = gen2[1], chr = chr2[1],
                        start = min(ord2), end = max(ord2)),
                  by = c("gen1","chr1","blkID")]
  gffo <- with(subset(gff, paste(genome, chr) %in% paste(chrd$genome, chrd$chr)), data.table(
    genome = genome, chr = chr, start = ord, end = ord, ofID = ofID))
  setkey(gffo, genome, chr, start, end)
  setkey(chro, genome, chr, start, end)
  fo <- foverlaps(gffo, chro)
  fo <- subset(fo, complete.cases(fo))
  fo <- subset(fo, paste(gen1, chr1) %in% paste(chrd$genome, chrd$chr))
  fo <- subset(fo, paste(genome, chr) %in% paste(chrd$genome, chrd$chr))
  fo <- with(fo, data.table(refChr = chr1, ofID2 = ofID))
  setkey(fo, refChr, ofID2)
  fo[,n := .N, by = "ofID2"]
  fo[,refChro := match(refChr, chrl[[refGenome]])]
  setkey(fo, refChro)
  fo <- subset(fo, !duplicated(ofID2))
  fo[,`:=`(n = NULL, refChro = NULL)]

  hitsRip <- merge(hitsRip, fo, by = "ofID2", all.x = T)
  hitsRip[,nu := uniqueN(refChr, na.rm = T), by = "blkID"]
  ho <- subset(hitsRip, nu > 1)

  if(nrow(ho) > 0){
    setkey(ho, ord1)
    spl <- split(ho, by = "blkID")
    hitsRip2 <- rbindlist(lapply(spl, function(x){
      x <- subset(x, !is.na(refChr))
      m <- find_modalValue(x$refChr)
      p <- sum(x$refChr == m) / nrow(x)
      if(p > minprp){
        x[,refChr := m]
      }else{
        x[,rl := add_rle(refChr)]
        x$refChr[x$rl < minrl] <- m
        x[,rl := add_rle(refChr,  which = "id")]
        x[,blkID := paste0(blkID, "_", rl)]
        x[,rl := NULL]
      }
      return(x)
    }))
    hitsRip2[,blkID := paste(blkID, refChr)]
    hitsRip <- rbind(subset(hitsRip, nu == 1), hitsRip2)
  }
  hitsRip[,nu := NULL]
  return(hitsRip)
}

#' @title reorder_gff
#' @description
#' \code{reorder_gff} reorder_gff
#' @rdname plot_utils
#' @export
reorder_gff <- function(gff, genomeIDs, minGenesOnChr, refGenome){
  genome <- ord <- n <- chr <- chrn <- medPos <- chrNameOrd <- og <- NULL
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)
  nGenes <- gff[,list(n = .N, medPos = median(ord)), by = c("genome","chr")]
  nGenes <- subset(nGenes, n > minGenesOnChr)
  nGenes[,chrn := as.numeric(gsub("[^0-9]", "", chr))]
  setorder(nGenes, genome, chrn, chr, medPos, -n)
  refChrs <- nGenes$chr[nGenes$genome == refGenome]
  ugc <- with(nGenes, paste(genome, chr))
  gff[,chrNameOrd := as.numeric(factor(paste(genome, chr), levels = ugc))]
  if(!"og" %in% colnames(gff))
    gff[,og := NA]
  gff <- subset(gff, !is.na(chrNameOrd))[,c("genome","chr","chrNameOrd","ofID","ord","start","end","og")]
  setkey(gff, genome, chrNameOrd, ord)

  gff[,ord := 1:.N, by = "genome"]
  return(gff)
}

#' @title load_refHits
#' @description
#' \code{load_refHits} load_refHits
#' @rdname plot_utils
#' @export
load_refHits <- function(gsParam, genomeIDs, refGenome, plotRegions){
  regID <- og <- blkAnchor <- blkID <- regAnchor <- NULL
  nonRefGen <- genomeIDs[genomeIDs != refGenome]
  fs <- file.path(gsParam$paths$results,
                  sprintf("%s_%s_synHits.txt.gz",
                          c(rep(refGenome, length(nonRefGen)), nonRefGen),
                          c(nonRefGen, rep(refGenome, length(nonRefGen)))))
  fs <- fs[file.exists(fs)]
  hitsRef <- rbindlist(lapply(fs, function(i){
    if(plotRegions){
      x <- subset(fread(
        i, select = c("gen1","ofID1","ofID2","regAnchor","regID","og"),
        na.strings = c("NA","")),
        regAnchor & !is.na(regID) & !is.na(og))
      setnames(x, c("regID", "regAnchor"), c("blkID", "blkAnchor"))
    }else{
      x <- subset(fread(
        i, select = c("gen1","ofID1","ofID2","blkAnchor","blkID","og"),
        na.strings = c("NA","")),
        blkAnchor & !is.na(blkID) & !is.na(og))
    }
    if(x$gen1[1] != refGenome)
      setnames(x, c("ofID1","ofID2"), c("ofID2","ofID1"))
    x <- x[,c("ofID1","ofID2","blkID","og")]

    return(x)
  }))
  return(hitsRef)
}

#' @title cosine curve source data
#' @description
#' \code{add_synChr2gff} vector of points for polygons based on cosine curves
#' @rdname plot_utils
#' @export
add_synChr2gff <- function(gff, refHits, genomeIDs, gapProp, refGenome){
  ord1 <- ord <- gen2 <- genome <- refOrd <- chrOrd <- chr <- end <- gap <- NULL
  start <- linOrd <- linBp <- NULL

  synChrs <- refHits[,list(refOrd = median(ord1, na.rm = T)), by = c("gen2", "chr2")]
  synChrs <- subset(synChrs, complete.cases(synChrs))
  synChrs[,`:=`(genome = factor(gen2, levels = genomeIDs), gen2 = NULL)]
  setkey(synChrs, genome, refOrd)
  synChrs[,chrOrd := 1:.N, by = "genome"]
  ugcSyn <- synChrs$chrOrd; names(ugcSyn) <- with(synChrs, paste(genome, chr2))
  tmp <- subset(subset(gff, genome == refGenome), !duplicated(chr))
  rv <- tmp$chrNameOrd; names(rv) <- with(tmp, paste(refGenome, chr))
  ugcSyn <- c(rv, ugcSyn)
  gff[,chrOrd := ugcSyn[paste(genome, chr)]]
  setkey(gff, genome, chrOrd, ord)
  gff[,ord := 1:.N, by = "genome"]

  nbp <- gff[,list(nbp = max(end)), by = c("chr","genome","chrOrd")]
  setkey(nbp, genome, chrOrd)
  gps <- with(nbp, tapply(nbp, genome, sum))
  gps <- max(gps) * ((max(gps) * gapProp)/gps)
  nbp[,gap := gps[genome]]
  nbp[,start := c(0, cumsum(nbp[-.N] + gap[-.N])), by = "genome"]
  lbpv <- nbp$start; names(lbpv) <- with(nbp, paste(genome, chr))
  gff[,`:=`(linBp = start + lbpv[paste(genome, chr)])]

  gps <- with(gff, tapply(ord, genome, max))
  gps <- max(gps) * ((max(gps) * gapProp)/gps)
  gff[,linOrd := ord + (gps[genome] * (chrOrd - 1))]
  gff[,`:=`(linBp = linBp - mean(range(linBp, na.rm = T)),
            linOrd = linOrd - mean(range(linOrd, na.rm = T))),
      by = "genome"]
  return(gff)
}

#' @title read_ripHits
#' @description
#' \code{read_ripHits} read_ripHits
#' @rdname plot_utils
#' @export
read_ripHits <- function(gsParam, genomeIDs, useBlks){
  og <- regAnchor <- regID <- blkAnchor <- blkID <- NULL
  genomeOrd <- data.table(
    gen1 = genomeIDs[-length(genomeIDs)],
    gen2 = genomeIDs[-1],
    y = 1:length(genomeIDs[-1]))

  fs <- with(genomeOrd, file.path(gsParam$paths$results,
                                  sprintf("%s_%s_synHits.txt.gz",
                                          c(gen1, gen2),
                                          c(gen2, gen1))))
  fs <- fs[file.exists(fs)]
  hitsRip <- rbindlist(lapply(fs, function(i){
    if(useBlks){
      x <- subset(fread(
        i, select = c("ofID1","ofID2","gen1","gen2","blkID","blkAnchor","og"),
        na.strings = c("NA","")),
        blkAnchor & !is.na(blkID) & !is.na(og))
    }else{
      x <- subset(fread(
        i, select = c("ofID1","ofID2","gen1","gen2","regID","regAnchor","og"),
        na.strings = c("NA","")),
        regAnchor & !is.na(regID) & !is.na(og))
      setnames(x, c("regID", "regAnchor"), c("blkID", "blkAnchor"))
    }

    if(!paste(x$gen1[1], x$gen2[1]) %in% paste(genomeOrd$gen1, genomeOrd$gen2))
      setnames(x, c(1,2), c("ofID2","ofID1"))

    return(x[,c("ofID1","ofID2","blkID")])
  }))
  return(hitsRip)
}

#' @title calc_refChrByGene
#' @description
#' \code{calc_refChrByGene} calc_refChrByGene
#' @rdname plot_utils
#' @export
calc_refChrByGene <- function(gff,
                              refGenome,
                              blkSize,
                              nCores){
  refChr <- chr <- genome <- rl <- ofID <- NULL
  # -- if a  has a ref chr in it, populate gff
  gff2 <- data.table(gff)
  linOrd <- NULL
  setkey(gff2, linOrd)
  gff2[,refChr := chr[genome == refGenome][1], by = "og"]
  gr <- subset(gff2, !is.na(refChr))

  # -- for these seed ref chrs, ensure that runs > blkSize
  for(i in 2:blkSize){
    gr[,rl := add_rle(refChr), by = c("genome","chr")]
    gr <- subset(gr, rl >= i)
  }
  grv <- gr$refChr;  names(grv) <- gr$ofID
  gff2[,`:=`(refChr = grv[ofID])]

  # -- populate NAs by chr and genome
  spl <- split(gff2, by = c("genome","chr"))
  gff2 <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){

    # -- pull genes with no refChr specificatin
    whna <- which(is.na(x$refChr))
    if(length(whna) > 0){

      # -- pull genes without NA ref Chrs
      whnona <- which(!is.na(x$refChr))

      # -- find gene too left and right
      suppressWarnings(nal <- sapply(whna, function(y)
        whnona[which(whnona < y & whnona == max(whnona[whnona < y]))[1]]))
      suppressWarnings(nar <- sapply(whna, function(y)
        whnona[which(whnona > y & whnona == min(whnona[whnona > y]))[1]]))

      # -- fill split gff with these values
      x[,`:=`(chrl = refChr, chrr = refChr)]
      x$chrl[whna] <- x$chrl[nal]
      x$chrr[whna] <- x$chrr[nar]

      # -- fill refChr with agreeing right and left values
      whmv <- with(x, which(!is.na(chrl) & !is.na(chrr) & chrl == chrr & is.na(refChr)))
      x$refChr[whmv] <- x$chrl[whmv]
      x[,`:=`(chrl = NULL, chrr = NULL)]
    }
    return(x)
  }))

  # -- return the vector of reference chrs
  rcv <- gff2$refChr; names(rcv) <- gff2$ofID
  return(rcv)
}


#' @title order_synChrs plot
#' @description
#' \code{order_synChrs} Make a riparian plot
#' @rdname plot_utils
#' @import data.table
#' @export
order_synChrs <- function(hits,
                          minHits2plot = 50){

  ofID1 <- ofID2 <- ord1 <- ord2 <- chr1 <- chr2 <- nHits1 <- nHits2 <- NULL
  chrn1 <- meanPos1 <- chrOrder1 <- meanPos2 <- chrn2 <- chrOrder2 <- NULL
  nHits2tot <- propTot <- NULL

  h <- hits[,list(nu = uniqueN(c(ofID1, ofID2)),
                  nHits1 = uniqueN(ofID1),
                  nHits2 = uniqueN(ofID2),
                  meanPos1 = mean(ord1)),
            by = c("gen1", "gen2", "chr1", "chr2")]
  h[,`:=`(chrn1 = as.numeric(gsub("[^0-9]", "", chr1)),
          chrn2 = as.numeric(gsub("[^0-9]", "", chr2)))]
  h <- subset(h, nHits1 >= minHits2plot & nHits2 >= minHits2plot)
  setorder(h, chrn1, -nHits1, meanPos1)
  h[,chrOrder1 := as.numeric(factor(chr1, levels = unique(chr1)))]
  setorder(h, chrOrder1)
  c1 <- subset(h, !duplicated(chr1))

  h2 <- hits[,list(nHits2tot = uniqueN(ofID2),
                   meanPos2 = mean(ord2)),
             by = "chr2"]
  h2 <- merge(h, h2, by = "chr2")
  h2[,propTot := nHits2 / nHits2tot]
  setorder(h2, chr2, -propTot, meanPos2)
  c2 <- subset(h2, !duplicated(chr2))
  setorder(c2, chrOrder1, chrn2, -propTot, meanPos2)
  c2[,chrOrder2 := 1:.N]
  return(data.table(
    genome = c(c1$gen1, c2$gen2),
    chr = c(c1$chr1, c2$chr2),
    chrOrder = c(c1$chrOrder1, c2$chrOrder2)))
}

#' @title get colors for chromosome backgrounds
#' @description
#' \code{color_chrBounds}  get colors for chromosome backgrounds
#' @rdname plot_utils
#' @import data.table
#' @importFrom grDevices colorRampPalette
#' @export
color_chrBounds <- function(hits, lightChrFill, darkChrFill, emptyChrFill){
  x <- y <- n <- chr1 <- chr2 <- max1 <- percMax1 <- NULL
  setkey(hits, x)
  hits[,`:=`(chr1 = factor(chr1, levels = unique(chr1)),
             chr2 = factor(chr2, levels = unique(chr2)))]
  eg <- hits[,list(CJ(levels(chr1), levels(chr2)))]
  setnames(eg, c("chr1", "chr2"))
  chrBx <- hits[, list(x0 = min(x), x1 = max(x), x12 = mean(range(x))),
                by = c("chr1")]
  chrBy <- hits[, list(y0 = min(y), y1 = max(y), y12 = mean(range(y))),
                by = c("chr2")]
  chrBnd <- merge(chrBx, merge(chrBy, eg, by = "chr2"), by = "chr1")
  chrN <- hits[,list(n = .N),
               by = c("chr1", "chr2")]
  chrBnd <- merge(chrBnd, chrN, by = c("chr1", "chr2"), all.x = T)
  chrBnd$n[is.na(chrBnd$n)] <- 0
  chrBnd[,n := n - min(n), by = "chr1"]
  chrBnd[,max1 := max(n), by = "chr1"]
  cb0 <- subset(chrBnd, n == 0)
  cb0[,col := emptyChrFill]
  cb1 <- subset(chrBnd, n > 0)
  cb1[,percMax1 := ceiling((n / max1)*100)]
  cscl <- colorRampPalette(c(lightChrFill, darkChrFill))(100)
  cb1[,col := cscl[percMax1]]
  cb <- rbind(cb0, cb1, fill = T)[,c(colnames(chrBnd), "col"), with = F]
  return(cb)
}

#' @title draw_scaleBar
#' @description
#' \code{draw_scaleBar} draw_scaleBar
#' @rdname plot_utils
#' @importFrom graphics segments
#' @export
draw_scaleBar <- function(x, y, yspan, xspan, label, lwd, cex){
  xstart <- x - (xspan / 2)
  xend <- x + (xspan / 2)
  ytop <- y + (yspan / 2)
  ybottom <- y - (yspan / 2)
  segments(x0 = xstart, x1 = xend, y0 = y, y1 = y, lwd = lwd)
  segments(x0 = xstart, x1 = xstart, y0 = ytop, y1 = ybottom, lwd = lwd)
  segments(x0 = xend, x1 = xend, y0 = ytop, y1 = ybottom, lwd = lwd)
  text(x = xstart + (xspan / 2), y = ytop + (yspan / 2), labels = label, adj = c(.5,0), cex = cex)
}


#' @title cosine curve source data
#' @description
#' \code{cosine_points} vector of points for polygons based on cosine curves
#' @rdname plot_utils
#' @export
cosine_points <- function(){
  npts = 1e4 # initial number of points
  keepat = round(npts / 20) # grid to keep always
  grid <- seq(from = 0, to = pi, length.out = npts) # grid
  x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
  y <- grid / max(grid) # scaled grid
  # calculate slope for each point
  x1 <- x[-1];  y1 <- y[-1]
  x2 <- x[-length(x)];  y2 <- y[-length(y)]
  s <-  (y1 - y2) / (x1 - x2)
  # choose points that capture changes in slope
  ds <- cumsum(abs(diff(s)))*5
  wh <- c(1,which(!duplicated(round(ds))), length(x))
  wh2 <- c(wh, seq(from = 0, to = length(x), by = round(keepat)))
  wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
  wh <- wh[order(wh)]
  return(cbind(x[wh], y[wh]))
}

#' @title convert cosine points to polygon
#' @description
#' \code{calc_curvePolygon} from 2d coordinates, make a curve
#' @rdname plot_utils
#' @export
calc_curvePolygon <- function(start1,
                              end1 = NULL,
                              start2,
                              end2 = NULL,
                              y1,
                              y2){
  scaledCurve <- cosine_points()
  if (!is.null(end1) | !is.null(end2)) {
    tp <- rbind(
      start1 = data.table(
        x = start1, y = y1),
      poly1 = data.table(
        x = scale_between(x = scaledCurve[,1], min = start1, max = start2),
        y = scale_between(x = scaledCurve[,2], min = y1, max = y2)),
      start2 = data.table(x = start2, y = y2),
      end2 = data.table(
        x = end2, y = y2),
      poly2 = data.table(
        x = scale_between(x = scaledCurve[,1], min = end2, max = end1),
        y = scale_between(x = scaledCurve[,2], min = y2, max = y1)),
      end1 = data.table(
        x = end1, y = y1))
  }else{
    tp <- data.table(
      x = scale_between(x = scaledCurve[,1], min = start1, max = start2),
      y = scale_between(x = scaledCurve[,2], min = y1, max = y2))
  }
  return(tp)
}

#' @title calculate coordinates for rounded rectange polygons
#' @description
#' \code{round_rect} from x-y coordinates, make a rounded rectangle
#' @rdname plot_utils
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft, ybottom, xright, ytop){

  if (ytop <= ybottom)
    stop("ytop must be > ybottom")
  if (xleft >= xright)
    stop("xleft must be < xright")

  # measure graphics device
  asp <- diff(par("usr")[3:4]) / diff(par("usr")[1:2])
  dev <- dev.size()[1] / dev.size()[2]

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = 200)
  circX <- cos(theta)
  circY <- sin(theta)
  leftC <- which(circX <= 0)
  rightC <- which(circX >= 0)

  xR <- circX[rightC]
  yR <- circY[rightC]
  ordYR <- rev(order(yR))
  xR <- xR[ordYR]
  yR <- yR[ordYR]

  xL <- circX[leftC]
  yL <- circY[leftC]
  ordYL <- order(yL)
  xL <- xL[ordYL]
  yL <- yL[ordYL]

  # project onto graphics device and scale
  xRightS <- xright - (radius / asp / dev)
  xLeftS <- xleft + (radius / asp / dev)
  if (centerX < xLeftS)
    xLeftS <- centerX
  if (centerX > xRightS)
    xRightS <- centerX
  xLS <- scale_between(xL, xleft, xLeftS)
  xRS <- scale_between(xR, xRightS, xright)
  yLS <- scale_between(yL, ybottom, ytop)
  yRS <- scale_between(yR, ybottom, ytop)
  return(data.table(x = c(xRS,xLS), y = c(yRS, yLS)))
}


#' @title check if a vector is coercible to R colors
#' @description
#' \code{are_colors} check if a vector is coercible to R colors
#' @rdname plot_utils
#' @importFrom grDevices col2rgb
#' @export
are_colors <- function(col) {
  sapply(col, function(X) {
    tryCatch(is.matrix(col2rgb(X)),
             error = function(e) FALSE)
  })
}

#' @title add transparency to a color
#' @description
#' \code{add_alpha} add transparency to a color
#' @rdname plot_utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col,
                      alpha = 1){

  if (missing(col))
    stop("Please provide a vector of colors.")
  if (!all(are_colors(col)))
    stop("Please provide a vector of colors.")

  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1],
              x[2],
              x[3],
              alpha = alpha))
}

#' @title calculate coordinates for rounded rectange polygons
#' @description
#' \code{round_rect} from x-y coordinates, make a rounded rectangle
#' @rdname plot_utils
#' @import data.table
#' @export
interp_linear <- function(x,
                          y,
                          interpTails = TRUE){
  rl <- ip <- NULL
  if(length(x) != length(y) || !is.numeric(x) || !is.numeric(y)){
    warning("x and y must be numeric/integer vectors of equal length\n")
  }else{

    # -- convert to data table
    z <- subset(data.table(x = x, y = y, i = 1:length(x)), !is.na(x))
    if(nrow(z) < 1 || all(is.na(y))){
      warning("no non-missing values in x or y\n")
    }else{

      # -- subset to complete cases in x and order by x
      z <- subset(z, !is.na(x))
      setkey(z, x)

      # -- find runs of NAs in y
      z[,rl := add_rle(is.na(y), which = "id")]

      # -- pull runs to infer (not first and last if they are NAs)
      if(interpTails){
        if(is.na(z$y[1])){
          z[,rl := rl + 1]
          z <- rbind(data.table(
            x = min(z$x, na.rm = T) - .5,
            y = min(z$y, na.rm = T) - .5,
            i = 0, rl = 1),
            z)
        }
        if(is.na(z$y[nrow(z)])){
          print(z)
          z <- rbind(z, data.table(
            x = max(z$x, na.rm = T) + .5,
            y = max(z$y, na.rm = T) + .5,
            i = max(z$i, na.rm = T) + 1,
            rl = max(z$rl, na.rm = T) + 1))
        }
        toinf <- subset(z, is.na(y))
      }else{
        toinf <- subset(z, is.na(y) & !rl %in% c(1, max(rl)))
      }

      if(nrow(toinf) < 1){
        warning("no missing values of y to interpolate")
      }else{
        # -- get max right and min left values for each non-missing run
        minr <- with(subset(z, !is.na(y)), tapply(y, rl, min))
        maxl <- with(subset(z, !is.na(y)), tapply(y, rl, max))

        # -- linear interpolation of runs of NAs from bounding values
        toinf[,ip := seq(from = maxl[as.character(rl-1)],
                         to = minr[as.character(rl+1)],
                         length.out = .N+2)[-c(1, .N+2)],
              by = "rl"]

        # -- fill NAs and return
        y[toinf$i] <- toinf$ip
      }
    }
  }

  return(y)
}
