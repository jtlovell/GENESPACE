#' @title Genespace plotting routines
#' @description
#' \code{plot_hits} Genespace plotting routines
#'
#' @param hits data.table of hits
#' @param plotRegions logical, should regions be plotted (instead of blocks)?
#' @param bufferOnly logical, should only hits in buffers be plotted?
#' @param anchorOnly logical, should only hits that are anchors be plotted?
#' @param round2 integer, specifying the rounding of gene rank order positions
#' to reduce the total number of points. 1 = don't round
#' @param alpha numeric (0-1) specifying transparency of the points
#' @param chrLabFun function to parse chr IDs to make them more readible
#' @param gapProp numeric (0-1) specifying the proportional size of gaps
#' relative to the length of the largest genome
#' @param useOrder logical, should gene rank order be used in lieu of physical
#' positions?
#' @param axisTitleCex character expansion for the axes and genome labels
#' @param chrLabCex character expansion for the chromosome labels
#' @param gff annotated gff-like data.table
#' @param darkChrFill color of the most dense chr backgrounds
#' @param lightChrFill color of the least populated chr backgrounds
#' @param emptyChrFill color of the empty chr backgrounds
#' @param minGenes2plot integer specifying the minimum number of genes on a
#' chr to plot
#' @param onlyOg logical, should only og hits be plotted?
#' @param cols vector of colors to use for points
#' @param plotTitle character string specifying the title of the plot
#'
#' @details ...
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @import data.table
#' @importFrom graphics title
#' @export
plot_hits <- function(hits,
                      plotType = c("allHits", "allOG", "regAnchor", "regBuffer", "blkAnchor", "blkBuffer"),
                      reorderChrs = TRUE,
                      minGenes2plot = 10,
                      gapProp = .005,
                      useOrder = TRUE,
                      round2 = ifelse(useOrder, 10, 5e4),
                      cols = NULL,
                      alpha = 1,
                      axisTitleCex = .6,
                      darkChrFill = "grey60",
                      lightChrFill = "grey85",
                      emptyChrFill = "grey97",
                      chrLabCex = .4,
                      returnSourceData = F,
                      chrLabFun = function(x)
                        gsub("^0","",gsub("^chr|^scaffold|^lg|_","",tolower(x)))){

  ##############################################################################
  # -- ad hoc function to find and code chromosome bounds by hit number
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

  ##############################################################################
  # -- ad hoc function to add linear coordinates to the hits
  add_linCoords2hits <- function(hits,
                                 gapProp,
                                 useOrder,
                                 reorderChrs){
    h <- data.table(hits)
    ho <- subset(h, isOg & blkAnchor)
    cov1 <- rank(with(ho, tapply(ord1, chr1, median)))
    if(reorderChrs){
      cov2 <- rank(with(ho, tapply(ord1, chr2, median)))
    }else{
      cov2 <- rank(with(ho, tapply(ord2, chr2, median)))
    }
    h[,`:=`(chrOrd1 = cov1[chr1],
            chrOrd2 = cov2[chr2])]

    if(useOrder){
      gap <- gapProp * max(c(h$ord1, h$ord2), na.rm = T)
      h[,`:=`(sclGap1 = gap * (chrOrd1 - 1),
              sclGap2 = gap * (chrOrd2 - 1))]
      setkey(h, chrOrd1, ord1)
      h[,x := 1:.N + sclGap1]
      setkey(h, chrOrd2, ord2)
      h[,y := 1:.N + sclGap2]
    }else{
      gap <- gapProp * max(c(h$end1, h$end2), na.rm = T)
      mxv1 <- cumsum(with(h, tapply(end1, chrOrd1, max))) + gap
      mxv2 <- cumsum(with(h, tapply(end2, chrOrd2, max))) + gap
      mxv1 <- c(0, mxv1[1:(length(mxv1) - 1)])
      names(mxv1) <- 1:length(mxv1)
      mxv2 <- c(0, mxv2[1:(length(mxv2) - 1)])
      names(mxv2) <- 1:length(mxv2)
      h[,`:=`(sclGap1 = mxv1[as.character(chrOrd1)],
              sclGap2 = mxv2[as.character(chrOrd2)])]
      h[,x := start1 + sclGap1]
      h[,y := start2 + sclGap2]
    }
    return(h)
  }

  ##############################################################################
  # -- ad hoc function to condense hits rounded to positions and count them
  condense_hits <- function(hits, round2){
    h <- data.table(hits)
    if(round2 > 0){
      h[,`:=`(x = round_toInteger(x, round2),
               y = round_toInteger(y, round2))]
      h <- h[,list(n = .N),
               by = c("chr1","chr2","x","y","col")]
      h <- subset(h, n > 0)
      h[,n := frank(n, ties.method = "dense")]
      h[,n := frank(round_toInteger(n, 5), ties.method = "dense")]
      h[,n := scale_between(n, min = .25, max = 1)]
      col <- n <- NULL
      h[,col := sapply(1:nrow(h), function(i)
        add_alpha(h$col[i], alpha = h$n[i]))]
      setkey(h, n)
    }else{
      h <- h[,c("chr1","chr2","x","y","col")]
    }
    h <- subset(h, !duplicated(paste(x, y)))
    return(h)
  }

  ofID1 <- ofID2 <- x <- y <- isOg <- n <- genome <- gen2 <- ref <- chr1 <- NULL
  regBuffer <- blkBuffer <- regAnchor <- blkAnchor <- og <- colGrp <- NULL
  blkID <- regID <- chr2 <- bSize <- NULL

  # -- subset hits to right stuff
  pltTypes <- c("allHits", "allOG", "regAnchor", "regBuffer", "blkAnchor", "blkBuffer")
  if(!plotType %in% pltTypes)
    plotType <- "allOG"
  p <- plotType[1]
  if(p == "allOG"){
    tp <- subset(hits, isOg)
  }else{
    tp <- data.table(hits)
  }
  plotTitle <- ifelse(
    p == "allHits", "all hits",
    ifelse(p == "allOG", "Orthogroup constained hits",
           ifelse(p == "blkAnchor", "syntenic block anchors",
                  ifelse(p == "regAnchor", "syntenic region anchors",
                         ifelse(p == "regBuffer", "syntenic region buffers",
                                "syntenic block buffers")))))

  if(p %in% c("allHits", "allOG") & length(cols) == 0){
    cols <- "blue2"
  }else{
    if(p %in% c("allHits", "allOG") & length(cols) >= 1){
      cols <- cols[1]
    }else{
      if(is.null(cols) || length(cols) != 1)
        cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E", "#E6AB02",
                  "#A6761D", "#666666", "darkred","darkblue")
    }
  }

  if(p == "regAnchor") tp <- subset(tp, regAnchor)
  if(p == "regBuffer") tp <- subset(tp, regBuffer)
  if(p == "blkAnchor") tp <- subset(tp, blkAnchor)
  if(p == "blkBuffer") tp <- subset(tp, blkBuffer)
  tp[,colGrp := 1]
  if(p %in% c("regAnchor", "regBuffer")) tp[,colGrp := regID]
  if(p %in% c("blkAnchor", "blkBuffer")) tp[,colGrp := blkID]

  tp[,nu1 := uniqueN(ofID1), by = "chr1"]
  tp[,nu2 := uniqueN(ofID2), by = "chr2"]
  tp <- subset(tp, nu1 >= minGenes2plot & nu2 >= minGenes2plot)

  # -- color hits
  if(length(cols) == 1){
    tp[,col := cols]
  }else{
    colScale <- colorRampPalette(cols)
    cols <- sample(colScale(uniqueN(tp$colGrp)))
    names(cols) <- as.character(unique(tp$colGrp))
    tp[,col := cols[as.character(colGrp)]]
  }

  # -- order chromosomes and add linear x/y positions
  tp <- add_linCoords2hits(
    hits = tp, gapProp = gapProp,
    useOrder = useOrder, reorderChrs = reorderChrs)

  # -- round to nearest position
  tp <- condense_hits(hits = tp, round2 = round2)
  tp <- subset(tp, complete.cases(tp))
  # calculate chr bounds / membership
  cb <- color_chrBounds(
    hits = tp, lightChrFill = lightChrFill,
    darkChrFill = darkChrFill, emptyChrFill = emptyChrFill)

  # make plot window
  par(mar = c(2,2,1,1))
  xoffset <- min(tp$x) - (diff(range(tp$x))/20)
  yoffset <- min(tp$y) - (diff(range(tp$y))/20)
  plot(
    NA, NA,
    xlim = c(xoffset, max(tp$x)),
    ylim = c(yoffset, max(tp$y)),
    type = "n", axes = F,
    xlab = "", ylab = "",
    asp = 1)

  if(useOrder){
    title(
      xlab = paste(hits$gen1[1], "chromosomes (gene rank order)"),
      ylab = paste(hits$gen2[1], "chromosomes (gene rank order)"),
      line = 0, cex.lab = axisTitleCex,
      main = plotTitle)
  }else{
    title(
      xlab = paste(hits$gen1[1], "chromosomes (physical gene position)"),
      ylab = paste(hits$gen2[1], "chromosomes (physical gene position)"),
      line = 0, cex.lab = axisTitleCex,
      main = plotTitle)
  }

  # plot chrs and label
  with(cb, rect(
    xleft = x0, xright = x1,
    ybottom = y0, ytop = y1,
    col = col, border = NA))

  with(subset(cb, !duplicated(chr1)),
       text(
         x = x12, y = yoffset*.925, labels = chrLabFun(as.character(chr1)),
         cex = chrLabCex, adj = c(1.05,.5), srt = 90))
  with(subset(cb, !duplicated(chr2)),
       text(
         x = xoffset*.925, cex = chrLabCex, labels = chrLabFun(as.character(chr2)),
         y = y12,  adj = c(1.05,.5)))

  # -- plot points
  setkey(tp, col, n)
  if(uniqueN(cols) > 1){
    with(tp, points(x, y, col = "white", pch = 16, cex= .4))
    with(tp, points(x, y, col = col, pch = 16, cex= .25))
  }else{
    if(nrow(tp) > 10e3){
      with(tp, points(x, y, col = col, pch = "."))
    }else{
      with(tp, points(x, y, col = col, pch = 16, cex= .25))
    }
  }
  if(returnSourceData)
    return(tp)
}
