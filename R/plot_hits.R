#' @title Genespace plotting routines
#' @description
#' \code{plot_hits} Genespace plotting routines
#'
#' @param hits data.table of hits
#' @param plotType character string specifying the plot type
#' @param reorderChrs logical, should the chromosomes be re-ordered based on
#' synteny?
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
#' @param darkChrFill color of the most dense chr backgrounds
#' @param lightChrFill color of the least populated chr backgrounds
#' @param emptyChrFill color of the empty chr backgrounds
#' @param minGenes2plot integer specifying the minimum number of genes on a
#' chr to plot
#' @param cols vector of colors to use for points
#' @param returnSourceData logical, should the source data to build the plot
#' be returned?
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
#' @importFrom grDevices colorRampPalette
#' @export
plot_hits <- function(hits,
                      onlyOG = TRUE,
                      onlyAnchors = TRUE,
                      useBlks = TRUE,
                      onlyArrayReps = TRUE,
                      onlyBuffer = FALSE,
                      reorderChrs = TRUE,
                      minGenes2plot = 10,
                      gapProp = .01,
                      cols = NULL,
                      alpha = 1,
                      axisTitleCex = .6,
                      darkChrFill = "grey60",
                      lightChrFill = "grey85",
                      emptyChrFill = "grey97",
                      missingHitCol = "grey40",
                      chrLabCex = .4,
                      returnSourceData = F,
                      chrLabFun = function(x)
                        gsub("^0","",gsub("^chr|^scaffold|^lg|_","",tolower(x)))){

  ##############################################################################
  # -- ad hoc function to find and code chromosome bounds by hit number
  color_chrBounds <- function(hits, lightChrFill, darkChrFill, emptyChrFill){
    setDTthreads(1)
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


  setDTthreads(1)

  # -- subset hits to right stuff
  tp <- data.table(hits)
  if(onlyOg)
    tp <- subset(tp, isOg)
  if(onlyAnchors)
    tp <- subset(tp, isAnchor)
  if(!useBlks)
    tp[,blkID := regID]
  if(onlyBuffer)
    tp <- subset(tp, inBuffer)
  if(onlyArrayReps){
    tp <- subset(tp, isRep1 & isRep2)
  }
  tp[,ord1 := frank(ord1, ties.method = "dense")]
  tp[,ord2 := frank(ord2, ties.method = "dense")]

  tp[,nu1 := uniqueN(ofID1), by = "chr1"]
  tp[,nu2 := uniqueN(ofID2), by = "chr2"]
  tp <- subset(tp, nu1 >= minGenes2plot & nu2 >= minGenes2plot)

  cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E", "#E6AB02",
            "#A6761D", "#666666", "darkred","darkblue")
  uniqBlks <- unique(tp$blkID)
  colScale <- colorRampPalette(cols)
  cols <- sample(colScale(uniqueN(uniqBlks)))
  names(cols) <- as.character(uniqBlks)
  tp[,col := cols[as.character(blkID)]]
  tp$col[is.na(tp$col)] <- missingHitCol

  plotTitle <- sprintf(
    "%s: %s",
    ifelse(onlyOG & onlyArrayReps, "only array rep OGs",
           ifelse(onlyOG, "only OGs",
                  ifelse(onlyArrayReps, "only array reps","all hits"))),
    ifelse(onlyAnchors, "constrained to anchors",
           ifelse(onlyBuffer, "in synteny buffer", "not synteny constrained")))

  # -- order chromosomes and add linear x/y positions
  gps <- max(unique(tp[,c("ord1", "ord2")]))
  setkey(tp, ord1)
  tp[,n := (as.numeric(factor(chr1, levels = unique(chr1)))-1)]
  tp[,x := ord1 + (gps * gapProp * n)]
  setkey(tp, ord2)
  tp[,n := (as.numeric(factor(chr2, levels = unique(chr2)))-1)]
  tp[,y := ord2 + (gps * gapProp * n)]


  tp <- subset(tp, complete.cases(tp[,c("chr1", "chr2", "x", "y")]))

    # calculate chr bounds / membership
  cb <- color_chrBounds(
    hits = tp, lightChrFill = lightChrFill,
    darkChrFill = darkChrFill, emptyChrFill = emptyChrFill)

  # make plot window
  # par(mar = c(2,2,1,1))
  xoffset <- min(tp$x) - (diff(range(tp$x))/20)
  yoffset <- min(tp$y) - (diff(range(tp$y))/20)

  plot(
    NA, NA,
    xlim = c(xoffset, max(tp$x)),
    ylim = c(yoffset, max(tp$y)),
    type = "n", axes = F,
    xlab = "", ylab = "",
    asp = 1)

  title(
    xlab = paste(hits$gen1[1], "chromosomes (gene rank order)"),
    ylab = paste(hits$gen2[1], "chromosomes (gene rank order)"),
    line = 0, cex.lab = axisTitleCex,
    main = plotTitle)

  # plot chrs and label
  x1 <- x0 <- y0 <- y1 <- col <- x12 <- chr1<- chr2 <- y12 <- NULL
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
  col <- n <- x <- y <- NULL
  setkey(tp, col, n)
  if(nrow(tp) > 20e3){
    with(tp, points(x, y, col = "white", pch = 16, cex = .2))
    with(tp, points(x, y, col = add_alpha(col, alpha = alpha), pch = "."))
  }else{
    with(tp, points(x, y, col = "white", pch = 16, cex = .4))
    with(tp, points(x, y, col = add_alpha(col, alpha = alpha), pch = 16, cex = .25))
  }
  if(returnSourceData)
    return(tp)
}
