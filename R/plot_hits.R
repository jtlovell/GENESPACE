#' @title Genespace plotting routines
#' @description
#' \code{plot_hits} Genespace plotting routines
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param hits data.table of hits
#' @param plotRegions logical, should regions be plotted (instead of blocks)?
#' @param bufferOnly logical, should only hits in buffers be plotted?
#' @param anchorOnly logical, should only hits that are anchors be plotted?
#' @param noBlkCol color for points that are not in blocks
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
#' @param fixedAspRat logical, should the plots have fixed aspect ratios?
#' @param plotTitle character string specifying the title of the plot
#' @param minHits2plot integer, the minimum number of hits to include
#' @details ...
#'
#' @import data.table
#' @importFrom graphics title
#' @export
plot_hits <- function(hits,
                      gff,
                      plotRegions = FALSE,
                      bufferOnly = FALSE,
                      anchorOnly = FALSE,
                      noBlkCol = add_alpha("gold",.2),
                      gsParam,
                      minGenes2plot = 50,
                      onlyOg = TRUE,
                      gapProp = .005,
                      useOrder = T,
                      round2 = ifelse(useOrder, 10, 5e4),
                      cols = NULL,
                      alpha = 1,
                      fixedAspRat = TRUE,
                      axisTitleCex = .6,
                      darkChrFill = "grey60",
                      lightChrFill = "grey85",
                      emptyChrFill = "grey97",
                      chrLabCex = .4,
                      chrLabFun = function(x)
                        gsub("^0","",gsub("^chr|^scaffold|^lg|_","",tolower(x))),
                      plotTitle = "diamond hit dotplot"){
  ofID1 <- ofID2 <- x <- y <- isOg <- n <- genome <- gen2 <- ref <- chr1 <- NULL
  regBuffer <- blkBuffer <- regAnchor <- blkAnchor <- og <- colGrp <- NULL
  blkID <- regID <- chr2 <- bSize <- NULL

  if(is.null(cols) || length(cols) != 1)
    cols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E", "#E6AB02",
              "#A6761D", "#666666", "darkred","darkblue")

  # -- subset hits to right stuff
  refh <- data.table(hits)
  if(bufferOnly){
    if(plotRegions){
      refh <- subset(refh, regBuffer)
    }else{
      refh <- subset(refh, blkBuffer)
    }
  }
  if(anchorOnly){
    if(plotRegions){
      refh <- subset(refh, regAnchor)
    }else{
      refh <- subset(refh, blkAnchor)
    }
  }
  if(onlyOg)
    refh <- subset(refh, !is.na(og))
  refh[,colGrp := 1]
  if(plotRegions)
    refh[,colGrp := regID]
  if(!plotRegions)
    refh[,colGrp := blkID]

  # -- check if self hits and, if so, rename paralogs if required
  g1 <- refh$gen1[1]; g2 <- refh$gen2[1]
  gf <- subset(gff, genome %in% c(g1, g2))
  if(g1 == g2){
    pselfAnch <- with(subset(refh, blkAnchor), mean(ofID1 == ofID2))
    if(pselfAnch > 0.05){
      gf1 <- data.table(gf)
      gf1[,genome := sprintf("%s_para", g1)]
      g2 <- sprintf("%s_para", g1)
      gf <- rbind(gf, gf1)
      refh[,gen2 := g2]
    }
  }

  # -- reorder gff and add new coordinates
  gf <- reorder_gff(
    gf,
    refGenome = g1,
    minGenesOnChr = minGenes2plot,
    genomeIDs = c(g1, g2))
  gv <- gf$genome; cv <- gf$chr; ov <- gf$ord; sv <- gf$start; ev <- gf$end
  names(gv) <- names(cv) <- names(ov) <- names(sv) <- names(ev) <- gf$ofID
  refh[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2],
             ord1 = ov[ofID1], ord2 = ov[ofID2])]
  refh <- subset(refh, complete.cases(refh[,c("gen1","gen2","chr1","chr2","ord1","ord2")]))

  # add syntenic orthogroup info
  if(!"synOg" %in% colnames(gff))
    gf <- add_synOg2gff(
      hits = refh,
      gff = gf,
      gsParam = gsParam,
      genomeIDs = c(g1, g2),
      useBlks = !plotRegions,
      allowRBHinOg = T)


  # -- get syntenic chrs and linear positions
  gf <- add_synChr2gff(
    gff = gf,
    refHits = subset(refh, !is.na(og)),
    refGenome = g1,
    genomeIDs = c(g1, g2),
    gapProp = gapProp)
  ov <- gf$linOrd; sv <- gf$linBp;  names(ov) <- names(sv) <-  gf$ofID

  # -- add linear positions to hits
  if(useOrder){
    refh[,`:=`(x = ov[ofID1], y = ov[ofID2])]
  }else{
    refh[,`:=`(x = sv[ofID1], y = sv[ofID2])]
  }

  # -- round to nearest position
  if(round2 > 0){
    refh[,`:=`(x = round_toInteger(x, round2), y = round_toInteger(y, round2))]
    tp <- refh[,list(n = .N),
                 by = c("chr1","chr2","x","y","colGrp")]
    tp <- subset(tp, n > 0)
    setkey(tp, n)
  }else{
    tp <- ref[,c("chr1","chr2","x","y","colGrp")]
  }
  tp <- subset(tp, !duplicated(paste(x, y)))
  # add colors
  if(length(cols) == 1){
    tp[,col := add_alpha(cols, alpha = alpha)]
  }else{
    if(any(is.na(tp$colGrp)))
      tp$colGrp[is.na(tp$colGrp)] <- "XXXX_noReg_XXXX"
    colScale <- colorRampPalette(cols)
    colPal <- sample(colScale(uniqueN(tp$colGrp)))
    names(colPal) <- unique(tp$colGrp)
    tp[,col := add_alpha(colPal[colGrp], alpha)]
    tp$col[tp$colGrp == "XXXX_noReg_XXXX"] <- noBlkCol
  }

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
      xlab = paste(g1,"chromosomes (gene rank order)"),
      ylab = paste(g2,"chromosomes (gene rank order)"),
      line = 0, cex.lab = axisTitleCex,
      main = plotTitle)
  }else{
    title(
      xlab = paste(g1,"chromosomes (physical gene position)"),
      ylab = paste(g2,"chromosomes (physical gene position)"),
      line = 0, cex.lab = axisTitleCex,
      main = plotTitle)
  }

  # plot chrs and label
  cb <- color_chrBounds(
    hits = subset(refh, !is.na(og)), lightChrFill = lightChrFill,
    darkChrFill = darkChrFill, emptyChrFill = emptyChrFill)

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
  if(!any(c(plotRegions, bufferOnly, anchorOnly))){
    setorder(tp, colGrp)
    with(tp, points(x, y, col = col, pch = "."))
  }else{
    tp[,bSize := .N, by = "colGrp"]
    setorder(tp, bSize)
    with(tp, points(x, y, col = "white", pch = 16, cex= .4))
    with(tp, points(x, y, col = col, pch = 16, cex= .25))
  }
  return(tp)
}

