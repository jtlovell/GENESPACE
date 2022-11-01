#' @title Make riparian plot using hits, not OGs
#' @description
#' \code{plot_riparian} Updated and more accurate version of plot_riparian
#'
#' @name plot_riparian
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param useBlks logical, should blocks be plot (and not regions)?
#' @param labelTheseGenomes character string of genomes that have labeled chrs
#' @param invertTheseChrs data.table with two columns, genome and chr containing
#' the lists of genomes and chromosomes that should be inverted in the plot
#' @param genomeLabCex chracter expansion of the genome labels
#' @param braidAlpha numeric (0-1) specifying the transparency of the braids
#' @param braidBorderLwd numeric specifying the weight of borders on the braids
#' @param genomeIDs character vector at least partially matching the genomeIDs
#' in gsParam
#' @param blackBg logical, should the background be dark?
#' @param refGenome single character string specifying which genome is the ref
#' @param reorderChrs logical, should chromosomes be re-ordered by synteny?
#' @param gapProp numeric (0-1) specifying the proportional size of gaps
#' relative to the length of the largest genome
#' @param useOrder logical, should gene rank order be used in lieu of physical
#' positions?
#' @param chrLabCex character expansion for the chromosome labels
#' @param chrBorder character or integer coercible to an R color for the border
#' of the chromosome regions
#' @param chrFill character or integer coercible to an R color for the interior
#' fill of the chromosome regions
#' @param refChrCols if braids should be colored by reference chromosome,
#' a vector of colors or integers that is coerced into a color palette. If
#' no coloring by chr, specify NULL.
#' @param labelChrBiggerThan integer specifying the minimum number of genes or
#' bp to label a chromosome
#' @param highlightRef color to highlight the reference chromosomes.
#' @param chrLabFun function to parse chr IDs to make them more readible

#'
#' @details ...
#'
#'
#'
#' @examples
#' \dontrun{
#' # see vignette dedicated to plot_riparian
#' }
#'
#' @title plot_riparian
#' @description
#' \code{plot_riparian} plot_riparian
#' @rdname plot_riparian
#' @export
#' @import data.table
#' @import ggplot2
#' @export
plot_riparian <- function(
    gsParam,
    reorderChrsBySynteny = TRUE,
    refGenome = NULL,
    labelTheseGenomes = NULL,
    genomeIDs = NULL,
    useOrder = TRUE,
    chrLengths = NULL,
    faidir = NULL,
    blockCoords = NULL,
    gapProp = 0.01,
    minBlkSize = 10,
    pdfFile = NULL,
    chrFill = "white",
    highlightTheseRegions = NULL,
    braidAlpha = .75,
    scalePlotHeight = 1,
    scalePlotWidth = 0.3,
    chrLabBorderLwd = .1,
    refChrCols = NULL,
    ylabBuff = NULL,
    refChrOrdFun = function(x)
      frank(list(as.numeric(gsub('\\D+','', x)), x), ties.method = "random"),
    chrLabFun = function(x)
      gsub("^0", "",
           gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x))),
    xlab = sprintf(
      "Chromosomes scaled by %s",
      ifelse(useOrder, "gene rank order", "physical position")),
    verbose = TRUE){


  check_highlightRegParam <- function(x){
    if(!is.null(x)){
      # if there are regions to highlight, check and make sure the data look ok
      if(!is.data.frame(x))
        stop("regions given to highlight, but not a data.table/data.frame\n")
      if(!all(c("genome", "chr", "start", "end") %in% colnames(x)))
        stop("highlightTheseRegions must be a data.table or data.frame with columns: geome, chr, start, end\n")
      x <- subset(x, complete.cases(x))
      if(nrow(x) < 1)
        stop("no complete (non-NA) rows in highlightTheseRegions\n")

      if("id" %in% colnames(x)){
        if(any(duplicated(x$id))){
          warning("found duplicate region ids, renaming")
          x[,id := NULL]
        }
      }
      if(!"id" %in% colnames(x))
        x[,id := sprintf("reg_%s", 1:.N)]

      if("col" %in% colnames(x)){
        if(!all(are_colors(x$col))){
          warning("some colors in highlightTheseRegions are not valid, re-coloring\n")
          x[,col := NULL]
        }
      }
      if(!"col" %in% colnames(x))
        x[,col := gs_colors(n = .N)]
    }
    return(x)
  }

  check_blockCoordsParam <- function(x){
    blkCols <- c(
      "genome1", "genome2", "chr1", "chr2", "start1", "end1", "start2",
      "end2", "blkID", "color")
    if(!is.null(x)){
      if(!is.data.frame(x))
        x <- NULL
      if(!all(blkCols %in% colnames(x)))
        x <- NULL
    }
    if(!is.null(x)){
      if(!"nhits1" %in% colnames(x))
        x[,nhits1 := Inf]
      if(!"nhits2" %in% colnames(x))
        x[,nhits2 := Inf]

      if("color" %in% colnames(x)){
        blk[,color := check_character(color)]
        blk$color[is.na(blk$color)] <- "xxx"
        if(!all(are_colors(x$color))){
          x[,color := NULL]
          warning("some values in blockCoords color column could not be converted to character (HEX or R named) color values. Coloring by genome1 chromosomes.")
        }
      }

      if(!"color" %in% colnames(x)){
        u <- unique(c(x$chr1, x$chr2))
        cols <- gs_colors(n = length(u))
        names(cols) <- as.character(u)
        x[,color := cols[as.character(chr1)]]
      }
    }
    return(x)
  }

  check_labelTheseGenomesParam <- function(x, genomeIDs){
    ltg <- x
    if(!is.null(x))
      x <- x[x %in% genomeIDs]
    if(length(x) == 0)
      x <- NULL
    if(is.null(x))
      x <- genomeIDs
    if(!is.null(ltg))
      if(!identical(ltg, x))
        warning("some labelTheseGenomes are not in genomeIDs.")
    return(x)
  }

  set_blkColors <- function(blk, runType, bg, cols, alpha, refChrOrdFun){
    if(runType == "highlight"){
      if(is.null(cols) || is.na(cols) || !are_colors(cols[1])){
        blk[,color := bg]
      }else{
        blk[,color := cols[1]]
      }
    }

    # if default get this figured out
    if(runType == "default"){
      u <- unique(blk$refChr)
      uord <- refChrOrdFun(u)
      names(uord) <- u
      uord <- uord[order(uord)]
      u <- names(uord)
      if(is.null(cols))
        cols <- gs_colors(length(u))
      if(length(cols) == 1)
        cols <- rep(cols, length(u))
      # if there is still a problem with the colors (or null) remake using gspal
      if(length(cols) != length(u) && !is.null(cols))
        cols <- colorRampPalette(cols)(length(u))
      names(cols) <- u
      blk[,color := cols[refChr]]
    }
    return(blk)
  }

  subset_toChainedGenomes <- function(genomeIDs,
                                      blk){
    genomeOrd <- data.table(
      genome1 = genomeIDs[-length(genomeIDs)], genome2 = genomeIDs[-1])
    genomeOrd[,`:=`(y1 = match(genome1, genomeIDs),
                    y2 = match(genome2, genomeIDs))]
    u <- with(genomeOrd, paste(genome1, genome2))
    if(!all(u %in% paste(blk$genome1, blk$genome2)))
      stop("problem with the block coordinates, some chained combinations of genomeIDs are not in the blocks\n")
    blk <- merge(genomeOrd, blk, by = c("genome1", "genome2"))

    # drop any rows with missing data
    blk <- subset(blk, complete.cases(blk))
    return(blk)
  }

  simplify_blkCoords <- function(blk, useOrder, runType){
    if(runType %in% c("highlight", "default")){
      if(useOrder){
        tpb <- with(blk, data.table(
          genome1 = genome1, genome2 = genome2, chr1 = chr1, chr2 = chr2,
          y1 = y1, y2 = y2, start1 = startOrd1, start2 = startOrd2,
          end1 = endOrd1, end2 = endOrd2, color = color, blkID = blkID))
        tpb[,minStart1 := min(c(start1, end1)) - 1, by = c("genome1", "chr1")]
        tpb[,minStart2 := min(c(start2, end2)) - 1, by = c("genome2", "chr2")]
        tpb[,`:=`(start1 = start1 - minStart1, end1 = end1 - minStart1,
                  start2 = start2 - minStart2, end2 = end2 - minStart2)]
        tpb[,`:=`(minStart1 = NULL, minStart2 = NULL)]

      }else{
        tpb <- with(blk, data.table(
          genome1 = genome1, genome2 = genome2, chr1 = chr1, chr2 = chr2,
          y1 = y1, y2 = y2, start1 = startBp1, start2 = startBp2,
          end1 = endBp1, end2 = endBp2, color = color, blkID = blkID))
      }
    }
    return(tpb)
  }

  get_chrLens <- function(chrLengths, faidir, blk){
    u <- unique(unlist(blk[,c("genome1", "genome2")]))
    if(is.data.frame(chrLengths)){
      if(!all(c("genome", "chr", "length") %in% colnames(chrLengths))){
        chrLengths <- NULL
      }else{
        if(!all(u %in% chrLengths$genome)){
          chrLengths <- NULL
        }
      }
    }else{
      chrLengths <- NULL
    }

    if(!is.null(faidir)){
      if(!dir.exists(faidir)){
        faif <- file.path(faidir, sprintf("%s.fa.fai", u))
        names(faif) <- u

        if(!all(file.exists(faif))){
          warning("cannot find all fa.fai files in faidir\n")
        }else{
          chrLengths <- rbindlist(names(faif), function(i){
            tmp <- fread(
              faif[i], select = c(1, 2), col.names = c("chr", "length"))
            tmp[,genome := i]
            return(tmp)
          })
          chrLengths <- chrLengths[,c("genome", "chr", "length"), with = F]
        }
      }
    }

    if(is.null(chrLengths)){
      if(verbose)
        cat("Cannot find valid fais or chrLengths using chromosome coordinates from the blocks.\n")
      chrLengths <- with(blk, data.table(
        genome = c(genome1, genome2, genome1, genome2),
        chr = c(chr1, chr2, chr1, chr2),
        end = c(start1, start2, end1, end2)))
      chrLengths <- chrLengths[,list(length = max(end)), by = c("genome", "chr")]
    }

    u <- with(blk, unique(paste(c(genome1, genome2), c(chr1, chr2))))
    if(!all(u %in% paste(chrLengths$genome, chrLengths$chr))){
      warning("could not parse chrLengths, some genome/chr combinations in block coordinates are not in chrLengths\n")
      chrLengths <- NULL
    }
    cout <- chrLengths$length
    names(cout) <- with(chrLengths, paste(genome, chr))
    return(chrLengths)
  }

  make_scaleBar <- function(chrLengths, ylabBuff){
    s <- min(with(chrLengths, tapply(length, genome, sum))) * .3
    if(useOrder){
      s <- ifelse(s > 1e5, round(s, -5),
                  ifelse(s > 1e4, round(s, -4),
                         ifelse(s > 1e3, round(s, -3),
                                ifelse(s > 100, round(s, -2), round(s, -1)))))
    }else{
      s <- ifelse(s > 1e9, round(s, -9),
                  ifelse(s > 1e8, round(s, -8),
                         ifelse(s > 1e7, round(s, -7),
                                ifelse(s > 1e6, round(s, -6), round(s, -5)))))
    }

    top <- max(chrLengths$y) + (ylabBuff * 4)
    bot <- max(top) - ylabBuff
    left <- max(with(chrLengths, tapply(chrstart, genome, min)))
    right <- left + s
    mid <- (bot + top)/2
    scbar <- data.table(
      line = c("left", "right","mid"),
      x = c(left, right, left),
      xend = c(left, right, right),
      y = c(bot, bot, mid),
      yend = c(top, top, mid))

    sclab <- sprintf(
      "%s %s", ifelse(useOrder, s, s/1e6), ifelse(useOrder, "genes", "Mbp"))

    return(list(label = sclab, segments = scbar))
  }

  ##############################################################################
  ##############################################################################
  # 1. Get the input data together
  # ref genome and genomeIDs
  rg <- refGenome
  gids <- genomeIDs
  if(is.null(genomeIDs))
    gids <- gsParam$genomeIDs
  if(is.null(rg))
    rg <- gids[1]

  # plot labels
  labgen <- check_labelTheseGenomesParam(labelTheseGenomes, genomeIDs = gids)

  # plot defaults

  if(is.null(ylabBuff))
    ylabBuff <- length(gids)/40
  plotHeight <- (sqrt(length(gids)) * scalePlotHeight) + 2

  if(!is.null(pdfFile))
    if(!dir.exists(dirname(pdfFile)))
      stop("if pdfFile is given, must have a valid parent directory\n")
  if(is.null(pdfFile))
    pdfFile <- file.path(gsParam$paths$riparian, sprintf(
      "%s_%s.riparian.pdf", rg, ifelse(useOrder, "geneOrder", "bp")))

  ##############################################################################
  # -- 1.1 determine runtype from parameters
  # -- make sure the highlighted regions conform, or set to NULL
  hlreg <- check_highlightRegParam(highlightTheseRegions)

  # -- make sure the blockCoords conform, or set to null
  blk <- check_blockCoordsParam(blockCoords)

  if(is.null(gsParam) && is.null(blk))
    stop("either gsParam or blockCoords needs to be specified\n")

  runType <- ifelse(!is.null(hlreg), "highlight",
                    ifelse(!is.null(blk), "custom", "default"))

  ##############################################################################
  # -- 1.2 read in the blocks
  # phased blocks if default
  if(runType == "default"){
    blk <- subset(
      fread(file.path(gsParam$paths$riparian, "refPhasedBlkCoords.txt"),
            na.strings = c("", "NA")), refGenome == rg)
    blk[,blkID := paste(blkID, refChr)]
  }

  # raw blocks if highlighting
  if(runType == "highlight"){
    blk <- fread(file.path(gsParam$paths$results, "blkCoords.txt"),
                 na.strings = c("", "NA"))
  }
  blk[,blkID := paste(genome1, genome2, blkID)]

  # subset to just daisy chained hits
  blk <- subset(blk, nHits1 >= minBlkSize & nHits2 >= minBlkSize)
  tp <- subset_toChainedGenomes(blk = blk, genomeIDs = gids)
  if(verbose)
    cat(sprintf("riparian plot type: %s, with %s collinear blocks\n",
                runType, nrow(tp)))

  if(gsParam$ploidy[rg] > 1 & verbose)
    cat(strwrap(sprintf(
      "**NOTE** The reference genome %s has ploidy > 1. If possible, try to use
      a haploid reference genome. Polyploid references will produce many
      overlapping blocks and increase the file size.", rg),
      indent = 0, exdent = 8), sep = "\n")

  ##############################################################################
  # -- 1.3 add colors to the blocks
  # if highlight, the background is a single color
  tp <- set_blkColors(
    blk = tp, runType = runType, alpha = braidAlpha,
    bg = "grey25", cols = refChrCols, refChrOrdFun = refChrOrdFun)

  ##############################################################################
  # -- 1.4 simplify the block coordinates, depending on run type
  tp <- simplify_blkCoords(blk = tp, useOrder = useOrder, runType = runType)

  ##############################################################################
  # -- 1.5 check chromosome length information
  clens <- get_chrLens(chrLengths = chrLengths, faidir = faidir, blk = tp)

  ##############################################################################
  # -- 1.6 if highlight these regions, add these in
  if(runType == "highlight"){
    print(highlightTheseRegions)
    highlightTheseRegions <<- highlightTheseRegions
    blkReg <- with(highlightTheseRegions, calc_regBlkCoords(
      gsParam = gsParam,
      regGenome = genome,
      regChr = chr,
      blk = tp,
      regStart = start,
      regEnd = end,
      regCol = color,
      regID = id,
      minBlkSize = minBlkSize))
    # blkReg <- subset_toChainedGenomes(blk = blkReg, genomeIDs = genomeIDs)
    blkReg <- subset_toChainedGenomes(blk = blkReg, genomeIDs = gids)
    blkReg <- simplify_blkCoords(
      blk = blkReg, useOrder = useOrder, runType = runType)
    tp <- rbind(tp, blkReg[,colnames(tp), with = F])
  }

  ##############################################################################
  ##############################################################################
  # 2. get coordinate system together

  ##############################################################################
  # -- 2.1 get chromosome orders relative to synteny with reference
  u <- with(tp, unique(paste(c(genome1, genome2), c(chr1, chr2))))
  blk <- subset(blk, paste(genome1, chr1) %in% u & paste(genome2, chr2) %in% u)

  chrOrd <- get_chrOrder(
    blk = blk, refGenome = rg, refChrOrdFun = refChrOrdFun,
    reorderChrsBySynteny = reorderChrsBySynteny)

  ##############################################################################
  # -- 2.2 re-calculate positions relative to chromosome order
  clens[,ord := chrOrd[paste(genome, chr)]]
  setkey(clens, genome, ord)

  gapSize <- round(gapProp * sum(clens$length[clens$genome == rg]))
  clens[,gap := c(0, rep(gapSize, (.N)-1)), by = "genome"]
  clens[,chrstart := cumsum(c(1, (length[-.N]+1)) + gap), by = "genome"]
  clens[,med := sum(length + gap)/2, by = "genome"]
  clens[,`:=`(chrend = chrstart + length, gap = NULL)]
  clens[,`:=`(chrstart = chrstart - med, chrend = chrend - med)]
  clens[,med := NULL]
  ##############################################################################
  # -- 2.3 project chr start positions onto block coordintates
  a2 <- clens$chrstart - 1
  names(a2) <- with(clens, paste(genome, chr))

  tp[,`:=`(xs1 = a2[paste(genome1, chr1)] + start1,
           xe1 = a2[paste(genome1, chr1)] + end1,
           xs2 = a2[paste(genome2, chr2)] + start2,
           xe2 = a2[paste(genome2, chr2)] + end2)]
  tpb <- data.table(tp)
  u <- with(tp, unique(paste(c(genome1, genome2), c(chr1, chr2))))
  clens <- subset(clens, paste(genome, chr) %in% u )
  clens[,y := match(genome, gids)]
  ##############################################################################
  ##############################################################################
  # 3. Make the plotting vectors

  # 3.1 get parameters together that are needed for coordinates
  clens[, chrLab := chrLabFun(chr)]
  maxChrChars <- with(clens, tapply(chrLab, genome, function(x) sum(nchar(x))))
  plotWidth <- (sqrt(max(maxChrChars)) * scalePlotWidth) + 6

  # -- 4.1 polygon matrix
  plist <- rbindlist(lapply(1:nrow(tpb), function(i){
    x <- with(tpb[i, ], calc_curvePolygon(
      start1 = xs1, end1 = xe1, start2 = xs2,
      end2 = xe2, y1 = y1, y2 = y2))
    x[,`:=`(blkID = tpb$blkID[i], col = tpb$col[i])]
    return(x)
  }))
  plist[,u := paste(blkID, col)]
  coll <- subset(plist, !duplicated(plist[,4:5, with = F]))[,4:5, with = F]
  plist[,u := factor(u, levels = coll$u)]

  # -- 4.2 rounded rectangle polygons
  clist <- rbindlist(lapply(1:nrow(clens), function(i){
    z <- clens[i,]
    wid <- ifelse(clens$genome[i] %in% labelTheseGenomes,
                  ylabBuff/1, ylabBuff/2)
    out <- data.table(z[,c("genome","chr")], with(z, round_rect(
      xleft = chrstart, xright = chrend, ybottom = y - wid,
      ytop = y + wid, yrange = range(clens$y),
      xrange = range(c(clens$chrstart, clens$chrend)),
      plotWidth = plotWidth, plotHeight = plotHeight)))
    return(out)
  }))

  # -- 4.3 scale bar
  sbData <- make_scaleBar(chrLengths = clens, ylabBuff = ylabBuff)

  glab <- subset(clens, !duplicated(genome))
  setkey(glab, y)

  ##############################################################################
  # -- more plotting parameter calculations

  ##############################################################################

  ##############################################################################
  # 5. make the plot


  p <- ggplot()+

    # -- all blocks
    geom_polygon(
      data = plist,
      aes(x = x, y = y, group = u, fill = u),
      alpha = braidAlpha)+

    # -- chromosome label backgrounds
    geom_polygon(
      data = clist,
      aes(x = x, y = y, group = paste(genome, chr)),
      fill = chrFill, lwd = chrLabBorderLwd)+

    # -- chr label text
    geom_text(
      data = clens,
      aes(x = (chrstart+chrend)/2, y = y, label = chrLab),
      col = "black", family = "Helvetica", size = 5*0.36, hjust = .5)+

    # -- scale bar
    geom_segment(
      data = sbData$segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      lwd = chrLabBorderLwd, col = chrFill)+
    geom_text(
      data = subset(sbData$segments, line == "mid"),
      aes(x = (x + xend)/2, y = y),
      label = sbData$label,
      col = chrFill, family = "Helvetica", size = 5*0.36, hjust = .5, vjust = 1.5)+

    # -- braid colors
    scale_fill_manual(values = coll$col, guide = "none")+

    # -- genome labels
    scale_y_continuous(breaks = glab$y, labels = glab$genome, expand = c(0.01, 0.01), name = NULL)+
    scale_x_continuous(expand = c(0.005, 0.005))+

    # -- themes
    theme(panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_text(family = "Helvetica", size = 6),
          axis.text.y = element_text(family = "Helvetica", size = 7))+
    labs(x = xlab)

  # -- write the plot
  if(verbose)
    cat(strwrap(sprintf("writing plot to %s", pdfFile), indent = 8, exdent = 8), sep = "\n")
  pdf(pdfFile, height = plotHeight, width = plotWidth)
  print(p)
  dev.off()

  return(list(
    sourceData = list(
      braidPolygon = plist,
      braidCoords = tpb,
      chrPolygon = clist,
      chrCoords = clens,
      genomeLabels = glab,
      scaleBar = sbData,
      cols = coll,
      braidAlpha = braidAlpha,
      chrLabBorderLwd = chrLabBorderLwd,
      chrFill = chrFill,
      plotHeight = plotHeight,
      plotWidth = plotWidth,
      pdfFile =  pdfFile),
    ggplotObj = p))


}

#' @title plot_riparianRegs
#' @description
#' \code{plot_riparianRegs} plot_riparianRegs
#' @rdname plot_riparianRegs
#' @export
#' @import data.table
#' @import ggplot2
#' @export
plot_riparianRegs <- function(...){
  cat("regionPlotting coming in the next genespace release\n")
}

#' @title convert cosine points to polygon
#' @description
#' \code{calc_curvePolygon} from 2d coordinates, make a curve
#' @rdname plot_riparian
#' @export
calc_curvePolygon <- function(start1,
                              end1 = NULL,
                              start2,
                              end2 = NULL,
                              y1,
                              y2,
                              npts = 250,
                              keepat = round(npts / 20)){
  cosine_points <- function(npts, keepat){
    # initial number of points
    # grid to keep always
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

  scaledCurve <- cosine_points(npts = npts, keepat = keepat)
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
#' @rdname plot_riparian
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft, ybottom, xright, ytop, plotWidth, plotHeight, xrange, yrange, npts = 50){

  if (ytop <= ybottom){
    tmp <- ytop
    ytop <- ybottom
    ybottom <- tmp
  }
  if (xright <= xleft){
    tmp <- xleft
    xleft <- xright
    xright <- tmp
  }

  # measure graphics device
  pars <- c(xrange, yrange)
  asp <- diff(pars[3:4]) / diff(pars[1:2])
  dev <- plotWidth / plotHeight

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = npts)
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


#' @title calc_regBlkCoords
#' @description
#' \code{calc_regBlkCoords} calc_regBlkCoords
#' @rdname plot_riparian
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
calc_regBlkCoords <- function(gsParam,
                              blk,
                              regGenome = NULL,
                              regChr = NULL,
                              regStart = NULL,
                              regEnd = NULL,
                              regID = NULL,
                              regCol = NULL,
                              minBlkSize){

  if(is.null(regGenome))
    stop("regGenome must be specified to calculate region block coords")
  if(is.null(regChr))
    stop("regChr must be specified to calculate region block coords")
  if(is.null(regStart))
    regStart <- 0
  if(is.null(regEnd))
    regEnd <- Inf
  if(is.null(regID))
    regID <- 1:length(regChr)
  if(is.null(regCol))
    regCol <- gs_colors(n = length(regChr))
  spFiles <- file.path(gsParam$paths$pangenome, sprintf(
    "%sintegratedSynPos.txt", gsParam$genomeIDs))
  u <- paste(regGenome, regChr)
  interp <- rbindlist(lapply(spFiles, function(x)
    subset(fread(x), isAnchor & paste(interpGenome, interpChr) %in% u & genome != interpGenome)))

  if(!all(regChr %in% interp$interpChr))
    warning("some regChr are not in the reference genome")
  if(!length(regStart) %in% c(1, length(regChr)))
    stop("regStart not of the same length as regChr")
  if(!length(regStart) %in% c(1, length(regChr)))
    stop("regStart not of the same length as regChr")
  if(!length(regEnd) %in% c(1, length(regChr)))
    stop("regEnd not of the same length as regChr")
  if(any(duplicated(regID)) || length(regID) != length(regChr))
    stop("regIDs must be unique and same length as regChrs")

  wh <- which(regChr %in% interp$interpChr)
  regs <- data.table(
    regGenome = regGenome[wh],
    regChr = regChr[wh],
    regStart = regStart[wh],
    regEnd = regEnd[wh],
    id = regID[wh],
    color = regCol)

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  regList <- rbindlist(lapply(1:nrow(regs), function(i){
    x <- regs[i,]
    bedsub <- subset(bed, genome == x$regGenome & chr == x$regChr &
                       start >= x$regStart & end <= x$regEnd & !noAnchor)
    os <- range(bedsub$ord)
    intsub <- subset(
      interp, interpGenome == x$regGenome & genome != x$regGenome &
        interpChr == x$regChr & interpOrd >= os[1] & interpOrd <= os[2])


    out <- rbind(bedsub[,c("genome", "ofID")], intsub[,c("genome", "ofID")])
    out[,`:=`(regID = x$id, color = x$color)]
    return(out)
  }))

  md <- gsParam$annotBlastMd
  ug <- with(blk, unique(paste(c(genome1, genome2), c(genome2, genome1))))
  uc <- with(blk, unique(paste(genome1, genome2)))
  mdo <- subset(md, paste(query, target) %in% ug)
  u <- unique(regList$ofID)
  rhits <- rbindlist(lapply(mdo$annotBlastFile, function(i)
    subset(fread(
      i, na.strings = c("", "NA"),
      select = c("ofID1", "chr1", "start1", "end1", "ord1", "genome1",
                 "ofID2", "chr2", "start2", "end2", "ord2", "genome2",
                 "isAnchor", "lgBlkID")),
      isAnchor & !is.na(lgBlkID) & ofID1 %in% u & ofID2 %in% u)))
  rg1 <- with(regList, data.table(ofID1 = ofID, regID1 = regID, color = color))
  rg2 <- with(regList, data.table(ofID2 = ofID, regID2 = regID))
  rh <- merge(rg1, merge(rg2, rhits, by = "ofID2", allow.cartesian = T),
              by = "ofID1", allow.cartesian = T)
  rh <- subset(rh, regID1 == regID2)
  rh[,`:=`(hlRegID = regID1, blkID = sprintf("%s_hlreg%s",lgBlkID, regID1),
           regID1 = NULL, regID2 = NULL)]

  out <- calc_blkCoords(rh, mirror = T)
  out <- subset(out, paste(genome1, genome2) %in% uc)
  tmp <- subset(rh, !duplicated(blkID))
  colv <- tmp$col; names(colv) <- tmp$blkID
  out[,color := colv[blkID]]
  return(out)
}

#' @title get_chrOrder
#' @description
#' \code{get_chrOrder} get_chrOrder
#' @rdname plot_riparian
#' @export
#' @import data.table
#' @export
get_chrOrder <- function(blk, refGenome, reorderChrsBySynteny, refChrOrdFun){

  if(reorderChrsBySynteny){
    b <- subset(blk, genome1 == refGenome & genome2 != refGenome)
    if(!"start1" %in% colnames(blk)){
      b[,`:=`(start1 = startOrd1, start2 = startOrd2,
              end1 = endOrd1, end2 = endOrd2)]
    }

    # -- reformat blocks
    b <- subset(b, !duplicated(b))

    # -- get unique ref chromosomes ordered
    rco <- refChrOrdFun(unique(b$chr1))
    names(rco) <- unique(b$chr1)
    b[,refOrder := rco[chr1]]

    # -- get proportion of hits to ref chromosomes
    bo <- b[,list(n = sum(abs(start2 - end2))), by = c("genome2", "chr2","refOrder")]
    bo[,prop := n/sum(n), by = c("genome2", "chr2")]
    bo[,wtOrd := prop * refOrder]
    bo[,keep := all(prop < 0.25), by = c("genome2", "chr2")]
    bo$keep[bo$prop > 0.25] <- TRUE
    bo <- subset(bo, keep)
    boc <- bo[,list(wtOrd = round(weighted.mean(refOrder, prop))),
              by = c("genome2", "chr2")]
    boc[,refOrder := refChrOrdFun(chr2), by = "genome2"]
    setkey(boc, genome2, wtOrd, refOrder)
    boc[,chrOrd := 1:.N, by = "genome2"]
    bor <- data.table(genome = refGenome, chr = names(rco), chrOrd = rco)
    boc <- with(boc, data.table(genome = genome2, chr = chr2, chrOrd = chrOrd))
    bocr <- rbind(boc, bor)
    chrOrds <- bocr$chrOrd; names(chrOrds) <- with(bocr, paste(genome, chr))
  }else{
    u <- with(blk, data.table(genome = c(genome1, genome2), chr = c(chr1, chr2)))
    u <- subset(u, !duplicated(u))
    u[,chrOrd := refChrOrdFun(chr), by = "genome"]
    setorder(u, genome, chr, chrOrd)
    u[,ord := 1:.N, by = "genome"]
    chrOrds <- u$ord; names(chrOrds) <- with(u, paste(genome, chr))
  }
  return(chrOrds)
}
