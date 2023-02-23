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
#' @param highlightBed data.frame (or coercible to a data.frame, e.g.
#' data.table) with four required columns: genome, chr, start, end ... these
#' are the bp coordinates of the regions of interest. Can also include a column
#' "color" with the desired color to plot for that particular regions. This
#' can be used to color chromosomes in polyploid references.
#' @param backgroundColor character or integer coercible to an R color. Only
#' used if highlightBed is also specified. This is background to show all blocks
#' across references. If the full background is not desired, can be set to NULL
#' or NA. In which case the plot will ONLY contain the regions specified in
#' highlightBed.
#' @param labelTheseGenomes character string of genomes that have labeled chrs
#' @param invertTheseChrs data.table with two columns, genome and chr containing
#' the lists of genomes and chromosomes that should be inverted in the plot.
#' Not currently implemented.
#' @param braidAlpha numeric (0-1) specifying the transparency of the braids
#' @param braidBorderLwd numeric specifying the weight of borders on the braids.
#' Not currently in use.
#' @param genomeIDs character vector at least partially matching the genomeIDs
#' in gsParam
#' @param refGenome single character string specifying which genome is the ref
#' @param reorderBySynteny logical, should chromosomes be re-ordered by synteny?
#' @param gapProp numeric (0-1) specifying the proportional size of gaps
#' relative to the length of the largest genome
#' @param useOrder logical, should gene rank order be used in lieu of physical
#' positions?
#' @param chrLabCex character expansion for the chromosome labels. Not currently
#' in use.
#' @param chrFill character or integer coercible to an R color for the interior
#' fill of the chromosome regions
#' @param chrLabFun function to parse chr IDs to make them more readible
#' @param pdfFile file.path to the pdf device to store the file
#' @param scalePlotWidth numeric, scaling factor for plot width. Larger values
#' make wider plots
#' @param scalePlotHeight numeric, scaling factor for plot height Larger values
#' make taller plots
#' @param chrLabBorderLwd character or integer coercible to an R color for the border
#' of the chromosome regions.
#' @param refChrOrdFun function, to re-order reference chromosomes
#' @param start1 numeric, start postion of the first block
#' @param end1 numeric, end postion of the first block
#' @param start2 numeric, start postion of the second block
#' @param end2 numeric, end postion of the second block
#' @param y1 numeric, bottom postion of the block
#' @param y2 numeric, top postion of the block
#' @param npts integer, the number of points to use in the curve.
#' @param keepat integer, the grid size of points to keep
#' @param xleft numeric, left position of rounded rectangle
#' @param xright numeric, right position of rounded rectangle
#' @param ytop numeric, top position of rounded rectangle
#' @param ybottom numeric, bottom position of rounded rectangle
#' @param plotWidth numeric, width of the plot
#' @param plotHeight numeric, height of the plot
#' @param xrange numeric, x range of data
#' @param yrange numeric, y range of the data
#' @param blk data.table of block coordinates
#' @param addvert funny stuff here
#' @param scaleGapSize numeric 0-1, specifying the gap scaling level, where 0 means
#' all gaps are the same size regardless of genome size. 1 means that the
#' genomes sizes are all nearly the same and smaller genomes have larger gaps.
#' @param addThemes ggplot2 themes to add to the riparian plot
#' @param verbose logical, should updates be printed to the console?
#' @param useRegions logical, should regions or smaller blocks be plotted?
#' @param minChrLen2plot replaces previous minSize2plot. Minimum size in the
#' scale plotted (order if useOrder = TRUE, bp if useOrder = FALSE)
#' @param scaleBraidGap integer specifying how much of a gap should exist
#' between the center of the chromosome labels (y position) and the begining
#' of the braid polygons. 0 = no gap between braids, 1 = gap equal to the
#' vertical size of the chromosome polygons so that the braids start and the
#' bottom/top of the chromosome polygons. Scales linearly form there.
#' @param howSquare numeric 0-inf, specifying how square the rounded rectangles
#' should be. 0 = no straight sections. Inf = no rounded sections
#' @param customRefChrOrder character vector with the order of chromosomes in
#' the reference genome
#' @param palette function coercible to an R color palette
#' @param xlabel label for the x axis of the plot
#' @param bed data.table with combined bed file
#' @param theseBlocksFirst internal parameter specifying which blocks should be
#' plotted first
#'
#' @details By default, acts directly on the reference-phased syntenic block
#' coordinates generated by integrate_synteny(). Uses these positions to
#' generate linear coordinates for block breakpoints and colors these by
#' the reference genome chromosomes.
#'
#'
#' @examples
#' \dontrun{
#' # see vignette dedicated to plot_riparian, coming soon.
#' }
#'
#' @title Make riparian plots
#' @description
#' \code{plot_riparian} The main GENESPACE plotting routine, which generate
#' braided river or 'riparian' plots.
#' @rdname plot_riparian
#' @import data.table
#' @import ggplot2
#' @export
plot_riparian <- function(gsParam,
                          highlightBed = NULL,
                          blk = NULL,
                          backgroundColor = add_alpha("grey20", .5),
                          pdfFile = NULL,
                          refGenome = gsParam$genomeIDs[1],
                          useRegions = TRUE,
                          genomeIDs = gsParam$genomeIDs,
                          labelTheseGenomes = gsParam$genomeIDs,
                          useOrder = TRUE,
                          gapProp = 0.005,
                          chrFill = "white",
                          chrLabBorderLwd = .5,
                          scalePlotHeight = 1,
                          scalePlotWidth = 1,
                          minChrLen2plot = 100,
                          braidAlpha = .5,
                          scaleBraidGap = 0,
                          reorderBySynteny = TRUE,
                          howSquare = 0,
                          customRefChrOrder = NULL,
                          palette = gs_colors,
                          chrLabCex = 1,
                          braidBorderLwd = 0,
                          invertTheseChrs = NULL,
                          chrLabFun = function(x)
                            gsub("^0", "",
                                 gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x))),
                          xlabel = sprintf(
                            "Chromosomes scaled by %s",
                            ifelse(useOrder, "gene rank order", "physical position")),
                          scaleGapSize = 0.25,
                          addThemes = NULL,
                          verbose = FALSE,
                          refChrOrdFun = function(x)
                            frank(list(as.numeric(gsub('\\D+','', x)), x), ties.method = "random")){

  color <- blkID <- pass <- genome <- chr <- start <- end <- NULL
  # -- there are two methods of plotting.

  # 1. Default phase the blocks by a reference genome chromosomes
  # 2. A background of one color and then highlighted colors overlapping it.

  if(!"synteny" %in% names(gsParam))
    gsParam <- set_syntenyParams(gsParam)

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  if(is.null(highlightBed)){

    hapGs <- names(gsParam$ploidy)[gsParam$ploidy == 1]
    if(length(hapGs) == 0)
      stop("riparian plots must have a haploid reference as an anchor. Either re-run GENESPACE with a haploid outgroup or use custom parameters in `riparian_engine`")
    if(is.null(refGenome))
      refGenome <- hapGs[1]
    if(!refGenome %in% hapGs)
      stop(sprintf("specified reference genome %s is not haploid, but haploid genomes %s are in this GENESPACE run. Either use one of the haploid references or use custom parameters in `riparian_engine`",
                   refGenome, paste(hapGs, collapse = ", ")))

    blksTp <- phase_blks(
      gsParam = gsParam,
      refGenome = refGenome,
      useRegions = useRegions,
      blkSize = gsParam$params$blkSize,
      synBuff = gsParam$params$synBuff,
      refChr = NULL,
      refStartBp = NULL,
      refEndBp = NULL)

    pltDat <- riparian_engine(
      blk = blksTp, bed = bed, pdfFile = pdfFile, refGenome = refGenome,
      genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
      useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
      chrLabBorderLwd = chrLabBorderLwd, scalePlotHeight = scalePlotHeight,
      scalePlotWidth = scalePlotWidth, minChrLen2plot = minChrLen2plot,
      braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap, verbose = verbose,
      reorderBySynteny = reorderBySynteny, howSquare = howSquare,
      customRefChrOrder = customRefChrOrder, addThemes = addThemes,
      invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
      scaleGapSize = scaleGapSize, palette = palette,refChrOrdFun = refChrOrdFun)
    if(is.null(pdfFile))
      print(pltDat$ggplotObj)
  }else{
    highlightBed <- data.table(highlightBed)
    if(nrow(highlightBed) == 0)
      stop("highlightBed must be a data.table or data.frame with >= 1 rows\n")
    if(!all(c("chr", "genome") %in% names(highlightBed)))
      stop("highlightBed must be a data.table or data.frame with columns named chr and genome\n")
    u <- unique(paste(bed$genome, bed$chr))
    highlightBed[,pass := paste(genome, chr) %in% u]
    if(!all(highlightBed$pass))
      stop("the following genome/chr combinations are not in the run:", subset(bed, !pass)[,1:2])
    highlightBed[,pass := NULL]
    if(!"start" %in% names(highlightBed))
      highlightBed[,start := 0]
    if(!"end" %in% names(highlightBed))
      highlightBed[,end := Inf]

    if(is.null(backgroundColor) || is.na(backgroundColor)){
      allBlks <- NULL
    }else{
      hapGs <- names(gsParam$ploidy)[gsParam$ploidy == 1]
      if(is.null(refGenome))
        refGenome <- hapGs[1]

      # -- get background all blocks
      allBlks <- nophase_blks(
        gsParam = gsParam,
        useRegions = useRegions,
        blkSize = gsParam$params$blkSize)
      allBlks[,color := backgroundColor]
      allBlks[,blkID := sprintf("all_%s", blkID)]
    }

    # -- get blocks in the intervals of interest

    if(!"color" %in%  names(highlightBed))
      highlightBed[,color := gs_colors(.N)]
    blksBed <- rbindlist(lapply(1:nrow(highlightBed), function(i){
      blksTp <- with(highlightBed, phase_blks(
        gsParam = gsParam,
        refGenome = genome[i],
        useRegions = useRegions,
        blkSize = gsParam$params$blkSize,
        synBuff = gsParam$params$synBuff,
        refChr = chr[i],
        refStartBp = start[i],
        refEndBp = end[i]))
      blksTp[,color := highlightBed$color[i]]
      return(blksTp)
    }))

    blksTp <- rbind(allBlks, blksBed, fill = T)

    pltDat <- riparian_engine(
      blk = blksTp, bed = bed, pdfFile = pdfFile, refGenome = refGenome,
      theseBlocksFirst = unique(allBlks$blkID), verbose = verbose,
      genomeIDs = genomeIDs, labelTheseGenomes = labelTheseGenomes,
      useOrder = useOrder, gapProp = gapProp, chrFill = chrFill,
      chrLabBorderLwd = chrLabBorderLwd, scalePlotHeight = scalePlotHeight,
      scalePlotWidth = scalePlotWidth, minChrLen2plot = minChrLen2plot,
      braidAlpha = braidAlpha, scaleBraidGap = scaleBraidGap,
      reorderBySynteny = reorderBySynteny, howSquare = howSquare,
      customRefChrOrder = customRefChrOrder, addThemes = addThemes,
      invertTheseChrs = invertTheseChrs, chrLabFun = chrLabFun, xlabel = xlabel,
      scaleGapSize = scaleGapSize, refChrOrdFun = refChrOrdFun,
      palette = colorRampPalette(backgroundColor))
    if(is.null(pdfFile))
      print(pltDat$ggplotObj)
  }
  return(list(blks = blksTp, plotData = pltDat))
}


#' @title engine to create the riparian plot
#' @description
#' \code{riparian_engine} underlying routines for plot_riparian, related to
#' plotting and data parsing.
#' @rdname plot_riparian
#' @import data.table
#' @import ggplot2
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
riparian_engine <- function(blk,
                            bed,
                            pdfFile = NULL,
                            refGenome = NULL,
                            genomeIDs = NULL,
                            theseBlocksFirst = NULL,
                            labelTheseGenomes = NULL,
                            useOrder = TRUE,
                            gapProp = 0.005,
                            chrFill = "white",
                            chrLabBorderLwd = .5,
                            scalePlotHeight = 1,
                            scalePlotWidth = 1,
                            minChrLen2plot = 100,
                            braidAlpha = .5,
                            scaleBraidGap = 0,
                            reorderBySynteny = TRUE,
                            howSquare = 0,
                            chrLabCex = 1,
                            customRefChrOrder = NULL,
                            palette = gs_colors,
                            invertTheseChrs = NULL,
                            braidBorderLwd = 0,
                            chrLabFun = function(x)
                              gsub("^0", "",
                                   gsub("chr|scaf|chromosome|scaffold|^lg|_", "", tolower(x))),
                            xlabel = sprintf(
                              "Chromosomes scaled by %s",
                              ifelse(useOrder, "gene rank order", "physical position")),
                            scaleGapSize = 0.25,
                            addThemes = NULL,
                            verbose = FALSE,
                            refChrOrdFun = function(x)
                              frank(list(as.numeric(gsub('\\D+','', x)), x), ties.method = "random")){

  ##############################################################################
  subset_toChainedGenomes <- function(genomeIDs, blk){
    genome1 <- genome2 <- NULL
    genomeOrd <- data.table(
      genome1 = genomeIDs[-length(genomeIDs)], genome2 = genomeIDs[-1])
    genomeOrd[,`:=`(y1 = match(genome1, genomeIDs),
                    y2 = match(genome2, genomeIDs))]
    u <- with(genomeOrd, paste(genome1, genome2))
    if(!all(u %in% paste(blk$genome1, blk$genome2)))
      warning("problem with the block coordinates, some chained combinations of genomeIDs are not in the blocks\n")
    blk <- merge(genomeOrd, blk, by = c("genome1", "genome2"))
    return(blk)
  }

  ##############################################################################
  get_chrLens <- function(blk, minChrLen2plot, useOrder){
    genome1 <- genome2 <- chr1 <- chr2 <- startBp1 <- startBp2 <- endBp1 <-
      endBp2 <- startOrd1 <- startOrd2 <- endOrd1 <- endOrd2 <- end <- length <-
      genome <- chr <- NULL
    if(useOrder){
      chrLengths <- with(blk, data.table(
        genome = c(genome1, genome2, genome1, genome2),
        chr = c(chr1, chr2, chr1, chr2),
        end = c(startOrd1, startOrd2, endOrd1, endOrd2)))
    }else{
      chrLengths <- with(blk, data.table(
        genome = c(genome1, genome2, genome1, genome2),
        chr = c(chr1, chr2, chr1, chr2),
        end = c(startBp1, startBp2, endBp1, endBp2)))
    }

    chrLengths <- chrLengths[,list(length = max(end)), by = c("genome", "chr")]
    u <- with(blk, unique(paste(c(genome1, genome2), c(chr1, chr2))))
    if(!all(u %in% paste(chrLengths$genome, chrLengths$chr))){
      warning("could not parse chrLengths, some genome/chr combinations in block coordinates are not in chrLengths\n")
      chrLengths <- NULL
    }

    out <- subset(chrLengths, length >= minChrLen2plot)
    return(out)
  }

  ##############################################################################
  make_scaleBarData <- function(chrLengths, ylabBuff, useOrder){
    length <- genome <- s <- chrstart <- genome <- bot <- top <- mid <- NULL
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

    top <- max(chrLengths$y2) + (ylabBuff * 4)
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
  pull_synChrOrd <- function(refGenome, bed, clens){

    ordByFun <- genome <- noAnchor <- refchr <- reford <- ord <- median <-
      plotOrd <- meanPos <- NULL

    setkey(clens, ordByFun)
    refChrOrder <- clens$chr[clens$genome == refGenome]

    bedr <- with(subset(bed, genome == refGenome & !noAnchor), data.table(
      refgenome = refGenome, refchr = chr, reford = ord, og = og))
    beda <- subset(bed, genome != refGenome & !noAnchor)[,c("genome", "chr", "ord", "og")]
    bedm <- merge(
      subset(bedr, !duplicated(bedr)),
      subset(beda, !duplicated(beda)), by = "og", allow.cartesian = T)
    bedm[,refChrOrder := match(refchr, refChrOrder)]
    bedm[,reford := frank(reford, ties.method = "dense"), by = "refchr"]
    bedm[,ord := (1:.N)/.N, by = "refchr"]

    chro <- bedm[,list(meanPos = round(median(refChrOrder + ord), 2)),
                 by = c("genome", "chr")]

    outa <- merge(clens, chro, by = c("genome", "chr"))
    outr <- subset(clens, genome == refGenome)
    outr[,plotOrd := ordByFun]

    setkey(outa, meanPos, ordByFun, genome)
    outa[,plotOrd := 1:.N, by = "genome"]

    # outa <<- outa
    # outr <<- outr

    outa <- outa[,colnames(outr), with = F]
    return(rbind(outr, outa))
  }

  ##############################################################################
  reorder_inChrs <- function(blk){
    genome1 <- genome2 <- chr1 <- chr2 <- startOrd1 <- startOrd2 <- start <-
      genome <- chr <- endOrd1 <- endOrd2 <- chr1st <- chr2st <- mpv <- NULL
    minPos <- with(blk, data.table(
      genome = c(genome1, genome2),
      chr = c(chr1, chr2),
      start = c(startOrd1, startOrd2)))
    minPos <- minPos[,list(chrStart = min(start)), by = c("genome", "chr")]
    mpv <- minPos$chrStart; names(mpv) <- with(minPos, paste(genome, chr))
    blk[,`:=`(chr1st = mpv[paste(genome1, chr1)],
              chr2st = mpv[paste(genome2, chr2)])]
    blk[,`:=`(startOrd1 = (startOrd1 - chr1st) + 1,
              startOrd2 = (startOrd2 - chr2st) + 1,
              endOrd1 = (endOrd1 - chr1st) + 1,
              endOrd2 = (endOrd2 - chr2st) + 1)]
    return(blk)
  }

  ##############################################################################
  calc_gapSize <- function(chrLens, scaleGapSize, gapProp, rg){
    maxGapSize <- genomeSize <- totLen <- ngaps <- totLen <- genome <-
      diffMax <- gapSize <- NULL
    maxGapSize <- round(gapProp * sum(chrLens$length[chrLens$genome == rg]))
    genomeSize <- chrLens[,list(totLen = sum(as.numeric(length)), ngaps = .N),
                          by = "genome"]
    genomeSize[,diffMax := (max(totLen) - totLen)/ngaps]
    genomeSize[,gapSize := maxGapSize + (diffMax * scaleGapSize)]
    gapSize <- genomeSize$gapSize; names(gapSize) <- genomeSize$genome
    return(gapSize)
  }


 blkID <- refChr <- genome1 <- genome2 <- genome <- ordByFun <- chr <-
   plotOrd <- gapSize <- chrstart <- med <- ht <- y1 <- y2 <- chrLab <- u1 <-
   xs1 <- cstart <- xe1 <- cend <- difstart <- difend <- u2 <- xs2 <- xe2 <-
   x <- y <- color <- x1 <- x2 <- xend <- yend <- line <- NULL

  ##############################################################################
  ##############################################################################
  refG <- refGenome
  gids <- genomeIDs

  blk <- data.table(blk)
  blksFileOrd <- unique(blk$blkID)
  refGenome <- NULL
  genomeIDs <- NULL

  if(is.null(gids))
    gids <- unique(unlist(blk[,c("genome1", "genome2")]))
  fullblk <- data.table(blk)
  if(is.null(refG))
    refG <- gids[1]
  if(is.null(labelTheseGenomes))
    labelTheseGenomes <- gids

  if(!is.null(pdfFile)){
    if(!dir.exists(dirname(pdfFile))){
      cat(strwrap(sprintf(
        "NOTE: parent directory for %s does not exist. Just writing the plot to the current graphics device. 'Square-ness' of the chromosome polygons may be affected", pdfFile),
        indent = 8, exdent = 8), sep = "\n")
      pdfFile <- NULL
    }else{
      if(verbose)
        cat(strwrap(sprintf(
          "writing plot to %s", pdfFile),
          indent = 8, exdent = 8), sep = "\n")
    }
  }

  if(!is.null(invertTheseChrs)){
    if(!is.data.frame(invertTheseChrs))
      stop("invertTheseChrs given, but this needs to be a data.frame\n")
    if(!all(c("genome", "chr") %in% colnames(invertTheseChrs)))
      stop("invertTheseChrs data.frame given, but does not have col names 'chr' and 'genome' \n")
    invchr <- with(invertTheseChrs, paste(genome, chr))
  }else{
    invchr <- NULL
  }

  if(!is.null(customRefChrOrder)){
    customRefChrOrder <- as.character(customRefChrOrder)
    if(any(is.null(customRefChrOrder)) || any(is.na(customRefChrOrder)))
      stop("customRefChrOrder is given but can't be coerced to character\n")
  }

  plotHeight <- (uniqueN(gids) * scalePlotHeight) + 1
  plotWidth <- 8 * scalePlotWidth
  ylabBuff <- length(gids) / 80

  ##############################################################################
  # 1 Get chromosome positions nailed down
  # -- 1.1 if useorder, re-scale order so that each chr starts at 1
  if("refChr" %in% colnames(blk))
    blk[,blkID := paste(blkID, refChr)]
  blk <- subset(blk, genome1 %in% gids & genome2 %in% gids)
  if(useOrder)
    blk <- reorder_inChrs(blk)

  # -- 1.2 calculate chromosome lengths
  clens <- get_chrLens(
    blk = blk,
    useOrder = useOrder,
    minChrLen2plot = minChrLen2plot)

  # -- 1.3 get chromosome order by chromosome name
  if(!is.null(customRefChrOrder)){
    rlen <- subset(clens, genome == refG)
    if(any(!rlen$chr %in% customRefChrOrder))
      stop("customRefChrOrder but not all refchrs are in this list\n")
    rlen[,ordByFun := as.integer(factor(chr, levels = customRefChrOrder))]
    nlen <- subset(clens, genome != refG)
    nlen[,ordByFun := refChrOrdFun(chr), by = "genome"]
    clens <- rbind(rlen, nlen)
  }else{
    clens[,ordByFun := refChrOrdFun(chr), by = "genome"]
  }
  setkey(clens, genome, ordByFun)

  # -- 1.4 [optional] replace non-reference genome chr order to max synteny
  if(reorderBySynteny){
    clens <- pull_synChrOrd(refGenome = refG, bed = bed, clens = clens)
  }else{
    clens[,plotOrd := ordByFun]
  }
  clens[,ordByFun := NULL]

  # -- 1.5 get genome-specific gap sizes
  setkey(clens, genome, plotOrd)
  gapSizes <- calc_gapSize(
    chrLens = clens, scaleGapSize = scaleGapSize,
    gapProp = gapProp, rg = refG)

  # -- 1.6 add gap to the chrlengths
  clens[,gapSize := gapSizes[genome]]

  # -- 1.7 add cumulative position to the chromosomes
  clens[,chrstart := cumsum(c(1, (length[-.N]+1)) + gapSize)-gapSize,
        by = "genome"]

  # -- 1.8 get median position for each genome
  clens[,med := max(chrstart + length)/2,
        by = "genome"]

  # -- 1.9 add start and end positions of each chromosome in linear coords
  clens[,`:=`(x1 = chrstart - med, x2 = (chrstart + length) - med)]

  # -- 1.10 add thickness of the chr polygons
  clens[,ht := ifelse(genome %in% labelTheseGenomes, ylabBuff/1, ylabBuff/2)]
  clens[,`:=`(y1 = match(genome, gids),
              y2 = match(genome, gids))]

  # -- 1.11 if desired, invert specified chromosomes
  # if(!is.null(invertTheseChrs)){
  #   cleni <- subset(clens, paste(genome, chr) %in% invchr)
  #   if(nrow(cleni) == 0)
  #     stop("invertTheseChrs given, but no matches found\n")
  #   setnames(cleni, c("x1","x2"), c("x2", "x1"))
  #   cleno <- subset(clens, !paste(genome, chr) %in% invchr)
  #   clens <- rbind(cleni, cleno)
  # }

  # -- 1.12 pull vector of chromosome start positions
  chrstartv <- clens$x1
  chrytop <- clens$y1 + (ylabBuff * scaleBraidGap)
  chrybot <- clens$y2 - (ylabBuff * scaleBraidGap)

  names(chrstartv) <- names(chrybot) <-
    names(chrytop) <- with(clens, paste(genome, chr))

  clens[,`:=`(y1 = y1 + ht,
              y2 = y2 - ht,
              ht = NULL)]

  # -- 1.11 make the polygons
  chrPolygons <- rbindlist(lapply(1:nrow(clens), function(i){
    z <- clens[i,]
    out <- data.table(z[,c("genome","chr")], with(z, round_rect(
      xleft = x1, xright = x2, ybottom = y2,
      ytop = y1, yrange = range(c(clens$y1, clens$y2)),
      xrange = range(c(clens$x1, clens$x2)),
      plotWidth = plotWidth + (plotWidth * howSquare),
      plotHeight = plotHeight)))
    return(out)
  }))


  # -- 1.13 add chromosome labels
  clens[,chrLab := chrLabFun(chr)]
  if(!is.null(invchr)){
    wh <- which(paste(clens$genome, clens$chr) %in% invchr)
    clens$chrLab[wh] <- sprintf("%s*", clens$chrLab[wh])
  }

  ##############################################################################
  # 2 Get block positions nailed down

  # -- 2.1 subset blocks to daisy chains
  dsychn <- subset_toChainedGenomes(genomeIDs = gids, blk = blk)

  # -- 2.2 add color column
  if("color" %in% colnames(dsychn)){
    if(any(!are_colors(dsychn$color)))
      stop("colors are given in the blocks but can't be coerced to R colors\n")
  }else{
    ulevs <- subset(clens, genome == dsychn$refGenome[1])[,c("plotOrd", "chr")]
    setkey(ulevs, plotOrd)
    ulevs[,`:=`(color = palette(.N), refChr = chr, plotOrd = NULL, chr = NULL)]
    dsychn <- merge(dsychn, ulevs, by = "refChr")
  }


  # -- 2.4 simplify and convert to projections depending on useOrder
  if(useOrder){
    tp <- with(dsychn, data.table(
      xs1 = startOrd1 + chrstartv[paste(genome1, chr1)],
      xs2 = startOrd2 + chrstartv[paste(genome2, chr2)],
      xe1 = endOrd1 + chrstartv[paste(genome1, chr1)],
      xe2 = endOrd2 + chrstartv[paste(genome2, chr2)],
      y1 = chrytop[paste(genome1, chr1)],
      y2 = chrybot[paste(genome2, chr2)],
      u1 = paste(genome1, chr1),
      u2 = paste(genome2, chr2),
      color = color, blkID = blkID))
  }else{
    tp <- with(dsychn, data.table(
      xs1 = startBp1 + chrstartv[paste(genome1, chr1)],
      xs2 = startBp2 + chrstartv[paste(genome2, chr2)],
      xe1 = endBp1 + chrstartv[paste(genome1, chr1)],
      xe2 = endBp2 + chrstartv[paste(genome2, chr2)],
      y1 = chrytop[paste(genome1, chr1)],
      y2 = chrybot[paste(genome2, chr2)],
      u1 = paste(genome1, chr1),
      u2 = paste(genome2, chr2),
      color = color, blkID = blkID))
  }

  # -- 2.5 if necessary invert some chromosomes
  if(!is.null(invchr)){
    tmp <- subset(clens, paste(genome, chr) %in% invchr)
    chrStarts <- tmp$x1; names(chrStarts) <- with(tmp, paste(genome, chr))
    chrEnds <- tmp$x2; names(chrEnds) <- with(tmp, paste(genome, chr))

    tp1 <- subset(tp, u1 %in% invchr)
    if(nrow(tp1) > 0){
      tp0 <- subset(tp, !u1 %in% invchr)
      tp1[,`:=`(cstart = chrStarts[u1], cend = chrEnds[u1])]
      tp1[,`:=`(difstart = xs1 - cstart, difend = xe1 - cstart)]
      tp1[,`:=`(xs1 = cend - difstart, xe1 = cend - difend)]
      tp1[,`:=`(cstart = NULL, cend = NULL, difstart = NULL, difend = NULL)]
      tp <- rbind(tp0, tp1)
    }

    tp2 <- subset(tp, u2 %in% invchr)
    if(nrow(tp2) > 0){
      tp0 <- subset(tp, !u2 %in% invchr)
      tp2[,`:=`(cstart = chrStarts[u2], cend = chrEnds[u2])]
      tp2[,`:=`(difstart = xs2 - cstart, difend = xe2 - cstart)]
      tp2[,`:=`(xs2 = cend - difstart, xe2 = cend - difend)]
      tp2[,`:=`(cstart = NULL, cend = NULL, difstart = NULL, difend = NULL)]
      tp <- rbind(tp0, tp2)
    }
  }

  # -- 2.6 make the braid polygons
  braidPolygons <- rbindlist(lapply(1:nrow(tp), function(i){
    x <- with(tp[i, ], calc_curvePolygon(
      start1 = xs1, end1 = xe1, start2 = xs2,
      end2 = xe2, y1 = y1, y2 = y2, addvert = 0))
    x[,`:=`(blkID = tp$blkID[i], color = tp$color[i])]
    return(x)
  }))

  # -- 6. Make the scalebar
  sbData <- make_scaleBarData(
    chrLengths = clens, ylabBuff = ylabBuff, useOrder = useOrder)

  glab <- subset(clens, !duplicated(genome))
  setkey(glab, y1)

  if(is.null(theseBlocksFirst)){
    p <- ggplot()+

      # -- all blocks
      geom_polygon(
        data = braidPolygons,
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = braidBorderLwd)
  }else{
    p <- ggplot()+

      # -- all blocks
      geom_polygon(
        data = subset(braidPolygons, blkID %in% theseBlocksFirst),
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = braidBorderLwd)+

      # -- all blocks
      geom_polygon(
        data = subset(braidPolygons, !blkID %in% theseBlocksFirst),
        aes(x = x, y = y, group = blkID, fill = color),
        alpha = braidAlpha, lwd = braidBorderLwd)
  }
  p <- p +

    # -- braid colors
    scale_fill_identity(guide = "none")+

    # -- chromosome label backgrounds
    geom_polygon(
      data = chrPolygons,
      aes(x = x, y = y, group = paste(genome, chr)),
      fill = chrFill, lwd = chrLabBorderLwd)+

    # -- chr label text
    geom_text(
      data = subset(clens, genome %in% labelTheseGenomes),
      aes(x = (x1+x2)/2, y = (y1+y2)/2, label = chrLab),
      col = "black", size = 5*0.36*chrLabCex, hjust = .5)+

    # -- scale bar
    geom_segment(
      data = sbData$segments,
      aes(x = x, xend = xend, y = y, yend = yend),
      lwd = chrLabBorderLwd, col = chrFill)+
    geom_text(
      data = subset(sbData$segments, line == "mid"),
      aes(x = (x + xend)/2, y = y),
      label = sbData$label,
      col = chrFill,size = 5*0.36*chrLabCex, hjust = .5, vjust = 1.5)+

    # -- genome labels
    scale_y_continuous(breaks = (glab$y1+glab$y2)/2, labels = glab$genome, expand = c(0.01, 0.01), name = NULL)+
    scale_x_continuous(expand = c(0.005, 0.005), name = xlabel)+

    # -- themes
    theme(panel.background = element_rect(fill = "black"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())

  if(!is.null(addThemes))
    p <- p + addThemes

  # -- write the plot
  if(!is.null(pdfFile)){
    pdf(pdfFile, height = plotHeight, width = plotWidth)
    print(p)
    dev.off()
  }

  return(list(ggplotObj = p,
              sourceData = list(blocks = dsychn, chromosomes = clens)))
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
                              addvert = 0,
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
  # print(scaledCurve)
  if (!is.null(end1) | !is.null(end2)) {
    sc1 <- scaledCurve[,1]
    sc2 <- scaledCurve[,2]
    if(addvert > 0){
      sc1 <- c(rep(sc1[1], addvert), sc1, rep(sc1[length(sc1)], addvert))
      sc12 <- c(rep(sc2[1], addvert), sc2, rep(sc2[length(sc2)], addvert))
    }
    tp <- rbind(
      start1 = data.table(
        x = start1, y = y1),
      poly1 = data.table(
        x = scale_between(x = sc1, min = start1, max = start2),
        y = scale_between(x = sc2, min = y1, max = y2)),
      start2 = data.table(x = start2, y = y2),
      end2 = data.table(
        x = end2, y = y2),
      poly2 = data.table(
        x = scale_between(x = sc1, min = end2, max = end1),
        y = scale_between(x = sc2, min = y2, max = y1)),
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
