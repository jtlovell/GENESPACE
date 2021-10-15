#' @title Genespace plotting routines
#' @description
#' \code{plot_riparian} Genespace plotting routines
#' @name plot_riparian
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param plotRegions logical, should regions be plot (and not blocks)?
#' @param onlyTheseChrs character vector specifying the reference chrs that
#' should be plot
#' @param labelTheseGenomes character string of genomes that have labeled chrs
#' @param genomeLabCex chracter expansion of the genome labels
#' @param minGenesOnChr integer, specifying the min number of genes a chromosome
#' that is plotted cna contain
#' @param refHits data.table of hits against the reference
#' @param useBlks logical, should blocks be used instead of regions?
#' @param blkSize integer specifying minimum block size
#' @param nCores integer specifying the number of parallel processes to run
#' @param braidAlpha numeric (0-1) specifying the transparency of the braids
#' @param braidBorderLwd numeric specifying the weight of borders on the braids
#' @param genomeIDs character vector at least partially matching the genomeIDs
#' in gsParam
#' @param blackBg logical, should the background be dark?
#' @param refGenome single character string specifying which genome is the ref
#' @param reorderChrs logical, should chromosomes be re-ordered by synteny?
#' @param minGenes integer specifying the minimum number of genes on a chr to
#' plot
#' @param gapProp numeric (0-1) specifying the proportional size of gaps
#' relative to the length of the largest genome
#' @param useOrder logical, should gene rank order be used in lieu of physical
#' positions?
#' @param chrLabCex character expansion for the chromosome labels
#' @param chrRectBuffer number, the amount of buffer around the center of chr
#' regions.
#' @param chrBorder character or integer coercible to an R color for the border
#' of the chromosome regions
#' @param chrFill character or integer coercible to an R color for the interior
#' fill of the chromosome regions
#' @param colByChrs if braids should be colored by reference chromosome,
#' a vector of colors or integers that is coerced into a color palette. If
#' no coloring by chr, specify NULL.
#' @param labelChrBiggerThan integer specifying the minimum number of genes or
#' bp to label a chromosome
#' @param chrLabFun function to parse chr IDs to make them more readible
#' @param gff annotated gff-like data.table
#' @param hitsRef data.table containing all syn og hits against the reference
#' @param hitsRip data.table containing all syn og hits among riparian genomes
#' @param chrd data.table converted chrList
#' @param chrl see chrList
#' @param minrl minimum run-length of consecutive refChr hits in a block to
#' warrant splitting.
#' @param minprp minimum proportion of hits to a refChr to get a block break
#' @param minGenes2plot integer specifying the minimum number of genes on a
#' chr to plot
#' @param onlyTheseRegions data.table with genome, chr, start and end columns
#' @param excludeChrOutOfRegion logical, should chromosome representations
#' be constrained to just those in synteny with the rest of the graph?
#' @param findRegHitsRecursive logical, should regional hit discovery be
#' recursive?
#' @param nGenomeLabChar number of characters for genome labes
#' @details ...
#'
#' @note \code{plot_riparian} is a generic name for the functions documented.
#' \cr
#' If called, \code{plot_riparian} returns its own arguments.
#'
#' @title Riparian plot
#' @description
#' \code{plot_riparian} Make a riparian plot
#' @rdname plot_riparian
#' @import data.table
#' @importFrom graphics strheight polygon
#' @importFrom stats quantile
#' @export
plot_riparian <- function(gsParam,
                          plotRegions = TRUE,
                          refGenome = NULL,
                          genomeIDs = NULL,
                          onlyTheseChrs = NULL,
                          onlyTheseRegions = NULL,
                          excludeChrOutOfRegion = FALSE,
                          findRegHitsRecursive = FALSE,
                          highlightRef = "white",
                          minGenes2plot = 50,
                          braidAlpha = .8,
                          braidBorderLwd = NULL,
                          reorderChrs = TRUE,
                          minGenes = 5,
                          gapProp = .005,
                          useOrder = TRUE,
                          chrLabCex = .5,
                          nGenomeLabChar = 20,
                          genomeLabCex = .75,
                          chrBorder = "black",
                          chrFill = "white",
                          chrRectBuffer = 1.5,
                          colByChrs = NULL,
                          labelChrBiggerThan = NULL,
                          labelTheseGenomes = NULL,
                          blackBg = TRUE,
                          chrLabFun = function(x)
                            gsub("^0","",gsub("^chr|^scaffold|^lg|_","",tolower(x)))){

  arrayID <- og <- synOG <- globOG <- inBlkOG <- synOg <- NULL
  genome <- ofID1 <- ofID2 <- chr1 <- chr <- gen1 <- ord1 <- ofID <- NULL
  rl <- refChr <- blkID <- gen2 <- startOrd1 <- endOrd1 <- end <- n <- NULL
  startOrd2 <- endOrd2 <- startBp1 <- endBp1 <- startBp2 <- NULL
  endBp2 <- firstGene1 <- x <- linBp <- linOrd <- y <- start <- NULL
  ##############################################################################
  # 1. rename a few things, check parameters, read in hits/gff
  ##############################################################################
  # -- specify all the genomes
  highlightRef <- highlightRef[1]
  if(!are_colors(highlightRef) || is.null(highlightRef))
    highlightRef <- "white"

   if(!all(are_colors(colByChrs)) || any(is.null(colByChrs))){
    if(blackBg)
      colByChrs <- c(
        "#BC4F43", "#F67243", "#FFA856", "#FFD469", "#F3E97F", "#C4EAE5",
        "#9DF4FF", "#9CC9FF", "#7D94FF", "#5758D6", "#8253C2",
        "#CF84FF", "#EDCDFF")
    if(!blackBg)
      colByChrs <- c(
        "#8B1D00", "#D05100", "#ED9004", "#F9C70E", "#EAE075", "#BAE0DB",
        "#8BEDF9", "#74B8FC", "#4871F9", "#040DC9","#0E004C","#5E09A3",
        "#C054F9","#E6BDFC")
  }

  verbose <- gsParam$params$verbose
  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
    genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]
  }else{
    tmp <- gsParam$genomes$genomeIDs
    tmp <- tmp[!tmp %in% gsParam$genomes$outgroup]
    genomeIDs <- genomeIDs[genomeIDs %in% tmp]
  }
  nGenomeLabChar <- min(max(nchar(genomeIDs)), nGenomeLabChar)

  # -- specify the ref genome
  if(is.null(refGenome) || !refGenome %in% genomeIDs || length(refGenome) > 1)
    refGenome <- genomeIDs[1]

  if(is.null(labelTheseGenomes))
    labelTheseGenomes <- genomeIDs
  labelTheseGenomes <- labelTheseGenomes[!duplicated(labelTheseGenomes)]
  labelTheseGenomes <- labelTheseGenomes[labelTheseGenomes %in% genomeIDs]
  if(length(genomeIDs) < 1)
    labelTheseGenomes <- genomeIDs
  # -- read the gff
  if(verbose)
    cat("\tLoading the gff ... ")
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, showProgress = F, na.strings = c("", "NA"))

  if("refGenome" %in% colnames(gff))
    gff[,refGenome:=NULL]

  if(!any(is.na(gff$inBlkOG))){
    gff[,og := inBlkOG]
  }else{
    if(!any(is.na(gff$synOG))){
      gff[,og := synOG]
    }else{
      gff[,og := globOG]
    }
  }

  gff[,synOg := og]
  gv <- gff$genome; cv <- gff$chr
  names(gv) <- names(cv) <- gff$ofID

  #-- read ref hits
  if(verbose)
    cat(sprintf("Done!\n\tMapping genes against %s chromosomes ... ", refGenome))
  refh <- load_refHits(
    gsParam = gsParam,
    plotRegions = plotRegions,
    genomeIDs = genomeIDs,
    refGenome = refGenome)
  refh[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2])]
  refh[,n := .N, by = c("gen1","gen2","chr1","chr2")]
  refh <- subset(refh, n >= minGenes2plot)

  if(!is.null(onlyTheseChrs)){
    if(all(onlyTheseChrs %in% refh$chr1)){
      refh <- subset(refh, chr1 %in% onlyTheseChrs)
      gcu <- with(refh, unique(c(paste(gen1, chr1), paste(gen2, chr2))))
      gff <- subset(gff, paste(genome, chr) %in% gcu)
    }
  }

  if(!is.null(onlyTheseRegions)){
    setkey(gff, genome, chr, start, end)
    setkey(onlyTheseRegions, genome, chr, start, end)
    fo <- foverlaps(gff, onlyTheseRegions)
    fo <- subset(fo, complete.cases(fo[,c("ofID", "start", "end")]))
    genesInReg <- gff$ofID[gff$og %in% unique(fo$og)]

    if(excludeChrOutOfRegion){
      gff <- subset(gff, ofID %in% genesInReg)
    }else{
      gcu <- with(refh, unique(c(paste(gen1, chr1), paste(gen2, chr2))))
      gff <- subset(gff, paste(genome, chr) %in% gcu)
    }
  }else{
    genesInReg <- unique(gff$ofID)
  }

  # -- get the colors
  rg <- reorder_gff(
    gff = subset(gff, genome == refGenome),
    minGenesOnChr = minGenes2plot,
    genomeIDs = genomeIDs,
    refGenome = refGenome)
  if(length(colByChrs) != uniqueN(rg$chr)){
    cols <- colorRampPalette(colByChrs)(uniqueN(rg$chr))
  }else{
    cols <- colByChrs
  }

  names(cols) <- unique(rg$chr)

  if(verbose)
    cat("Done!\n\tProjecting linear coordinate system ... ")
  # -- load the gff
  gff <- reorder_gff(
    gff = gff,
    minGenesOnChr = minGenes2plot,
    genomeIDs = genomeIDs,
    refGenome = refGenome)
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  refh[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]

  gff <- add_synChr2gff(
    gff = data.table(gff),
    refHits = refh,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    gapProp = gapProp)
  gv <- gff$genome; cv <- gff$chr; ov <- gff$linOrd; sv <- gff$linBp
  names(gv) <- names(cv) <- names(ov) <- names(sv) <- gff$ofID

  # -- get vector of inferred (or actual) reference genome chrs
  rcv <- calc_refChrByGene(
    gff = gff,
    refGenome = refGenome,
    blkSize = 5,
    nCores = gsParam$params$nCores)

  # -- load hits in the riparian path
  if(verbose)
    cat("Done!\n\tGenerating block coordinates ... ")
  riph <- read_ripHits(
    gsParam = gsParam,
    genomeIDs = genomeIDs,
    useBlks = !plotRegions)
  riph <- subset(riph, ofID1 %in% genesInReg & ofID2 %in% genesInReg)

  # -- add the reference chr and make new blocks therein
  riph[,`:=`(
    gen1 = gv[ofID1], gen2 = gv[ofID2],chr1 = cv[ofID1], chr2 = cv[ofID2],
    refChr = rcv[ofID1], ord1 = ov[ofID1], ord2 = ov[ofID2],
    start1 = sv[ofID1], start2 = sv[ofID2], end1 = sv[ofID1], end2 = sv[ofID2])]
  setkey(riph, gen1, chr1, ord1)
  riph <- subset(riph, complete.cases(riph))
  riph[,rl := add_rle(refChr, which = "id"),
       by = c("gen1", "chr1")]
  riph$refChr[riph$gen1 == refGenome] <- riph$chr1[riph$gen1 == refGenome]
  riph$refChr[riph$gen2 == refGenome] <- riph$chr1[riph$gen2 == refGenome]
  riph$rl[riph$gen1 == refGenome] <- 1
  riph$rl[riph$gen2 == refGenome] <- 1
  riph[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, blkID, refChr, rl)]

  # -- get the block coordinates
  bc <- calc_blkCoords(riph)
  bc[,`:=`(y1 = match(gen1, genomeIDs), y2 = match(gen2, genomeIDs))]
  if(useOrder){
    bc[,`:=`(xs1 = startOrd1, xe1 = endOrd1, xs2 = startOrd2, xe2 = endOrd2)]
  }else{
    bc[,`:=`(xs1 = startBp1, xe1 = endBp1, xs2 = startBp2, xe2 = endBp2)]
  }

  # -- add in the referenc chr and get the colors
  bc[,`:=`(refChr = rcv[firstGene1])]
  bc$refChr[bc$gen1 == refGenome] <- bc$chr1[bc$gen1 == refGenome]
  bc$refChr[bc$gen2 == refGenome] <- bc$chr2[bc$gen2 == refGenome]
  bc[,col := cols[refChr]]

  # -- calculate chromosome start/end coords
  gff[,x := linBp]
  if(useOrder)
    gff[,x := linOrd]
  gff[,y := match(genome, genomeIDs)]
  gff <- subset(gff, complete.cases(gff))
  chrPos <- gff[,list(start = min(x, na.rm = T), end = max(x, na.rm = T)),
                by = c("genome", "chr", "y")]

  ##############################################################################
  # 3. make the plot
  ##############################################################################
  # -- determine plotting scale
  if(verbose)
    cat("Done!\n\tRendering plot ... ")
  pmar <- par()["mar"]
  par(mar = c(1, 1, 1, 1))

  sevenChrWidth <- 0.1111654*nGenomeLabChar
  XinHt <- 0.1193034
  genomeLabs <- substr(genomeIDs, 1, 7)
  pwid <- par("din")[1]
  pht <- par("din")[2]
  wg <- genomeIDs[which.max(nchar(genomeIDs))]
  yPlotBuffProp <- XinHt / pht
  xPlotBuffProp <- sevenChrWidth / pwid
  yBuff <- yPlotBuffProp * length(genomeIDs) * 2
  xBuff <- xPlotBuffProp * diff(range(gff$x, na.rm = T))

  # -- make the plot
  chrPos <- subset(chrPos, complete.cases(chrPos))
  with(chrPos,
       plot(
         1, 1, type = "n", axes = F, bty = "n",
         xlim = c(min(start) - xBuff,
                  max(end)),
         ylim = c(min(y)-yBuff, max(y) + yBuff),
         ylab = "",
         xlab = "",
         main = ""))

  # -- make the labels
  chrRectWidth <- strheight("chr", cex = chrLabCex) * chrRectBuffer
  subchrText <- sprintf(
    "Chromosomes scaled by %s",
    ifelse(useOrder, "gene rank order", "physical position"))
  text(
    x = 0, y = 1 - (chrRectWidth*2),
    label = subchrText, adj = c(.5,.5), cex = genomeLabCex)

  if(blackBg)
    with(chrPos, rect(
      xleft = min(start) - xBuff, xright = max(end),
      ybottom = min(y)-(chrRectWidth/2),  ytop = max(y)+(chrRectWidth/2),
      border = NULL, col = "grey15"))

  # -- make the scale bar
  mo <- with(subset(gff, genome == refGenome), max(ord))
  mo <- ifelse(mo > 1e5, 5e4, ifelse(mo > 1e4, 5e3, ifelse(mo > 1e3, 500, 50)))
  mb <- sum(with(subset(gff, genome == refGenome), tapply(end, chr, max)))
  mb <- ifelse(mb > 1e9, 5e8, ifelse(mb > 1e8, 5e7, ifelse(mb > 1e7, 5e6, ifelse(mb > 1e6, 5e5, ifelse(mb > 1e5, 5e4, 5e3)))))
  if(mo >= 1e3){
    mol <- sprintf("%sk genes", mo/1e3)
  }else{
    mol <- sprintf("%s genes", mo)
  }
  if(mb >= 1e6){
    mbl <- sprintf("%s Mb", mb/1e6)
  }else{
    mbl <- sprintf("%s kb", mb/1e3)
  }

  if(useOrder){
    draw_scaleBar(
      x = quantile(chrPos$start,.5, na.rm = T),
      y = length(genomeIDs) + (chrRectWidth),
      xspan = mo,
      yspan = chrRectWidth/2,
      label = mol,
      lwd = .5, cex = chrLabCex)
    if(is.null(labelChrBiggerThan))
      labelChrBiggerThan <- mo/5
  }else{
    draw_scaleBar(
      x = quantile(chrPos$start, .5, na.rm = T),
      y = length(genomeIDs) + (chrRectWidth * 1.5),
      xspan = mb,
      yspan = chrRectWidth,
      label = mbl,
      lwd = .5, cex = chrLabCex)
    if(is.null(labelChrBiggerThan))
      labelChrBiggerThan <- mb/5
  }

  # -- draw the braid polygons
  polygonList <- lapply(1:nrow(bc), function(i){
    pg <-  with(bc[i,], calc_curvePolygon(
      start1 = xs1, start2 = xs2,
      end1 = xe1, end2 = xe2,
      y1 = y1 + chrRectWidth / 4, y2 = y2 - chrRectWidth / 4))
    if(!is.null(braidBorderLwd)){
      polygon(
        pg,
        border = ifelse(is.null(braidBorderLwd), NA, bc$col[i]),
        lwd = ifelse(is.null(braidBorderLwd), NA, braidBorderLwd),
        col = add_alpha(bc$col[i], braidAlpha))
    }else{
      polygon(
        pg,
        border = NA,
        col = add_alpha(bc$col[i], braidAlpha))
    }
    return(pg)
  })

  # -- draw and lable the chrs
  for (i in 1:nrow(chrPos)) {
    wid <- ifelse(chrPos$genome[i] %in% labelTheseGenomes,
                  chrRectWidth/2, chrRectWidth/4)
    polygon(with(chrPos[i,], round_rect(
      xleft = start, xright = end,
      ybottom = y - wid,
      ytop = y + wid)),
      border = chrBorder,
      col = ifelse(chrPos$genome[i] == refGenome, highlightRef, chrFill),
      lwd = .5)
    if(chrPos$genome[i] %in% labelTheseGenomes)
      if(with(chrPos[i,], end - start) > labelChrBiggerThan)
        with(chrPos[i,], text(
          (start + end)/2, y,
          labels = chrLabFun(chr), cex = chrLabCex))
  }

  # label the genomes
  cp <- chrPos[,list(x = min(start)), by = c("genome","y")]
  if(blackBg){
    with(cp, text(
      x = x - (max(chrPos$end) / 20), y = y, col = "white",
      label = substr(genome, 1, nGenomeLabChar), adj = c(1.2, .5), cex = genomeLabCex))

  }else{
    with(cp, text(
      x = x - (max(chrPos$end) / 20), y = y, col = "black",
      label = substr(genome, 1, nGenomeLabChar), adj = c(1.2, .5), cex = genomeLabCex))
  }

  par(mar = pmar)
  if(verbose)
    cat("Done!\n")
  return(list(
    chrPos = chrPos,
    blockCoord = bc,
    polygonList = polygonList,
    useOrder = useOrder,
    genomeIDs = genomeIDs,
    ripHits = riph,
    gff = gff,
    chrFill = chrFill,
    chrBorder = chrBorder,
    labelChrBiggerThan = labelChrBiggerThan,
    genomeLabCex = genomeLabCex,
    chrLabFun = chrLabFun,
    chrLabCex = chrLabCex,
    chrRectBuffer = chrRectBuffer))
}

#' @title split_blksByRef
#' @description
#' \code{split_blksByRef} split_blksByRef
#' @rdname plot_riparian
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
#' @rdname plot_riparian
#' @export
reorder_gff <- function(gff, genomeIDs, minGenesOnChr, refGenome){
  genome <- ord <- n <- chr <- chrn <- medPos <- chrNameOrd <- synOg <- NULL
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
#' @rdname plot_riparian
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
#' @rdname plot_riparian
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
#' @rdname plot_riparian
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
#' @rdname plot_riparian
#' @export
calc_refChrByGene <- function(gff,
                              refGenome,
                              blkSize,
                              nCores){
  refChr <- chr <- genome <- rl <- ofID <- NULL
  # -- if a synOg has a ref chr in it, populate gff
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


