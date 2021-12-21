#' @title Genespace plotting routines
#' @description
#' \code{plot_riparian} Genespace plotting routines
#'
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
#' @param highlightRef color to highlight the reference chromosomes.
#' @param chrLabFun function to parse chr IDs to make them more readible
#' @param minGenes2plot integer specifying the minimum number of genes on a
#' chr to plot
#' @param onlyTheseRegions data.table with genome, chr, start and end columns
#' @param excludeChrOutOfRegion logical, should chromosome representations
#' be constrained to just those in synteny with the rest of the graph?
#' @param findRegHitsRecursive logical, should regional hit discovery be
#' recursive?
#' @param nGenomeLabChar number of characters for genome label
#' @param annotatePlot logical, should names be printed on the plot
#' @param add2plot logical, should plotting use existing graphics device?
#' @param start1 numeric of length 1, specifying the start coordinate on genome1
#' @param end1 numeric of length 1, specifying the end coordinate on genome1
#' @param start2 numeric of length 1, specifying the start coordinate on genome2
#' @param end2 numeric of length 1, specifying the end coordinate on genome2
#' @param y1 numeric of length 1, specifying the y position of genome1
#' @param y2 numeric of length 1, specifying the y position of genome2
#' @param xleft numeric of length 1, specifying the left x coordinate
#' @param xright numeric of length 1, specifying the right x coordinate
#' @param ybottom numeric of length 1, specifying the bottom y coordinate
#' @param ytop numeric of length 1, specifying the top y coordinate
#' @param verbose should progress updates be printed to the console? If not
#' specified, taken from gsParam.
#' @param returnSourceData logical, should the source data to build the plot
#' be returned?
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
#' @import data.table
#' @importFrom graphics strheight polygon text segments
#' @importFrom stats quantile
#' @importFrom dbscan dbscan frNN
#' @export
plot_riparian <- function(gsParam,
                          plotRegions = TRUE,
                          refGenome = NULL,
                          genomeIDs = NULL,
                          onlyTheseChrs = NULL,
                          onlyTheseRegions = NULL,
                          excludeChrOutOfRegion = FALSE,
                          findRegHitsRecursive = FALSE,
                          invertTheseChrs = NULL,
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
                          returnSourceData = F,
                          verbose = NULL,
                          chrLabFun = function(x)
                            gsub("^0","",gsub("chr|scaf","", gsub("chr|chromosome|scaffold|^lg|_","",tolower(x)))),
                          annotatePlot = TRUE,
                          add2plot = FALSE){

  ##############################################################################
  # -- ad hoc function to draw scale bar on the riparian plot
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

  ##############################################################################
  # -- ad hoc function to
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

  ##############################################################################
  # -- ad hoc function to read riparian path hits
  read_ripHits <- function(gsParam, genomeIDs, useBlks, gff){
    og <-regID <- blkID <- isAnchor <- NULL
    genomeOrd <- data.table(
      gen1 = genomeIDs[-length(genomeIDs)],
      gen2 = genomeIDs[-1],
      y = 1:length(genomeIDs[-1]))

    fs <- with(genomeOrd, file.path(gsParam$paths$results, sprintf(
      "%s_%s_synHits.txt.gz",
      c(gen1, gen2), c(gen2, gen1))))

    fs <- fs[file.exists(fs)]
    hitsRip <- rbindlist(lapply(fs, function(i){
      if(useBlks){
        x <- subset(fread(
          i, select = c("ofID1","ofID2","gen1","gen2","blkID","isAnchor"),
          na.strings = c("NA","")),
          isAnchor & !is.na(blkID))
      }else{
        x <- subset(fread(
          i, select = c("ofID1","ofID2","gen1","gen2","regID","isAnchor"),
          na.strings = c("NA","")),
          isAnchor & !is.na(regID))
        setnames(x, "regID", "blkID")
      }

      if(!paste(x$gen1[1], x$gen2[1]) %in% paste(genomeOrd$gen1, genomeOrd$gen2))
        setnames(x, c(1,2), c("ofID2","ofID1"))

      return(x[,c("ofID1","ofID2","blkID")])
    }))
    ogv <- gff$og; names(ogv) <- gff$ofID
    hitsRip <- subset(hitsRip, ogv[ofID1] == ogv[ofID2])
    hitsRip[,og := ogv[ofID1]]
    return(hitsRip)
  }

  ##############################################################################
  # -- ad hoc function to add syntenic chrs to gff
  add_synChr2gff <- function(gff, refHits, genomeIDs, gapProp, refGenome){
    ord1 <- ord <- gen2 <- genome <- refOrd <- chrOrd <- chr <- end <- gap <- NULL
    start <- linOrd <- linBp <- NULL

    synChrs <- refHits[,list(refOrd = mean(ord1, na.rm = T)), by = c("gen2", "chr2")]
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

  ##############################################################################
  # -- ad hoc function to load reference hits
  load_refHits <- function(gsParam, genomeIDs, refGenome, plotRegions, gff){
    regID <- og <- blkID <- isAnchor <- NULL
    nonRefGen <- genomeIDs[genomeIDs != refGenome]
    fs <- file.path(gsParam$paths$results,
                    sprintf("%s_%s_synHits.txt.gz",
                            c(rep(refGenome, length(nonRefGen)), nonRefGen),
                            c(nonRefGen, rep(refGenome, length(nonRefGen)))))
    fs <- fs[file.exists(fs)]
    hitsRef <- rbindlist(lapply(fs, function(i){
      if(plotRegions){
        x <- subset(fread(
          i, select = c("gen1","ofID1","ofID2","isAnchor","regID"),
          na.strings = c("NA","")),
          isAnchor & !is.na(regID))
        setnames(x, "regID", "blkID")
      }else{
        x <- subset(fread(
          i, select = c("gen1","ofID1","ofID2","isAnchor","blkID"),
          na.strings = c("NA","")),
          isAnchor & !is.na(blkID))
      }
      if(x$gen1[1] != refGenome)
        setnames(x, c("ofID1","ofID2"), c("ofID2","ofID1"))
      x <- x[,c("ofID1","ofID2","blkID")]

      return(x)
    }))
    ogv <- gff$og; names(ogv) <- gff$ofID
    hitsRef <- subset(hitsRef, ogv[ofID1] == ogv[ofID2])
    hitsRef[,og := ogv[ofID1]]
    return(hitsRef)
  }

  ##############################################################################
  # -- ad hoc function to reorder the gff
  reorder_gff <- function(gff, genomeIDs, minGenesOnChr, refGenome){
    genome <- ord <- n <- chr <- chrn <- medPos <- chrNameOrd <- og <- NULL
    gff <- subset(gff, genome %in% genomeIDs)
    gff[,genome := factor(genome, levels = genomeIDs)]
    setkey(gff, genome, ord)
    nGenes <- gff[,list(n = .N, medPos = mean(ord)), by = c("genome","chr")]
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

  arrayID <- og <- synOG <- globOG <- inBlkOG <- ord <- ord2 <-
    genome <- ofID1 <- ofID2 <- chr1 <- chr <- gen1 <- ord1 <- ofID <-
    rl <- refChr <- blkID <- gen2 <- startOrd1 <- endOrd1 <- end <- n <-
    startOrd2 <- endOrd2 <- startBp1 <- endBp1 <- startBp2 <- regID <-
    endBp2 <- firstGene1 <- x <- linBp <- linOrd <- y <- start <- NULL


  ##############################################################################
  # 1. rename a few things, check parameters, read in hits/gff
  ##############################################################################
  # -- specify the ref genome
  if(is.null(refGenome) || !refGenome %in% genomeIDs || length(refGenome) > 1)
    refGenome <- genomeIDs[1]

  # -- genomeID checking
  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
    genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]
  }else{
    tmp <- gsParam$genomes$genomeIDs
    tmp <- tmp[!tmp %in% gsParam$genomes$outgroup]
    genomeIDs <- genomeIDs[genomeIDs %in% tmp]
  }

  # -- check that the refgenome and genomeIDs are OK
  gparIDs <- gsParam$genomes$genomeIDs
  if(!refGenome %in% gparIDs)
    stop(sprintf("The reference genome %s is not in the gsParam genomes\n",
                 refGenome))
  if(!all(genomeIDs %in% gparIDs))
    stop(sprintf("Specified genomeIDs: %s are not in the gsParam genomes\n",
                 paste(genomeIDs[!genomeIDs %in% gparIDs], collapse = ", ")))

  # -- check that all genomeIDs have synHits
  allHits <- sapply(genomeIDs, function(i){
    hi1 <- file.path(gsParam$paths$results,
                    sprintf("%s_%s_synHits.txt.gz", i, genomeIDs))
    hi2 <- file.path(gsParam$paths$results,
                     sprintf("%s_%s_synHits.txt.gz", genomeIDs, i))
    return(all(apply(cbind(hi1, hi2), 1, function(x) any(file.exists(x)))))
  })

  if(!all(allHits))
    stop("Some genomes do not have syntenic hits\n\tdid you use all specified genomeIDs in synteny?\n")


  # -- choose the colors
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

  # -- logical verbose checking
  if(is.null(verbose) || !is.logical(verbose[1])){
    verbose <- gsParam$params$verbose
  }else{
    verbose <- verbose[1]
  }

  # -- get the number of characters to use as the genome labels
  nGenomeLabChar <- min(max(nchar(genomeIDs)), nGenomeLabChar)

  # -- choose which genomes chrs will be labeled
  if(is.null(labelTheseGenomes))
    labelTheseGenomes <- genomeIDs
  labelTheseGenomes <- labelTheseGenomes[!duplicated(labelTheseGenomes)]
  labelTheseGenomes <- labelTheseGenomes[labelTheseGenomes %in% genomeIDs]
  if(length(genomeIDs) < 1)
    labelTheseGenomes <- genomeIDs

  # -- choose which chromosomes will be inverted
  invChrs <- invertTheseChrs
  if(!is.data.frame(invChrs))
    invChrs <- NULL
  if(!is.null(invChrs)){
    if(!all(colnames(invChrs) %in% c("genome", "chr"))){
      invChrs <- NULL
    }else{
      invChrs <- with(invChrs, paste(genome, chr))
    }
  }
  # -- read the gff
  if(verbose)
    cat("\tLoading the gff ... ")
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, showProgress = F, na.strings = c("", "NA"))

  if("refGenome" %in% colnames(gff))
    gff[,refGenome := NULL]

  # invert chromosomes if necessary
  if(!is.null(invChrs)){
    u <- with(gff, unique(paste(genome, chr)))
    invChrs <- invChrs[invChrs %in% u]
    gi <- subset(gff, paste(genome, chr) %in% invChrs)
    if(nrow(gi) > 0){
      ga <- subset(gff, !paste(genome, chr) %in% invChrs)
      setorder(gi, genome, chr, -ord)
      gi[,`:=`(ord = 1:.N,
               start = max(start) - start,
               end = max(end - end)), by = c("genome", "chr")]
      gff <- rbind(gi, ga)
      setkey(gff, genome, chr, ord)
      gff[,ord := 1:.N, by = "genome"]
    }else{
      invChrs <- NULL
    }
  }

  # -- get chromosome and genome vectors from gff
  gv <- gff$genome; cv <- gff$chr
  names(gv) <- names(cv) <- gff$ofID

  #-- read ref hits
  if(verbose)
    cat(sprintf("Done!\n\tMapping genes against %s chromosomes ... ", refGenome))
  refh <- load_refHits(
    gsParam = gsParam,
    gff = gff,
    plotRegions = plotRegions,
    genomeIDs = genomeIDs,
    refGenome = refGenome)
  refh[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
             chr1 = cv[ofID1], chr2 = cv[ofID2])]
  refh[,n := .N, by = c("gen1","gen2","chr1","chr2")]
  refh <- subset(refh, n >= minGenes2plot)

  # -- subset to only the chromosomes that will be plotted
  if(!is.null(onlyTheseChrs)){
    if(all(onlyTheseChrs %in% refh$chr1)){
      refh <- subset(refh, chr1 %in% onlyTheseChrs)
      gcu <- with(refh, unique(c(paste(gen1, chr1), paste(gen2, chr2))))
      gff <- subset(gff, paste(genome, chr) %in% gcu)
    }
  }

  # -- extend graph to include regions syntenic to onlyTheseChrs
  genesInReg <- unique(gff$ofID)
  if(!is.null(onlyTheseRegions) & excludeChrOutOfRegion){
    setkey(gff, genome, chr, start, end)
    setkey(onlyTheseRegions, genome, chr, start, end)
    fo <- foverlaps(gff, onlyTheseRegions)
    fo <- subset(fo, complete.cases(fo[,c("ofID", "start", "end")]))
    genesInReg <- gff$ofID[gff$og %in% unique(fo$og)]
    gff <- subset(gff, ofID %in% genesInReg)
  }

  # -- reorder the gff by the chromosomes/genomes to be plotted
  # -- this is just to get the colors nailed down
  rg <- reorder_gff(
    gff = subset(gff, genome == refGenome),
    minGenesOnChr = minGenes2plot,
    genomeIDs = genomeIDs,
    refGenome = refGenome)

  # -- get the colors
  if(length(colByChrs) != uniqueN(rg$chr)){
    cols <- colorRampPalette(colByChrs)(uniqueN(rg$chr))
  }else{
    cols <- colByChrs
  }
  names(cols) <- unique(rg$chr)

  # -- reorder the gff again
  if(verbose)
    cat("Done!\n\tProjecting linear coordinate system ... ")
  gff <- reorder_gff(
    gff = gff,
    minGenesOnChr = minGenes2plot,
    genomeIDs = genomeIDs,
    refGenome = refGenome)

  # -- pull vectors for plotting
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- add vectors to reference hits
  refh[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]

  # -- get syntenic chromosomes to reference figured out
  gff <- add_synChr2gff(
    gff = data.table(gff),
    refHits = refh,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    gapProp = gapProp)

  # -- get vectors of syntenic positions
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
    gff = gff,
    genomeIDs = genomeIDs,
    useBlks = !plotRegions)

  # -- subset those to the hits in the regions of interest
  riph <- subset(riph, ofID1 %in% genesInReg & ofID2 %in% genesInReg)

  # -- add position vectors
  riph[,`:=`(
    gen1 = gv[ofID1], gen2 = gv[ofID2],chr1 = cv[ofID1], chr2 = cv[ofID2],
    refChr = rcv[ofID1], ord1 = ov[ofID1], ord2 = ov[ofID2],
    start1 = sv[ofID1], start2 = sv[ofID2], end1 = sv[ofID1], end2 = sv[ofID2])]
  setkey(riph, gen1, chr1, ord1)
  riph <- subset(riph, complete.cases(riph))

  # -- if only regions, give new colors and get coords
  if(!is.null(onlyTheseRegions)){
    # -- extend the network across genomes
    setkey(gff, genome, chr, start, end)
    setkey(onlyTheseRegions, genome, chr, start, end)
    onlyTheseRegions[,regID := 1:.N]
    fo <- foverlaps(gff, onlyTheseRegions)
    fo <- subset(fo, complete.cases(fo[,c("ofID", "start", "end")]))
    genesInReg <- gff$ofID[gff$og %in% unique(fo$og)]
    ogInReg <- lapply(split(fo, by = "regID"), function(x) unique(x$og))
    genesInRegList <- lapply(ogInReg, function(x) subset(gff, og %in% x))
    if(is.null(gsParam$params$synteny$synBuff)){
      warning("can't find synteny parameters in the gsParam object\n\t using radius = 100 and blkSize = 5")
      radius <- 100
      blkSize <- 5
    }else{
      radius <- max(gsParam$params$synteny$synBuff)
      blkSize <- max(gsParam$params$synteny$blkSize)
    }

    # -- assign regIDs to these
    riphl <- lapply(genesInRegList, function(x)
      subset(riph, ofID1 %in% x$ofID & riph$ofID2 %in% x$ofID))
    riph <- rbindlist(lapply(1:length(riphl), function(i){
      riphl[[i]][,refChr := i]
      return(riphl[[i]])
    }))

    # -- secondary clustering if needed
    riph[,rl := dbscan(frNN(cbind(ord1, ord2), eps = radius),
                        minPts = blkSize)$cluster,
         by = c("gen1", "gen2","chr1", "chr2", "refChr")]
    riph <- subset(riph, rl > 0)
  }else{
    riph[,rl := add_rle(refChr, which = "id"),
         by = c("gen1", "chr1")]
    riph$refChr[riph$gen1 == refGenome] <- riph$chr1[riph$gen1 == refGenome]
    riph$refChr[riph$gen2 == refGenome] <- riph$chr1[riph$gen2 == refGenome]
    riph$rl[riph$gen1 == refGenome] <- 1
    riph$rl[riph$gen2 == refGenome] <- 1
  }

  # -- add the reference chr and make new blocks therein
  riph[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, blkID, refChr, rl)]
  blkv <- riph$refChr; names(blkv) <- riph$blkID
  blkv <- blkv[!duplicated(names(blkv))]

  # -- get the block coordinates
  bc <- calc_blkCoords(riph)
  bc[,`:=`(y1 = match(gen1, genomeIDs), y2 = match(gen2, genomeIDs))]
  if(useOrder){
    bc[,`:=`(xs1 = startOrd1, xe1 = endOrd1, xs2 = startOrd2, xe2 = endOrd2)]
  }else{
    bc[,`:=`(xs1 = startBp1, xe1 = endBp1, xs2 = startBp2, xe2 = endBp2)]
  }

  # -- add in the referenc chr and get the colors
  bc[,`:=`(refChr = blkv[blkID])]
  if(is.null(onlyTheseRegions)){
    bc$refChr[bc$gen1 == refGenome] <- bc$chr1[bc$gen1 == refGenome]
    bc$refChr[bc$gen2 == refGenome] <- bc$chr2[bc$gen2 == refGenome]
  }
  if(!is.null(onlyTheseRegions)){
    cols <- onlyTheseRegions$col
    names(cols) <- as.character(onlyTheseRegions$regID)
  }
  bc[,col := cols[as.character(refChr)]]
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
  if(!add2plot){
    pmar <- par()["mar"]
    par(mar = c(1, 1, 1, 1))
  }

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
  if(!add2plot){
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
  }

  # -- make the labels
  chrRectWidth <- strheight("chr", cex = chrLabCex) * chrRectBuffer
  subchrText <- sprintf(
    "Chromosomes scaled by %s",
    ifelse(useOrder, "gene rank order", "physical position"))
  if(annotatePlot)
    text(
    x = 0, y = 1 - (chrRectWidth*2),
    label = subchrText, adj = c(.5,.5), cex = genomeLabCex)

  if(!add2plot){
    if(blackBg)
      with(chrPos, rect(
        xleft = min(start) - xBuff, xright = max(end),
        ybottom = min(y)-(chrRectWidth/2),  ytop = max(y)+(chrRectWidth/2),
        border = NULL, col = "grey15"))
  }


  # -- make the scale bar
  if(annotatePlot){
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
        labelChrBiggerThan <- mo/20
    }else{
      draw_scaleBar(
        x = quantile(chrPos$start, .5, na.rm = T),
        y = length(genomeIDs) + (chrRectWidth * 1.5),
        xspan = mb,
        yspan = chrRectWidth,
        label = mbl,
        lwd = .5, cex = chrLabCex)
      if(is.null(labelChrBiggerThan))
        labelChrBiggerThan <- mb/20
    }
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
  if(annotatePlot){
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
            labels = ifelse(paste(genome, chr) %in% invChrs,
                            paste0(chrLabFun(chr), "*"), chrLabFun(chr)),
            cex = chrLabCex))
    }

    # label the genomes
    cp <- chrPos[,list(x = min(start)), by = c("genome","y")]
    if(blackBg){
      with(cp, text(
        x = x - (max(chrPos$end) / 100), y = y, col = "white",
        label = substr(genome, 1, nGenomeLabChar), adj = c(1, .5), cex = genomeLabCex))

    }else{
      with(cp, text(
        x = x - (max(chrPos$end) / 100), y = y, col = "black",
        label = substr(genome, 1, nGenomeLabChar), adj = c(1, .5), cex = genomeLabCex))
    }
  }

  # par(mar = pmar)
  if(verbose)
    cat("Done!\n")
  if(returnSourceData)
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
                              y2){
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
#' @rdname plot_riparian
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


