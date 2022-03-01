#' @title Make riparian plot using hits, not OGs
#' @description
#' \code{plot_riparianHits} Updated and more accurate version of plot_riparian
#'
#' @name plot_riparianHits
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
#' @param chrRectBuffer number, the amount of buffer around the center of chr
#' regions.
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
#' @param minGenes2plot integer specifying the minimum number of genes on a
#' chr to plot
#' @param onlySameChrs logical, should only chromosomes with the same names
#' be shown?
#' @param maxSyntenyFun function, how synteny is calculated relative to the
#' reference chromosome order.
#' @param onlyTheseRegions data.table with genome, chr, start and end columns
#' @param excludeNoRegChr logical, should chromosome representations
#' be constrained to just those in synteny with the rest of the graph?
#' @param excludeNoRegPos logical, should chromosome representations
#' be constrained to just syntenic positions with the rest of the graph?
#' @param nGenomeLabChar number of characters for genome label
#' @param annotatePlot logical, should names be printed on the plot
#' @param add2plot logical, should plotting use existing graphics device?
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
plot_riparianHits <- function(gsParam,
                              useBlks = TRUE,
                              useOrder = TRUE,
                              refGenome = NULL,
                              genomeIDs = NULL,

                              onlyTheseRegions = NULL,
                              excludeNoRegChr = FALSE,
                              excludeNoRegPos = FALSE,
                              onlySameChrs = FALSE,
                              invertTheseChrs = NULL,
                              reorderChrs = TRUE,
                              maxSyntenyFun = function(x) median(x),
                              refChrCols = NULL,

                              minGenes2plot = 50,
                              gapProp = ifelse(useOrder, 0.01, 0.05),

                              braidAlpha = .8,
                              braidBorderLwd = NULL,
                              chrLabCex = .5,
                              nGenomeLabChar = 20,
                              genomeLabCex = .75,
                              chrBorder = "black",
                              chrFill = "white",
                              highlightRef = "white",
                              chrRectBuffer = 1.5,
                              labelChrBiggerThan = NULL,
                              labelTheseGenomes = NULL,
                              chrLabFun = function(x)
                                gsub("^0","",gsub("chr|scaf","", gsub("chr|chromosome|scaffold|^lg|_","",tolower(x)))),

                              blackBg = TRUE,
                              returnSourceData = F,
                              verbose = NULL,
                              annotatePlot = TRUE,
                              add2plot = FALSE){

  read_refHits <- function(synParamsDt,
                           refGenome,
                           genomeIDs,
                           gff,
                           nCores,
                           minGenes,
                           plotRegions){

    hitsFile <- genome1 <- genome2 <- regID <- isRep1 <- isRep2 <- isAnchor <-
      ofID1 <- ofID2 <- regID <- blkID <- n1 <- ofID1 <- n2 <- ofID2 <- NULL

    # -- check that the hits files are in the synParam obj and keep only those
    # that exist.
    if(!"hitsFile" %in% colnames(synParamsDt))
      stop("need to add hits file path to synParams\n")
    synParamsDt <- subset(synParamsDt, file.exists(hitsFile))
    synParamsDt <- subset(synParamsDt, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
    if(nrow(synParamsDt) < 1)
      stop("no synHits files exist\n")

    # -- get vector of orthogroup IDs
    ov <- gff$og; names(ov) <- gff$ofID

    # -- split into two groups, one with the ref as genome1
    altSyn1 <- subset(synParamsDt, genome1 == refGenome)
    altSyn2 <- subset(synParamsDt, genome2 == refGenome & genome1 != refGenome)
    cols <- c("gen1", "gen2", "ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2",
              "blkID", "regID")

    # -- for each group ...
    if(nrow(altSyn1) > 0){
      alt1 <- rbindlist(mclapply(altSyn1$hitsFile, mc.cores = nCores, function(i){

        # -- read in the synHits file
        tmp <- fread(
          i, showProgress = FALSE, na.strings = c("", "NA"),
          select = c(cols, "isRep1", "isRep2", "isAnchor"))

        # -- keep only syntenic orthogroups and array reps
        tmp <- subset(tmp, !is.na(regID) & isRep1 & isRep2 & isAnchor)
        tmp <- subset(tmp, ov[ofID1] == ov[ofID2])

        return(tmp[, cols, with = F])
      }))
    }else{
      alt1 <- NULL
    }
    if(nrow(altSyn2) > 0){
      tmpCols <- c("gen2", "gen1", "ofID2", "ofID1", "chr2", "chr1", "ord2",
                   "ord1", "blkID", "regID")
      alt2 <- rbindlist(mclapply(altSyn2$hitsFile, mc.cores = nCores, function(i){
        tmp <- fread(i, showProgress = FALSE, na.strings = c("", "NA"),
                     select = c(cols, "isRep1", "isRep2", "isAnchor"))
        tmp <- subset(tmp, !is.na(regID) & isRep1 & isRep2 & isAnchor)
        tmp <- subset(tmp, ov[ofID1] == ov[ofID2])
        tmp <- tmp[,tmpCols, with = F]
        setnames(tmp, cols)
        return(tmp)
      }))
    }else{
      alt2 <- NULL
    }

    # -- combine the hits
    hits <- rbind(alt1, alt2)
    if(plotRegions){
      hits <- subset(hits, !is.na(regID))
      hits[,`:=`(blkID = regID, regID = NULL)]
    }else{
      hits <- subset(hits, !is.na(blkID))
      hits[,`:=`(regID = NULL)]
    }

    # -- drop chromosomes without enough hits
    hits[,n1 := uniqueN(ofID1), by = c("gen1", "chr1")]
    hits[,n2 := uniqueN(ofID2), by = c("gen2", "chr2")]
    hits <- subset(hits, n1 >= minGenes & n2 >= minGenes)
    hits[,`:=`(n1 = NULL, n2 = NULL)]
    return(hits)
  }

  prep_gff4rip <- function(gff,
                           useOrder,
                           genomeIDs,
                           chrOrd,
                           gapProp){


    isArrayRep <- pos <- linPos <- minChr <- maxChr <- genome <- chr <-
      gapSize <- isFirstGeneChr <- gap <- x <- maxPos <- y <- NULL

    # -- merge with the gff
    if(useOrder){
      gff <- with(subset(gff, isArrayRep), data.table(
        genome = genome, chr = chr, ofID = ofID, pos = ord))
    }else{
      gff <- with(subset(gff, isArrayRep), data.table(
        genome = genome, chr = chr, ofID = ofID, pos = start))
    }


    # -- convert positions to linear position within chromosome
    gff[,`:=`(minChr = 1, maxChr = 1 + (max(pos) - min(pos))),
        by = c("genome", "chr")]
    gff[,linPos := scale_between(pos, minChr[1], maxChr[1]),
        by = c("genome", "chr")]
    gff[,`:=`(minChr = NULL, maxChr = NULL)]
    gff <- merge(gff, chrOrd, by = c("genome", "chr"))
    gff[,genome := factor(genome, levels = genomeIDs)]
    setkey(gff, genome, chrOrd, pos)

    # -- add gap sizes
    gff[,`:=`(genome = as.character(genome), chr = as.character(chr))]
    gff[,gapSize := gapProp * max(pos)]

    # -- flag first gene on each chr and add in gap size
    gff[,isFirstGeneChr := c(FALSE, chr[-length(chr)] != chr[-1]), by = "genome"]
    gff[,gap := ifelse(isFirstGeneChr, gapSize, 0)]

    # -- add in linear position of previous chr max
    chrmax <- gff[,list(maxPos = max(linPos)), by = c("genome", "chrOrd")]
    chrmax[,chrOrd := chrOrd + 1]
    gff <- merge(gff, chrmax, by = c("genome", "chrOrd"), all.x = T)
    gff$maxPos[is.na(gff$maxPos) | !gff$isFirstGeneChr] <- 0

    # -- sum to get the x position, add y and simplify
    gff[, x := cumsum(maxPos) + cumsum(gap) + linPos, by = "genome"]
    gff[, y := as.numeric(factor(genome, levels = genomeIDs))]
    gff[,`:=`(gapSize = NULL, isFirstGeneChr = NULL,
              gap = NULL, maxPos = NULL)]

    # -- center the x positions for each genome
    gff[,x := x - median(x), by = "genome"]
    setkey(gff, y, x)
    return(gff)
  }

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

  pull_regHits <- function(onlyTheseRegions,
                           gff,
                           synParamsDt,
                           plotRegions,
                           minGenes2plot,
                           nCores){

    gen1 <- gen2 <- chr1 <- ofID1 <- regID <- NULL

    u <- unique(paste(onlyTheseRegions$genome, onlyTheseRegions$chr))
    ug <- unique(onlyTheseRegions$genome)
    sv <- gff$start; names(sv) <- gff$ofID
    ev <- gff$start; names(ev) <- gff$ofID
    regHits <- sapply(ug, USE.NAMES = T, simplify = F, function(i){
      x <- subset(read_refHits(
        synParamsDt = synParamsDt,
        refGenome = i,
        genomeIDs = genomeIDs,
        gff = gff,
        plotRegions = plotRegions,
        minGenes = minGenes2plot,
        nCores = nCores), gen1 != gen2)
      return(subset(x, paste(gen1, chr1) %in% u))
    })
    regHits <- rbindlist(lapply(1:nrow(onlyTheseRegions), function(i){
      x <- data.table(regHits[[onlyTheseRegions$genome[i]]])
      x <- subset(x, chr1 == onlyTheseRegions$chr[i])
      x <- subset(x, ev[ofID1] >= onlyTheseRegions$start[i])
      x <- subset(x, sv[ofID1] <= onlyTheseRegions$end[i])
      x[,regID := onlyTheseRegions$regID[i]]
      return(x)
    }))
    return(regHits)
  }

  ##############################################################################
  # 1. check the basic parameters
  # -- genomeIDs

  hitsFile <- genome1 <- genome2 <- regID <- genome <- chr <- ord <- gen1 <-
    gen2 <- gen1 <- gen2 <- ofID1 <- ofID2 <- ofID <- chrn <- medPos <- ord1 <-
    ord2 <- medOrd <- x1 <- x2 <- regID <- ofID1 <- ofID2 <- refChr1 <-
    refChr2 <- rl <- blkID <- chr1 <- chr2 <- x <- xstart <- NULL

  nCores <- gsParam$params$nCores

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(!all(genomeIDs %in% gsParam$genomes$genomeIDs))
    stop(sprintf(
      "genomeIDs (%s) must all be present in the gsParam genomeIDs (%s)",
      paste(genomeIDs, collapse = ","),
      paste(gsParam$genomes$genomeIDs, collapse = ",")))

  # -- refGenome
  if(is.null(refGenome))
    refGenome <- genomeIDs[1]
  if(!refGenome %in% genomeIDs)
    stop(sprintf(
      "refGenome %s not one of the genomeIDs %s",
      refGenome, paste(genomeIDs, collapse = ",")))

  # -- check that synteny params exist
  synp <- gsParam$params$synteny
  if(!is.data.table(synp))
    stop("Must run set_syntenyParams first!\n")
  synp <- data.table(synp)
  blkSize <- max(synp$blkSize)
  synBuff <- max(synp$selfRegionMask)

  # -- check the synhits exist
  synp[, hitsFile := file.path(
    gsParam$paths$results,
    sprintf("%s_%s_synHits.txt.gz", genome1, genome2))]
  synp <- subset(synp, file.exists(hitsFile))
  if(nrow(synp) < 1)
    stop("can't find synHits accompanying synParams ... has synteny been run?\n")

  # -- get the orthofinder directory if needed
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)
  verbose <- gsParam$params$verbose

  # -- check the gff file
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(!file.exists(gffFile))
    stop("can't find the annotated gff-like text file\t\n ... have you run annotate_gff yet?\n")
  gf <- fread(gffFile, showProgress = F, na.strings = c("NA", ""))
  genomeIDs <- genomeIDs[genomeIDs %in% gf$genome]
  # -- check that the reference is in the gff
  if(!refGenome %in% gf$genome)
    stop(sprintf("%s (specified refGenome) not in the gff. Available genomes are: \n\t%s\n",
                 refGenome, paste(unique(gf$genome), collapse = ",")))

  # -- check and get the regions in order
  if(!is.null(onlyTheseRegions)){
    regs <- data.table(onlyTheseRegions)
    regs[,regID := paste0("reg", 1:.N)]
    if(!"genome" %in% colnames(regs))
      regs[,genome := NA]
    if(!"chr" %in% colnames(regs))
      regs[,chr := NA]
    if(!"start" %in% colnames(regs))
      regs[,start := 0]
    if(!"end" %in% colnames(regs))
      regs[,end := 1e10]
    if(!"col" %in% colnames(regs))
      regs[,col := NA]
    if(any(!are_colors(regs$col)))
      regs[,col := NA]
    if(any(is.na(regs$col))){
      if(blackBg)
        tmp <- c("#C4645C", "#F5915C", "#FFC765", "#FCF8D5","#BEF3F9",
                 "#66B8FF", "#6666FF", "#9C63E1", "#F4BDFF")
      if(!blackBg)
        tmp <- c("#62322E", "#C60000", "#FF7500", "#FEDF99", "#BEF3F9",
                 "#43B8FF", "#204DBF", "#9C63E1", "#F4BDFF")
      regs[,col := colorRampPalette(tmp)(nrow(regs))]
    }
    regs <- subset(regs, complete.cases(regs[,c("genome", "chr", "start", "end")]))
    uchr <- with(gf, unique(paste(genome, chr)))
    if(!all(paste(regs$genome, regs$chr) %in% uchr))
      warning("some onlyTheseRegions not in gff .. entry rows:",
              which(!paste(regs$genome, regs$chr) %in% uchr))

    regs <- subset(regs, end > start & paste(genome, chr) %in% uchr)
    if(nrow(regs) < 1)
      stop("onlyTheseRegions specified but not correct format - must be a data.table where columns:\n\tgenome (exact matches to genomeIDs)\n\tchr (chr in that genome)\n\tstart (if specified, numeric within 0:max bp position on chr)\n\tend(if specified, same as start, but must be > start)\n\tcol (if specified a vector of colors, if this is missing or miss-specified, uses default palette)")
  }else{
    excludeNoRegChr <- FALSE
    excludeNoRegPos <- FALSE
    regs <- NULL
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

  # -- choose which chrs can be labeled
  if(is.null(labelChrBiggerThan))
    labelChrBiggerThan <- NA
  labelChrBiggerThan <- as.numeric(labelChrBiggerThan[1])
  if(is.null(labelChrBiggerThan) || is.na(labelChrBiggerThan)){
    if(useOrder)
      labelChrBiggerThan <- 50
    if(!useOrder)
      labelChrBiggerThan <- 5e6
  }

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

  # -- invert the chromosomes if necessary
  if(!is.null(invChrs)){
    u <- with(gf, unique(paste(genome, chr)))
    invChrs <- invChrs[invChrs %in% u]
    gi <- subset(gf, paste(genome, chr) %in% invChrs)
    if(nrow(gi) > 0){
      ga <- subset(gf, !paste(genome, chr) %in% invChrs)
      setorder(gi, genome, chr, -ord)
      gi[,`:=`(ord = 1:.N,
               start = max(start) - start,
               end = max(end - end)), by = c("genome", "chr")]
      gf <- rbind(gi, ga)
      setkey(gf, genome, chr, ord)
      gf[,ord := 1:.N, by = "genome"]
    }else{
      invChrs <- NULL
    }
  }

  # -- make index of genomes
  genomeOrd <- data.table(
    gen1 = genomeIDs[-length(genomeIDs)],
    gen2 = genomeIDs[-1],
    y = 1:length(genomeIDs[-1]))


  # -- get reference genome colors in order
  if(!is.null(refChrCols) & is.null(regs)){
    refChrCols <- refChrCols[are_colors(refChrCols)]
    if(length(refChrCols) == 0)
      refChrCols <- NULL
  }
  if(is.null(refChrCols) & is.null(regs)){
    if(blackBg)
      refChrCols <- c("#C4645C", "#F5915C", "#FFC765", "#FCF8D5","#BEF3F9",
                      "#66B8FF", "#6666FF", "#9C63E1", "#F4BDFF")
    if(!blackBg)
      refChrCols <- c("#62322E", "#C60000", "#FF7500", "#FEDF99", "#BEF3F9",
                      "#43B8FF", "#204DBF", "#9C63E1", "#F4BDFF")
  }

  # -- determine if we are coloring by a reference
  colorByRefChr <- is.null(regs) & length(refChrCols) > 1

  ##############################################################################
  # 1. Read in the hits
  # -- read reference hits
  if(verbose)
    cat("\tReading hits ... ")
  refHits <- subset(read_refHits(
    synParamsDt = synp,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    gff = gf,
    plotRegions = !useBlks,
    minGenes = minGenes2plot,
    nCores = nCores), gen1 != gen2)

  # -- read hits for each set of riparian links
  ripHits <- rbindlist(lapply(1:nrow(genomeOrd), function(i)
    subset(read_refHits(
      synParamsDt = synp,
      refGenome = genomeOrd$gen1[i],
      genomeIDs = c(genomeOrd$gen1[i], genomeOrd$gen2[i]),
      gff = gf,
      minGenes = minGenes2plot,
      plotRegions = !useBlks,
      nCores = nCores), gen1 != gen2)))
  if(verbose)
    cat("Done!\n")

  ##############################################################################
  # 2. If using specific regions, subset here
  if(!is.null(regs)){
    if(verbose)
      cat(sprintf("\tBuilding database of hits in %s regions ... ", nrow(regs)))
    regHits <- pull_regHits(
      onlyTheseRegions = regs,
      gff = gf,
      synParamsDt = synp,
      plotRegions = !useBlks,
      minGenes2plot = minGenes2plot,
      nCores = nCores)
    u <- with(regHits, unique(c(ofID1, ofID2)))
    if(verbose)
      cat("Done!\n")
  }else{
    regHits <- NULL
  }

  if(verbose)
    cat("\tGenerating plot coordinates ... ")

  # -- if excluding chrs out of the region, do that here
  if(excludeNoRegChr){
    u <- with(regHits, unique(c(ofID1, ofID2)))
    ripHits <- subset(ripHits, ofID1 %in% u & ofID2 %in% u)
    refHits <- subset(refHits, ofID1 %in% u & ofID2 %in% u)
    u <- with(ripHits, unique(paste(c(gen1, gen2), c(chr1, chr2))))
    gf <- subset(gf, paste(genome, chr) %in% u)
  }

  # -- if excluding positions out of the region, do that here
  if(excludeNoRegPos){
    u <- with(regHits, unique(c(ofID1, ofID2)))
    ripHits <- subset(ripHits, ofID1 %in% u & ofID2 %in% u)
    refHits <- subset(refHits, ofID1 %in% u & ofID2 %in% u)
    u <- with(ripHits, unique(c(ofID1, ofID2)))
    gf <- subset(gf, ofID %in% u)
  }

  ##############################################################################
  # 2. get linear positions of all genes
  # -- median position of each genome/chr against the reference
  gf[,genome := factor(genome, levels = genomeIDs)]
  setkey(gf, genome, ord)
  nGenes <- gf[,list(n = .N, medPos = mean(ord)), by = c("genome","chr")]
  nGenes[,chrn := as.numeric(gsub("[^0-9]", "", chr))]
  nGenes$chrn[is.na(nGenes$chrn)] <- 0
  setorder(nGenes, genome, chrn, -n, chr, medPos)
  chrord <- with(nGenes, paste(genome, chr))
  gf[,u := factor(paste(genome, chr), levels = chrord)]
  setkey(gf, genome, u, ord)
  gf[,ord := 1:.N, by = "genome"]
  ov <- gf$ord; names(ov) <- gf$ofID
  refHits[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]
  ripHits[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]
  if(reorderChrs){
    medChr <- refHits[,list(med = maxSyntenyFun(as.numeric(ord1))), by = c("gen2", "chr2")]
  }else{
    medChr <- refHits[,list(med = maxSyntenyFun(as.numeric(ord2))), by = c("gen2", "chr2")]
  }
  setnames(medChr, c("genome", "chr", "medOrd"))
  refChr <- refHits[,list(med = median(ord1)), by = c("gen1", "chr1")]
  setnames(refChr, c("genome", "chr", "medOrd"))
  chrOrd <- rbind(refChr, medChr)
  setkey(chrOrd, genome, medOrd)
  chrOrd[,chrOrd := 1:.N, by = "genome"]
  chrOrd[,medOrd := NULL]

  # -- convert gff to linear positions
  gf <- prep_gff4rip(
    gff = gf,
    chrOrd = chrOrd,
    useOrder = useOrder,
    genomeIDs = genomeIDs,
    gapProp = gapProp)
  pv <- gf$x; names(pv) <- gf$ofID

  ripHits[,`:=`(x1 = pv[ofID1], x2 = pv[ofID2], ord1 = NULL, ord2 = NULL)]
  ripHits[,`:=`(start1 = x1, start2 = x2, end1 = x1, end2 = x2, ord1 = x1,
                ord2 = x2, isAnchor = TRUE, arrayOrd1 = x1, arrayOrd2 = x2)]

  if(!is.null(regs)){
    bc <- rbindlist(lapply(1:nrow(regs), function(i){
      tmp <- subset(regHits, regID == regs$regID[i])
      tmp <- unique(c(tmp$ofID1, tmp$ofID2))
      ripreg <- subset(ripHits, ofID1 %in% tmp & ofID2 %in% tmp)
      xbc <- calc_blkCoords(ripreg)
      xbc[,regID := regs$regID[i]]
      xbc[,col := regs$col[i]]
      return(xbc)
    }))
  }else{
    tmp1 <- with(refHits, data.table(
      ofID1 = c(ofID1, ofID2), refChr1 = c(chr1, chr1)))
    tmp2 <- with(refHits, data.table(
      ofID2 = c(ofID1, ofID2), refChr2 = c(chr1, chr1)))
    if(colorByRefChr){
      tmp1 <- subset(tmp1, !duplicated(tmp1))
      tmp2 <- subset(tmp2, !duplicated(tmp2))
      ripHits <- merge(ripHits, tmp1, by = "ofID1", allow.cartesian = T)
      ripHits <- subset(ripHits, !duplicated(ripHits))
      ripHits <- merge(ripHits, tmp2, by = "ofID2", allow.cartesian = T)
      ripHits <- subset(ripHits, !duplicated(ripHits))
      ripHits <- subset(ripHits, refChr1 == refChr2)
      ripHits[,`:=`(refChr = refChr1, refChr1 = NULL, refChr2 = NULL)]
      ripHits <- subset(ripHits, complete.cases(ripHits))
      setkey(ripHits, gen1, ord1)
      ripHits[,rl := add_rle(refChr, which = "id"), by = c("gen1", "chr1", "blkID")]
      ripHits[,blkID := as.numeric(as.factor(paste(gen1, gen2, chr1, chr2, rl, refChr, blkID)))]
      ripHits[,rl := NULL]
      ripHits[,blkID := sprintf("%s_%s", refChr, blkID)]

      tmp <- subset(gf, genome == refGenome & chr %in% ripHits$refChr)
      refChrs <- unique(tmp$chr)

      print(refChrs)
      urchr <- unique(gf$chr[gf$genome == refGenome])
      if(length(refChrCols) != length(refChrs)){
        cols <- colorRampPalette(refChrCols)(length(refChrs))
      }else{
        cols <- refChrCols
      }
      names(cols) <- refChrs
    }else{
      ripHits[,refChr := 1]
      ripHits[,blkID := as.numeric(as.factor(paste(gen1, gen2, chr1, chr2, refChr, blkID)))]
      ripHits[,blkID := sprintf("%s_%s", refChr, blkID)]
      cols <- refChrCols[1]
      names(cols) <- "1"
    }

    bc <- calc_blkCoords(ripHits)
    bc[,c("refChr", "blkID") := tstrsplit(blkID, "_")]
    bc[,col := cols[refChr]]
  }

  bc <- with(bc, data.table(
    xstart1 = startBp1, xend1 = endBp1, xstart2 = startBp2, xend2 = endBp2,
    gen1 = gen1, gen2 = gen2, chr1 = chr1, chr2 = chr2,
    blkID = blkID, col = col))

  bc[,`:=`(y1 = as.numeric(factor(gen1, levels = genomeIDs)),
           y2 = as.numeric(factor(gen2, levels = genomeIDs)))]

  ##############################################################################
  # 3. make plotting data
  # -- make block coordinates
  if(onlySameChrs){
    ripHits <- subset(ripHits, chr1 == chr2)
    warning(sprintf("\nsubsetting to only the same chrs resulted in %s hits\n\t",
                    nrow(ripHits)))
  }

  # -- get chromosome coordinates
  chrPos <- gf[,list(xstart = min(x), xend = max(x)),
                by = c("genome", "chr", "y")]

  # -- get the colors
  highlightRef <- highlightRef[1]
  if(!are_colors(highlightRef) || is.null(highlightRef))
    highlightRef <- "white"

  bc <- subset(bc, complete.cases(bc))
  bc <- subset(bc, !duplicated(bc))
  ##############################################################################
  # 4. make the plot
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
  xBuff <- xPlotBuffProp * diff(range(gf$x, na.rm = T))

  # -- make the plot
  if(!add2plot){
    with(gf,
         plot(
           1, 1, type = "n", axes = F, bty = "n",
           xlim = c(min(x) - xBuff,
                    max(x)),
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
        xleft = min(xstart) - xBuff, xright = max(xend),
        ybottom = min(y)-(chrRectWidth/2),  ytop = max(y)+(chrRectWidth/2),
        border = NULL, col = "grey15"))
  }

  # -- make the scale bar
  if(annotatePlot){
    if(useOrder){
      n <- max(gf$pos)
      mo <- ifelse(n > 1e5, 20e3,
                   ifelse(n > 50e3, 10e3,
                          ifelse(n > 20e3, 5e3,
                                 ifelse(n > 8e3, 2e3,
                                        ifelse(n > 2e3, 500,
                                               ifelse(n > 500, 100,
                                                      ifelse(n > 100, 50, 10)))))))
    }else{
      n <- max(gf$pos)
      mo <- ifelse(n > 1e9, 5e8,
                   ifelse(n > 1e8, 5e7,
                          ifelse(n > 1e7, 5e6,
                                 ifelse(n > 1e6, 5e5,
                                        ifelse(n > 1e5, 5e4, 5e3)))))
    }

    if(mo >= 1e3 & useOrder)
      mol <- sprintf("%sk genes", mo/1e3)
    if(mo < 1e3 & useOrder)
      mol <- sprintf("%s genes", mo)
    if(mo >= 1e6 & !useOrder)
      mol <- sprintf("%s Mb", mo/1e6)
    if(mo <= 1e6 & !useOrder)
      mol <- sprintf("%s kb", mo/1e3)

    draw_scaleBar(
      x = quantile(gf$x, .5, na.rm = T),
      y = length(genomeIDs) + (chrRectWidth),
      xspan = mo,
      yspan = chrRectWidth/2,
      label = mol,
      lwd = .5,
      cex = chrLabCex)
  }
  # -- draw the braid polygons
  # setkey(bc, col)

  polygonList <- lapply(1:nrow(bc), function(i){
    pg <-  with(bc[i,], calc_curvePolygon(
      start1 = xstart1, start2 = xstart2,
      end1 = xend1, end2 = xend2,
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
        xleft = xstart, xright = xend,
        ybottom = y - wid,
        ytop = y + wid)),
        border = chrBorder,
        col = ifelse(chrPos$genome[i] == refGenome, highlightRef, chrFill),
        lwd = .5)
      if(chrPos$genome[i] %in% labelTheseGenomes)
        if(with(chrPos[i,], xend - xstart) > labelChrBiggerThan)
          with(chrPos[i,], text(
            (xstart + xend)/2, y,
            labels = ifelse(paste(genome, chr) %in% invChrs,
                            paste0(chrLabFun(chr), "*"), chrLabFun(chr)),
            cex = chrLabCex))
    }

    # label the genomes
    cp <- chrPos[,list(x = min(xstart)), by = c("genome","y")]
    if(blackBg){
      with(cp, text(
        x = x - (max(chrPos$xend) / 100), y = y, col = "white",
        label = substr(genome, 1, nGenomeLabChar), adj = c(1, .5), cex = genomeLabCex))

    }else{
      with(cp, text(
        x = x - (max(chrPos$xend) / 100), y = y, col = "black",
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
      ripHits = ripHits,
      refHits = refHits,
      gff = gf,
      chrFill = chrFill,
      chrBorder = chrBorder,
      labelChrBiggerThan = labelChrBiggerThan,
      genomeLabCex = genomeLabCex,
      chrLabFun = chrLabFun,
      chrLabCex = chrLabCex,
      chrRectBuffer = chrRectBuffer))
}
