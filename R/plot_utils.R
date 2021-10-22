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
#' @param radius numeric of length 1 specifying the eps dbscan parameter; the
#' search radius within which to count clustered density-based xy points.
#' @param blkSize integer of length 1 specifying the minimum size for a syntenic
#' block and the -s 'size' MCScanX parameter
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#' @param minRbhScore integer of length 1, see set_syntenyParams
#' @param genome1 character string specifying first of two genomeIDs
#' @param genome2 character string specifying second of two genomeIDs
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param synBuff integer of length 1 specifying the maximum euclidean distance
#' from an 'anchor' so that it can be considered syntenic
#' @param selfOnly logical, should only self hits be considered
#' @param overwrite logical, should the results be overwrittem?
#' @param maskHits data.table of hits that should be excluded
#' @param synParam data.table with synteny parameters. See set_syntenyParams.
#' @param selfRegionMask integer, the radius around self hits that should be
#' masked
#' @param type either 'primary' or 'secondary' depending on the scale of
#' inference
#' @param dropInterleavesSmallerThan integer, the minimum block size to retain
#' after splitting overlapping blocks
#' @param minPropDup numeric (0-1) specifying the minimum proportion of
#' duplicated hits to allow two overlapping blocks to not be split
#' @param maxIter integer, the maximum number of block splitting interations
#' @param nhits integer, the number of hits to retain
#' @param blks data.table containing the block coordinates
#' @param allowRBHinOG logical, for cross compatibility with plot_hits
#' @param useBlks logical, for cross compatibility with plot_hits
#' @param verbose logical, should updates be printed to the console?
#' @param refGenome character string matching one of the genomeIDs in gsParam

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
