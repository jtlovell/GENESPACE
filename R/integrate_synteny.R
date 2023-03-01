#' @title Integrate syntenic positions across multiple genomes
#' @name integrate_synteny
#' @description
#' \code{integrate_synteny} Internal functions to connect syntenic positions for
#' riparian plotting and interpolated gene order calculation for placement of
#' syntenic pan-genes.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#' @param refGenome Character string specifying the reference genome to use
#' @param useRegions Logical, should regions be used instead of blocks?
#' @param blkSize see init_genesapce
#' @param synBuff see init_genesapce
#' @param refChr string matching chromosome IDs in refGenome
#' @param refStartBp physical (bp) region start coordinate on refChr
#' @param refEndBp physical (bp) region end coordinate on refChr
#' @param mirror logical, should the block coordinates be mirrored between
#' query and target genomes
#' @param hits data.table with hits (see read_synHits)
#' @param verbose logical, should updates be printed to the console?

#' @details Functions here should not be called directly by the user except
#' if customizing riparian plots.

#' @title linear interpolation of the entire genespace run
#' @description
#' \code{interp_synPos} pipes data to interp_approx and combines across genomes
#' @rdname integrate_synteny
#' @import data.table
#' @export
interp_synPos <- function(gsParam, verbose = TRUE){
  ##############################################################################
  ##############################################################################
  # -- add hoc function to get all genes in a block

  ofID2 <- ofID <- n <- blkID <- nGenes <- noAnchor <- isArrayRep <-
    isAnchor <- ord1 <- ord2 <- ofID1 <- ofID2 <- blkID <- refGenome <- NULL
  ##############################################################################
  ##############################################################################
  # -- add hoc function to summarize syntenic mapping positions
  summarize_synPos <- function(synPos, verbose){
    cnts <- subset(synPos, !is.na(ofID2))
    n0 <- subset(bed, !ofID %in% synPos$ofID2)
    cnts[,n := uniqueN(blkID), by = "ofID2"]
    tab <- cnts[,list(nGenes = uniqueN(ofID2)), by = c("genome2", "n")]
    n0 <- n0[,list(nGenes = .N), by = c("genome")]
    n0[,n := 0]
    setnames(tab, "genome2", "genome")
    tab <- rbind(tab, n0)
    tab[,n := factor(ifelse(n <= 2, n, 3), levels = 0:3)]

    tab <- dcast(
      tab, n ~ genome, value.var = "nGenes", fill = 0, drop = FALSE,
      fun.aggregate = sum)
    tab <- melt(
      tab, id.vars = "n", variable.name = "genome", value.name = "nGenes")
    labs <- align_charLeft(genomeIDs)
    names(labs) <- genomeIDs

    tab <- tab[,list(nGenes = sum(nGenes)), by = c("genome", "n")]
    out <- data.table(tab)

    if(verbose){
      tab[,nGenes := align_charRight(nGenes)]
      spl <- split(tab, by = "genome")
      for(j in names(spl)){
        with(spl[[j]], cat(sprintf(
          "\t\t%s: %s / %s / %s / %s\n",
          labs[j], nGenes[n == 0], nGenes[n == 1],
          nGenes[n == 2], nGenes[n == 3])))
      }
    }
    setnames(out, "n", "nSynPos")
    return(out)
  }

  ##############################################################################

  ##############################################################################
  ##############################################################################

  # ... read in bedfile
  genomeIDs <- gsParam$genomeIDs
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  bed <- subset(bed, !noAnchor & isArrayRep)
  ##############################################################################
  # ... for each genome ID
  cnts <- rbindlist(lapply(genomeIDs, function(i){
    if(verbose)
      cat(sprintf("\t%s:  (0 / 1 / 2 / >2 syntenic positions)\n", i))
    ############################################################################
    # linear interpolation of positions
    # -- 1.1 read in syntenic hits
    hits <- read_refGenomeSynHits(
      gsParam = gsParam, refGenome = i)

    # -- 1.2 subset to syntenic anchor hits
    anchors <- subset(hits, isAnchor)

    # -- 1.3 get hits within block coordinates
    bedInBlk <- get_bedInBlk(hits = anchors, bed = bed)

    # -- 1.4 Do interpolation for each block
    x1 <- subset(bedInBlk, !is.na(ord1))
    setorder(x1, ord1, na.last = T)
    x1[,`:=`(ord2 = interp_approx(x = ord1, y = ord2)), by = "blkID"]

    x2 <- subset(bedInBlk, !is.na(ord2))
    setorder(x2, ord1, na.last = T)
    x2[,`:=`(ord1 = interp_approx(x = ord2, y = ord1)), by = "blkID"]

    # -- 1.5 combine
    interpInBlk <- rbind(x1, x2)
    interpInBlk <- subset(interpInBlk, !duplicated(paste(ofID1, ofID2, blkID)))

    # -- 1.6 write to file
    outf <- file.path(gsParam$paths$pangenes,
                      sprintf("%s_integratedSynPos.txt", i))
    write_intSynPos(interpInBlk, filepath = outf)

    # -- 1.7 summarize interpolated positions
    md <- summarize_synPos(interpInBlk, verbose = verbose)
    md[,refGenome := i]
    return(md)
  }))

  return(cnts)
}

#' @title phase syntenic blocks by a reference
#' @description
#' \code{phase_blks} splits blocks so that each is anchored to a single
#' reference genome chromosome, which allows the construction of default
#' riparian plots
#' @rdname integrate_synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
phase_blks <- function(gsParam,
                       refGenome,
                       useRegions,
                       blkSize = 5,
                       synBuff = 100,
                       refChr = NULL,
                       refStartBp = NULL,
                       refEndBp = NULL){

  hogID <- OG <- blkID <- regID <- chr1 <- end1 <- start1 <- isAnchor <-
    blkID <- regID <- refChr1 <- refChr2 <- n <- ord1 <- ord2 <- clus <-
    tord1 <- tord2 <- NULL
  # ... read in bedfile
  genomeIDs <- gsParam$genomeIDs
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  nCores <- gsParam$params$nCores

  if(length(refChr) > 1 || length(refEndBp) > 1 || length(refStartBp) > 1)
    stop("if specified, refChr, refStartBp & refEndBp must all have length 1\n")
  if(sum(is.null(refChr), is.null(refStartBp), is.null(refEndBp)) == 3){
    inReg <- FALSE
  }else{
    if(sum(is.null(refChr), is.null(refStartBp), is.null(refEndBp)) == 0){
      inReg <- TRUE
    }else{
      stop("if specified, refChr, refStartBp & refEndBp must all have length 1\n")
    }
  }

  refHits <- read_refGenomeSynHits(
    gsParam = gsParam, refGenome = refGenome)

  if(!all(is.na(refHits$regID)) && useRegions)
    refHits[,blkID := regID]

  if(inReg){
    refHits <- subset(
      refHits, chr1 == refChr & end1 >= refStartBp & start1 <= refEndBp)
  }

  rh <- with(subset(refHits, !is.na(blkID)), data.table(
    ofID = c(ofID1, ofID2), refChr = c(chr1, chr1)))
  rh <- subset(rh, complete.cases(rh))
  rh <- subset(rh, !duplicated(rh))

  rh1 <- data.table(rh); rh2 <- data.table(rh)
  setnames(rh1, c("ofID1", "refChr1"))
  setnames(rh2, c("ofID2", "refChr2"))

  synhitFiles <- gsParam$synteny$blast$synHits
  hcols <- c("ofID1", "ofID2", "genome1", "genome2", "chr1", "chr2", "ord1", "ord2", "start1", "end1", "start2", "end2","blkID")
  refBlks <- rbindlist(mclapply(synhitFiles, mc.cores = nCores, function(j){

    # -- read in the hits
    h <- subset(read_synHits(j), isAnchor)
    if(!all(is.na(h$regID)) && useRegions){
      h[,blkID := regID]
    }

    h <- h[,hcols, with = F]
    h <- subset(h, complete.cases(h))

    # -- merge with reference hit chromosome info
    h <- merge(rh2, h, by = "ofID2", allow.cartesian = T, all.y = T)
    h <- merge(rh1, h, by = "ofID1", allow.cartesian = T, all.y = T)

    # -- re-calculate syntenic blocks
    anchs <- subset(h, refChr1 == refChr2)
    anchs <- subset(anchs, complete.cases(anchs))
    anchs[,n := .N, by = c("chr1", "chr2", "blkID", "refChr1")]
    anchs[,`:=`(tord1 = frank(ord1, ties.method = "dense"),
                tord2 = frank(ord2, ties.method = "dense")),
          by = "blkID"]
    if(nrow(anchs) < blkSize){
      bc <- NULL
    }else{
      anchs <- subset(anchs, n >= blkSize)
      anchs[,clus := dbscan(frNN(
        x = cbind(tord1, tord2),
        eps = synBuff),
        minPts = blkSize)$cluster,
        by = c("chr1", "chr2", "blkID", "refChr1")]
      anchs[,blkID := sprintf("%s=refChr_XXX_%s_%s",refChr1, blkID, clus)]

      bc <- calc_blkCoords(anchs, mirror = T)
      bc[,c("refChr", "blkID") := tstrsplit(blkID, "=refChr_XXX_")]
    }
    return(bc)
  }))
  refBlks[,refGenome := refGenome]
  refBlks$refGenome[is.na(refBlks$refChr)] <- NA

  return(refBlks)
}

#' @title calculate block coordinates
#' @description
#' \code{nophase_blks} simple wrapper for calc_blkCoords so that block
#' coordinates are calculated across all combination of genomes.
#' @rdname integrate_synteny
#' @import data.table
#' @importFrom parallel mclapply
#' @export
nophase_blks <- function(gsParam,
                         useRegions,
                         blkSize = 5){

  isAnchor <- regID <- blkID <- n <- NULL
  synhitFiles <- gsParam$synteny$blast$synHits
  nCores <- gsParam$params$nCores
  hcols <- c("ofID1", "ofID2", "genome1", "genome2", "chr1", "chr2", "ord1", "ord2", "start1", "end1", "start2", "end2","blkID")
  allBlks <- rbindlist(mclapply(synhitFiles, mc.cores = nCores, function(j){

    # -- read in the hits
    h <- subset(read_synHits(j), isAnchor)
    if(!all(is.na(h$regID)) && useRegions){
      h[,blkID := regID]
    }

    h <- h[,hcols, with = F]
    anchs <- subset(h, complete.cases(h))
    anchs[,n := .N, by = c("chr1", "chr2", "blkID")]
    anchs <- subset(anchs, n >= blkSize)

    bc <- calc_blkCoords(anchs, mirror = T)
    return(bc)
  }))

  return(allBlks)
}

#' @title calculate syntenic block coordinates
#' @description
#' \code{calc_blkCoords} from a hits object, determine block coordinates,
#' orientation and membership
#' @rdname integrate_synteny
#' @import data.table
#' @importFrom stats cor
#' @export
calc_blkCoords <- function(hits, mirror = FALSE){
  setDTthreads(1)

  # -- get the columns and complete observations for these
  hcols <- c("blkID", "start1", "start2", "end1", "end2", "ord1", "ord2",
             "chr1", "chr2", "genome1", "genome2", "ofID1", "ofID2")
  bhits <- subset(hits, complete.cases(hits[,hcols, with = F]))

  if(mirror){
    tmp <- data.table(bhits)
    setnames(tmp, gsub("2$", "3", colnames(tmp)))
    setnames(tmp, gsub("1$", "2", colnames(tmp)))
    setnames(tmp, gsub("3$", "1", colnames(tmp)))
    bhits <- rbind(bhits, tmp[,colnames(bhits), with = F])
    bhits <- subset(bhits, !duplicated(paste(ofID1, ofID2, blkID)))
  }

  # -- get the genome1 coordinates
  ofID1 <- start1 <- end1 <- ofID1 <- ord1 <- blkID <- NULL
  setkey(bhits, ord1)
  blks1 <- bhits[,list(
    startBp1 = min(start1), endBp1 = max(end1),
    startOrd1 = min(ord1), endOrd1 = max(ord1),
    firstGene1 = first(ofID1), lastGene1 = last(ofID1),
    nHits1 = uniqueN(ofID1)),
    by = c("blkID", "genome1","genome2", "chr1", "chr2")]

  # -- get the genome2 coordinates
  ofID2 <- start2 <- end2 <- ofID2 <- ord2 <- NULL
  setkey(bhits, ord2)
  blks2 <- bhits[,list(
    minBp2 = min(start2), maxBp2 = max(end2),
    minOrd2 = min(ord2), maxOrd2 = max(ord2),
    minGene2 = first(ofID2), maxGene2 = last(ofID2),
    nHits2 = uniqueN(ofID2),
    orient = ifelse(length(ord1) <= 1, "+",
                    ifelse(cor(jitter(ord1),
                               jitter(ord2)) > 0,"+", "-"))),
    by = c("blkID", "genome1","genome2", "chr1", "chr2")]

  # -- merge the two coordinates
  blks <- merge(blks1, blks2, by = c("genome1","genome2","chr1","chr2","blkID"))

  # -- fix the coordinates for inverted blocks
  orient <- NULL
  bgfor <- subset(blks, orient == "+")
  bgrev <- subset(blks, orient == "-")

  maxBp2 <- minBp2 <- maxOrd2 <- minOrd2 <- maxGene2 <- minGene2 <- NULL
  bgrev[,`:=`(startBp2 = maxBp2, endBp2 = minBp2,
              startOrd2 = maxOrd2, endOrd2 = minOrd2,
              firstGene2 = maxGene2, lastGene2 = minGene2)]
  bgfor[,`:=`(startBp2 = minBp2, endBp2 = maxBp2,
              startOrd2 = minOrd2, endOrd2 = maxOrd2,
              firstGene2 = minGene2, lastGene2 = maxGene2)]
  blks <- rbind(bgfor, bgrev)
  return(blks)
}
