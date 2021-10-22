#' @title Parsing of pairwise hits into synteny
#' @description
#' \code{synteny} The main GENESPACE engine to constrain genome-wide homology
#' hits to synteny.
#'
#' @name synteny
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
#' @param verbose logical, should updates be printed to the console?

#' @details The main engine for GENESPACE synteny searching. This
#' finds syntenic 'anchors' that are high-confidence synteny- and homology-
#' constrained hits, then pulls nearby hits within a specified buffer around
#' 'anchor' hits between two genomes. Combined, this provides a framework to
#' both analyze syntenic duplicates (e.g. tandem arrays) and have high
#' confidence that only the desired types of hits (orthologs, homoeologs, etc.)
#' are considered.
#'
#' The basic premise is that we can find synteny in haploid genome comparisons (most
#' diploid species have haploid genome representations) best by first
#' constraining the initial search to the single best scoring pairwise hits for
#' each gene. Then, if desired further subsetting this set to only gene pairs
#' that fall into the same orthogroups. This approach effectively removes
#' relic syntenic regions from ancient whole genome duplications and produces
#' a set of high-confidence synteny "anchors" which serve as known positions of
#' syntenic (ortho/para/homeolo)gous regions. We then search in a fixed-radius
#' for nearest neighbors within a gene-rank order buffer around the anchors. For
#' intra-genomic hits, the self hits are the anchors and a buffer is calculated
#' via euclidean distance. For intra-genomic hits in polyploids, the self-hit
#' regions are masked and the synteny search is re-run to more effectively find
#' homeologs.
#'
#' It is important to note that this does NOT produce finalized block
#' coordinates, but only large-scale regions that are syntenic. These results
#' are usually passed to an additional orthofinder run, either globally, or for
#' polyploids or searches with paralogs, within-block pairwise runs. See
#' rerun_orthofinder.
#'
#' Parameterization of this function is complex and varies by the type of
#' contrast desired. To simplify parameterization, we have build a convenience
#' function that infers various parameters based on basic
#' genespace parameters and ploidy. Set set_syntenyParams for more information
#' about the input synParam data.table.
#'
#' @return a 'hits' data.table for a pairwise combination of genomes.
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @title synteny
#' @description
#' \code{synteny} the pipeline to form sytenic regions, blocks and coordiantes
#' @rdname synteny
#' @import data.table
#' @importFrom grDevices pdf dev.off
#' @export
synteny <- function(gsParam, genomeIDs = NULL, overwrite = F){
  ##############################################################################
  # 0. Set up the environment and parse the synteny data.table
  ##############################################################################
  genome1 <- genome2 <- gnum1 <- gnum2 <- synOg <- genome <- og <- nChr1 <- NULL
  nChr2 <- blkID <- ofID1 <- ofID2 <- score <- scrRank1 <- scrRank2 <- NULL
  isAnchor <- gen1 <- gen2 <- chr1 <- chr2 <- regBuffer <- globOG <- NULL
  ord1 <- ord2 <- NULL

  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
  }else{
    if(!any(genomeIDs %in% gsParam$genomes$genomeIDs))
      stop("specified genomeIDs dont look right\n")
  }
  genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]

  verbose <- gsParam$params$verbose
  writeTo <- gsParam$paths$results
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  blksFile <- file.path(gsParam$paths$results, "syntenicBlocks.txt.gz")

  # set the synteny parameters
  if(is.data.table(gsParam$params$synteny))
    if(!all(genomeIDs %in% gsParam$params$synteny$genome1))
      gsParam$params$synteny <- NULL
  if(!is.data.table(gsParam$params$synteny)){
    cat("Synteny Parameters have not been set! Setting to defaults\n")
    gsParam <- set_syntenyParams(gsParam)
  }

  # -- get an annotated gff file
  gsParam <- annotate_gff(
    gsParam = gsParam,
    genomeIDs = genomeIDs,
    overwrite = overwrite)

  # -- enforce genome order from genomeIDs, keep only uni-directional blasts
  syn <- data.table(gsParam$params$synteny)
  syn[,`:=`(gnum1 = match(genome1, genomeIDs),
            gnum2 = match(genome2, genomeIDs))]
  syn <- subset(syn, gnum1 <= gnum2)
  setkey(syn, gnum1, gnum2)

  # -- split by intra and intergenomic hits
  selfSyn <- subset(syn, genome1 == genome2)
  nonSelfSyn <- subset(syn, genome1 != genome2)

  ##############################################################################
  # 1. load and parse the gff
  ##############################################################################
  # -- find the location of orthofinder results files
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)

  # -- annotate gff with arrays, syntenic orthogroups, etc.
  gff <- fread(gffFile, na.strings = c("NA",""), showProgress = F)
  arrayReps <- gff$ofID[gff$isArrayRep]

  # -- flag syntenic arrays and choose optimal representative gene
  isArrayRep <- ord <- NULL
  if(verbose)
    with(gff, cat(sprintf(
      "Found %s genes, %s orthogroups and %s arrays with %s genes\n",
      nrow(gff), uniqueN(globOG),
      uniqueN(arrayID, na.rm = T), sum(!is.na(arrayID)))))

  ##############################################################################
  # 2. self culling
  # -- run the following steps for each genome and only within genome hits
  ##############################################################################
  if(verbose)
    cat("Pulling within-genome synteny ... \n\tGenome:",
        "n raw hits / hits in (regions) / hits in (blks)\n")
  intraBcs <- rbindlist(lapply(1:nrow(selfSyn), function(i){

    # -- check if analysis should be run
    genome1 <- selfSyn$genome1[i]
    if(verbose)
      cat(sprintf("\t%s (selfhit): ", pull_strWidth(genome1, 7)))
    hitsFile <- file.path(
      writeTo,
      sprintf("%s_%s_synHits.txt.gz", genome1, genome1))
    if(file.exists(hitsFile) & !overwrite){
      hits <- fread(hitsFile, showProgress = F, na.strings = c("", "NA"))
      outBlks <- calc_blkCoords(subset(hits, !is.na(blkID)))
      if(verbose)
        cat("found existing hits and !overwrite, so reusing\n")
    }else{
      noParalog <- with(selfSyn[i,],
                        ploidy1 + ploidy2 + nSecondHits1 + nSecondHits2 == 2)

      # -- read in the blast
      ahits <- parse_blast4synteny(
        gsParam = gsParam,
        genome1 = genome1,
        genome2 = genome1,
        gff = gff,
        selfOnly = FALSE)

      # -- subset to the same chromosome
      phits <- subset(ahits, chr1 == chr2)
      cat(sprintf("%s / ", nrow(phits)))
      shits <- subset(
        phits,
        nChr1 >= selfSyn$blkSize[i] & nChr2 >= selfSyn$blkSize[i])

      # -- self synteny
      sr <- pull_selfRegion(
        hits = shits,
        synBuff = selfSyn$synBuff[i])
      blkv <- sr$blkID; names(blkv) <- with(sr, paste(ofID1, ofID2))
      phits[,`:=`(blkID = blkv[paste(ofID1, ofID2)],
                  blkAnchor = ofID1 == ofID2,
                  blkBuffer = paste(ofID1, ofID2) %in% names(blkv),
                  regID = blkv[paste(ofID1, ofID2)],
                  regAnchor = ofID1 == ofID2,
                  regBuffer = paste(ofID1, ofID2) %in% names(blkv),
                  lgRegionID = blkv[paste(ofID1, ofID2)])]
      cat(sprintf("%s (%s) / %s (%s)\n",
                  sum(phits$regBuffer), uniqueN(phits$regID, na.rm = T),
                  sum(phits$blkBuffer), uniqueN(phits$blkID, na.rm = T)))
      if(!noParalog){
        mhits <- subset(phits, chr1 == chr2)
        mhits  <- subset(mhits, abs(ord1 - ord2) <= selfSyn$selfRegionMask[i] | regBuffer)
      }
      phits <- subset(phits, regBuffer)

      ############################################################################
      # 2.2 homeologs, if desired
      ############################################################################
      if(noParalog){
        u <- with(phits, unique(paste(ofID1, ofID2)))
        shits <- subset(ahits, !paste(ofID1, ofID2) %in% u)
      }else{
        if(verbose)
          cat(sprintf("\t%s (homeolg): ",  pull_strWidth(genome1, 7)))
        shits <- pipe_synteny(
          gsParam = gsParam,
          gff = gff,
          maskHits = mhits,
          synParam = selfSyn[i,],
          type = "secondary")
        shits$regID[!is.na(shits$regID)] <- paste0("sec_", shits$regID)[!is.na(shits$regID)]
        shits$blkID[!is.na(shits$blkID)] <- paste0("sec_", shits$blkID)[!is.na(shits$blkID)]
        shits$lgRegionID[!is.na(shits$lgRegionID)] <- paste0("sec_", shits$lgRegionID)[!is.na(shits$lgRegionID)]
      }

      outHits <- rbind(phits, shits, fill = T)
      outBlks <- calc_blkCoords(subset(outHits, !is.na(blkID)))

      # -- write the results
      fwrite(
        outHits, quote = F, showProgress = F, sep = "\t",
        file = file.path(writeTo, sprintf("%s_%s_synHits.txt.gz",
                                          genome1, genome1)))

      pm <- par()["mfrow"]
      par(mfrow = c(1,1))
      pmar <- par()["mar"]
      par(mar = c(3,2,2,2))

      # -- make the plots
      pdf(file.path(gsParam$paths$results,
                    sprintf("%s_%s_dotplots.pdf", genome1, genome1)),
          height = 6, width = 6)
      par(mfrow = c(1,1))

      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = subset(outHits, score > 100),
        alpha = .025,
        gff = gff,
        onlyOg = F,
        cols = "blue2",
        bufferOnly = F,
        anchorOnly = F,
        round2 = 5,
        plotTitle = "All diamond hits score > 100"))
      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = T,
        plotTitle = "hits in syntenic regions"))
      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = F, anchorOnly = T,
        plotTitle = "anchor hits in blocks"))
      dev.off()

      par(mfrow = pm)
      par(mar = pmar)
      # -- return block coordinates
    }
    return(outBlks)
  }))


  ##############################################################################
  # 3. intergenomic regions
  # -- run the following steps for each non-self pairwise genome combination
  ##############################################################################
  if(verbose)
    cat("Pulling intergenomic synteny ...\n")

  interBcs <- rbindlist(lapply(1:nrow(nonSelfSyn), function(i){
    gid1 <- nonSelfSyn$genome1[i]
    gid2 <- nonSelfSyn$genome2[i]
    syn <- nonSelfSyn[i,]
    if(verbose)
      cat(sprintf(
        "\t%s-%s (primary): ", pull_strWidth(gid1, 7), pull_strWidth(gid2, 7)))

    hitsFile <- file.path(
      writeTo,
      sprintf("%s_%s_synHits.txt.gz", gid1, gid2))
    if(file.exists(hitsFile) & !overwrite){
      hits <- fread(hitsFile, showProgress = F, na.strings = c("", "NA"))
      outBlks <- calc_blkCoords(subset(hits, !is.na(blkID)))
      if(verbose)
        cat("found existing hits and !overwrite, so reusing\n")
    }else{
      phits <- pipe_synteny(
        gsParam = gsParam,
        gff = gff,
        maskHits = NULL,
        synParam = syn)

      # -- optionally get secondary hits too.
      if(nonSelfSyn$nSecondHits[i] == 0){
        shits <- NULL
      }else{
        if(verbose)
          cat(sprintf(
            "\t%s-%s (secondy): ", pull_strWidth(gid1, 7), pull_strWidth(gid2, 7)))
        shits <- pipe_synteny(
          gsParam = gsParam,
          gff = gff,
          maskHits = subset(phits, regBuffer),
          synParam = syn,
          type = "secondary")
        shits$regID[!is.na(shits$regID)] <- paste0("sec_", shits$regID)[!is.na(shits$regID)]
        shits$blkID[!is.na(shits$blkID)] <- paste0("sec_", shits$blkID)[!is.na(shits$blkID)]
      }

      outHits <- rbind(phits, shits, fill = T)
      outBlks <- calc_blkCoords(subset(outHits, !is.na(blkID)))

      # -- write the results
      fwrite(
        outHits, quote = F, showProgress = F, sep = "\t",
        file = file.path(writeTo, sprintf("%s_%s_synHits.txt.gz",
                                          gid1, gid2)))

      nr1 <- sum(gff$genome == gid1)
      nr2 <- sum(gff$genome == gid2)
      plw <- 6
      plh <- 6
      if(nr1 < nr2)
        plw <- plw * (nr2/nr1)
      if(nr1 > nr2)
        plh <- plh * (nr1/nr2)
      # -- make the plots

      pdf(file.path(gsParam$paths$results,
                    sprintf("%s_%s_dotplots.pdf", gid1, gid2)),
          height = plw, width = plh)

      pm <- par()["mfrow"]
      par(mfrow = c(1,1))
      pmar <- par()["mar"]
      par(mar = c(3,2,2,2))

      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = subset(outHits, score > 100),
        gff = gff,
        alpha = .025,
        onlyOg = F,
        cols = "blue2",
        bufferOnly = F,
        anchorOnly = F,
        round2 = 5,
        plotTitle = "All diamond hits score > 100"))
      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = T,
        plotTitle = "hits in syntenic regions"))
      suppressWarnings(pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = F, anchorOnly = T,
        plotTitle = "anchor hits in blocks"))
      dev.off()

      par(mfrow = pm)
      par(mar = pmar)
    }
    # return block coordinates
    return(outBlks)
  }))

  # return complete block coordinate data.table
  blks1 <- rbind(intraBcs, interBcs)
  blks2 <- data.table(blks1)
  setnames(blks2, gsub("2$","3", colnames(blks2)))
  setnames(blks2, gsub("1$","2", colnames(blks2)))
  setnames(blks2, gsub("3$","1", colnames(blks2)))
  cn <- intersect(colnames(blks1), colnames(blks2))
  blks <- rbind(blks1[,cn, with = F], blks2[,cn, with = F])
  blks <- subset(blks, !duplicated(paste(gen1, gen2, blkID)))

  fwrite(blks, file = blksFile, sep = "\t", quote = F, showProgress = F)
  if(verbose)
    cat(sprintf(
      "\tSynteny constraints - Done!\n\tSyntenic block coordinates written to /results/%s\n",
        basename(blksFile)))

  # -- get syntenic orthologs
  gsParam <- pull_synOGs(gsParam = gsParam)
  if(verbose)
    cat("\tWrote gff to file: /results/gffWithOgs.txt.gz\n\tDone!\n")
  return(gsParam)
}

#' @title calc_blkCoords
#' @description
#' \code{calc_blkCoords} calc_blkCoords
#' @rdname synteny
#' @import data.table
#' @importFrom stats cor
#' @export
calc_blkCoords <- function(hits){
  ord1 <- ord2 <- start1 <- start2 <- end1 <- end2 <- ofID1 <- ofID2 <-minBp2 <- NULL
  blkID <- orient <- maxBp2 <- maxBp1 <- maxOrd2 <- minOrd2 <-maxGene2 <- minGene2 <- NULL

  bhits <- subset(hits, !is.na(blkID))
  setkey(bhits, ord1)
  blks1 <- bhits[,list(startBp1 = min(start1, na.rm = T), endBp1 = max(end1, na.rm = T),
                       startOrd1 = min(ord1, na.rm = T), endOrd1 = max(ord1, na.rm = T),
                       firstGene1 = first(ofID1), lastGene1 = last(ofID1),
                       nHits1 = uniqueN(ofID1)),
                 by = c("blkID", "gen1","gen2", "chr1", "chr2")]
  setkey(bhits, ord2)
  blks2 <- bhits[,list(minBp2 = min(start2, na.rm = T), maxBp2 = max(end2, na.rm = T),
                       minOrd2 = min(ord2, na.rm = T), maxOrd2 = max(ord2, na.rm = T),
                       minGene2 = first(ofID2), maxGene2 = last(ofID2, 1),
                       nHits2 = uniqueN(ofID2),
                       orient = ifelse(length(ord1) <= 1, "+",
                                       ifelse(cor(jitter(ord1),
                                                  jitter(ord2)) > 0,"+", "-"))),
                 by = c("blkID")]
  blks <- merge(blks1, blks2, by = "blkID")

  bgfor <- subset(blks, orient == "+")
  bgrev <- subset(blks, orient == "-")
  bgrev[,`:=`(startBp2 = maxBp2, endBp2 = minBp2,
              startOrd2 = maxOrd2, endOrd2 = minOrd2,
              firstGene2 = maxGene2, lastGene2 = minGene2)]
  bgfor[,`:=`(startBp2 = minBp2, endBp2 = maxBp2,
              startOrd2 = minOrd2, endOrd2 = maxOrd2,
              firstGene2 = minGene2, lastGene2 = maxGene2)]
  blks <- rbind(bgfor, bgrev)
  return(blks)
}


#' @title finalize_blocks
#' @description
#' \code{finalize_blocks} finalize_blocks
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
finalize_blocks <- function(hits,
                            synParam,
                            gsParam,
                            type = "primary",
                            verbose = FALSE,
                            minRbhScore = 50){

  isOg <- maxScr1 <- maxScr2 <- score <- og <- ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  isCollin <- blkID <- bord1 <- bord2 <- gen1 <- gen2 <- nblk <- regID <- NULL
  regAnchor <- blkIDn <- chr1 <- chr2 <- isAnchor <- NULL
  nhits1 <- nSecondHits1 <- nhits2 <- nSecondHits2 <- blkSizeSecond <- NULL
  synBuffSecond <- nGapsSecond <- NULL
  ##############################################################################
  # 1. Pull anchor hits and re-rank
  ##############################################################################
  sp <- data.table(synParam[1,])
  if(type == "Second")
    sp[,`:=`(nhits1 = (nhits1-1 + nSecondHits1),
             nhits2 = (nhits2-1 + nSecondHits2),
             blkSize = blkSizeSecond,
             nGaps = nGapsSecond,
             synBuff = synBuffSecond)]

  dropInterleavesSmallerThan <- sp$dropInterleavesSmallerThan
  blkSize <- sp$blkSize
  nGaps <- sp$nGaps
  synBuff <- sp$synBuff
  path2mcscanx <- gsParam$paths$mcscanxCall
  nCores <- gsParam$params$nCores

  anchs <- subset(hits, regAnchor & !is.na(regID))

  # -- calculate rank order within regions
  anchs[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "regID"]

  ##############################################################################
  # 2. remove non-colinear hits within regions
  ##############################################################################
  # -- choose representative hits within each region
  spl <- split(anchs, by = "regID")
  colAnch <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    suppressWarnings(x[,isCollin := !is.na(run_mcscanx(
      hits = x,
      gsParam = gsParam,
      blkSize = blkSize,
      nGaps = nGaps,
      path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)])])
    if(!"isCollin" %in% colnames(x))
      x[,isCollin := FALSE]
    return(x)
  }), use.names=TRUE)

  # -- drop non colinear hits and re-rank
  colAnch <- subset(colAnch, isCollin)
  colAnch[,isCollin := NULL]
  colAnch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense")),
          by = "regID"]

  ##############################################################################
  # 3. cluster within regions into blks
  ##############################################################################
  # -- initial clustering, allowing to drop some non-clustering points
  spl <- split(colAnch, by = "regID")
  iblks <- rbindlist(mclapply(spl, mc.cores = nCores, function(x){
    x[,blkID := dbscan(frNN(cbind(ord1, ord2),
                            eps = sqrt((blkSize^2)*2)),
                       minPts = blkSize)$cluster]
    return(x)
  }))
  iblks <- subset(iblks, blkID != 0)

  # -- clustering within blks, again dropping blocks too small
  iblks[,blkID := sprintf("%s_%s_%s_%s",gen1, gen2, regID, blkID)]
  iblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "blkID"]
  iblks[,nblk := dbscan(frNN(cbind(ord1, ord2),
                             eps = sqrt(2)+.1),
                        minPts = blkSize)$cluster,
        by = "blkID"]
  fblks <- subset(iblks, blkID != 0)

  # -- final clustering within blks
  fblks[,blkID := sprintf("%s_%s", blkID, nblk)]
  fblks[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")),
        by = "blkID"]
  fblks[,nblk := dbscan(frNN(cbind(ord1, ord2),
                             eps = sqrt(2)+.1),
                        minPts = 0)$cluster,
        by = "blkID"]

  # -- split interleaved blocks
  outBlks <- fblks[,colnames(hits), with = F]

  blv <- outBlks$blkID; names(blv) <- with(outBlks, paste(ofID1, ofID2))
  blAnc <- with(outBlks, paste(ofID1, ofID2))

  ##############################################################################
  # 4. format and return blk regs
  ##############################################################################

  # get block regions and hits therein
  hits[,blkID := blv[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(hits, !is.na(blkID)))
  inReg <- get_hitsInBlks(blks = bc, hits = hits, nCores = nCores)
  blv <- inReg$blkID; names(blv) <- with(inReg, paste(ofID1, ofID2))

  # -- get hits in buffer from block anchors
  inReg[,isAnchor := paste(ofID1, ofID2) %in% blAnc]
  inBuff <- find_hitsInBuff(
    hits = inReg, nCores = nCores, synBuff = synBuff)
  buffv <- with(inBuff, paste(ofID1, ofID2))

  # -- return a new hits obj
  hits[,`:=`(blkID = blv[paste(ofID1, ofID2)],
             blkAnchor = paste(ofID1, ofID2) %in% blAnc,
             blkBuffer = paste(ofID1, ofID2) %in% buffv)]
  return(hits)
}

#' @title parse_blast4synteny
#' @description
#' \code{parse_blast4synteny} parse_blast4synteny
#' @rdname synteny
#' @import data.table
#' @export
parse_blast4synteny <- function(gsParam,
                                genome1,
                                genome2,
                                gff,
                                selfOnly){

  score <- ofID1 <- ofID2 <- isOg <- og1 <- scrRank1 <- scrRank2 <- NULL
  chr1 <- chr2 <- ord1 <- ord2 <- NULL
  # -- get vectors from the gff
  gv <- gff$genome; cv <- gff$chr; ogv <- gff$globOG
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(gv) <- names(cv) <- names(ogv) <- gff$ofID
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- read in the blast hits
  ofsp <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  if(genome1 == genome2){
    bl <- read_blast(
      path = gsParam$paths$blastDir,
      ofID1 = ofsp[genome1],
      ofID2 = ofsp[genome1])
  }else{
    bl <- rbind(read_blast(
      path = gsParam$paths$blastDir,
      ofID1 = ofsp[genome1],
      ofID2 = ofsp[genome2]),
      with(read_blast(
        path = gsParam$paths$blastDir,
        ofID1 = ofsp[genome2],
        ofID2 = ofsp[genome1]),
        data.table(
          ofID1 = ofID2, ofID2 = ofID1, score = score)))
  }

  # -- choose only the best scoring non-duplicated hits
  setorder(bl, -score)
  bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))

  # -- add annotation information in
  bl[,`:=`(gen1 = gv[ofID1], gen2 = gv[ofID2],
           start1 = sv[ofID1], start2 = sv[ofID2],
           end1 = ev[ofID1], end2 = ev[ofID2],
           chr1 = cv[ofID1], chr2 = cv[ofID2],
           ord1 = ov[ofID1], ord2 = ov[ofID2],
           og1 = ogv[ofID1], og2 = ogv[ofID2],
           scrRank1 = 1, scrRank2 = 1,
           isOg = ogv[ofID1] == ogv[ofID2])]
  if(selfOnly)
    bl <- subset(bl, chr1 == chr2)
  bl <- subset(bl, !is.na(ord1) & !is.na(ord2))

  # -- get the number of unique hits by chr (for later culling)
  bl[,`:=`(nChr1 = uniqueN(ofID1[isOg]),
           nChr2 = uniqueN(ofID2[isOg])),
     by = c("gen1","gen2","chr1","chr2")]
  bl[,`:=`(og = ifelse(isOg, og1, NA), og1 = NULL, og2 = NULL)]

  # -- optionally add the rank scores by geneID
  if(!selfOnly){
    setorder(bl, ofID1, -score)
    bl[,scrRank1 := 1:.N, by = "ofID1"]
    setorder(bl, ofID2, -score)
    bl[,scrRank2 := 1:.N, by = "ofID2"]
  }
  return(bl)
}

#' @title find_hitsInBuff
#' @description
#' \code{find_hitsInBuff} find_hitsInBuff
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
find_hitsInBuff <- function(hits,
                            nCores,
                            synBuff){

  inBuffer <- ord1 <- ord2 <- isAnchor <- NULL

  if(!"isAnchor" %in% colnames(hits))
    stop("column name isAnchor must be in hits\n")
  splHit <- split(hits, by = c("gen1","gen2","chr1","chr2"))
  out <- rbindlist(mclapply(splHit, mc.cores = nCores, function(x){
    if(any(x$isAnchor)){
      nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synBuff))
      wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
      x[,inBuffer := 1:.N %in% wh]
      return(x)
    }
  }), fill = T)
  if(!"inBuffer" %in% colnames(out)){
    return(NULL)
  }else{
    return(subset(out, inBuffer))
  }
}

#' @title clus_dbscan
#' @description
#' \code{clus_dbscan} clus_dbscan
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
clus_dbscan <- function(hits,
                        radius,
                        blkSize,
                        nCores){
  blkID <- ord1 <- ord2 <- chr1 <- chr2 <- NULL

  x <- split(hits, by = c("gen1","gen2","chr1","chr2"))
  x <- rbindlist(mclapply(x, mc.cores = nCores, function(y){
    y[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = radius),
                       minPts = blkSize)$cluster]
    return(y)
  }))
  x <- subset(x, blkID > 0)
  x[,blkID := sprintf("blk_%s_%s_%s", chr1, chr2, blkID)]
  return(x)
}

#' @title find_globalAnchors
#' @description
#' \code{find_globalAnchors} find_globalAnchors
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
find_globalAnchors <- function(hits,
                               gsParam,
                               synParam,
                               type = "primary"){

  add_rbh2og <- function(hits, minRbhScore = 50){
    ofID1 <- ofID2 <- NULL
    maxScr1 <- maxScr2 <- score <- og <- NULL
    hits[,maxScr1 := max(score), by = c("blkID", "ofID1")]
    hits[,maxScr2 := max(score), by = c("blkID", "ofID2")]

    outrb <- rbindlist(lapply(split(hits, by = "blkID"), function(x){
      xog <- subset(x, !is.na(og))
      xrbh <- subset(x, maxScr1 == score & maxScr2 == score & score > minRbhScore &
                       !ofID1 %in% xog$ofID1 & !ofID2 %in% xog$ofID2)
      xrbh[,og := "RBH"]
      return(rbind(xrbh, xog))
    }))
    outrb$og[outrb$og == "RBH"] <- sprintf("RBH%s",1:sum(outrb$og == "RBH"))
    return(outrb)
  }

  nhits1 <- nSecondHits1 <- nhits2 <- nSecondHits2 <- blkSizeSecond <- NULL
  nGapsSecond <- synBuffSecond <- isMasked <- ord1 <- ord2 <- isArrayRep <- NULL
  og <- ofID1 <- score <- scrRank1 <- ofID2 <- scrRank2 <- blkID <- NULL
  isAnchor <- inBuffer <- onlyOgAnchorsSecond <- NULL
  # -- get synteny parameters in line
  sp <- data.table(synParam[1,])
  if(type != "primary")
    sp[,`:=`(nhits1 = (nhits1-1 + nSecondHits1),
             nhits2 = (nhits2-1 + nSecondHits2),
             blkSize = blkSizeSecond,
             nGaps = nGapsSecond,
             synBuff = synBuffSecond,
             onlyOgAnchors = onlyOgAnchorsSecond)]

  # -- drop masked hits and re-rank order
  hits <- subset(hits, !isMasked)
  hits[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]

  # -- drop non-array rep genes
  anch <- subset(hits, isArrayRep)

  # -- drop non-og hits, if necessary
  if(sp$onlyOgAnchors)
    anch <- subset(anch, !is.na(og))

  # -- subset to top n hits
  setorder(anch, ofID1, -score)
  anch[,scrRank1 := 1:.N, by = "ofID1"]
  setorder(anch, ofID2, -score)
  anch[,scrRank2 := 1:.N, by = "ofID2"]
  anch <- subset(anch, scrRank1 <= sp$nhits1 & scrRank2 <= sp$nhits2)

  # -- re-rank
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]

  # -- initial anchor building
  anchu <- run_mcscanx(
    hits = anch,
    gsParam = gsParam,
    blkSize = sp$blkSize,
    nGaps = sp$nGaps,
    path2mcscanx = gsParam$paths$mcscanxCall)
  anch <- subset(anch, paste(ofID1, ofID2) %in% names(anchu))

  # -- initial anchor clustering
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch <- clus_dbscan(
    hits = anch,
    radius = sp$synBuff,
    blkSize = sp$blkSize,
    nCores = gsParam$params$nCores)
  anchv <- anch$blkID; names(anchv) <- with(anch, paste(ofID1, ofID2))

  # -- find all hits in regions in the array rep data
  anch <- subset(hits, isArrayRep)
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch[,blkID := anchv[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(anch, !is.na(blkID)))
  anch <- get_hitsInBlks(
    blks = bc,
    hits = anch,
    nCores = gsParam$params$nCores)
  anch[,isAnchor := paste(ofID1, ofID2) %in% names(anchv)]
  anch <- subset(find_hitsInBuff(
    hits = anch,
    nCores = gsParam$params$nCores,
    synBuff = sp$synBuff), inBuffer)

  # -- add in RBHs and re cluster
  anch <- add_rbh2og(
    hits = subset(anch, !is.na(blkID) & inBuffer))
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anchu <- run_mcscanx(
    hits = anch,
    gsParam = gsParam,
    blkSize = sp$blkSize,
    nGaps = sp$nGaps,
    path2mcscanx = gsParam$paths$mcscanxCall)
  anch <- subset(anch, paste(ofID1, ofID2) %in% names(anchu))

  # -- final anchor clustering
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  anch <- clus_dbscan(
    hits = anch,
    radius = sp$synBuff,
    blkSize = sp$blkSize,
    nCores = gsParam$params$nCores)
  anchv <- anch$blkID; names(anchv) <- with(anch, paste(ofID1, ofID2))
  hits[,blkID := anchv[paste(ofID1, ofID2)]]

  # -- split overlapping blocks
  anch2 <- split_overlaps(
    hits = subset(hits, !is.na(blkID)),
    dropInterleavesSmallerThan = 2,
    minPropDup = .05,
    verbose = F)
  anchv2 <- anch2$blkID; names(anchv2) <- with(anch2, paste(ofID1, ofID2))

  # -- final assignment
  hits[,blkID := anchv2[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(hits, !is.na(blkID)))
  out <- get_hitsInBlks(
    blks = bc,
    hits = hits,
    nCores = gsParam$params$nCores)
  out[,isAnchor := paste(ofID1, ofID2) %in% names(anchv)]
  out <- find_hitsInBuff(
    hits = out,
    nCores = gsParam$params$nCores,
    synBuff = sp$synBuff)
  return(out)
}

#' @title pipe_synteny
#' @description
#' \code{pipe_synteny} pipe_synteny
#' @rdname synteny
#' @import data.table
#' @export
pipe_synteny <- function(gsParam,
                         gff,
                         maskHits,
                         synParam,
                         type = "primary"){

  gen1 <- gen2 <- isArrayRep <- ofID1 <- ofID2 <- isMasked <- regBuffer <- NULL
  blkIDn <- blkID <- chr1 <- chr2 <- inBuffer <- isAnchor <- lgRegionID <- NULL

  verbose <- gsParam$params$verbose
  # -- read in the blasts
  bl <- with(synParam, parse_blast4synteny(
    gsParam = gsParam,
    genome1 = genome1,
    genome2 = genome2,
    gff = gff,
    selfOnly = FALSE))
  bl <- subset(bl, !is.na(gen1) & !is.na(gen2))
  if(verbose)
    cat(sprintf("%s / ", nrow(bl)))

  # -- flag array reps
  arrayReps <- gff$ofID[gff$isArrayRep]
  bl[,isArrayRep := ofID1 %in% arrayReps & ofID2 %in% arrayReps]

  # -- optionally remove masked regions
  bl[,isMasked := FALSE]
  if(!is.null(maskHits)){
    u <- with(maskHits, paste(ofID1, ofID2))
    bl$isMasked[with(bl, paste(ofID1, ofID2)) %in% u] <- TRUE
  }

  # -- get the syntenic region hits
  regs <- find_globalAnchors(
    hits = bl,
    gsParam = gsParam,
    synParam = synParam,
    type = type)
  regs[,blkIDn := as.numeric(as.factor(blkID)), by = c("chr1", "chr2")]
  regs[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, blkIDn)]

  if(verbose)
    cat(sprintf("%s (%s) / ", sum(regs$inBuffer), uniqueN(regs$blkID, na.rm = T)))

  # -- get the syntenic blocks
  regs[,`:=`(regID = blkID, regBuffer = inBuffer, regAnchor = isAnchor)]
  blks <- finalize_blocks(
    hits = regs,
    synParam = synParam,
    gsParam = gsParam,
    verbose = T,
    type = type)

  out <- merge(
    bl,
    blks[,c("ofID1","ofID2","regID","blkID","regAnchor", "blkAnchor","blkBuffer","regBuffer")],
    by = c("ofID1", "ofID2"), all = T)

  h <- clus_dbscan(
    subset(out, regBuffer),
    radius = synParam$synBuff+1,
    blkSize = 1,
    nCores = gsParam$params$nCores)
  hb <- h$blkID;  names(hb) <- with(h, paste(ofID1, ofID2))
  out[,lgRegionID := hb[paste(ofID1, ofID2)]]
  # -- reformat and return
  if(verbose)
    cat(sprintf("%s (%s)\n", sum(blks$blkBuffer), uniqueN(blks$blkID, na.rm = T)))
  for(i in c("regAnchor", "blkAnchor","blkBuffer","regBuffer","lgRegionID"))
    out[[i]][is.na(out[[i]])] <- FALSE
  return(out)
}

#' @title pull_nonSelfReg
#' @description
#' \code{pull_nonSelfReg} pull_nonSelfReg
#' @rdname synteny
#' @import data.table
#' @export
pull_nonSelfReg <- function(hits,
                            selfRegionMask,
                            nhits){
  ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  u <- with(pull_selfRegion(hits, synBuff = selfRegionMask),
            unique(paste(ofID1, ofID2)))
  x <- subset(hits, !paste(ofID1, ofID2) %in% u)
  x[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
          ord2 = frank(ord2, ties.method = "dense"))]
  return(x)
}

#' @title pull_selfRegion
#' @description
#' \code{pull_selfRegion} pull_selfRegion
#' @rdname synteny
#' @import data.table
#' @export
pull_selfRegion <- function(hits,
                            synBuff){
  chr1 <- chr2 <- ord1 <- ord2 <- ofID1 <- ofID2 <- NULL
  buff <- sqrt(synBuff^2 * 2) + 1
  x <- subset(
    hits,
    chr1 == chr2 &
      abs(ord1 - ord2) <= synBuff*2)
  x[,`:=`(blkID = sprintf("blk_%s_%s_self", chr1, chr2),
          isAnchor = ofID1 == ofID2)]
  return(x)
}

#' @title get_hitsInBlks
#' @description
#' \code{get_hitsInBlks} get_hitsInBlks
#' @rdname synteny
#' @import data.table
#' @importFrom parallel mclapply
#' @export
get_hitsInBlks <- function(blks, hits, nCores){
  ord1 <- ord2 <- blkID <- NULL
  splh <- split(hits, by = c("chr1","chr2","gen1","gen2"))
  splb <- split(blks, by = c("chr1","chr2","gen1","gen2"))
  splh <- splh[names(splb)]
  hits <- rbindlist(mclapply(1:nrow(blks), mc.cores = 1, function(i){
    y <- blks[i,]
    x <- splh[[with(y, paste(chr1, chr2, gen1, gen2, sep = "."))]]

    x <- subset(x, ord1 >= y$startOrd1 & ord1 <= y$endOrd1 &
                  ord2 >= y$minOrd2 & ord2 <= y$maxOrd2)
    x[,blkID := y$blkID]
    return(x)
  }))
  return(hits)
}

#' @title split_overlaps
#' @description
#' \code{split_overlaps} split_overlaps
#' @rdname synteny
#' @export
split_overlaps <- function(hits,
                           dropInterleavesSmallerThan,
                           minPropDup = 0.05,
                           maxIter = 20,
                           verbose){
  ############################################################################
  ############################################################################

  ############################################################################
  # -- function to split blocks by RLEs
  split_ovlBlks <- function(blk1, blk2, hits, dropInterleavesSmallerThan){
    blkID <- ofID1 <- ofID2 <- ord1 <- rl1 <- ord2 <- rl2 <- n <- NULL
    y <- subset(hits, blkID %in% c(blk1, blk2))
    y <- subset(y, !duplicated(ofID1))
    y <- subset(y, !duplicated(ofID2))
    setkey(y, ord1)
    y[,rl1 := add_rle(blkID)]
    setkey(y, ord2)
    y[,rl2 := add_rle(blkID)]
    y[,blkID := sprintf("%s_%s%s", blkID, rl1, rl2)]
    y[,n := .N, by = blkID]
    y <- subset(y, n >= dropInterleavesSmallerThan)
    y[,`:=`(rl1 = NULL, rl2 = NULL, n = NULL)]
    return(y)
  }
  ############################################################################
  # -- function to count and find overlapping blocks
  count_ovlHits <- function(hits, minPropDup){
    blkID <- i.blkID <- blk1 <- blk2 <- r1 <- r2 <- propDup <- n <- hasPriority <- NULL
    # -- calculate block coordinates
    hitstmp <- data.table(hits)
    bc <- calc_blkCoords(hitstmp)
    bc1 <- with(bc, data.table(
      blkID = blkID, chr = chr1, start = startOrd1, end = endOrd1,
      key = c("chr","start","end")))
    bc2 <- with(bc, data.table(
      blkID = blkID, chr = chr2, start = minOrd2, end = maxOrd2,
      key = c("chr","start","end")))

    # -- find non-self overlapping regions
    fo1 <- subset(foverlaps(bc1, bc1), blkID != i.blkID)
    fo2 <- subset(foverlaps(bc2, bc2), blkID != i.blkID)

    # -- add overlapping info into an output object
    out <- data.table(
      blk1 = c(fo1$blkID, fo2$blkID),
      blk2 = c(fo1$i.blkID, fo2$i.blkID))

    # -- choose the non-repeated overlaps
    u <- unique(c(out$blk1, out$blk2))
    un <- order(u); names(un) <- u
    out[,`:=`(r1 = un[blk1], r2 = un[blk2])]
    out[,u := ifelse(r1 > r2,  paste(blk1, blk2), paste(blk2, blk1))]
    out <- subset(out, !duplicated(u))
    out[,`:=`(r1 = NULL, r2 = NULL, u = NULL)]

    # -- count the number of duplicated hits between each pair of blocks
    if(nrow(out) > 0){
      tmp <- split(subset(hitstmp, blkID %in% unlist(out)), by = "blkID")
      out[,propDup := sapply(1:nrow(out), function(j){
        tb1 <- tmp[[out$blk1[j]]]
        tb2 <- tmp[[out$blk2[j]]]
        tot <- uniqueN(c(tb1$ofID1, tb1$ofID2, tb2$ofID1, tb2$ofID2))
        dup <- sum(tb1$ofID1 %in% tb2$ofID1, tb1$ofID2 %in% tb2$ofID2,
                   tb2$ofID1 %in% tb1$ofID1, tb2$ofID2 %in% tb1$ofID2)
        return(dup/tot)
      })]
      out <- subset(out, propDup <= minPropDup)
      if(nrow(out) > 0){
        # -- get the total number of hits in each block
        out[,`:=`(n1 = sapply(tmp[blk1], nrow),
                  n2 = sapply(tmp[blk2], nrow))]

        # -- prioritize by the pairs with the largest small block
        out[,n := apply(.SD, 1, min), .SDcols = c("n1", "n2")]

        setorder(out, -n)
        out[,n := NULL]
        v <- matrix(duplicated(as.character(unlist(t(out[,1:2])))), ncol = 2, byrow = T)
        out[,hasPriority := !apply(v, 1, any)]
        return(out)
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
  }
  ############################################################################
  ############################################################################
  # -- iterate through trying to split overlapping blocks
  hasPriority <- blkID <- NULL
  anchs <- data.table(hits)
  nov <- 1
  iter <- 0
  while(nov > 0 && iter <= maxIter){
    iter <- iter + 1
    if(verbose)
      cat(sprintf("\titer %s: %s hits / %s blks\n",
                  iter, nrow(anchs), uniqueN(anchs$blkID)))
    ovlBlks <- count_ovlHits(anchs, minPropDup = minPropDup)
    if(!is.null(ovlBlks)){
      ovlBlks <- subset(ovlBlks, hasPriority)
      nov <- nrow(ovlBlks)
      if(is.null(nov))
        nov <- 0
      if(nov > 0){
        splAnch <- rbindlist(lapply(1:nrow(ovlBlks), function(j)
          split_ovlBlks(
            blk1 = ovlBlks$blk1[j],
            blk2 = ovlBlks$blk2[j],
            hits = anchs,
            dropInterleavesSmallerThan)))
        if(!is.null(splAnch)){
          anchs <- rbind(
            splAnch,
            subset(anchs, !blkID %in% unique(unlist(ovlBlks[,c(1:2)]))))
        }else{
          nov <- 0
        }
      }else{
        nov <- 0
      }
    }else{
      nov <- 0
    }
  }
  if(verbose)
    cat("\nfound no additional overlapping blocks\n")
  return(anchs)
}

#' @title add number of anchors to gff
#' @description
#' \code{annotate_gff} add number of anchors to gff
#' @rdname synteny
#' @importFrom parallel mclapply
#' @export
annotate_gff <- function(gsParam,
                         genomeIDs = NULL,
                         overwrite = FALSE){

  genome <- ord  <- arrayID <- isArrayRep <- med <- pepLen <- rnk <- ofID <- NULL
  medbp <- start <- dist2med <- dist2bp <- og <- n <- chr <- globOG <- NULL

  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(file.exists(gffFile) & !overwrite){
    tmp <- fread(gffFile, na.strings = c("-", "NA", ""), showProgress = F)
    if(all(c("pepLen", "globOG", "arrayID","isArrayRep") %in% colnames(tmp)))
      stop("annotated gff file exists and !overwrite, so not running ...\n")
  }
  if(!is.data.table(gsParam$params$synteny))
    stop("Must run set_syntenyParams first!\n")

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)
  if(verbose)
    cat("Loading annotations ...\n")
  # -- get paths to the orthofinder run
  if(is.na(gsParam$paths$orthogroupsDir)){
    if(verbose)
      cat("\tIndexing location of orthofinder results ... ")
    gsParam <- find_orthofinderResults(gsParam)
    if(verbose)
      cat("Done!\n")
  }

  ##############################################################################
  # 1. Load the gff and add global metadata
  ##############################################################################
  # -- read in the gff
  if(verbose)
    cat("\tReading the gffs ... ")
  gff <- read_gff(gsParam$paths$gff)
  gff <- add_ofID2gff(gff, gsParam$paths$blastDir)
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add peptide length
  if(verbose)
    cat("Done!\n\tPulling gene lengths ... ")
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)

  # -- add global orthogroupsto the gff
  if(verbose)
    cat("Done!\n\tParsing global orthogroups ... ")
  ogs <- parse_ogs(gsParam)
  gff <- merge(gff, ogs, by = c("genome","id"), all.x = T)
  gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
  setnames(gff, "ogID", "globOG")

  ##############################################################################
  # 2. Define collinear arrays
  ##############################################################################
  # -- build arrays from method in gsParam
  if(verbose)
    cat("Done!\nDefining collinear orthogroup arrays ... \n")
  gff <- add_arrays2gff(gsParam = gsParam, gff = gff)
  gff[,arrayID := paste(arrayID, og)]
  gff[,n := .N, by = "arrayID"]
  gff$arrayID[gff$n == 1] <- NA
  gff[,`:=`(n = NULL, isArrayRep = TRUE, og = NULL)]

  # -- choose the array reps
  if(verbose)
    cat("\tChoosing array representative genes ... ")
  gffa <- subset(gff, !is.na(arrayID))
  setkey(gffa, genome, ord)
  gffn <- subset(gff, is.na(arrayID))
  gffa[,arrayID := sprintf(
    "%s_%s_%s",
    genome, chr, as.numeric(factor(arrayID, levels = unique(arrayID))))]
  gffa[,med := as.numeric(median(ord)), by = "arrayID"]
  gffa[,medbp := as.numeric(median(start)), by = "arrayID"]
  gffa[,`:=`(dist2med = abs(med - ord),
             dist2bp = abs(medbp - start))]
  setorder(gffa, arrayID, dist2med, dist2bp, -pepLen)
  gffa[,rnk := 1:.N, by = "arrayID"]
  gffa[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL, dist2bp = NULL,
             med = NULL, medbp = NULL)]
  gff <- rbind(gffa, gffn)
  setkey(gff, genome, ord)

  if(verbose)
    cat("Done!\n")
  gff[,`:=`(synOG = NA, inBlkOG = NA, combOG = NA, og = globOG, refCoord = NA)]
  fwrite(gff, file = gffFile, sep = "\t", quote = F, showProgress = F)
  return(gsParam)
}

#' @title add array representative to a gff object
#' @description
#' \code{add_arrayRep2gff} choose most central gene by orthogroup
#' @rdname synteny
#' @importFrom parallel mclapply
#' @export
add_arrayRep2gff <- function(gff,
                             gsParam){
  synArr <- chr <- NULL
  # -- count peptides
  gff <- add_pepLen2gff(gff = gff, gsParam = gsParam)
  genomeIDs <- unique(gff$genome)
  nCores <- gsParam$params$nCores

  # -- count number of orthologs
  di <- dir.exists(gsParam$paths$orthologuesDir)
  dl <- length(list.files(gsParam$paths$orthologuesDir)) > 1

  nGenome <- nGenes <- gen2 <- id1 <- gen1 <- genome <- id <- og <- NULL
  if(di && dl){
    ogcnt <- rbindlist(mclapply(genomeIDs, mc.cores = nCores, function(i){
      ogs <- parse_orthologues(gsParam = gsParam, refGenome = i)
      ogn <- ogs[,list(nGenome = uniqueN(gen2),
                       nGenes = .N), by = c("gen1","id1")]
      return(ogn)
    }))
    setorder(ogcnt, -nGenome, -nGenes)
    ogcnt <- subset(ogcnt, !duplicated(paste(gen1, id1)))
    nog <- ogcnt$nGenome; ng <- ogcnt$nGenes
    names(ng) <- names(nog) <- with(ogcnt, paste(gen1, id1))
    gff[,`:=`(nGenomeOrthologs = nog[paste(genome, id)],
              nTotalOrthologs = ng[paste(genome, id)])]
    gff$nGenomeOrthologs[is.na(gff$nGenomeOrthologs)] <- 0
    gff$nTotalOrthologs[is.na(gff$nTotalOrthologs)] <- 0
  }else{
    gff[,`:=`(nGenomeOrthologs = 0,
              nTotalOrthologs = 0)]
  }

  # -- split into single and multiple member arrays
  gffi <- data.table(gff)
  d2h <- gsParam$params$maxDistBtwPgHits
  gff[,synArr := as.integer(as.factor(paste(genome, chr, og)))]
  gff[,nog := .N, by = "synArr"]
  g1 <- subset(gff, nog == 1)
  g1[,nog := NULL]
  g2 <- subset(gff, nog > 1)

  # -- calculate the maximum distance between genes in an array
  maxJump <- ord <- NULL
  g2[,maxJump := max(diff(ord[order(ord)])), by = "synArr"]
  g2r <- subset(g2, maxJump > d2h)
  g2 <- subset(g2, maxJump <= d2h)

  # -- cluster genes in arrays with big jumps
  clus <- NULL
  if(nrow(g2r) > 1){
    g2[,maxJump := NULL]
    g2r[,clus := dbscan(frNN(cbind(ord, ord), eps = d2h), minPts = 0)$cluster,
        by = "synArr"]
    g2r[,synArr := paste(synArr, clus)]
    g2 <- rbind(g2,  g2r[,colnames(g2),with = F])
  }

  # -- calculate distance to the median
  dist2median <- nGenomeOrthologs <- pepLen <- nTotalOrthologs <- ofID <- NULL
  g2[,dist2median := abs(as.numeric(median(ord, na.rm = T)) - ord),
     by = c("synArr","genome","chr")]

  # -- order and rank genes, choosing representatives for each array
  setorder(g2, genome, chr, synArr, -nGenomeOrthologs, -nTotalOrthologs,
           dist2median, -pepLen, ord)
  arep <- rbind(g1, g2[,colnames(g1), with = F])
  sar <- as.numeric(as.factor(arep$synArr)); names(sar) <- arep$ofID
  arep <- subset(arep, !duplicated(synArr))
  gffi[,`:=`(synArray = sar[ofID],
             isArrayRep = ofID %in% arep$ofID)]
  return(gffi)
}

#' @title add_arrays2gff
#' @description
#' \code{add_arrays2gff} add_arrays2gff
#' @rdname synteny
#' @import data.table
#' @importFrom Biostrings readAAStringSet
#' @export
add_arrays2gff <- function(gsParam,
                           gff){
  arrayID <- genome <- chr <- globOG <- n <- rng <- ord <- clus <- og <- NULL
  collinearOG <- NULL

  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  synBuff <- max(gsParam$params$synteny$synBuff)

  # -- make global arrays from orthogroups
  gff[,arrayID := sprintf("%s_%s_%s", genome, chr, globOG)]
  gff[,n := .N, by = "arrayID"]

  # -- combine 1x ogs with ogs in regions < synBuff
  g1 <- subset(gff, n == 1)
  g2 <- subset(gff, n > 1)
  g2[,rng := diff(range(ord)),  by = "arrayID"]
  g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
  g2 <- subset(g2, rng > synBuff)

  # -- combine above with ogs without max gap < synBuff
  g2[,rng := max(diff(ord[order(ord)])), by = "arrayID"]
  g1 <- rbind(g1, subset(g2, rng <= synBuff)[,colnames(g1), with = F])
  g2 <- subset(g2, rng > synBuff)

  # -- split ogs with gaps
  g2[,clus := dbscan(frNN(cbind(ord, ord), eps = synBuff), minPts = 1)$cluster,
     by = "arrayID"]
  g2[,arrayID := sprintf("%s_%s", arrayID, clus)]
  gff <- rbind(g1, g2[,colnames(g1), with = F])

  # -- NA out arrays with just one member
  gff[,n := .N, by = "arrayID"]
  gff[,og := arrayID]
  gff$arrayID[gff$n == 1] <- NA
  gff[,n := NULL]
  setkey(gff, genome, ord)

  # -- print updates and number of global orthogroups
  if(verbose){
    cat("\tUsing collinear orthogroups for array identity:\n")
    nu <- lapply(split(subset(gff, !is.na(arrayID)), by = "genome"), function(x)
      cat(sprintf("\t%s: %s genes in %s collinear arrays\n",
                  x$genome[1], nrow(x), uniqueN(x$arrayID))))
  }

  return(gff)
}

#' @title pull_synOGs
#' @description
#' \code{pull_synOGs} pull_synOGs
#' @rdname synteny
#' @import data.table
#' @importFrom grDevices pdf dev.off
#' @export
pull_synOGs <- function(gsParam,
                        genomeIDs = NULL){

  isArrayRep <- ogInblk <- ofID1 <- ofID2 <- inBlkOG <- ofID <- NULL
  og <- synOG <- globOG <- u <- arrayID <- NULL
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, na.strings = c("", "NA"), showProgress = F)
  verbose <- gsParam$params$verbose

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs

  # -- pull global syntenic orthogroups
  if(verbose)
    cat(sprintf(
      "Checking synteny-constrained global orthogroups for synOGs\n\tn. global OGs = %s\n",
      uniqueN(gff$globOG)))
  synGff <- add_synOg2gff(
    gff = subset(gff, isArrayRep),
    useBlks = F,
    gsParam = gsParam,
    genomeIDs = genomeIDs,
    allowRBHinOg = T)
  synGff[,og := globOG]
  # -- add the syntenic OGs to the original gff
  gff[,u := ifelse(is.na(arrayID), ofID, arrayID)]
  sog <- synGff$synOG
  names(sog) <- with(synGff, ifelse(is.na(arrayID), ofID, arrayID))
  gff[,`:=`(synOG = sog[u], u = NULL)]

  if(verbose)
    cat(sprintf("\tn. syntenic OGs = %s\n", uniqueN(gff$synOG)))

  # -- if desired, pull orthogroups within blocks
  if(gsParam$params$orthofinderInBlk){
    if(verbose)
      cat("Adding syntenic orthogroups from pairwise w/in-region hits\n")
    ofh <- blkwise_orthofinder(
      gsParam = gsParam,
      gff = synGff)
    ofh[, ogInblk := clus_igraph(ofID1, ofID2)]
    ofv <- ofh$ogInblk; names(ofv) <- ofh$ofID1
    synGff[,inBlkOG := ofv[ofID]]
    mol <- max(synGff$inBlkOG, na.rm = T)
    nmis <- sum(is.na(synGff$inBlkOG))
    synGff$inBlkOG[is.na(synGff$inBlkOG)] <- (mol + 1):(mol + nmis)
    synGff[,inBlkOG := as.integer(factor(inBlkOG, levels = unique(inBlkOG)))]

    gff[,u := ifelse(is.na(arrayID), ofID, arrayID)]
    bog <- synGff$inBlkOG
    names(bog) <- with(synGff, ifelse(is.na(arrayID), ofID, arrayID))
    gff[,`:=`(inBlkOG = bog[u], u = NULL)]
    gff[,og := inBlkOG]
  }else{
    gff[,inBlkOG := NA]
    gff[,og := synOG]
  }
  fwrite(gff, file = gffFile, sep = "\t", quote = F, showProgress = F)
  return(gsParam)
}

#' @title add_synOg2gff
#' @description
#' \code{add_synOg2gff} add_synOg2gff
#' @rdname synteny
#' @export
add_synOg2gff <- function(gff,
                          hits = NULL,
                          gsParam,
                          genomeIDs,
                          allowRBHinOg,
                          useBlks){
  blkBuffer <- isSelf <- ofID1 <- ofID2 <- og <- ofID <- synOG <- blkID <- NULL
  # -- find hits files
  if(is.null(hits)){
    eg <- CJ(genomeIDs, genomeIDs)
    fs <- file.path(gsParam$paths$results,
                    sprintf("%s_%s_synHits.txt.gz",eg[[1]], eg[[2]]))
    fs <- fs[file.exists(fs)]

    # -- read all hits
    nCores <- gsParam$params$nCores
    hts <- rbindlist(mclapply(fs, mc.cores = nCores, function(i){
      if(useBlks){
        x <- fread(
          i, select = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"),
          na.strings = c("","NA"))
      }else{
        x <- fread(
          i, select = c("ofID1","ofID2","og","regBuffer","regAnchor","regID"),
          na.strings = c("","NA"),
          col.names = c("ofID1","ofID2","og","blkBuffer","blkAnchor","blkID"))
      }
      x <- subset(x, blkBuffer)
      x[,isSelf := any(ofID1 == ofID2), by = "blkID"]
      x <- subset(x, !isSelf & !is.na(og))
      if(!allowRBHinOg)
        x <- subset(x, !grepl("RBH", og))
      return(x[,c("ofID1", "ofID2", "blkAnchor", "blkID")])
    }))
  }else{
    hts <- data.table(hits)
    hts <- subset(hts, !is.na(og))
    if(!allowRBHinOg)
      hts <- subset(hts, !grepl("RBH", og))
    if(useBlks){
      hts <- hts[,c("ofID1", "ofID2", "blkBuffer", "blkAnchor", "blkID")]
    }else{
      hts <- with(hts, data.table(
        ofID1 = ofID1, ofID2 = ofID2,
        blkBuffer = regBuffer, blkAnchor = regAnchor, blkID = regID))
    }
    hts <- subset(hts, blkBuffer)
    hts[,blkBuffer := NULL]
  }

  # -- convert to syntenic orthogroups
  ic <- with(subset(hts, !is.na(blkID)), clus_igraph(
    id1 = c(ofID1, ofID2), id2 = c(ofID2, ofID1)))
  gff[,synOG := ic[ofID]]
  nmis <- sum(is.na(gff$synOG))
  mol <- max(gff$synOG, na.rm = T)
  gff$synOG[is.na(gff$synOG)] <- (mol + 1):(mol + nmis)
  gff[,synOG := as.integer(factor(synOG, levels = unique(synOG)))]
  return(gff)
}
