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
#' @param dropInterleavesSmallerThan integer of length 1, see set_syntenyParams
#' @param nGaps integer of length 1 specifying the -m 'gaps' MCScanX paramerter
#' @param path2mcscanx character string file.path pointing to the install
#' directory for the MCScanX program. In particular, there must be an
#' executable in this directory for MCScanX_h.
#' @param minRbhScore integer of length 1, see set_syntenyParams
#' @param genome1 character string specifying first of two genomeIDs
#' @param genome2 character string specifying second of two genomeIDs
#' @param gff annotated gff with orthogroups included, see read_gff
#' @param synBuff integer of length 1 specifying the maximum euclidean distance
#' from an 'anchor' so that it can be considered syntenic
#' @param selfOnly logical, should only self hits be considered
#' @param overwrite logical, should the results be overwrittem?
#' @param onlyOgAnchors logical, should anchors be restricted to hits in the
#' same orthogroup? See set_syntenyParams.
#' @param nhits1 Integer, how many hits should be kept for each gene in genome1.
#' See set_syntenyParams.
#' @param nhits2 Integer, how many hits should be kept for each gene in genome2.
#' See set_syntenyParams.
#' @param maskHits data.table of hits that should be excluded
#' @param synParam data.table with synteny parameters. See set_syntenyParams.
#' @param selfRegionMask integer, the radius around self hits that should be
#' masked
#' @param nhits integer, the number of hits to retain
#' @param blks data.table containing the block coordinates
#' @param blastDir file.path to the location of the blast results.
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
  isAnchor <- gen1 <- gen2 <- chr1 <- chr2 <-  NULL

  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
  }else{
    if(!any(genomeIDs %in% gsParam$genomes$genomeIDs))
      stop("specified genomeIDs dont look right\n")
  }
  genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]

  verbose <- gsParam$params$verbose
  syn <- data.table(gsParam$params$synteny)
  writeTo <- gsParam$paths$results
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  blksFile <- file.path(gsParam$paths$results, "syntenicBlocks.txt.gz")

  # -- enforce genome order from genomeIDs, keep only uni-directional blasts
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
  if(verbose)
    cat("Setting up input data\n\tParsing annotations ... \n\t")

  # -- find the location of orthofinder results files
  gsParam <- find_orthofinderResults(gsParam)

  # -- annotate gff with arrays, syntenic orthogroups, etc.
  gff <- annotate_gff(gsParam = gsParam, genomeIDs = genomeIDs)

  # -- flag syntenic arrays and choose optimal representative gene
  isArrayRep <- ord <- NULL
  if(verbose)
    with(gff, cat(sprintf(
      "Found %s genes, %s orthogroups and %s arrays with %s genes\n",
      nrow(gff), uniqueN(gff$og),
      sum(!isArrayRep[!duplicated(synArray)]),
      sum(!isArrayRep[!duplicated(synArray)]) + sum(duplicated(synArray)))))
  gffAll <- data.table(gff)
  gff <- subset(gff, isArrayRep)
  gff[,ord := frank(ord, ties.method = "dense"), by = "genome"]
  ##############################################################################
  # 2. self culling
  # -- run the following steps for each genome and only within genome hits
  ##############################################################################
  if(verbose)
    cat("Pulling within-genome synteny\n\tGenome:",
        "n raw hits / hits in (regions) / anchors / hits in (blks)\n")
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
      phits <- subset(ahits, chr1 == chr2)
      cat(sprintf("%s / ", nrow(phits)))

      # -- self synteny
      sr <- pull_selfRegion(
        hits = subset(phits, nChr1 >= selfSyn$blkSize[i] & nChr2 >= selfSyn$blkSize[i]),
        synBuff = selfSyn$synBuff[i])
      blkv <- sr$blkID; names(blkv) <- with(sr, paste(ofID1, ofID2))
      phits[,`:=`(blkID = blkv[paste(ofID1, ofID2)],
                  blkAnchor = ofID1 == ofID2,
                  blkBuffer = paste(ofID1, ofID2) %in% names(blkv),
                  regID = blkv[paste(ofID1, ofID2)],
                  regAnchor = ofID1 == ofID2,
                  regBuffer = paste(ofID1, ofID2) %in% names(blkv))]
      cat(sprintf("%s (%s) / %s / %s (%s)\n",
                  nrow(phits), uniqueN(phits$regID),
                  sum(phits$ofID1 == phits$ofID2),
                  nrow(phits), uniqueN(phits$regID)))

      ############################################################################
      # 2.2 homeologs, if desired
      ############################################################################
      if(noParalog){
        u <- with(phits, unique(paste(ofID1, ofID2)))
        shits <- subset(ahits, !paste(ofID1, ofID2) %in% u)
      }else{
        if(verbose)
          cat(sprintf("\t%s (homeolg): ",  pull_strWidth(genome1, 7)))
        shits <- with(selfSyn[i,], pipe_synteny(
          gsParam = gsParam,
          gff = gff,
          maskHits = subset(phits, !is.na(regBuffer)),
          genome1 = genome1,
          genome2 = genome1,
          nhits1 = (nhits1-1 + nSecondHits1),
          nhits2 = (nhits2-1 + nSecondHits2),
          blkSize = blkSizeSecond,
          nGaps = nGapsSecond,
          synBuff = synBuffSecond,
          onlyOgAnchors = onlyOgAnchorsSecond,
          dropInterleavesSmallerThan = dropInterleavesSmallerThan,
          verbose = verbose))
        shits$regID[!is.na(shits$regID)] <- paste0("sec_", shits$regID)[!is.na(shits$regID)]
        shits$blkID[!is.na(shits$blkID)] <- paste0("sec_", shits$blkID)[!is.na(shits$blkID)]
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
      pd <- plot_hits(
        gsParam = gsParam,
        hits = subset(outHits, score > 100),
        alpha = .025,
        gff = gff,
        onlyOg = F,
        cols = "blue2",
        bufferOnly = F,
        anchorOnly = F,
        round2 = 5,
        plotTitle = "All diamond hits score > 100")
      pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = T,
        plotTitle = "hits in syntenic regions")
      pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = F, anchorOnly = T,
        plotTitle = "anchor hits in blocks")
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
    cat("Pulling intergenomic synteny\n")

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
      phits <- with(syn, pipe_synteny(
        gsParam = gsParam,
        gff = gff,
        maskHits = NULL,
        genome1 = genome1,
        genome2 = genome2,
        nhits1 = nhits1,
        nhits2 = nhits2,
        blkSize = blkSize,
        nGaps = nGaps,
        synBuff = synBuff,
        onlyOgAnchors = onlyOgAnchors,
        dropInterleavesSmallerThan = dropInterleavesSmallerThan,
        verbose = verbose))

      # -- optionally get secondary hits too.
      if(nonSelfSyn$nSecondHits[i] == 0){
        shits <- NULL
      }else{
        if(verbose)
          cat(sprintf(
            "\t\t%s-%s (secondy): ", pull_strWidth(genome1, 7), pull_strWidth(genome2, 7)))
        shits <- with(syn, pipe_synteny(
          gsParam = gsParam,
          gff = gff,
          maskHits = subset(phits, !is.na(regID)),
          genome1 = genome1,
          genome2 = genome2,
          nhits1 = nSecondHits1,
          nhits2 = nSecondHits2,
          blkSize = blkSizeSecondSize,
          nGaps = nGapsSecond,
          synBuff = synBuffSecond,
          onlyOgAnchors = onlyOgAnchorsSecond,
          dropInterleavesSmallerThan = dropInterleavesSmallerThan,
          verbose = verbose))
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

      pd <- plot_hits(
        gsParam = gsParam,
        hits = subset(outHits, score > 100),
        gff = gff,
        alpha = .025,
        onlyOg = F,
        cols = "blue2",
        bufferOnly = F,
        anchorOnly = F,
        round2 = 5,
        plotTitle = "All diamond hits score > 100")
      pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = T,
        plotTitle = "hits in syntenic regions")
      pd <- plot_hits(
        gsParam = gsParam,
        hits = outHits, gff = gff, onlyOg = T, round2 = 10,
        bufferOnly = T, plotRegions = F, anchorOnly = T,
        plotTitle = "anchor hits in blocks")
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
    cat("Adding syntenic orthogroups to gff ... ")
  gffAll <- add_synOg2gff(
    gff = gffAll,
    useBlks = F,
    gsParam = gsParam,
    genomeIDs = genomeIDs,
    allowRBHinOg = T)
  fwrite(gffAll, sep = "\t", quote = F, showProgress = F, file = gffFile)
  if(verbose)
    cat("Done\n")
  return(blks)
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
  blks1 <- bhits[,list(startBp1 = min(start1), endBp1 = max(end1),
                       startOrd1 = min(ord1), endOrd1 = max(ord1),
                       firstGene1 = first(ofID1), lastGene1 = last(ofID1),
                       nHits1 = uniqueN(ofID1)),
                 by = c("blkID", "gen1","gen2", "chr1", "chr2")]
  setkey(bhits, ord2)
  blks2 <- bhits[,list(minBp2 = min(start2), maxBp2 = max(end2),
                       minOrd2 = min(ord2), maxOrd2 = max(ord2),
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
                            dropInterleavesSmallerThan,
                            blkSize,
                            nGaps,
                            synBuff,
                            path2mcscanx,
                            nCores,
                            minRbhScore){

  isOg <- maxScr1 <- maxScr2 <- score <- og <- ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  isCollin <- blkID <- bord1 <- bord2 <- gen1 <- gen2 <- nblk <- regID <- NULL
  regAnchor <- blkIDn <- chr1 <- chr2 <- isAnchor <- NULL

  ##############################################################################
  # 1. Pull anchor hits and re-rank
  ##############################################################################
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
  fblks <- fblks[,colnames(hits), with = F]
  outBlks <- split_ovlpBlks(
    hits = fblks,
    dropInterleavesSmallerThan = dropInterleavesSmallerThan,
    verbose = F)
  outBlks[,blkIDn := as.numeric(as.factor(blkID)), by = c("chr1", "chr2")]
  outBlks[,blkID := sprintf("%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, blkIDn)]
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
  gv <- gff$genome; cv <- gff$chr; ogv <- gff$og
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

#' @title split_ovlpBlks
#' @description
#' \code{split_ovlpBlks} split_ovlpBlks
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
split_ovlpBlks <- function(hits,
                           dropInterleavesSmallerThan,
                           verbose){
  pull_ovlNoDupBlks <- function(hits,
                                chrCol,
                                genomeCol,
                                ordcol){
    ordx <- genomex <- chrx <- start <- end <- blkID <- i.blkID <- anyDup <- NULL
    u <- ovlStart <- i.start <- i.end <- ovlEnd <- NULL
    blks <- data.table(hits)
    blks[,`:=`(ordx = blks[[ordcol]], chrx = blks[[chrCol]],
               genomex = blks[[genomeCol]])]
    spl <- split(blks, by = "blkID")
    b <- blks[,list(start = min(ordx), end = max(ordx)),
              by = c("chrx","blkID","genomex")]
    setkey(b, genomex, chrx, start, end)
    bo <- subset(foverlaps(b, b), blkID != i.blkID)
    if(nrow(bo) > 0){
      bo[,anyDup := unlist(lapply(1:nrow(bo), function(i){
        x <- unique(spl[[as.character(bo$blkID[i])]]$ofID1)
        y <- unique(spl[[as.character(bo$i.blkID[i])]]$ofID1)
        x2 <- unique(spl[[as.character(bo$blkID[i])]]$ofID2)
        y2 <- unique(spl[[as.character(bo$i.blkID[i])]]$ofID2)
        return(any(x %in% y) || any(y %in% x) || any(x2 %in% y2) || any(y2 %in% x2))
      }))]
      chk <- subset(bo, !anyDup)
      chk[,u := paste(c(blkID, i.blkID)[order(c(blkID, i.blkID))], collapse = ","),
          by = c("blkID", "i.blkID")]
      chk <- subset(chk, !duplicated(u))
      chk[,ovlStart := c(start, end, i.start, i.end)[order(c(start, end, i.start, i.end))][2],
          by = c("blkID", "i.blkID")]
      chk[,ovlEnd := c(start, end, i.start, i.end)[order(c(start, end, i.start, i.end))][3],
          by = c("blkID", "i.blkID")]
      return(chk)
    }else{
      return(NULL)
    }
  }

  blks <- end1 <- end2 <- start1 <- start2 <- blk1 <- blk2 <- ord1 <- ord2 <- NULL
  nblk <- no1 <- no2 <- blkID <- ofID1 <- ofID2 <- NULL
  tmp <- data.table(hits)

  if(dropInterleavesSmallerThan > 1){
    minPts <- dropInterleavesSmallerThan-1
    eps <- sqrt(((dropInterleavesSmallerThan-1)^2) * 2)
  }else{
    minPts <- 0
    eps <- sqrt(2)
  }

  novlReg <- 1
  iter <- 0
  nb1 <- 0
  nb2 <- 1
  while(novlReg > 0 & nb2 > nb1){
    iter <- iter + 1
    nb1 <- uniqueN(tmp$blkID[!is.na(tmp$blkID)])
    if(verbose)
      cat(sprintf("\t\tIter %s: %s blks ... ",iter, nb1))
    ovl1 <- pull_ovlNoDupBlks(
      hits = tmp, chrCol = "chr1", genomeCol = "gen2", ordcol = "ord1")
    if(!is.null(ovl1)){
      ovl1 <- with(ovl1, data.table(blks = u, start1 = ovlStart, end1 = ovlEnd))
    }

    ovl2 <- pull_ovlNoDupBlks(
      hits = tmp, chrCol = "chr2", genomeCol = "gen2", ordcol = "ord2")
    if(!is.null(ovl2)){
      ovl2 <- with(ovl2, data.table(blks = u, start2 = ovlStart, end2 = ovlEnd))
    }
    if(!is.null(ovl1) & !is.null(ovl2)){
      if(nrow(ovl1) > 0 & nrow(ovl2) > 0){
        ovl <- merge(ovl1, ovl2, by = "blks", all = T)
        novlReg = nrow(ovl)
        if(novlReg > 0){
          ovl[,c("blk1", "blk2") := tstrsplit(blks, ",")]
          ovl[,width := ifelse(is.na(end1), end2 - start2,
                               ifelse(is.na(end2), end1 - start1,
                                      (end1 - start1) + (end2 - start2)))]
          setorder(ovl, -width)
          ovl <- subset(ovl, !duplicated(blk1))
          ovl <- subset(ovl, !duplicated(blk2))
          ovl <- subset(ovl, !blk2 %in% blk1)
          setorder(ovl, -width)
          if(verbose)
            cat(sprintf("%s ovlps", nrow(ovl)))

          spl <- split(tmp, by = "blkID")

          reblk <- rbindlist(lapply(1:nrow(ovl), function(i){
            x <- ovl[i,]
            y <- rbind(data.table(spl[[x$blk1]]), data.table(spl[[x$blk2]]))
            y[,`:=`(no1 = frank(ord1, ties.method = "dense"),
                    no2 = frank(ord2, ties.method = "dense"))]
            y[,nblk := dbscan(frNN(cbind(no1, no2), eps = eps), minPts = 2)$cluster]
            y[,blkID := sprintf("%s_%s", blkID, nblk)]
            y$blkID[y$nblk == 0] <- NA
            y[,`:=`(no1 = NULL, no2 = NULL, nblk = NULL)]
            return(y)
          }))

          tmp <- subset(tmp, !paste(ofID1, ofID2) %in% paste(reblk$ofID1, reblk$ofID2))
          tmp <- rbind(tmp, reblk)
          nb2 <- uniqueN(tmp$blkID[!is.na(tmp$blkID)])
          if(verbose)
            cat(sprintf(" --> %s split blks / %s hits dropped\n",
                        uniqueN(tmp$blkID[!is.na(tmp$blkID)]),
                        sum(is.na(tmp$blkID))))
        }else{
          nb2 <- uniqueN(tmp$blkID[!is.na(tmp$blkID)])
          if(verbose)
            cat('no more overlaps')
        }
      }
    }else{
      nb2 <- uniqueN(tmp$blkID[!is.na(tmp$blkID)])
      novlReg <- 0
      if(verbose)
        cat('no more overlaps')
    }
  }
  return(tmp)
}

#' @title pull_synReg
#' @description
#' \code{pull_synReg} pull_synReg
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
pull_synReg <- function(hits,
                        onlyOgAnchors,
                        nhits1,
                        nhits2,
                        blkSize,
                        nGaps,
                        path2mcscanx,
                        dropInterleavesSmallerThan,
                        nCores,
                        synBuff){

  add_rbh2og <- function(hits, minRbhScore = 50){
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

  og <- ofID1 <- score <- ofID2 <- nhits1 <- nhits2 <- ord1 <- ord2 <- NULL
  scrRank1 <- scrRank2 <- isAnchor <- NULL
  hits <- data.table(hits)

  # -- subset to just ogs if specified
  if(onlyOgAnchors){
    anch <- subset(hits, !is.na(og))
  }else{
    anch <- data.table(hits)
  }

  ##############################################################################
  # 1. Initial region building on culled topscore hits
  ##############################################################################
  # -- calculate score rank, subset to topn and re-order genes

  setorder(anch, ofID1, -score)
  anch[,scrRank1 := 1:.N, by = "ofID1"]
  setorder(anch, ofID2, -score)
  anch[,scrRank2 := 1:.N, by = "ofID2"]
  anch <- subset(anch, scrRank1 <= nhits1 & scrRank2 <= nhits2)
  anch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]

  # -- initial anchor building
  anch[,isAnchor := !is.na(run_mcscanx(
    hits = anch,
    blkSize = blkSize,
    nGaps = nGaps,
    path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)])]

  ##############################################################################
  # 2. Secondary culling with all hits in regions around initial anchors
  ##############################################################################
  # -- subset hits to those near the anchors
  chr1 <- chr2 <- NULL
  ua <- with(subset(anch, isAnchor), unique(paste(ofID1, ofID2)))
  u <- with(anch, paste(chr1, chr2))
  regs <- subset(hits, paste(chr1, chr2) %in% u)
  regs[,isAnchor := paste(ofID1, ofID2) %in% ua]
  regs[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
  regs <- find_hitsInBuff(
    hits = regs,
    nCores = nCores,
    synBuff = synBuff)
  bufv <- with(regs, paste(ofID1, ofID2))

  # -- initial large block clustering
  regBlk <- clus_dbscan(
    hits = regs,
    radius = synBuff,
    blkSize = blkSize,
    nCores = nCores)

  # -- find RBHs in the larger blocks, add these as orthogroups
  blkID <- inBuffer <- NULL
  regBlk <- add_rbh2og(
    hits = subset(regBlk, !is.na(blkID) & inBuffer))

  # -- re-rank and re-cluster anchors
  regAnch <- subset(regBlk, !is.na(og))
  regAnch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense"))]
  regAnch[,isAnchor := !is.na(run_mcscanx(
    hits = regAnch,
    blkSize = blkSize,
    nGaps = nGaps,
    path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)])]

  # -- re-rank and cluster anchors into blocks
  regAnch <- subset(regAnch, isAnchor)
  regAnch[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                ord2 = frank(ord2, ties.method = "dense"))]
  regsAnch <- clus_dbscan(
    hits = regAnch,
    radius = synBuff,
    blkSize = blkSize,
    nCores = nCores)
  ancv <- with(regsAnch, paste(ofID1, ofID2))

  # -- split blocks if there are interleaves
  regBlks <- split_ovlpBlks(
    hits = regsAnch, verbose = F,
    dropInterleavesSmallerThan = dropInterleavesSmallerThan)
  regBlks <- subset(regBlks, !is.na(blkID))
  blv <- regBlks$blkID; names(blv) <- with(regBlks, paste(ofID1, ofID2))

  ##############################################################################
  # 3. Finalize the region coordinates / hits
  ##############################################################################
  # -- find hits in the regs
  hits[,blkID := blv[paste(ofID1, ofID2)]]
  bc <- calc_blkCoords(subset(hits, !is.na(blkID)))
  inreg <- get_hitsInBlks(blks = bc, hits = hits, nCores = nCores)
  blv <- inreg$blkID; names(blv) <- with(inreg, paste(ofID1, ofID2))

  # -- get vector of new OGs including RBHs
  regOgv <- regBlk$og; names(regOgv) <- with(regBlk, paste(ofID1, ofID2))
  allOgv <- hits$og; names(allOgv) <- with(hits, paste(ofID1, ofID2))
  allOgv <- allOgv[!is.na(allOgv)]
  newOgv <- regOgv[!names(regOgv) %in% names(allOgv)]
  ogv <- c(allOgv, newOgv)

  # -- return a new hits obj
  hits[,`:=`(blkID = NA, blkAnchor = NA, blkBuffer = NA,
             og = ogv[paste(ofID1, ofID2)],
             regID = blv[paste(ofID1, ofID2)],
             regAnchor = paste(ofID1, ofID2) %in% ancv,
             regBuffer = paste(ofID1, ofID2) %in% bufv)]

  return(hits)
}

#' @title pipe_synteny
#' @description
#' \code{pipe_synteny} pipe_synteny
#' @rdname synteny
#' @import data.table
#' @export
pipe_synteny <- function(gsParam,
                         genome1,
                         genome2,
                         gff,
                         nhits1,
                         nhits2,
                         blkSize,
                         nGaps,
                         synBuff,
                         maskHits,
                         synParam,
                         onlyOgAnchors,
                         dropInterleavesSmallerThan,
                         verbose){

  # -- subset gff to array representatives and re-rank order
  genome <- gen1 <- gen2 <- isAnchor <- NULL
  gf1 <- subset(gff, genome == genome1)
  gf2 <- subset(gff, genome == genome2)

  bl <- parse_blast4synteny(
    gsParam = gsParam,
    genome1 = genome1,
    genome2 = genome2,
    gff = gff,
    selfOnly = FALSE)
  bl <- subset(bl, !is.na(gen1) & !is.na(gen2))

  # -- optionally remove masked regions
  chr1 <- chr2 <- ofID1 <- ofID2 <- NULL
  if(!is.null(maskHits)){
    uchr <- with(maskHits, paste(chr1, chr2))
    tmp <- subset(bl, paste(chr1, chr2) %in% uchr)
    nb <- subset(bl, !paste(chr1, chr2) %in% uchr)
    u <- with(maskHits, unique(paste(ofID1, ofID2)))
    tmp[,isAnchor := paste(ofID1, ofID2) %in% u]
    ib <- find_hitsInBuff(
      hits = tmp,
      nCores = gsParam$params$nCores,
      synBuff = synBuff)
    u <- with(ib, unique(paste(ofID1, ofID2)))
    bl <- rbind(nb, subset(bl, !paste(ofID1, ofID2) %in% u))
  }

  if(verbose)
    cat(sprintf("%s / ", nrow(bl)))

  # -- get the syntenic region hits
  regs <- pull_synReg(
    hits = data.table(bl),
    onlyOgAnchors = onlyOgAnchors,
    synBuff = synBuff,
    nhits1 = nhits1,
    nhits2 = nhits2,
    blkSize = blkSize,
    nGaps = nGaps,
    dropInterleavesSmallerThan = dropInterleavesSmallerThan,
    path2mcscanx = gsParam$paths$mcscanxCall,
    nCores = gsParam$params$nCores)
  if(verbose)
    cat(sprintf("%s (%s) / ", sum(regs$regBuffer), uniqueN(regs$regID, na.rm = T)))

  # -- get the syntenic blocks
  blks <- finalize_blocks(
    hits = regs,
    dropInterleavesSmallerThan = dropInterleavesSmallerThan,
    blkSize = blkSize,
    nGaps = nGaps,
    synBuff = synBuff,
    minRbhScore = 50,
    path2mcscanx = gsParam$paths$mcscanxCall,
    nCores = gsParam$params$nCores)

  # -- reformat and return
  if(verbose)
    cat(sprintf("%s (%s)\n", sum(blks$blkBuffer), uniqueN(blks$blkID, na.rm = T)))

  return(blks)
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
      abs(ord1 - ord2) <= buff)
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

#' @title add orthofinder ID to a gff object
#' @description
#' \code{add_ofID2gff} read the orthofinder species and gene IDs and merge
#' these with the gff-like data.table
#' @rdname synteny
#' @export
add_ofID2gff <- function(gff, blastDir){
  id <- ofID <- genomeNum <- genome <- NULL
  specIDs <- read_orthofinderSpeciesIDs(blastDir)
  gv <- names(specIDs); names(gv) <- as.character(specIDs)
  seqIDs <- read_orthofinderSequenceIDs(blastDir)
  seqIDs[,genome :=  gv[as.character(genomeNum)]]
  idv <- seqIDs$ofID; names(idv) <- with(seqIDs, paste(genome, id))
  gff[,ofID := idv[paste(genome, id)]]
  return(gff)
}

#' @title add peptide length to a gff object
#' @description
#' \code{add_pepLen2gff} read the peptide lengths and merge with the gff
#' @rdname synteny
#' @export
add_pepLen2gff <- function(gff, gsParam){
  pepLen <- id <- ofID <-  NULL
  spl <- split(gff, by = "genome")
  naa <- rbindlist(lapply(spl, function(x){
    x[,pepLen := get_nAA(gsParam$paths$peptide[x$genome[1]], raw = T)[id]]
    return(x[,c("ofID","pepLen")])
  }))
  nao <- naa$pepLen; names(nao) <- naa$ofID
  gff[,pepLen := nao[ofID]]
  return(gff)
}

#' @title add array representative to a gff object
#' @description
#' \code{add_arrayRep2gff} choose most central gene by orthogroup
#' @rdname synteny
#' @importFrom parallel mclapply
#' @export
add_arrayRep2gff <- function(gff, gsParam){
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

#' @title add number of anchors to gff
#' @description
#' \code{annotate_gff} add number of anchors to gff
#' @rdname synteny
#' @importFrom parallel mclapply
#' @export
annotate_gff <- function(gsParam, genomeIDs){
  allowRBHinOg <- any(gsParam$genomes$ploidy > 1)
  nCores <- gsParam$params$nCores
  verbose <- gsParam$params$verbose
  d2h <- gsParam$params$maxDistBtwPgHits

  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")

  # -- find the location of orthofinder results files
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)

  # -- load the orthogroups
  ogs <- parse_ogs(gsParam)

  # -- load the gff as a data.table and add the orthofinder gene IDs to it
  genome <- ord <- NULL
  gff <- add_ofID2gff(read_gff(gsParam$paths$gff), gsParam$paths$blastDir)
  gff <- subset(gff, genome %in% genomeIDs)
  gff[,genome := factor(genome, levels = genomeIDs)]
  setkey(gff, genome, ord)

  # -- add orthogroups to the gff
  gff <- merge(gff, ogs, by = c("genome","id"), all.x = T)

  # -- populate orthogroup column with ogIDs for genes without an OG
  gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
  setnames(gff, "ogID", "og")

  gff <- add_arrayRep2gff(gff = gff, gsParam = gsParam)

  return(gff)
}

#' @title add_synOg2gff
#' @description
#' \code{add_synOg2gff} add_synOg2gff
#' @rdname plot_riparian
#' @export
add_synOg2gff <- function(gff,
                          hits = NULL,
                          gsParam,
                          genomeIDs,
                          allowRBHinOg,
                          useBlks){
  blkBuffer <- isSelf <- ofID1 <- ofID2 <- og <- ofID <- synOg <- blkID <- NULL
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
      x <- subset(x, !isSelf)
      if(!allowRBHinOg)
        x <- subset(x, !grepl("RBH", og))
      return(x[,c("ofID1", "ofID2", "blkAnchor", "blkID")])
    }))
  }else{
    hts <- data.table(hits)
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
  gff[,synOg := ic[ofID]]
  nmis <- sum(is.na(gff$synOg))
  mol <- max(gff$synOg, na.rm = T)
  gff$synOg[is.na(gff$synOg)] <- (mol + 1):(mol + nmis)
  gff[,synOg := as.integer(factor(paste(og, synOg), levels = unique(paste(og, synOg))))]
  return(gff)
}

