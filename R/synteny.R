#' @title Parsing of pairwise hits into synteny
#' @description
#' \code{synteny} The main GENESPACE engine to constrain genome-wide homology
#' hits to synteny.
#'
#' @param gsParam a list containing all parameters for a GENESPACE run. See
#' init_genespace
#' @param genomeIDs an optional vector of genomeIDs to consider. If not
#' specified (default) taken from gsParam$genomeIDs$genomeIDs
#'
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
#' @return updated gsParam list
#'
#' @import R.utils
#' @import data.table
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
