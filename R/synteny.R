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
#' @param hits data.table of hits
#' @param nGaps integer of length 1, specifying the -m param to mcscanx
#' for the primary MCScanX run. This acts on the results from the initial
#' MCScanX run.
#' @param blkSize integer of length 1, specifying the -s param to mcscanx
#' @param nSecondHits integer of length 1, specifying the number of blast
#' hits to include after masking.
#' @param synBuff Numeric > 0, specifying the distance from an anchor
#' to consider a hit syntenic. This parameter is also used to limit the search
#' radius in dbscan-based blk calculation. Larger values will return larger
#' tandem arrays but also may permit inclusion of spurious non-syntenic networks
#' @param nGapsSecond see nGaps, but passed to secondary hits after masking
#' primary hits.
#' @param blkSizeSecond see blkSize, but passed to the secondary scan if
#' nSecondaryHits > 0.
#' @param synBuffSecond see syntenyBuffer. Applied only to synteny
#' construction of secondary hits.
#' @param onlyOgAnchors logical, should only hits in orthogroups be considered
#' for anchors?
#' @param onlyOgAnchorsSecond logical should only hits in orthogroups be
#' considered for anchors in secondary blocks?
#' @param selfRegionMask integer specifying the size of the region (radius, gene
#' rank order) surrounding self hits to mask for secondary/homeologous hits
#' @param dropInterleavesSmallerThan integer specifying the smallest block to
#' keep in the split interleaves synteny step
#' @param minRbhScore integer specifying the minimum blast bit score to allow
#' for a RBH to be included
#' @param path2mcscanx character string file.path pointing to the install
#' directory for the MCScanX program. In particular, there must be an
#' executable in this directory for MCScanX_h.
#' @param overwrite logical, should existing directories be overwritten?
#' @param ... additional arguments passed to set_syntenyParam().
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
#' \cr
#' If called, \code{synteny} returns its own arguments.
#'
#' @examples
#' \dontrun{
#'
#' runwd <- file.path(getwd(), "testGenespace")
#' make_exampleDataDir(writeDir = runwd)
#'
#' gpar <- init_genespace(
#'   genomeIDs = c("human","chimp","rhesus"),
#'   speciesIDs = c("human","chimp","rhesus"),
#'   versionIDs = c("human","chimp","rhesus"),
#'   ploidy = rep(1,3),
#'   diamondMode = "fast",
#'   orthofinderMethod = "fast",
#'   wd = runwd,
#'   nCores = 4,
#'   minPepLen = 50,
#'   gffString = "gff",
#'   pepString = "pep",
#'   path2orthofinder = "orthofinder",
#'   path2mcscanx = "~/MCScanX",
#'   rawGenomeDir = file.path(runwd, "rawGenomes"))
#'
#' parse_annotations(
#'   gsParam = gpar,
#'   gffEntryType = "gene",
#'   gffIdColumn = "locus",
#'   gffStripText = "locus=",
#'   headerEntryIndex = 1,
#'   headerSep = " ",
#'   headerStripText = "locus=")
#'
#' gpar <- run_orthofinder(gsParam = gpar, overwrite = F)
#'
#' gpar <- synteny(gsParam = gpar) # use defaults
#'
#' # -- run again (need to set overwrite = T) with blkSize = 10
#' gpar$params$synteny <- NULL
#' gpar <- synteny(gsParam = gpar, blkSize = 10, overwrite = T)
#'
#' # -- run again with custom specs
#' # **NOTE** if params$synteny is a data.table, synteny will respect it. So,
#' # this can be generated and modified before running. If it is NULL (as above)
#' #synteny will calculate parameters internally with additional arguments.
#'
#' # here, increase the synteny buffer and
#' # make the blocks need to be bigger between human and chimp
#' gpar <- set_syntenyParams(gpar, synBuff = 200)
#' wh <- with(gpar$params$synteny, which(
#'   genome1 %in% c("humnan", "chimp") & genome2 %in% c("humnan", "chimp")))
#' gpar$params$synteny$blkSize[wh] <- 10
#' gpar <- synteny(gsParam = gpar, overwrite = T)
#' }
#'
#'
#' @rdname synteny
#' @import R.utils
#' @import data.table
#' @export
synteny <- function(gsParam, genomeIDs = NULL, overwrite = F, ...){
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
    tmp <- set_syntenyParams(gsParam, ...)
    gsParam <- set_syntenyParams(gsParam)
    if(identical(tmp$params$synteny, gsParam$params$synteny)){
      cat("Synteny Parameters have not been set! Setting to defaults\n")
    }else{
      cat("Synteny Parameters have not been set! Using user-defined settings\n")
      gsParam$params$synteny <- tmp$params$synteny
    }
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
  gff <- pull_synOGs(gsParam = gsParam)
  gffo <- combine_inblkSynOG(
    genomeIDs = genomeIDs,
    gff = gff,
    gsParam = gsParam)
  fwrite(gff, file = gffFile, sep = "\t", quote = F, showProgress = F)
  if(verbose)
    cat("\tWrote gff to file: /results/gffWithOgs.txt.gz\n\tDone!\n")
  return(gsParam)
}

#' @title Set synteny parameters
#' @description
#' \code{set_syntenyParams} Calculate paramters for pairwise synteny search.
#' Other synteny functions require this as input.
#' @rdname synteny
#' @importFrom Biostrings readAAStringSet
#' @import data.table
#' @export
set_syntenyParams <- function(gsParam,
                              onlyOgAnchors = TRUE,
                              onlyOgAnchorsSecond = TRUE,
                              blkSize = 5,
                              blkSizeSecond = blkSize,
                              nGaps = 5,
                              nGapsSecond = nGaps,
                              nSecondHits = 0,
                              synBuff = 100,
                              synBuffSecond = synBuff,
                              selfRegionMask = synBuff * 2,
                              dropInterleavesSmallerThan = 2,
                              minRbhScore = 50){

  ##############################################################################
  # 1. setup environment
  ##############################################################################

  genome1 <- genome2 <- tiebreak <- nGenes1 <- nGenes2 <- query <- target <- NULL
  nSecondHits1 <- nSecondHits2 <- nhits1 <- nhits2 <- NULL

  # -- check the parameters
  onlyOgAnchors <- as.logical(onlyOgAnchors[1])
  if(is.null(onlyOgAnchors) || is.na(onlyOgAnchors) || length(onlyOgAnchors) == 0)
    stop("cannot coerce onlyOgAnchors to logical\n")
  onlyOgAnchorsSecond <- as.logical(onlyOgAnchorsSecond[1])
  if(is.null(onlyOgAnchorsSecond) || is.na(onlyOgAnchorsSecond) || length(onlyOgAnchorsSecond) == 0)
    stop("cannot coerce onlyOgAnchors to logical\n")

  blkSize <- as.integer(blkSize)[1]
  if(is.null(blkSize) || is.na(blkSize) || length(blkSize) == 0)
    stop("cannot coerce blkSize to integer\n")
  blkSizeSecond <- as.integer(blkSizeSecond)[1]
  if(is.null(blkSizeSecond) || length(blkSizeSecond) == 0)
    stop("cannot coerce blkSizeSecond to integer\n")

  nGaps <- as.integer(nGaps)[1]
  if(is.null(nGaps) || is.na(nGaps) || length(nGaps) == 0)
    stop("cannot coerce nGaps to integer\n")
  nGapsSecond <- as.integer(nGapsSecond)[1]
  if(is.null(nGapsSecond) || length(nGapsSecond) == 0)
    stop("cannot coerce nGapsSecond to integer\n")

  synBuff <- as.integer(synBuff)[1]
  if(is.null(synBuff) || is.na(synBuff) || length(synBuff) == 0)
    stop("cannot coerce synBuff to integer\n")
  synBuffSecond <- as.integer(synBuffSecond)[1]
  if(is.null(synBuffSecond) || is.na(synBuffSecond) || length(synBuffSecond) == 0)
    stop("cannot coerce synBuffSecond to integer\n")

  selfRegionMask <- as.integer(selfRegionMask)[1]
  if(is.null(selfRegionMask) || is.na(selfRegionMask) || length(selfRegionMask) == 0)
    stop("cannot coerce selfRegionMask to integer\n")
  nSecondHits <- as.integer(nSecondHits)[1]
  if(is.null(nSecondHits) || is.na(nSecondHits) || length(nSecondHits) == 0)
    stop("cannot coerce nSecondHits to integer\n")
  dropInterleavesSmallerThan <- as.integer(dropInterleavesSmallerThan)[1]
  if(is.null(dropInterleavesSmallerThan) || is.na(dropInterleavesSmallerThan) || length(dropInterleavesSmallerThan) == 0)
    stop("cannot coerce dropInterleavesSmallerThan to integer\n")
  minRbhScore <- as.integer(minRbhScore)[1]
  if(is.null(minRbhScore) || is.na(minRbhScore) || length(minRbhScore) == 0)
    stop("cannot coerce minRbhScore to integer\n")

  # -- shorten some names
  genomeIDs <- gsParam$genomes$genomeIDs
  ploidy <- gsParam$genomes$ploidy[genomeIDs]
  blastDir <- gsParam$paths$orthofinder
  pepFiles <- gsParam$paths$peptide

  # -- Check the peptide files
  if(!all(file.exists(pepFiles))){
    wh <- which(!file.exists(pepFiles))
    stop(sprintf("the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
                 paste(names(pepFiles)[wh], collapse = "\n\t")))
  }

  # -- Calculate the number of genes for each genome
  nGenes <- sapply(pepFiles, function(x) length(readAAStringSet(x)))
  if(min(nGenes) == 0){
    wh <- which(nGenes == 0)
    stop(sprintf("the following genome annotations have not been correctly parsed:\n\t%s\n\tRerun parse_annotations\n",
                 paste(names(pepFiles)[wh], collapse = "\n\t")))
  }

  ##############################################################################
  # 2. build the pairwise parameter data.table scaffold
  ##############################################################################

  # make the database of unique combinations of genomes
  cmb <- data.table(
    genome1 = names(ploidy),
    genome2 = names(ploidy))
  cmb <- cmb[,CJ(genome1, genome2)]

  # -- add ploidy, orthofinder IDs, n genes
  speciesIDs <- (1:length(genomeIDs))-1
  names(speciesIDs) <- genomeIDs[order(genomeIDs)]
  cmb[,`:=`(u = paste(genome1, genome2),
            ploidy1 = ploidy[genome1],
            ploidy2 = ploidy[genome2],
            ofID1 = speciesIDs[genome1],
            ofID2 = speciesIDs[genome2],
            nGenes1 = nGenes[genome1],
            nGenes2 = nGenes[genome2],
            nhits1 = ploidy[genome2],
            nhits2 = ploidy[genome1])]

  # -- Choose query and target genomes, based on the total number of genes
  cmb[,tiebreak := as.numeric(factor(genome1, levels = genomeIDs)) -
        as.numeric(factor(genome2, levels = genomeIDs))]
  cmb[,`:=`(query = ifelse(nGenes1 > nGenes2, genome1,
                           ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                  genome1, genome2)),
            target = ifelse(nGenes1 > nGenes2, genome2,
                            ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                   genome2, genome1))),
      by = "u"]
  cmb[,`:=`(runBlast = genome1 == query,
            mirrorBlast = genome1 == target & genome1 != query,
            tiebreak = NULL)]

  # -- choose the number of secondary hits
  cmb[,nSecondHits1 :=
        ifelse(nhits1 == 1 & nhits2 == 1,
               nSecondHits,
               ifelse(genome1 == genome2 & nhits2 > 1 & nSecondHits > 0,
                      (nSecondHits * nhits2/2) + (nhits2/2),
                      ifelse(nhits2 > 1 & nSecondHits > 0,
                             nSecondHits * nhits2,
                             nSecondHits)))]
  cmb[,nSecondHits2 :=
        ifelse(nhits1 == 1 & nhits2 == 1,
               nSecondHits,
               ifelse(genome1 == genome2 & nhits1 > 1 & nSecondHits > 0,
                      (nSecondHits * nhits1/2) + (nhits1/2),
                      ifelse(nhits1 > 1 & nSecondHits > 0,
                             nSecondHits * nhits1,
                             nSecondHits)))]

  ##############################################################################
  # 3. Add all the parameters into the output data.table.
  ##############################################################################
  params <- c("onlyOgAnchors","onlyOgAnchorsSecond","blkSize", "blkSizeSecond",
              "nGaps", "nGapsSecond",
              "synBuff", "synBuffSecond", "selfRegionMask",
              "nSecondHits", "dropInterleavesSmallerThan", "minRbhScore")
  for(i in params)
    cmb[[i]] <- get(i)

  cmb <- cmb[with(cmb, order(-(nGenes1 + nGenes2))),]
  gsParam$params$synteny <- cmb
  return(gsParam)
}


#' @title Run MCScanX from R
#' @description
#' \code{run_mcscanx} Internal GENESPACE function to run MCScanX on hits when
#' searching for synteny.
#' @rdname synteny
#' @import R.utils
#' @import data.table
#' @export
run_mcscanx <- function(hits,
                        blkSize,
                        nGaps,
                        gsParam,
                        path2mcscanx){
  ##############################################################################
  # parameter argument checking
  gn2 <- mcsID1 <- mcsID2 <- blkID <- ord2 <- ofID1 <- ofID2 <- ord1 <- NULL
  runBlast <- chr1 <- chr2 <- blkID <- NULL
  hitCols <- c("ofID1","ofID2","chr1","chr2","ord1","ord2","score")
  if(!all(hitCols %in% colnames(hits)))
    stop("Looks like the hits data.table is misformed. It must be a data.table",
         "comparing two genomes with the columns:", hitCols,"\n")
  if(length(blkSize) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(length(nGaps) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(!is.integer(blkSize)){
    blkSize <- as.integer(blkSize)
    if(is.na(blkSize))
      stop("Cannot coerce blkSize to an integer\n")
  }
  if(!is.integer(nGaps)){
    nGaps <- as.integer(nGaps)
    if(is.na(nGaps))
      stop("Cannot coerce nGaps to an integer\n")
  }
  MCScanX_h <- file.path(path2mcscanx, "MCScanX_h")
  if(!file.exists(MCScanX_h))
    stop("Cannot find MCScanX_h executable in", path2mcscanx,"\n")

  ##############################################################################
  # set tmp directory
  tmpDir <- file.path(
    gsParam$paths$results,
    paste0("tmp_",
           paste(sample(c(letters,LETTERS), 20, replace = T), collapse = "")))
  if (dir.exists(tmpDir))
    unlink(tmpDir, recursive = T)
  dir.create(tmpDir)
  on.exit(expr = unlink(tmpDir, recursive = T))

  ##############################################################################
  # convert gene locations to 'gff' like mcscanx file
  u1 <- subset(hits, !duplicated(ofID1))
  u2 <- subset(hits, !duplicated(ofID2))
  setkey(u1, chr1, ord1)
  setkey(u2, chr2, ord2)
  u1[,mcsID1 := paste0("aa", as.numeric(factor(chr1, levels = unique(chr1))))]
  u2[,mcsID2 := paste0("bb", as.numeric(factor(chr2, levels = unique(chr2))))]

  mcs1 <- u1[,c("mcsID1", "ofID1","ord1","ord1")]
  setnames(mcs1, c("chr", "id","start","end"))

  mcs2 <- u2[,c("mcsID2", "ofID2","ord2","ord2")]
  setnames(mcs2, c("chr", "id","start","end"))
  mcsGffIn <- rbind(mcs1, mcs2)

  ##############################################################################
  # convert hits to mcscanx_h blast input
  mcsBlsIn <- hits[,c("ofID1", "ofID2", "score")]
  # mcsBlsIn[,score := 1]

  mcsBlsIn[,ofID2 := paste0(ofID2,"xxxx")]
  mcsGffIn$id[grepl("bb", mcsGffIn$chr)] <- paste0(mcsGffIn$id[grepl("bb", mcsGffIn$chr)],"xxxx")

  blFile <- file.path(tmpDir, "mcs.homology")
  gfFile <- file.path(tmpDir, "mcs.gff")
  colFile <- file.path(tmpDir, "mcs.collinearity")

  fwrite(
    mcsGffIn,
    file = gfFile,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    showProgress = FALSE,
    verbose = FALSE)
  fwrite(
    mcsBlsIn,
    file = blFile,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    showProgress = FALSE,
    verbose = FALSE)

  ##############################################################################
  # run mcscanx_h
  mcsCom <- sprintf(
    "%s -a -b 2 -c 2 -m %s -s %s %s 1>/dev/null 2>&1",
    MCScanX_h, nGaps, blkSize, file.path(tmpDir,"mcs"))
  system(mcsCom)
  ##############################################################################
  idg <- strsplit(as.character(hits$ofID1[1]), "_")[[1]][1]
  # parse collinearity file
  suppressWarnings(collin <-  fread(
    cmd = sprintf("cat %s | grep %s_ | grep :", colFile, idg),
    col.names = c("blkID","gn1","gn2"),
    select = 1:3,
    header = F))
  if(nrow(collin) > 1){
    collin[,blkID := as.numeric(sapply(blkID, function(x)
      strsplit(x, "-")[[1]][1])) + 1]
    collin[,gn2 := gsub("xxxx$", "", gn2)]
    mcsb <- collin$blkID
    names(mcsb) <- with(collin, paste(gn1, gn2))
    return(mcsb)
  }
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
  return(gff)
}

