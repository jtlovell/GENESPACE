#' @title Build genespace pangenome
#'
#' @description
#' \code{pangenome} Convert orthogroup and synteny information into a
#' pangenome database. Predict locations of orthogroups that are missing a
#' node in the reference.
#'
#' @name pangenome
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param refGenome character string matching one of the genomeIDs in gsParam
#' @param genomeIDs character vector, specifying which genomes to use. Defaults
#' to all genomeIDs specification in gsParam.
#' @param allowRBHinOg logical, should RBHs be allowed in the construction of
#' syntenic orthologs
#' @param hitsRef hits data.table against the reference genome
#' @param blks data.table of block coordinates
#' @param gff gff with additional annotations columns (see synteny)
#' @param gffRaw raw gff from read_gff
#' @param infPos data.table with inferred gene positions
#' @param allowRBHinOg logical, should synOGs be constructed aware of RBHs?
#' @param hits data.table of hits
#' @param maxGap integer specifying the maxiumum allowable gap between synOgs
#' @param infGenePos data.table with inferred gene positions
#' @param nCores integer specifying the number of parallel processes to run
#'
#' @param verbose logical, specifying whether to print updates to the console.
#' Defaults to verbose specification in gsParam, unless specified here
#'
#' @details ...
#' \enumerate{
#' \item genome1/2: the genome IDs for each pairwise run
#' \item uGenes1/2: the number of unique genes for each genome
#' \item query/target: the identity of query/target genomes in each blast run
#' \item run/mirrorBlast: logical whether the blast was run or mirrored
#' \item gn1/2: orthofinder genomeID numbers for each genome
#' \item db/fa/1/2/blFile: diamond database, fasta files and blast files.
#' }
#'
#' @return ...
#'
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @note \code{pangenome} is a generic name for the functions documented.
#' \cr
#' If called, \code{pangenome} returns its own arguments.
#'
#' @title build pangenome database
#' @description
#' \code{pangenome} Predict locations in reference genome for all
#' genes and collapse into pangenome database
#' @rdname pangenome
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom dbscan dbscan frNN
#' @importFrom stats complete.cases median
#' @export
pangenome <- function(gsParam,
                      refGenome = NULL,
                      genomeIDs = NULL,
                      verbose = NULL,
                      allowRBHinOg = NULL){

  genome<- isArrayRep  <- ord  <- ofID  <- ofID1  <- ofID2  <- regAnchor <- NULL
  blkAnchor  <- regID  <- gen1  <- gen2  <- blkID  <- synArray  <- repOrd <- NULL
  synOg  <- repOfID  <- pgID  <- refChr  <- clus  <- id1  <- id2 <- NULL
  isOrtholog  <- isSyntenic <- NULL

  # - get the distance among genes to merge ref ogs
  if(is.null(gsParam$params$synteny))
    stop("must run set_syntenyParams first\n")
  di <- max(gsParam$params$synteny$synBuff)*2
  if(is.null(verbose))
    verbose <- gsParam$params$verbose
  if(length(verbose) > 1 || !is.logical(verbose))
    verbose <- TRUE
  nCores <- gsParam$params$nCores

  # -- specify all the genomes
  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
    genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]
  }else{
    tmp <- gsParam$genomes$genomeIDs
    tmp <- tmp[!tmp %in% gsParam$genomes$outgroup]
    genomeIDs <- genomeIDs[genomeIDs %in% tmp]
  }

  # -- specify the ref genome
  if(is.null(refGenome))
    refGenome <- genomeIDs[1]
  if(!refGenome %in% genomeIDs)
    stop(sprintf("refgenome %s is not in genome IDs", refGenome))
  if(verbose)
    cat(sprintf("Initiating pan-genome, using %s as reference\n", refGenome))

  # -- if any polyploids, allow RBH in net
  if(is.null(allowRBHinOg) || !is.logical(allowRBHinOg) || length(allowRBHinOg) != 1)
    allowRBHinOg <- TRUE
  if(all(gsParam$genomes$ploidy == 1))
    allowRBHinOg <- FALSE

  # -- specify the output files
  blksFile <- file.path(gsParam$paths$results, "syntenicBlocks.txt.gz")
  pgFile <- file.path(gsParam$paths$results, sprintf("%s_pangenomeDB.txt.gz", refGenome))

  #############################################################################
  # 1. Parse/read gff, ref hits and blks
  #############################################################################
  # -- read gff
  if(verbose)
    cat(sprintf("Preparing input data ... \n\tReading gff and adding syn OGs ... "))
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gffRaw <- subset(fread(gffFile), genome %in% genomeIDs)
  gff <- subset(gffRaw, isArrayRep)
  ov <- gffRaw$ord; sv <- gffRaw$start; ev <- gffRaw$end
  names(ov) <- names(sv) <- names(ev) <- gffRaw$ofID
  gff[,ord := ov[ofID]]

  # -- read all syn OGs hits
  if(verbose)
    cat(sprintf(
      "found %s genes and %s syn OGs\n\tReading syntenic hits ... ",
      nrow(gff), uniqueN(gff$synOg)))
  hts <- read_refSynOgHits(
    gsParam = gsParam,
    refGenome = refGenome,
    genomeIDs = genomeIDs,
    allowRBHinOg = allowRBHinOg)

  # -- count the number of anchor hits
  if(verbose)
    cat(sprintf(
      "found %s OG hits in %s syn regions and %s blks\n\tCounting syntenic anchors by gene ... ",
      nrow(hts), uniqueN(hts$regID), uniqueN(hts$blkID)))
  gff <- add_nAnchor2gff(gsParam = gsParam, gff = gff, hits = hts)

  # -- read in block coordinates
  if(verbose)
    cat(sprintf(
      "Done!\n\tCalculating block and region coordinates ... ",
      nrow(gff), uniqueN(gff$synOg)))
  hts[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2],
            start1 = sv[ofID1], start2 = sv[ofID2],
            end1 = ev[ofID1], end2 = ev[ofID2])]
  regHit <- subset(hts, regAnchor)
  blkHit <- subset(hts, blkAnchor)
  regHit[,`:=`(blkID = regID, regID = NULL, blkAnchor = regAnchor, regAnchor = NULL)]
  blkHit[,`:=`(regID = NULL, regAnchor = NULL)]
  blkCoord <- calc_blkCoords(blkHit)
  regCoord <- calc_blkCoords(regHit)

  # -- read in syn orthogroup hits against reference
  if(verbose)
    cat("Done!\n\tPulling hits against reference ... ")
  cn1 <-c("gen1","gen2","chr1","chr2","ofID1","ofID2","ord1","ord2","og","blkID","regID","blkAnchor","regAnchor")
  cn2 <- c("gen2","gen1","chr2","chr1","ofID2","ofID1","ord2","ord1","og","blkID","regID","blkAnchor","regAnchor")
  hts <- hts[,cn1, with = F]
  hitsRef <- subset(hts, gen1 == refGenome)
  tmp <- subset(hts, gen2 == refGenome & gen1 != refGenome)[,cn2, with = F]
  setnames(tmp, cn2)
  hitsRef <- subset(rbind(hitsRef, tmp), !is.na(regID) | !is.na(blkID))

  if(verbose)
    cat(sprintf("found %s syn anchors\n", nrow(hitsRef)))

  #############################################################################
  # 2. Infer gene positions against reference
  #############################################################################

  # -- initial inference
  if(verbose)
    cat(sprintf("Infering positions against %s ...\n\t(genome: missing / 1x / 2x / >2x placements)",
                refGenome))
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)

  infReg <- lift_genePos(
    hits = blkHit,
    refGenome = refGenome,
    blks = blkCoord,
    gff = gffRaw,
    verbose = verbose,
    nCores = nCores,
    maxGap = 10,
    genomeIDs = genomeIDs)

  # -- pull gapped regions
  infPos <- fill_missingGenePos(
    hits = blkHit,
    infGenePos = infReg,
    gff = gffRaw,
    genomeIDs = genomeIDs,
    maxGap = 10,
    verbose = verbose)

  #############################################################################
  # 3. Build pan-genome
  #############################################################################

  # -- build framework with the reference
  if(verbose)
    cat(sprintf("Building pan-genome\n\tScaffolding array reps ... "))
  refPg <- build_repPangenome(
    gsParam = gsParam,
    gff = gff,
    infPos = infPos,
    refGenome = refGenome)

  # -- add in array members
  if(verbose)
    with(refPg, cat(sprintf(
      "%s genes assigned to %s positions\n\tAdding array members ... ",
      uniqueN(ofID), uniqueN(paste(repOfID, refChr, clus)))))
  gffArray <- subset(gffRaw, !ofID %in% refPg$ofID)
  gffRep <- subset(gff, isArrayRep & synArray %in% gffArray$synArray)
  gffmer <- merge(
    gffRep[,c("ofID", "synArray")],
    data.table(arrayMem = gffArray$ofID, synArray = gffArray$synArray),
    by = "synArray", allow.cartesian = T)
  pgmer <- merge(refPg, gffmer, by = "ofID", allow.cartesian = T, all = T)
  if(verbose)
    with(gffmer, cat(sprintf(
      "found %s arrays with %s genes\n",
      uniqueN(synArray), uniqueN(c(ofID, arrayMem)))))

  # -- reformat and finalize
  setorder(pgmer, repOrd, synOg, repOfID, genome, na.last = T)
  pgmer[,pgID := as.integer(factor(
    paste(synOg, refChr, clus),
    levels = unique(paste(synOg, refChr, clus))))]
  pgmer[,`:=`(clus = NULL, synArray = NULL)]
  setcolorder(pgmer, c("refGenome", "pgID", "synOg", "refChr", "repOrd","repOfID","genome","ofID", "arrayMem"))

  #############################################################################
  # 4. annotate pangenome
  #############################################################################
  if(verbose)
    cat("Annotating pangenome\n\tFlagging syntenic orthologs ... ")
  if(is.na(gsParam$paths$orthologuesDir)){
    cat("Cannot find orthologs. Perhaps orthofinder -og was run?\n\tNot adding non-syntenic orthologs\n\t")
    pgo <- data.table(pgmer)
  }else{
    # -- pull genes that were in orthogroups without projected positions against ref
    olDt <- parse_orthologues(
      gsParam = gsParam,
      refGenome = refGenome,
      nCores = gsParam$params$nCores)
    ov <- gff$ofID; names(ov) <- with(gff, paste(genome, id))
    olDt[,`:=`(ofID1 = ov[paste(gen1, id1)], ofID2 = ov[paste(gen2, id2)])]
    ug <- with(olDt, c(paste(ofID1, ofID2), paste(ofID2, ofID1)))
    pgmer[,isOrtholog := paste(repOfID, ofID) %in% ug]

    if(verbose)
      cat(sprintf("found %s orthologs\n\tChecking for non-syntenic orthologs ... ",
                  sum(pgmer$isOrtholog)))
    pgout2 <- data.table(pgmer)
    pgout2[,isSyntenic := TRUE]

    #   # -- Pull non-syntenic orthlogs
    u <- unique(with(pgout2, paste(repOfID, ofID)))
    olDtMiss <- subset(olDt, !paste(ofID1, ofID2) %in% u)
    pgout3 <- subset(pgout2, repOfID %in% olDtMiss$ofID1)
    pgout3 <- subset(pgout3, !duplicated(pgID))
    pgout3 <- pgout3[,c("synOg", "refChr", "repOrd", "repOfID","pgID")]
    pgout3 <- merge(
      pgout3,
      with(olDtMiss, data.table(
        repOfID = ofID1, ofID = ofID2, genome = gen2, refGenome = refGenome)),
      by = "repOfID", allow.cartesian = T)
    pgout3[,`:=`(isOrtholog = TRUE, isSyntenic = FALSE)]
    if(verbose)
      cat(sprintf("added %s entries\n\tWriting pangenome to file ... ", nrow(pgout3)))

    pgo <- rbind(pgout2, pgout3, fill = T)

  }

  fwrite(pgo, sep = "\t", quote = F, showProgress = F, file = pgFile)
  if(verbose)
    cat(sprintf("Done!\n\tWritten to %s\n",pgFile))

  return(pgo)
}

#' @title infer_genePos
#' @description
#' \code{infer_genePos} infer_genePos
#' @rdname pangenome
#' @import data.table
#' @export
infer_genePos <- function(hitsRef,
                          blks,
                          gff,
                          gsParam,
                          genomeIDs,
                          refGenome){
  infer_pos <- function(ord1, # named order vector from gff
                        ord2,
                        maxGap){ # all alt ofIDs within blk coordinates
    o2 <- l <- r <- o1 <- runID <- NULL
    o12 <- subset(
      data.table(o1 = ord1, o2 = ord2, index = 1:length(ord1), key = "o2"),
      !is.na(o2))
    if(sum(stats::complete.cases(o12)) < 2){
      return(NA)
    }else{
      whNaSt <- which(is.na(o12$o1) & o12$o2 == min(o12$o2, na.rm = T))
      whNaEn <- which(is.na(o12$o1) & o12$o2 == max(o12$o2, na.rm = T))
      o12$o2[whNaSt] <-  o12$o2[whNaSt] + .001
      o12$o2[whNaEn] <-  o12$o2[whNaEn] - .001
      setkey(o12, o2)
      o12[,runID := add_rle(is.na(o1), which = "id")]
      if(is.na(o12$o1[1]))
        o12 <- subset(o12, runID > 1)
      if(is.na(o12$o1[nrow(o12)]))
        o12 <- subset(o12, runID != o12$runID[nrow(o12)])

      oi <- subset(o12, !is.na(o1))
      om <- subset(o12, is.na(o1))
      si <- oi[,list(r = o1[1], l = o1[length(o1)]), by = "runID"]
      rv <- si$r; names(rv) <- as.character(si$runID-1)
      lv <- si$l; names(lv) <- as.character(si$runID+1)
      om[,l := lv[as.character(runID)]]
      om[,r := rv[as.character(runID)]]
      om[,o1 := as.numeric(o1)]
      om[,o1 := as.numeric(seq(
        from = l[1], to = r[1], length.out = .N+2)[2:(.N+1)]),
        by = "runID"]
      om$o1[with(om, r - l) > maxGap] <- NA
      ord1[om$index] <- om$o1
      return(ord1)
    }
  }
  isArrayRep <- genome <- ord <- blkID <- gen1 <- genome <- ofID2 <- ord2 <- NULL
  ord1i <- ord1 <- ofID <- n <- refOrd <- rng <- refChr <- clus <- ord <- NULL
  maxGap <- min(gsParam$params$synteny$blkSize)*2
  verbose <- gsParam$params$verbose
  if(is.null(verbose))
    verbose <- gsParam$params$verbose
  if(length(verbose) > 1 || !is.logical(verbose))
    verbose <- TRUE
  nCores <- gsParam$params$nCores

  g <- subset(gff, isArrayRep)
  setkey(g, genome, ord)
  spl <- split(hitsRef, by = "blkID")
  wh2keep <- sapply(spl, function(x) sum(x$blkAnchor) > 1)
  spli <- spl[wh2keep]
  spln <- spl[!wh2keep]
  blki <- subset(blks, blkID %in% names(spli) & gen1 == refGenome)
  spli <- spli[blki$blkID]

  # -- parse the gff so to easily get position information for ref and alts
  gf1 <- with(subset(g, genome == refGenome), data.table(
    ofID1 = ofID, ord1 = ord))
  gf2 <- sapply(genomeIDs, simplify = F, USE.NAMES = T, function(x)
    with(subset(g, genome == x), data.table(
      ofID2 = ofID, ord2 = ord)))
  splbi <- split(blki, by = "gen2")
  d2h <- gsParam$params$maxDistBtwPgHits
  infPos <- rbindlist(lapply(names(splbi), function(i){
    # -- prep genome i
    if(verbose)
      cat(sprintf("\n\t%s: ", i))
    gfi2 <- gf2[[i]]
    bi <- splbi[[i]]
    ip <- rbindlist(mclapply(1:nrow(bi), mc.cores = 1, function(j){
      # -- pull genes and anchor in blk i
      b <- bi[j,]
      g1 <- gf1[with(gf1, which(ofID1 == b$firstGene1):which(ofID1 == b$lastGene1)),]
      g2 <- gfi2[with(gfi2[[b$gen2]], which(ofID2 == b$firstGene2):which(ofID2 == b$lastGene2)),]
      h <- data.table(spli[[b$blkID]][,c("ofID1","ofID2","og")])
      if("2_2446" %in% h$ofID2)
        print(j)
      # -- merge genes by anchor og
      g1 <- merge(g1, h[,c("ofID1","og")], by = "ofID1", all.x = T)
      g2 <- merge(g2, h[,c("ofID2","og")], by = "ofID2", all.x = T)
      g12 <- merge(subset(g1, complete.cases(g1)),
                   subset(g2, complete.cases(g2)),
                   by = "og", allow.cartesian = T, all = T)
      g12 <- rbind(g12,
                   subset(g2, !ofID2 %in% g12$ofID2),
                   fill = T)
      g12 <- subset(g12, !duplicated(ofID2))

      # -- prep input and infer the positions
      setkey(g12, ord2)
      x1 <- g12$ord1; names(x1) <- g12$ofID1
      x2 <- g12$ord2; names(x2) <- g12$ofID2

      if(any(is.na(g12$ord1))){
        g12[,ord1i := infer_pos(ord1 = x1, ord2 = x2, maxGap = maxGap)]
      }else{
        g12[,ord1i := ord1]
      }

      o12 <- with(g12, data.table(
        ofID = ofID2, genome = b$gen2, refGenome = b$gen1,
        chr = b$chr2, refChr = b$chr1, ord = ord2, refOrd = ord1i))
      return(o12)
    }))

    # -- add genes without a hit
    ip0 <- subset(gff, genome == i & !ofID %in% unique(ip$ofID))
    if(nrow(ip0) > 0){
      ip0 <- with(ip0, data.table(
        ofID = ofID, genome = genome, refGenome = refGenome, chr = chr,
        refChr = NA, ord = ord, refOrd = NA))
      ip <- rbind(ip, ip0)
    }

    # -- separate genes with a single location / chr from those with multiple
    ip[,n := uniqueN(refOrd[!is.na(refOrd)]), by = c("refChr", "ofID")]
    ip2 <- subset(ip, n >= 2)
    ip <- subset(ip, n < 2)
    ip[,n := NULL]

    # -- pull out genes that with a range of locations < min dist
    if(nrow(ip2) > 1){
      ip2[,rng := diff(range(refOrd)),  by = c("refChr", "ofID")]
      ip2a <- subset(ip2, rng < d2h)
      ip2 <- subset(ip2, rng >= d2h)
      ip2a[,refOrd := median(refOrd),  by = c("refChr", "ofID")]
      ip2a <- subset(ip2a, !duplicated(paste(refChr, ofID)))[,colnames(ip), with = F]
      ip <- rbind(ip, ip2a)

      # -- check and split genes that hit multiple places on a single chr
      ip2[,clus := dbscan(frNN(
        cbind(refOrd, refOrd),
        eps = d2h),
        minPts = 1)$cluster,
        by = c("refChr", "ofID")]
      ip2[,refOrd := median(refOrd),  by = c("refChr", "ofID","clus")]
      ip2b <- subset(ip2, !duplicated(paste(refChr, ofID)))[,colnames(ip), with = F]
      ip <- rbind(ip, ip2b)
    }

    setkey(ip, ord)
    ip <- subset(ip, !duplicated(paste(ofID, refChr, refOrd)))
    ipn <-  ip[,list(n = sum(!is.na(refOrd))), by = "ofID"]
    if(verbose)
      with(ipn, cat(sprintf("%s / %s / %s / %s", sum(n == 0), sum(n == 1), sum(n == 2), sum(n > 2))))
    return(ip)
  }))

  return(infPos)
}

#' @title build_repPangenome
#' @description
#' \code{build_repPangenome} build_repPangenome
#' @rdname pangenome
#' @import data.table
#' @export
build_repPangenome <- function(gsParam, gff, gffRaw, infPos, refGenome){

  isArrayRep <- genome <- synOg <- posType <- ofID <- refChr <- chr <- NULL
  refOrd <- ord <- clus <- n <- maxJump <- anyAnnot <- nanch <- NULL
  nog <- pepLen <- repOfID <- repOrd <- ofIDs <- medOrd <- NULL
  dist2median <- nBlkAnchors <- nRegAnchors <- nGenomeOrthologs  <- NULL
  nTotalOrthologs <- NULL

  d2h <- gsParam$params$maxDistBtwPgHits

  # -- build the scaffold against genes in the reference

  nogv <- with(gff, nGenomeOrthologs + nTotalOrthologs)
  sogv <- gff$synOg; nahv <- gff$nBlkAnchors; narv <- gff$nRegAnchors; pepv <- gff$pepLen
  names(pepv) <- names(nahv) <- names(nogv) <- names(sogv) <- names(narv) <- gff$ofID

  ##############################################################################
  # 1. Build the scaffold with representatives from the refgenome
  ##############################################################################
  g <- subset(gff, isArrayRep)
  refSynOgs <- subset(g, genome == refGenome)$synOg
  gref <- subset(g, synOg %in% refSynOgs)
  pgRef <- subset(gref, genome == refGenome)[,c("chr","ord","ofID","synOg")]
  pgRef[,posType := "annotation"]

  # -- add in alternate positions for reference scaffold syn ogs
  posRef <- subset(infPos, ofID %in% pgRef$ofID & !is.na(refChr))
  posRef <- subset(posRef, chr != refChr | abs(refOrd - ord) > d2h)

  # -- combine annotation and inferred positions
  g1 <- rbind(
    pgRef,
    with(posRef, data.table(
      chr = refChr, ord = refOrd, ofID = ofID,
      synOg = sogv[ofID], posType = "inferred")))
  g1[,clus := 1]
  g1[,n := uniqueN(paste(ord, chr)), by = c("chr","synOg")]

  # -- split out those with just one position and multi-mappers
  g2 <- subset(g1, n > 1)
  g1 <- subset(g1, n == 1)
  g1[,n := NULL]

  # -- split out those that may need to be clustered
  setkey(g2, synOg, ord)
  g2[,maxJump := max(diff(ord)), by = c("chr", "synOg")]
  g2r <- subset(g2, maxJump > d2h * 2)
  g2 <- subset(g2, maxJump <= d2h * 2)

  # -- cluster and add back in those with broad hits
  if(nrow(g2r) > 1){
    g2r[,clus := dbscan(frNN(cbind(ord, ord), eps = d2h*2), minPts = 0)$cluster,
        by = c("synOg", "chr")]
    g2 <- rbind(g2, g2r)
  }

  # -- if any of the hits are there because of annotation position, override pos
  g2[,anyAnnot := any(posType == "annotation"), by = c("chr", "synOg", "clus")]
  g2a <- subset(g2, anyAnnot)
  g2a[,ord := median(ord[posType == "annotation"]), by = c("chr", "synOg", "clus")]
  g2a[,`:=`(maxJump = NULL, anyAnnot = NULL, n = NULL)]

  # -- for those without annot positions, add weights and choose the best one
  g2n <- subset(g2, !anyAnnot)
  g2n[,`:=`(pepLen = pepv[ofID], nanch = nahv[ofID], nog = nogv[ofID])]
  setorder(g2n, chr, synOg, -nanch, -nog, -pepLen)
  g2n[, ord := ord[1], by = c("chr", "synOg", "clus")]
  g2 <- rbind(g2n[,colnames(g2a), with = F], g2a)

  # -- build the reference scaffold
  pg <- rbind(g1, g2)
  setkey(pg, synOg, posType, ord)
  setnames(pg, 1:2, c("refChr", "refOrd"))

  # -- assign the representative gene models for each position
  pg[,`:=`(repOfID = ofID[1], repOrd = refOrd[1]),
     by = c("refChr", "synOg", "clus")]
  pgr <- subset(pg, ofID == repOfID)[,c("synOg","refChr","refOrd","clus","repOfID")]
  pgo <- pg[,list(ofIDs = list(unique(ofID[ofID != repOfID]))),
            by = c("refChr", "repOrd", "synOg", "clus","repOfID")]
  setkey(pgo, synOg, repOrd)
  pgo[,refGenome := refGenome]

  ##############################################################################
  # 2. add in genes from other genomes
  ##############################################################################
  pgi <- merge(
    pgo[,c("refGenome","refChr", "repOrd", "synOg", "clus", "repOfID")],
    subset(g, genome != refGenome)[,c("synOg", "ofID", "genome")],
    by = "synOg", allow.cartesian = T)

  # -- reformat reference pangenome
  pgr <- pgo[,list(ofID = c(repOfID, unlist(ofIDs)),
                   genome = refGenome),
             by = c("refGenome","refChr", "repOrd", "synOg", "clus", "repOfID")]
  setkey(pgr, synOg, repOrd)

  # -- join together
  pgOut <- rbind(pgi, pgr)
  setkey(pgOut, synOg, repOrd)

  ##############################################################################
  # 3. construct pangenome for genes that don't have refgenome in synog
  ##############################################################################

  # subset the gff and inferred positions to missing genes
  pgAlt <- subset(g, genome != refGenome)[,c("ofID","synOg")]
  pgAlt <- subset(pgAlt, !synOg %in% unique(pgOut$synOg))
  pgAlt[,posType := "annotation"]

  # -- add in alternate positions for reference scaffold syn ogs
  posAlt <- subset(infPos, ofID %in% pgAlt$ofID & !is.na(refChr) & !is.na(refOrd))

  # -- combine annotation and inferred positions
  g1 <- merge(pgAlt, posAlt, by = "ofID", allow.cartesian = T)
  g1[,clus := 1]
  g1[,n := uniqueN(refOrd), by = c("refChr","synOg")]

  # -- split out those with just one position and multi-mappers
  g2 <- subset(g1, n > 1)
  g1 <- subset(g1, n == 1)
  g1[,n := NULL]

  # -- split out those that may need to be clustered
  setkey(g2, synOg, refOrd)
  g2[,maxJump := max(diff(refOrd)), by = c("refChr", "synOg")]
  g2r <- subset(g2, maxJump > d2h)
  g2 <- subset(g2, maxJump <= d2h)

  # -- cluster and add back in those with broad hits
  if(nrow(g2r) > 1){
    g2r[,clus := dbscan(frNN(cbind(refOrd, refOrd), eps = d2h), minPts = 0)$cluster,
        by = c("synOg", "refChr")]
    g2 <- rbind(g2, g2r)
  }
  g2[,`:=`(n = NULL, maxJump = NULL)]
  g3 <- rbind(g1, g2)
  g3[,medOrd := as.numeric(median(refOrd, na.rm = T)), by = c("refChr","clus","synOg")]
  g3[,dist2median := abs(medOrd - refOrd), by = c("refChr","clus","synOg")]
  g3[,`:=`(pepLen = pepv[ofID], nanch = nahv[ofID], nog = nogv[ofID])]
  setorder(g3, refChr, synOg, -nanch, dist2median, -nog, -pepLen)
  g3[, `:=`(repOfID = ofID[1],
            repOrd = refOrd[1],
            refGenome = refGenome),
     by = c("refChr", "synOg", "clus")]
  g3o <- g3[,colnames(pgOut), with = F]
  pgOut <- rbind(pgOut, g3o)

  ##############################################################################
  # 4. add in unanchored synOgs
  ##############################################################################
  galt <- subset(g, !synOg %in% unique(pgOut$synOg))
  setorder(galt, synOg, -nBlkAnchors, -nRegAnchors, -nGenomeOrthologs, -nTotalOrthologs, -pepLen)
  galt[, `:=`(repOfID = ofID[1]), by = c("synOg")]
  galto <- with(galt, data.table(
    synOg = synOg, refGenome = refGenome, refChr = NA, repOrd = NA, clus = 1,
    repOfID = repOfID, ofID = ofID, genome = genome))
  pgOut <- rbind(pgOut, galto)
  setorder(pgOut, repOrd, synOg, genome, na.last = T)
  return(pgOut)
}

#' @title read_refSynOgHits
#' @description
#' \code{read_refSynOgHits} read_refSynOgHits
#' @rdname pangenome
#' @import data.table
#' @export
read_refSynOgHits <- function(gsParam, refGenome, genomeIDs, allowRBHinOg){
  gen1 <- gen2 <- regID <- og <- isRBH <- chr1 <- chr2 <- blkID <- regID <- NULL
  eg <- subset(
    CJ(gen1 = genomeIDs, gen2 = genomeIDs),
    gen1 == refGenome | gen2 == refGenome)
  fs <- file.path(gsParam$paths$results,
                  sprintf("%s_%s_synHits.txt.gz",eg$gen1, eg$gen2))
  fs <- fs[file.exists(fs)]

  hts <- rbindlist(mclapply(fs, mc.cores = gsParam$params$nCores, function(i){
    x <- fread(
      i, select = c("gen1", "gen2", "chr1", "chr2", "ofID1", "ofID2",
                    "og", "regID", "blkID", "blkAnchor", "regAnchor"),
      na.strings = c("","NA"))
    x <- subset(x, !is.na(regID) | !is.na(og))
    x[,`:=`(isRBH = grepl("^RBH", og))]
    if(!allowRBHinOg)
      x <- subset(x, !isRBH)
    return(x[,c("gen1", "gen2", "chr1", "chr2", "og",
                "ofID1", "ofID2", "blkID", "regID", "blkAnchor", "regAnchor")])
  }))
  whnb <- which(is.na(hts$blkID))
  whnr <- which(is.na(hts$regID))
  hts[,`:=`(
    blkID = sprintf(
      "%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, as.numeric(as.factor(blkID))),
    regID = sprintf(
      "%s_%s_%s_%s_%s", gen1, gen2, chr1, chr2, as.numeric(as.factor(regID))))]
  hts$blkID[whnb] <- NA
  hts$regID[whnr] <- NA
  return(hts)
}

#' @title add_nAnchor2gff
#' @description
#' \code{add_nAnchor2gff} add_nAnchor2gff
#' @rdname pangenome
#' @import data.table
#' @export
add_nAnchor2gff <- function(gsParam, gff, hits){
  regAnchor <- blkAnchor <- blkID <- regID <- ofID <- NULL
  tmp <- with(subset(hits, regAnchor | blkAnchor), data.table(
    ofID = c(ofID1, ofID2),
    regAnchor = c(regAnchor, regAnchor),
    blkAnchor = c(blkAnchor, blkAnchor),
    regID = c(regID, regID),
    blkID = c(blkID, blkID)))
  htcnt <- tmp[,list(nBlkAnch = uniqueN(blkID[blkAnchor]),
                     nRegAnch = uniqueN(regID[regAnchor])),
               by = "ofID"]
  ploidy <- sum(gsParam$genomes$ploidy)*2
  htcnt$nBlkAnch[htcnt$nBlkAnch >= ploidy] <- ploidy-1
  htcnt$nRegAnch[htcnt$nRegAnch >= ploidy] <- ploidy-1

  a <- htcnt$nBlkAnch; r <- htcnt$nRegAnch; names(a) <- names(r) <- htcnt$ofID
  gff[,`:=`(nBlkAnchors = a[ofID], nRegAnchors = r[ofID])]
  gff$nBlkAnchors[is.na(gff$nBlkAnchors)] <- 0
  gff$nRegAnchors[is.na(gff$nRegAnchors)] <- 0
  return(gff)
}

#' @title lift_genePos
#' @description
#' \code{lift_genePos} lift_genePos
#' @rdname pangenome
#' @import data.table
#' @export
lift_genePos <- function(gff,
                         hits,
                         blks,
                         maxGap,
                         verbose,
                         nCores,
                         refGenome,
                         genomeIDs){
  refOrd <- ofID <- ord <- ord2 <- genome <- d <- ord1 <- rl <- ord2 <- NULL
  ofID2 <- blkID <- gen1 <- gen2 <- ord <- genome <- NULL
  setkey(gff, genome, ord)
  splgff <- lapply(split(gff, by = c("genome")), function(x) split(x, by = "chr"))
  infout <- rbindlist(lapply(genomeIDs, function(thisGen){
    genBlks <- subset(blks, gen2 == thisGen & gen1 == refGenome)
    genHits <- subset(hits, gen2 == thisGen & gen1 == refGenome)
    genBlks <- subset(genBlks, blkID %in% genHits$blkID)
    genHits <- split(genHits, by = "blkID")[genBlks$blkID]
    if(verbose)
      cat(sprintf("\n\t\t%s: ", thisGen))
    outgen <- rbindlist(mclapply(names(genHits), mc.cores = 1, function(thisBlk){
      x <- genHits[[thisBlk]]
      y <- subset(genBlks, blkID == thisBlk)
      z2 <- data.table(splgff[[y$gen2]][[y$chr2]])
      z2 <- z2[which(z2$ofID == y$firstGene2):which(z2$ofID == y$lastGene2),c("ofID","ord")]
      setnames(z2, c("ofID2","ord2"))
      z2 <- subset(z2, !ofID2 %in% x$ofID2)
      z <- rbind(x[,c("ofID1","ord1","ofID2","ord2")], z2, fill = T)
      setkey(z, ord2)
      if(all(!is.na(z$ord1))){
        return(z)
      }else{
        z[,rl := add_rle(is.na(ord1), which = "id")]
        z <- subset(z, !(is.na(ord1) & rl %in% c(1, max(rl))))

        # -- get closest two obs from each run
        uv <- subset(z, !is.na(ord1))
        setkey(uv, ord1)
        uv <- with(uv, split(ord1, rl))
        iv <- lapply(2:length(uv), function(i){
          il <- uv[[i-1]]
          ir <- uv[[i]]
          eg <- CJ(il, ir)[,d := abs(ir - il)]
          eg <- eg[which.min(d),]
          return(c(eg[[1]], eg[[2]]))
        })

        ov <- split(subset(z, is.na(ord1)), by = "rl")
        names(iv) <- names(ov)
        for(i in names(ov)){
          l <- iv[[i]][1]
          r <- iv[[i]][2]
          if(abs(l - r) <= maxGap){
            ov[[i]]$ord1 <- seq(from = l, to = r, length.out = nrow(ov[[i]]))
          }
        }

        out <- rbind(subset(z, !is.na(ord1)), rbindlist(ov))
        out[,`:=`(rl = NULL, refChr = y$chr1, blkID = thisBlk)]
        return(out)
      }
    }), fill = T)
    setkey(outgen, ord1)
    outgen <- subset(outgen, !is.na(ord1))
    tab <- table(table(outgen$ofID2))
    gout <- subset(gff, genome == thisGen)[,c("genome","ofID","ord","chr")]
    outgen[,ord2 := NULL]
    setnames(outgen, 1:3, c("refOfID","refOrd","ofID"))
    outgff <- merge(gout, outgen, by = "ofID", all.x = T, allow.cartesian = T)
    setorder(outgff, ord, na.last = T)
    outgff <- subset(outgff, !duplicated(paste(ofID, blkID)))
    ipn <- outgff[,list(n = sum(!is.na(refOrd))), by = "ofID"]
    if(verbose)
      with(ipn, cat(sprintf("%s / %s / %s / %s", sum(n == 0), sum(n == 1), sum(n == 2), sum(n > 2))))
    return(outgff)
  }))
  return(infout)
}

#' @title fill_missingGenePos
#' @description
#' \code{fill_missingGenePos} fill_missingGenePos
#' @rdname pangenome
#' @import data.table
#' @export
fill_missingGenePos <- function(hits,
                                infGenePos,
                                genomeIDs,
                                gff,
                                maxGap,
                                verbose){
  gen2 <- chr2 <- genome <- chr <- refOrd <- ofID <- closeGene <- ord <- NULL
  n <- closeOrd <- grp <- NULL
  uchr <- with(hits, unique(paste(gen2, chr2)))
  gffr <- subset(gff, paste(genome, chr) %in% uchr)
  inreg <- subset(infGenePos, !is.na(refOrd))
  gmiss <- subset(gffr, !ofID %in% inreg$ofID & !is.na(genome) & !is.na(chr))
  if(verbose)
    cat(sprintf("\n\tfound %s genes with no inferred position\n", nrow(gmiss)))
  splg <- split(gmiss, by = c("genome","chr"))
  splr <- split(inreg, by = c("genome", "chr"))[names(splg)]
  gfill <- rbindlist(lapply(names(splg), function(i){
    x <- splg[[i]]
    y <- splr[[i]]
    if(is.null(y)){
      x[,closeGene := NA]
    }else{
      x[,closeGene := sapply(ord, function(j)
        y$ofID[which.min(abs(j - y$ord))][1])]
    }
    x[,closeGene := as.character(unlist(closeGene))]
    return(x)
  }))
  ov <- gff$ord; names(ov) <- gff$ofID
  gfill[,closeOrd := ov[closeGene]]
  gfill$closeGene[abs(gfill$closeOrd - gfill$ord) > maxGap] <- NA
  if(verbose)
    cat(sprintf("\t%s genes are > %s genes from a block - no inferred pos. for these\n",
                sum(is.na(gfill$closeGene)), maxGap))
  gfill <- subset(gfill, !is.na(closeGene))
  iclose <- subset(inreg, ofID %in% gfill$closeGene)
  itmp <- with(iclose, data.table(
    closeGene = ofID, refOfID = refOfID, refOrd = refOrd,
    refChr = refChr, blkID = blkID))
  infGap <- merge(gfill[,c("ofID","genome","ord","chr","closeGene")],
                  itmp, by = "closeGene", allow.cartesian = T)
  infGap[,`:=`(blkID = NA, closeGene = NULL)]
  if(verbose)
    cat(sprintf("\tadded %s inferred positions for %s genes outside of blocks\n",
                nrow(infGap), uniqueN(infGap$ofID)))

  infPos <- rbind(infGenePos, infGap)
  gffHit <- subset(gff, paste(genome, chr) %in% uchr)
  icnt <- subset(infPos, !is.na(refOrd))[,list(n = .N),
                                         by = c("genome","ofID")]
  ocnt <- rbind(icnt, with(subset(gffHit, !ofID %in% icnt$ofID), data.table(
    genome = genome, ofID = ofID, n = 0)))
  ocnt[,grp := ifelse(n > 2, "2+", as.character(n))]
  ocnt <- ocnt[,list(cnt = .N), by = c("genome","grp")]
  ocnt[,genome := factor(genome, levels = genomeIDs)]
  setkey(ocnt, genome, grp)
  if(verbose)
    cat("\tFinal counts (genome: missing / 1x / 2x / >2x placements)\n")
  if(verbose){
    for(i in genomeIDs){
      cat(sprintf("\t\t%s: %s / %s / %s / %s \n",
                  i,
                  ocnt$cnt[ocnt$genome == i & ocnt$grp == "0"],
                  ocnt$cnt[ocnt$genome == i & ocnt$grp == "1"],
                  ocnt$cnt[ocnt$genome == i & ocnt$grp == "2"],
                  ocnt$cnt[ocnt$genome == i & ocnt$grp == "2+"]))
    }
  }
  return(infPos)
}

