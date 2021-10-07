#' @title syntenicOGs
#' @description
#' \code{syntenicOGs} syntenicOGs
#'
#' @name syntenicOGs
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
#' @title pull_synOGs
#' @description
#' \code{pull_synOGs} pull_synOGs
#' @rdname syntenicOGs
#' @import data.table
#' @importFrom grDevices pdf dev.off
#' @export
pull_synOGs <- function(gsParam, genomeIDs = NULL){
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile, na.strings = c("", "-", "NA"), showProgress = F)
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
  if(verbose)
    cat(sprintf("\tn. syntenic OGs = %s\n", uniqueN(synGff$synOG)))

  # -- if desired, pull orthogroups within blocks
  if(gsParam$params$orthofinderInBlk){
    if(verbose)
      cat("Adding syntenic orthogroups from pairwise w/in-region hits\n")
    ofh <- blkwise_orthofinder(
      gsParam = gsParam,
      gff = synGff)
    ofh[, ogInblk := clus_igraph(ofID1, ofID2)]
    ofv <- ofh$ogInblk; names(ofv) <- ofh$ofID1
    synGff[,inblkOG := ofv[ofID]]
    mol <- max(synGff$inblkOG, na.rm = T)
    nmis <- sum(is.na(synGff$inblkOG))
    synGff$inblkOG[is.na(synGff$inblkOG)] <- (mol + 1):(mol + nmis)
    synGff[,inblkOG := as.integer(factor(inblkOG, levels = unique(inblkOG)))]
  }else{
    synGff[,inblkOG := NA]
  }
  out <- synGff[,c("genome", "chr", "start", "end", "ord", "id", "ofID", "globOG", "synOG", "inblkOG")]
  fwrite(out, file = file.path(gsParam$paths$results, "gffWithSynOgs.txt.gz"))
  return(out)
}

#' @title add_synOg2gff
#' @description
#' \code{add_synOg2gff} add_synOg2gff
#' @rdname syntenicOGs
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


#' @title blkwise_orthofinder
#' @description
#' \code{blkwise_orthofinder} blkwise_orthofinder
#' @rdname syntenicOGs
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom Biostrings readAAStringSet
#' @export
blkwise_orthofinder <- function(gsParam,
                                gff,
                                genomeIDs = NULL,
                                minGenes4of = 40,
                                verbose = NULL){

  #############################################################################
  # 1.Checking
  #############################################################################
  # -- get various parameters
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(is.null(verbose))
    verbose <- gsParam$params$verbose
  nCores <- gsParam$params$nCores
  # -- make sure orthofinder has been run and blast results exist
  if(is.na(gsParam$paths$blastDir))
    gsParam <- find_orthofinderResults(gsParam)

  # -- read in gff
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- get the synParams in
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, gnum1 <= gnum2)
  setkey(synp, gnum1, gnum2)

  # -- get orthofinder species IDs
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)

  # -- read in and split up the peptide files
  pepspl <- sapply(genomeIDs, USE.NAMES = T, simplify = F, function(i)
    readAAStringSet(file.path(gsParam$paths$blastDir,
                              sprintf("Species%s.fa", ofSpId[i]))))

  # -- run orthofinder within each region
  p5 <- function(x)
    paste(c(x, rep(" ", max(0, 7 - nchar(x)))), collapse = "")
  if(verbose)
    cat(sprintf(
      "\tRunning orthofinder by region ...\n\tGenome1-Genome2|Copy Number: %s%s%s%s\n",
      p5("absent"),p5("1x"),p5("2x"),p5("2+x")))

  blnk <- paste(rep(" ", 16), collapse = "")
  synOgInBlkHits <- rbindlist(lapply(1:nrow(synp), function(i){
    geno1 <- synp$genome1[i]
    geno2 <- synp$genome2[i]
    if(verbose)
      cat(sprintf("\t%s-%s|", pull_strWidth(geno1, 7),pull_strWidth(geno2, 7)))

    # -- load the syntenic hits
    fs <- file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",geno1, geno2))
    hits <- subset(fread(
      fs,
      na.strings = c("NA", ""),
      showProgress = F),
      !is.na(lgRegionID) & regBuffer)

    # -- calculate block coordinates from lgRegions (aggregated regions)
    hits[,`:=`(blkID = lgRegionID, blkAnchor = regAnchor, blkBuffer = regBuffer)]
    bv <- with(hits, c(blkID, blkID))
    names(bv) <- with(hits, c(paste(ofID1, ofID2), paste(ofID2, ofID1)))
    bv <- bv[!duplicated(names(bv))]

    expn <- count_expectn(hits = hits, gff = gff)
    tb <- subset(expn, !duplicated(paste(ofID, regID)))
    tb <- tb[,list(n = sum(!is.na(regID))), by = "ofID"]
    if(verbose)
      cat(sprintf("expected CN: %s%s%s%s\n",
                  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))
    oggl <- with(subset(hits, !is.na(og)), unique(paste(
      c(regID, regID), c(ofID1, ofID2))))
    expn[,hasOgGlob := paste(regID, ofID) %in% oggl]
    tb <- expn[,list(n = sum(!is.na(regID) & hasOgGlob)), by = "ofID"]
    if(verbose)
      cat(sprintf("\t%sglob. OG CN: %s%s%s%s\n",
                  blnk,  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))

    # -- read raw blast hits
    h <- read_hits4of(
      gsParam = gsParam,
      genome1 = geno1,
      genome2 = geno2)
    h[,blkID := bv[paste(V1, V2)]]
    h <- subset(h, !is.na(blkID))
    hc <- h[,list(n1 = uniqueN(V1), n2 = uniqueN(V2)), by = "blkID"]
    hc <- subset(hc, n1 >= minGenes4of & n2 >= minGenes4of)
    setorder(hc, -n1, -n2)
    # -- split hits by blocks with enuf genes
    splh <- sapply(split(h, by = "blkID"), simplify = F, USE.NAMES = T, function(x) x[,1:12])

    # -- run of in each block
    synOgHits <- rbindlist(mclapply(hc$blkID, mc.cores = nCores, function(j){
      hj <- splh[[j]]
      rj <- subset(hits, lgRegionID == j)
      tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogdt <- run_ofFromObj(
        blast00 = subset(h, V1 %in% rj$ofID1 & V2 %in% rj$ofID1),
        blast01 = subset(h, V1 %in% rj$ofID1 & V2 %in% rj$ofID2),
        blast10 = subset(h, V1 %in% rj$ofID2 & V2 %in% rj$ofID1),
        blast11 = subset(h, V1 %in% rj$ofID2 & V2 %in% rj$ofID2),
        pep0 = pepspl[[geno1]],
        pep1 = pepspl[[geno2]],
        writeDir = tmpDir)
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      ogv <- ogdt$og; names(ogv) <- ogdt$ofID
      rj[,`:=`(synOg1 = ogv[ofID1], synOg2 = ogv[ofID2])]
      return(subset(rj, !is.na(og) | synOg1 == synOg2))
    }))

    ogpw <- with(subset(synOgHits, synOg1 == synOg2), unique(paste(
      c(regID, regID), c(ofID1, ofID2))))

    expn[,hasOgPw := paste(regID, ofID) %in% ogpw]
    tb <- expn[,list(n = sum(!is.na(regID) & hasOgPw)), by = "ofID"]
    if(verbose)
      cat(sprintf("\t%sinblk OG CN: %s%s%s%s\n",
                  blnk,  p5(sum(tb$n == 0)), p5(sum(tb$n == 1)), p5(sum(tb$n == 2)), p5(sum(tb$n > 2))))

    return(with(subset(synOgHits, synOg1 == synOg2), data.table(
      ofID1 = ofID1, ofID2 = ofID2, synOg = synOg1)))
  }))

  return(synOgInBlkHits)
}

#' @title count_nBlks
#' @description
#' \code{count_nBlks} Gcount_nBlks
#' @rdname syntenicOGs
#' @import data.table
#' @export
count_nBlks <- function(gsParam, genomeIDs, minGenes4of){
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, gnum1 <= gnum2)
  setkey(synp, gnum1, gnum2)
  blks <- rbindlist(mclapply(1:nrow(synp), mc.cores = gsParam$params$nCores, function(i){
    fs <- with(synp[i,], file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",genome1, genome2)))
    h <- fread(
      fs, na.strings = c("NA", ""), showProgress = F,
      select = c("lgRegionID", "regBuffer","regAnchor","ofID1","ofID2","gen1","gen2","chr1","chr2","start1","start2","end1","end2","ord1","ord2"))
    h <- subset(h, complete.cases(h))
    h <- subset(h, regBuffer)
    h[,`:=`(blkID = lgRegionID, lgRegionID = NULL)]
    b <- calc_blkCoords(subset(h, regAnchor))
    return(b)
  }))
  if(all(blks$nHits1 >= minGenes4of) && all(blks$nHits2 >= minGenes4of)){
    cat(sprintf("\tAll %s regions have > %s genes\n", nrow(blks), minGenes4of))
  }else{
    bd <- subset(blks, nHits1 < minGenes4of | nHits2 < minGenes4of)
    bk <- subset(blks, nHits1 >= minGenes4of | nHits2 >= minGenes4of)
    cat(sprintf("\tUsing synteny-constrained global OGs for %s/%s regions/hits, with < %s hits\n",
                nrow(bd), sum(bd$nHits1), minGenes4of))
    cat(sprintf("\tRe-building syntenic orthogroups for %s regions with %s hits\n",
                nrow(bk), sum(bk$nHits1)))
  }
  return(blks)
}

#' @title read_hits4of
#' @description
#' \code{read_hits4of} read_hits4of
#' @rdname syntenicOGs
#' @import data.table
#' @export
read_hits4of <- function(gsParam, genome1, genome2){
  read_invertBlast <- function(gsParam, genome1, genome2, ofSpId, invert = T){
    h <- read_blast(
      path = gsParam$paths$blastDir, onlyIDScore = F,
      ofID1 = ofSpId[genome1],
      ofID2 = ofSpId[genome2])
    if(invert){
      h1 <- h[,c(2,1,3:6,8,7,10,9,11,12)]
      setnames(h1, colnames(h))
      h <- rbind(h, h1)
      setorder(h, -V12)
      h <- subset(h, !duplicated(paste(V1, V2)))
    }
    return(h)
  }
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  if(genome1 == genome2){
    h <- read_invertBlast(
      gsParam = gsParam, genome1 = genome1, genome2 = genome2, ofSpId = ofSpId)
  }else{
    h <- rbind(
      read_invertBlast(
        gsParam = gsParam, genome1 = genome1, genome2 = genome1,
        ofSpId = ofSpId, invert = F),
      read_invertBlast(
        gsParam = gsParam, genome1 = genome2, genome2 = genome2,
        ofSpId = ofSpId, invert = F),
      read_invertBlast(
        gsParam = gsParam, genome1 = genome1, genome2 = genome2,
        ofSpId = ofSpId, invert = F),
      read_invertBlast(
        gsParam = gsParam, genome1 = genome2, genome2 = genome1,
        ofSpId = ofSpId, invert = F))
  }
  return(h)
}

#' @title run_ofInReg
#' @description
#' \code{run_ofInReg} run_ofInReg
#' @rdname syntenicOGs
#' @import data.table
#' @export
run_ofInReg <- function(gsParam, genome1, genome2, blks, gff, pepspl){

  nCores <- gsParam$params$nCores
  geno1 <- genome1[1]
  geno2 <- genome2[1]
  genome1 <- genome2 <- NULL

  # -- subset the blks object to the genomes of interst
  b <- subset(blks, gen1 == geno1 & gen2 == geno2)
  setorder(b, -nHits1)

  # -- read in the hits
  fs <- file.path(
    gsParam$paths$results,
    sprintf("%s_%s_synHits.txt.gz",geno1, geno2))
  hits <- fread(
    fs,
    na.strings = c("NA", ""),
    showProgress = F)

  # -- read the raw blasts
  h <- read_hits4of(
    gsParam = gsParam,
    genome1 = geno1,
    genome2 = geno2)

  # -- parse the gff
  gs1 <- split(subset(gff, genome == geno1), by = "chr")
  gs2 <- split(subset(gff, genome == geno2), by = "chr")

  # -- for each block ...
  synOgHits <- rbindlist(mclapply(1:nrow(b), mc.cores = nCores, function(j){
    bj <- b[j,]
    hitj <- subset(hits, lgRegionID == bj$blkID)
    gj1 <- subset(
      gs1[[as.character(bj$chr1)]],
      ord >= bj$startOrd1 & ord <= bj$endOrd1)
    gj2 <- subset(
      gs2[[as.character(bj$chr2)]],
      ord >= bj$minOrd2 & ord <= bj$maxOrd2)
    uj <- unique(c(gj1$ofID, gj2$ofID))
    hj <- subset(h, V1 %in% uj & V2 %in% uj)
    hij <- subset(hits, ofID1 %in% uj & ofID2 %in% uj)
    tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    ogdt <- run_ofFromObj(
      blast00 = subset(hj, V1 %in% gj1$ofID & V2 %in% gj1$ofID),
      blast01 = subset(hj, V1 %in% gj1$ofID & V2 %in% gj2$ofID),
      blast10 = subset(hj, V1 %in% gj2$ofID & V2 %in% gj1$ofID),
      blast11 = subset(hj, V1 %in% gj2$ofID & V2 %in% gj2$ofID),
      pep0 = pepspl[[geno1]],
      pep1 = pepspl[[geno2]],
      writeDir = tmpDir)
    if(dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    ogv <- ogdt$og; names(ogv) <- ogdt$ofID
    hij[,`:=`(synOg1 = ogv[ofID1], synOg2 = ogv[ofID2])]
    return(subset(hij, synOg1 == synOg2 | !(is.na(og) | is.na(regID))))
  }))
  synOgHits <- subset(synOgHits, !duplicated(paste(ofID1, ofID2)))
  synOgHits[,cls := ifelse(is.na(og), "onlyInBlk", ifelse(synOg1 != synOg2, "onlyGlob", "both"))]
  return(synOgHits)
}

#' @title count_ogtypes
#' @description
#' \code{count_ogtypes} count_ogtypes
#' @rdname syntenicOGs
#' @import data.table
#' @export
count_ogtypes <- function(gff, hits){
  allGenes <- subset(gff, genome %in% c(hits$gen1[1], hits$gen2[1]))
  sogr <- subset(hits, synOg1 == synOg2 & !is.na(regID))
  sogr[,synOg := clus_igraph(ofID1, ofID2)]

  sogv <- with(sogr, c(synOg, synOg)); names(sogv) <- with(sogr, c(ofID1, ofID2))
  sogv <- sogv[!duplicated(names(sogv))]
  tabs <- table(table(sogv))
  nls <- unlist(list(
    no = sum(!allGenes$ofID %in% names(sogv)),
    na = sum(allGenes$ofID %in% names(sogv)),
    n1 = as.numeric(tabs["2"]*2),
    n2 = as.numeric(tabs["4"]*4)))

  sogr <- subset(hits, !is.na(og) & !is.na(regID))
  sogv <- with(sogr, c(og, og)); names(sogv) <- with(sogr, c(ofID1, ofID2))
  sogv <- sogv[!duplicated(names(sogv))]
  tabs <- table(table(sogv))
  nlg <- unlist(list(
    no = sum(!allGenes$ofID %in% names(sogv)),
    na = sum(allGenes$ofID %in% names(sogv)),
    n1 = as.numeric(tabs["2"]*2),
    n2 = as.numeric(tabs["4"]*4)))

  sogr <- subset(hits, !is.na(og))
  sogv <- with(sogr, c(og, og)); names(sogv) <- with(sogr, c(ofID1, ofID2))
  sogv <- sogv[!duplicated(names(sogv))]
  tabs <- table(table(sogv))
  nla <- unlist(list(
    no = sum(!allGenes$ofID %in% names(sogv)),
    na = sum(allGenes$ofID %in% names(sogv)),
    n1 = as.numeric(tabs["2"]*2),
    n2 = as.numeric(tabs["4"]*4)))
  outs <-  rbind(nls, nlg, nla)
  rownames(outs) <- c("inBlk", "syn", "glob")
  return(outs)
}


#' @title count_expectn
#' @description
#' \code{count_expectn} count_expectn
#' @rdname syntenicOGs
#' @import data.table
#' @export
count_expectn <- function(hits, gff){
  b <- hits[,list(start = min(ord1), end = max(ord1)),
            by = c("gen1","chr1","regID")]
  if(hits$gen1[1] != hits$gen2[1] ){
    b2 <- hits[,list(start = min(ord2), end = max(ord2)),
               by = c("gen2","chr2","regID")]
    setnames(b2, colnames(b))
    b <- rbind(b, b2)
  }
  b <- with(b, data.table(
    genome = gen1, chr = chr1, regID = regID, start = start, end = end))
  g <- with(subset(gff, genome %in% c(hits$gen1[1], hits$gen2[1])), data.table(
    genome = genome, chr = chr, start = ord, end = ord, ofID = ofID))
  setkey(g, genome, chr, start, end)
  setkey(b, genome, chr, start, end)
  f <- subset(foverlaps(g, b), !is.na(ofID))[,c("genome","regID","ofID")]
  f[,nReg := uniqueN(regID[!is.na(regID)]), by = "ofID"]
  return(f)
}

#' @title run_ofFromObj
#' @description
#' \code{run_ofFromObj} run_ofFromObj
#' @rdname syntenicOGs
#' @import data.table
#' @importFrom Biostrings writeXStringSet
#' @export
run_ofFromObj <- function(blast00,
                          blast01,
                          blast10,
                          blast11,
                          pep0,
                          pep1,
                          writeDir){

  if(dir.exists(writeDir))
    stop(sprintf("%s exists. Specify non-existing directory\n",
                 writeDir))
  dir.create(writeDir)

  if(colnames(blast00)[1] != "ofID1"){
    setnames(blast00, 1:2, c("ofID1", "ofID2"))
    setnames(blast01, 1:2, c("ofID1", "ofID2"))
    setnames(blast10, 1:2, c("ofID1", "ofID2"))
    setnames(blast11, 1:2, c("ofID1", "ofID2"))
  }

  # -- make gene ID dictionaries
  id0 <- unique(c(blast00$ofID1, blast00$ofID2, blast01$ofID1, blast10$ofID2))
  id1 <- unique(c(blast11$ofID1, blast11$ofID2, blast01$ofID2, blast10$ofID1))
  id0 <- id0[order(id0)]
  id1 <- id1[order(id1)]
  names(id0) <- sprintf("0_%s", (1:length(id0))-1)
  names(id1) <- sprintf("1_%s", (1:length(id1))-1)

  # -- ensure that all ids are in the peptide files
  id0 <- id0[id0 %in% names(pep0)]
  id1 <- id1[id1 %in% names(pep1)]

  # -- rename peptides and invert dictionary
  p0 <- pep0[id0]; names(p0) <- names(id0)
  p1 <- pep1[id1]; names(p1) <- names(id1)
  di0 <- names(id0); di1 <- names(id1); names(di0) <- id0; names(di1) <- id1

  # -- write the peptide files / fake diamond dbs
  writeXStringSet(p0, filepath = file.path(writeDir, "Species0.fa"))
  writeXStringSet(p1, filepath = file.path(writeDir, "Species1.fa"))
  cat("NA", file = file.path(writeDir, "diamondDBSpecies0.dmnd"))
  cat("NA", file = file.path(writeDir, "diamondDBSpecies1.dmnd"))

  # -- write the species and sequence IDs
  sid <- data.table(
    of = paste0(c(names(id0), names(id1)), ":"),
    id = c(id0, id1))
  fwrite(
    sid, file = file.path(writeDir, "SequenceIDs.txt"),
    sep = " ", quote = F, row.names = F, col.names = F)
  cat(
    c("0: species1.fa", "1: species2.fa"),
    sep = "\n", file = file.path(writeDir, "SpeciesIDs.txt"))

  # -- rename the blast files
  bl00 <- subset(blast00, ofID1 %in% names(di0) & ofID2 %in% names(di0))
  bl01 <- subset(blast01, ofID1 %in% names(di0) & ofID2 %in% names(di1))
  bl10 <- subset(blast10, ofID1 %in% names(di1) & ofID2 %in% names(di0))
  bl11 <- subset(blast11, ofID1 %in% names(di1) & ofID2 %in% names(di1))
  bl00[,`:=`(ofID1 = di0[ofID1], ofID2 = di0[ofID2])]
  bl01[,`:=`(ofID1 = di0[ofID1], ofID2 = di1[ofID2])]
  bl10[,`:=`(ofID1 = di1[ofID1], ofID2 = di0[ofID2])]
  bl11[,`:=`(ofID1 = di1[ofID1], ofID2 = di1[ofID2])]
  bl00 <- subset(bl00, complete.cases(bl00[,1:12]))[,1:12]
  bl01 <- subset(bl01, complete.cases(bl01[,1:12]))[,1:12]
  bl10 <- subset(bl10, complete.cases(bl10[,1:12]))[,1:12]
  bl11 <- subset(bl11, complete.cases(bl11[,1:12]))[,1:12]

  # -- write the blasts
  fwrite(
    bl00, file = file.path(writeDir, "Blast0_0.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl01, file = file.path(writeDir, "Blast0_1.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl10, file = file.path(writeDir, "Blast1_0.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl11, file = file.path(writeDir, "Blast1_1.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)

  # -- run orthofinder
  com <- sprintf("orthofinder -b %s -og -a 1 -t 1  1>/dev/null 2>&1", writeDir)
  system(com)

  # -- find the files
  ogf <- order_filesByMtime(
    path = writeDir,
    recursive = T,
    pattern = "Orthogroups.tsv")[1]

  # -- read the orthogroups.tsv file and process
  ogdt <- fread(ogf, showProgress = F, verbose = F)
  ogdt <- melt(
    ogdt, id.vars = "Orthogroup", variable.name = "genome", value.name = "id")
  ogdt <- ogdt[,list(id = strsplit(id, ",")[[1]]), by = c("Orthogroup", "genome")]
  ogdt[,`:=`(genome = NULL, ofID = trimws(id), id = NULL,
             og = trimws(Orthogroup), Orthogroup = NULL)]
  ogdt <- subset(ogdt, !duplicated(ogdt))
  hasDup <- subset(ogdt, ofID %in% subset(ogdt, duplicated(ofID))$ofID)
  if(nrow(hasDup) > 1){
    m <- merge(hasDup, hasDup, by = "ofID", all = T, allow.cartesian = T)
    ci <- clus_igraph(m$og.x, m$og.y)
    ci <- ci[!duplicated(names(ci))]
    ogdt[,og := ifelse(og %in% names(ci), ci[og], og)]
    ogdt <- subset(ogdt, !duplicated(ogdt))
  }
  ogdt[,og := as.integer(factor(og, unique(og)))]

  nog <- c(id0, id0)
  nog <- nog[!nog %in% ogdt$ofID]

  # -- return data.table of ogs
  if(length(nog) > 0){
    ogdt <- rbind(ogdt, data.table(
      ofID = nog, og = sprintf("NOG_%s",1:length(nog))))
  }

  return(ogdt)
}


