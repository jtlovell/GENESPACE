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
                      genomeIDs = NULL,
                      refGenome = NULL){
  calc_geneBlkReford <- function(anchCoords, gff, anchHits, gsParam){
    blkc <- data.table(anchCoords)
    splg <- split(gff, by = c("genome", "chr"))
    spla <- split(ba$anchors, by = "blkID")
    nCores <- gsParam$params$nCores
    nonRefPos <- rbindlist(mclapply(1:nrow(blkc), mc.cores = nCores, function(i){
      x <- blkc[i, ]
      g1 <- splg[[sprintf("%s.%s", refGenome, x$chr1)]]
      g2 <- splg[[sprintf("%s.%s", x$gen2, x$chr2)]]
      a <- spla[[x$blkID]]
      h <- merge(
        data.table(ofID1 = g1$ofID, og = g1$combOG, ord1 = g1$ord),
        data.table(ofID2 = g2$ofID, og = g2$combOG, ord2 = g2$ord),
        by = "og", allow.cartesian = T, all = T)
      h <- subset(h, !is.na(ofID2))
      setkey(h, ord2)
      h[,isAnchor := paste(ofID1, ofID2) %in% paste(a$ofID1, a$ofID2)]
      h$ord1[!h$isAnchor] <- NA
      h[,ord1 := interp_linear(x = ord2, y = ord1, interpTails = F)]
      h <- subset(h, !is.na(ord1))[,list(ord = median(ord1, na.rm = T)),
                                   by = "ofID2"]
      setnames(h, c("ofID", "interpRefOrd"))
      h[,`:=`(blkID = x$blkID, refChr = x$chr1)]
      return(h)
    }))
    gffRef <- merge(gff, nonRefPos, by = "ofID", allow.cartesian = T, all.x = T)
    gffs <- subset(gff, genome == refGenome)
    gffs[,`:=`(blkID = sprintf("blk_%s_%s_self", chr, chr),
               interpRefOrd = ord,
               refChr = chr)]
    gffRef <- subset(rbind(gffRef, gffs), !is.na(refChr))
    return(gffRef)
  }
  assign_og2reford <- function(gffRef, synBuff){
    gffRef[,rngOrd := diff(range(interpRefOrd)), by = c("combOG", "refChr")]
    g1 <- subset(gffRef, rngOrd <= synBuff)
    g1[,clus := 1]
    g2 <- subset(gffRef, rngOrd > synBuff)
    g2[,clus := dbscan(frNN(cbind(interpRefOrd, interpRefOrd), eps = synBuff),
                       minPts = 0)$cluster,
       by = c("combOG", "refChr")]
    gffRef <- rbind(g1, g2)
    ogPos <- gffRef[,list(ord = median(interpRefOrd, na.rm = T)),
                    by = c("combOG", "refChr", "clus")]
    return(ogPos)
  }

  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(is.null(refGenome))
    refGenome <- genomeIDs[1]
  synBuff <- max(gsParam$params$synteny$synBuff)
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)
  verbose <- gsParam$params$verbose
  # -- load synteny-constrained orthogroups
  gff <- combine_inblkSynOG(
    genomeIDs = genomeIDs,
    refGenome = refGenome,
    gsParam = gsParam)

  # -- pull genes that are missing
  nsOrtho <- pull_nonSynOrthologs(
    gff = gff,
    gsParam = gsParam)
  if(verbose)
    cat(sprintf("\tFlagged %s non-syntenic orthologs\n", nrow(nsOrtho)))

  # -- get the block coordinates
  if(verbose)
    cat(sprintf("\tPulling block coordinates against %s ...\n", refGenome))
  ba <- pull_blkAnchors(
    gsParam = gsParam,
    gff = gff,
    refGenome = refGenome,
    genomeIDs = genomeIDs)
  if(verbose)
    cat(sprintf("\t\tFound %s anchors for %s blocks\n",
                nrow(ba$anchors), nrow(ba$coords)))
  blkc <- data.table(ba$coords)

  # -- get expected blocks and reference order for every gene
  if(verbose)
    cat("Determining reference position for combined orthogroups ...\n\tLinear interpolation of gene order in reference ... ")
  gffRef <- calc_geneBlkReford(
    anchCoords = ba$coords,
    anchHits = ba$anchors,
    gff = gff,
    gsParam = gsParam)
  if(verbose)
    cat("Done!\n\tOrthogroup clustering and median position determination ... ")

  # -- summarize by orthogroup
  ogPos <- assign_og2reford(
    gffRef = gffRef,
    synBuff = synBuff)
  ogn <- ogPos[,list(n = .N), by = "combOG"]
  if(verbose)
    cat(sprintf("Done!\n\tPos. count: 0x = %s, 1x = %s, 2x = %s, 3x = %s 4x = %s, 4+x = %s",
        sum(!gff$combOG %in% ogn$combOG),
        sum(ogn$n == 1), sum(ogn$n == 2), sum(ogn$n == 3),
        sum(ogn$n == 4), sum(ogn$n > 4)))

  # -- build pg for ogs with placements on reference
  if(verbose)
    cat("\nBuilding the pan-genome annotation ...\n")
  pg <- merge(
    with(gff, data.table(genome = genome, ofID = ofID, og = combOG)),
    with(ogPos, data.table(og = combOG, chr = refChr, clus = clus, ord = ord)),
    by = "og", all= T, allow.cartesian = T)
  pg[,u:=paste(chr,clus,ord,og)]
  pg[,`:=`(pgID = as.numeric(factor(u,levels= unique(u))),u=NULL)]
  if(verbose)
    with(pg, cat(sprintf(
      "\tInitial build with ... \n\t\t%s genes, %s OGs and %s placements and %s unplaced OGs",
      uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))

  # -- adding array members
  if(verbose)
    cat("\n\tAdding in array members ... \n")
  gffFilea <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gffa <- fread(gffFilea, na.strings = c("", "-", "NA"), showProgress = F)
  gffa[,n := .N, by = "arrayID"]
  arrRep <- subset(gffa, n > 1)[,list(
    mem = ofID[!isArrayRep], rep = ofID[isArrayRep]),
    by = c("arrayID", "genome", "chr")]
  arrRep[,`:=`(inPg = rep %in% pg$ofID, ofID = rep)]
  pga <- merge(
    subset(pg, ofID %in% arrRep$rep),
    with(arrRep, data.table(ofID = rep, mem = mem)),
    by  = "ofID", allow.cartesian = T)
  pga[,`:=`(ofID = mem, mem = NULL)]
  with(pga, cat(sprintf(
    "\t\tFound %s genes, %s OGs and %s/%s placed/unplaced OGs",
    uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))
  pg <- rbind(pg, pga)

  # -- integrate orthologs
  if(verbose)
    cat("\n\tAdding in non-syntenic orthologs ... \n")
  gv <- gffa$genome; names(gv) <- gffa$ofID
  pgs <- merge(
    subset(pg, ofID %in% nsOrtho$ofID),
    with(nsOrtho, data.table(ofID = ofID, mem = orthIDs)),
    by  = "ofID", allow.cartesian = T)
  pgs[,`:=`(ofID = mem, mem = NULL, nonSynOrtho = TRUE)]
  pgs[,genome := gv[ofID]]
  pgs <- subset(pgs, !duplicated(pgs))
  pg[,nonSynOrtho := FALSE]
  with(pgs, cat(sprintf(
    "\t\tFound %s genes, %s OGs and new %s entries\n",
    uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))
  pg <- rbind(pg, pgs)
  setorder(pg, ord, na.last = T)
  pg[,`:=`(pgID = as.numeric(factor(pgID, levels = unique(pgID))),
           og = as.numeric(factor(og, levels = unique(og))))]
  if(verbose)
    cat("\tFormating and writing the pangenome ... ")
  # -- add non-syn orthos
  pg[,tmp := gv[ofID]]
  pgout <- data.table(pg)
  iv  <- gffa$id; names(iv) <- gffa$ofID
  pgout[,id := iv[ofID]]
  pgout$id[pgout$nonSynOrtho] <- paste0(pgout$id[pgout$nonSynOrtho], "*")
  pgout <- dcast(pgout, pgID + og + chr + ord ~ genome,
                 value.var = "id", fun.aggregate = function(x) list(x))
  setorder(pgout, ord, na.last = T)
  pgf <- file.path(
    gsParam$paths$results,
    sprintf("%s_pangenomeDB.txt.gz", refGenome))
  fwrite(pg, file = pgf, sep = "\t")
  if(verbose)
    cat(sprintf("\nPangenome written to results/%s_pangenomeDB.txt.gz", refGenome))
  return(pgout)
}

#' @title combine_inblkSynOG
#' @description
#' \code{combine_inblkSynOG} combine_inblkSynOG
#' @rdname pangenome
#' @import data.table
#' @export
combine_inblkSynOG <- function(refGenome,
                                genomeIDs,
                                gsParam){

  if(gsParam$params$verbose)
    cat("Combining synteny-constrained and inblock orthogroups ...\n")
  gffFile <- file.path(gsParam$paths$results, "gffWithSynOgs.txt.gz")
  gff <- fread(gffFile, showProgress = F, na.strings = c("-", "NA", ""))

  genomeIDs <- c(refGenome, genomeIDs[genomeIDs != refGenome])
  if(gsParam$params$verbose)
    cat(sprintf("\tsyn OGs: %s, inblk OGs: %s",
                uniqueN(gff$synOG), uniqueN(gff$inblkOG)))
  inblk <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "inblkOG"]
  syn <- gff[,list(ofID1 = ofID[-.N], ofID2 = ofID[-1]), by = "synOG"]
  u <- with(inblk, paste(ofID1, ofID2))
  syn <- subset(syn, !paste(ofID1, ofID2) %in% u & ofID1 != ofID2)
  tmp <- rbind(syn[,c("ofID1", "ofID2")], inblk[,c("ofID1", "ofID2")])
  tmp[,clus := clus_igraph(ofID1, ofID2)]
  ov <- with(tmp, c(clus, clus)); names(ov) <- with(tmp, c(ofID1, ofID2))
  gff[,combOG := ov[ofID]]
  nmis <- sum(is.na(gff$combOG))
  mol <- max(gff$combOG, na.rm = T)
  gff$combOG[is.na(gff$combOG)] <- (mol + 1):(mol + nmis)
  gff[,combOG := as.integer(factor(combOG, levels = unique(combOG)))]

  if(gsParam$params$verbose)
    cat(sprintf(", combined OGs: %s\n",
                uniqueN(gff$combOG)))
  return(gff)
}

#' @title pull_nonSynOrthologs
#' @description
#' \code{pull_nonSynOrthologs} pull_nonSynOrthologs
#' @rdname pangenome
#' @import data.table
#' @export
pull_nonSynOrthologs <- function(gsParam,
                                gff){
  idv <- gff$ofID; names(idv) <- with(gff, paste(genome, id))
  ogv <- gff$combOG; names(ogv) <- gff$ofID

  orths <- rbindlist(lapply(unique(gff$genome), function(i){
    x <- parse_orthologues(
      gsParam = gsParam,
      refGenome = i,
      nCores = gsParam$params$nCores)
    x[,`:=`(ofID = idv[paste(gen1, id1)], orthIDs = idv[paste(gen2, id2)])]
    x[,`:=`(og1 = ogv[ofID], og2 = ogv[orthIDs])]
    x <- subset(x, og1 != og2)
    return(x)
  }))
  return(orths[,c("ofID", "orthIDs", "orthID")])
}

#' @title pull_blkAnchors
#' @description
#' \code{pull_blkAnchors} pull_blkAnchors
#' @rdname pangenome
#' @import data.table
#' @export
pull_blkAnchors <- function(gsParam,
                            gff,
                            refGenome,
                            genomeIDs){
  # -- get the hit files
  fs <- list.files(path = gsParam$paths$results, pattern = "_synHits.txt.gz$")
  fs <- fs[grepl(sprintf("^%s_", refGenome), fs) |
             grepl(sprintf("%s_synHits.txt.gz", refGenome), fs)]
  fs <- fs[sapply(fs, function(x)
    any(sapply(genomeIDs, function(y)
      grepl(y, x))))]

  # -- block coords from hits
  out <- rbindlist(mclapply(fs, mc.cores = gsParam$params$nCores, function(i){
    x <- fread(file.path(gsParam$paths$results, i),
               select = c("ofID1", "ofID2","blkID","blkAnchor","gen1"),
               showProgress = F,
               na.strings = c("NA", "-", ""))
    if(x$gen1 != refGenome){
      setnames(x, c("ofID2", "ofID1","blkID","blkAnchor","gen1"))
      x <- x[,c("ofID1", "ofID2","blkID","blkAnchor","gen1")]
    }
    x[,isSelf := any(ofID1 %in% ofID2), by = "blkID"]
    return(subset(x, blkAnchor & !isSelf)[,1:3])
  }))

  gv <- gff$genome; ov <- gff$ord; cv <- gff$chr
  names(gv) <- names(ov) <- names(cv) <- gff$ofID
  out[,`:=`(gen2 = gv[ofID2], chr1 = cv[ofID1], chr2 = cv[ofID2],
            ord1 = ov[ofID1], ord2 = ov[ofID2])]
  bc <- out[,list(start1 = min(ord1, na.rm = T), end1 = max(ord1, na.rm = T),
                  start2 = min(ord2, na.rm = T), end2 = max(ord2, na.rm = T)),
            by = c("gen2", "chr1", "chr2", "blkID")]
  return(list(anchors = out, coords = bc))
}


