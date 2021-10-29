#' @title Build genespace pangenome
#'
#' @description
#' \code{pangenome} Convert orthogroup and synteny information into a
#' pangenome database. Predict locations of orthogroups that are missing a
#' node in the reference.
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param refGenome character string matching one of the genomeIDs in gsParam
#' @param genomeIDs character vector, specifying which genomes to use. Defaults
#' to all genomeIDs specification in gsParam.
#' @param maxGapsBetweenEntries numeric of length 1, giving the maximum number
#' of pangenome entries between two unique orthogroup-prediction positions to
#' merge those positions. This cannot be smaller than the synteny buffer
#' in gsParam$params$synteny$synBuff and will be set to 2x that number if not specified
#' (default) or smaller than the synBuff maximum value.
#'
#' @details The pangenome annotation is a projection of syntenic orthogroups
#' on the physical coordinate system of a reference genome. The pangenome
#' function runs the following pipeline.
#' \enumerate{
#' \item within-block orthogroups and synteny-constrained global orthogroups
#' are merged.
#' \item the reference position is projected for all genes within syntenic
#' block bounds.
#' \item physical projected positions for each syntenic orthgroup is clustered,
#' permitting multiple reference locations (e.g. in a polyploid)
#' \item orthogroups with placements on the reference are populated and un-
#' placed orthgroups are added within reference position information
#' \item array members are added to the orthogroups
#' \item non-syntenic orthologs are added and flagged
#' }
#'
#' @return a data.table with lists of gene ids for each genome. Each row
#' corresponds to a unique combination of orthogroup and physical position
#' on the reference. Where multiple positions are inferred on a single
#' chromosome, the positions are broken out by the column 'clus'.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom dbscan dbscan frNN
#' @importFrom stats complete.cases median
#' @export
pangenome <- function(gsParam,
                      genomeIDs = NULL,
                      refGenome = NULL,
                      maxGapsBetweenEntries = NULL){
  ##############################################################################
  ##############################################################################
  # -- internal function to calculate gene positions by ref block coords
  calc_geneBlkReford <- function(anchCoords, gff, anchHits, gsParam, refGenome){
    ofID2 <- refChr <- ord <- chr <- genome <- ord1 <- ofID1 <- isAnchor <- ord2 <- NULL
    blkc <- subset(data.table(anchCoords), gen1 == refGenome)
    splg <- split(gff, by = c("genome", "chr"))
    spla <- split(subset(ba$anchors, gen1 == refGenome), by = "blkID")
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
      u <- with(a, unique(paste(ofID1, ofID2), paste(ofID2, ofID1)))
      h[,isAnchor := paste(ofID1, ofID2) %in% u]
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
  ##############################################################################
  ##############################################################################
  # -- internal function to assign and cluster orthogroups by ref position
  assign_og2reford <- function(gffRef, synBuff){
    rngOrd <- interpRefOrd <- clus <- NULL
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

  ##############################################################################
  # -- check the basic parameters
  u <- id <- pgID <- nonSynOrtho <- genome <- mem <- isArrayRep <- ofID <- NULL
  n <- og <- ord <- clus <- chr <- gen1 <- arrayID <- NULL
  # -- genomeIDs
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

  if(!is.data.table(gsParam$params$synteny))
    stop("Must run set_syntenyParams first!\n")

  # -- other needed parameters
  if(is.null(maxGapsBetweenEntries)){
    maxGapsBetweenEntries <- max(gsParam$params$synteny$synBuff)*2
  }else{
    maxGapsBetweenEntries <- as.numeric(maxGapsBetweenEntries[1])
    if(is.na(maxGapsBetweenEntries))
      stop("count not coerce maxGapsBetweenEntries to numeric\n")
  }
  mxBuff <- max(gsParam$params$synteny$synBuff)
  if(maxGapsBetweenEntries < mxBuff)
    maxGapsBetweenEntries <- mxBuff
  synBuff <- maxGapsBetweenEntries
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)
  verbose <- gsParam$params$verbose

  ##############################################################################
  # -- load synteny-constrained orthogroups
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(!file.exists(gffFile))
    stop("can't find the annotated gff-like text file\t\n ... have you run annotate_gff yet?\n")
  if(verbose)
    cat("Building pangenome source data ...\n\tReading gff ...")
  gffa <- fread(gffFile, showProgress = F, na.strings = c("NA", ""))
  gffa <- subset(gffa, genome %in% genomeIDs)
  gff <- subset(gffa, isArrayRep)
  if(verbose)
    cat(sprintf("\n\t... found %s genes %s array reps in %s merged OGs\n",
                nrow(gffa), nrow(gff), uniqueN(gff$combOG)))
  ##############################################################################
  # -- pull genes that are missing
  if(is.na(gsParam$paths$orthologuesDir)){
    cat("\tNo orthologues found ... will not be added to the pangenome\n")
    nsOrtho <- NA
  }else{
    nsOrtho <- pull_nonSynOrthologs(
      gff = gff,
      gsParam = gsParam)
    if(verbose)
      cat(sprintf("\tFlagged %s non-syntenic orthologs\n", nrow(nsOrtho)))
  }


  ##############################################################################
  # -- get the block coordinates
  if(verbose)
    cat(sprintf("\tPulling block coordinates against %s ...\n", refGenome))
  ba <- pull_blkAnchors(
    gsParam = gsParam,
    gff = gff,
    refGenome = refGenome)
  if(verbose)
    cat(sprintf("\t... found %s anchors for %s blocks\n",
                nrow(ba$anchors), nrow(ba$coords)))
  blkc <- data.table(ba$coords)

  ##############################################################################
  # -- get expected blocks and reference order for every gene
  if(verbose)
    cat("\tLinear interpolation of gene order in reference ... ")
  gffRef <- calc_geneBlkReford(
    refGenome = refGenome,
    anchCoords = ba$coords,
    anchHits = ba$anchors,
    gff = gff,
    gsParam = gsParam)

  ##############################################################################
  # -- summarize by orthogroup
  if(verbose)
    cat("Done!\n\tOrthogroup clustering and median position determination ... ")
  ogPos <- assign_og2reford(
    gffRef = gffRef,
    synBuff = synBuff)
  ogn <- ogPos[,list(n = .N), by = "combOG"]
  if(verbose)
    cat(sprintf("Done!\n\tPos. count: 0x = %s, 1x = %s, 2x = %s, 3x = %s 4x = %s, 4+x = %s",
        sum(!gff$combOG %in% ogn$combOG),
        sum(ogn$n == 1), sum(ogn$n == 2), sum(ogn$n == 3),
        sum(ogn$n == 4), sum(ogn$n > 4)))

  ##############################################################################
  # -- build pg
  if(verbose)
    cat("\nBuilding the pan-genome annotation ...\n")

  # -- combine the array rep gff with inferred positions of all OGs
  pg <- merge(
    with(gff, data.table(genome = genome, ofID = ofID, og = combOG)),
    with(ogPos, data.table(og = combOG, chr = refChr, clus = clus, ord = ord)),
    by = "og", all= T, allow.cartesian = T)

  # -- rename the pangenome IDs
  pg[,u := paste(chr, clus, ord, og)]
  pg[,`:=`(pgID = as.numeric(factor(u, levels = unique(u))), u = NULL)]
  if(verbose)
    with(pg, cat(sprintf(
      "\tInitial build with ... \n\t... %s genes, %s OGs and %s placements and %s unplaced OGs",
      uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))

  ##############################################################################
  # -- adding array members
  # -- load the raw gff
  if(verbose)
    cat("\n\tAdding in array members ... \n")

    # -- summarize by array, pulling the representative from the members
  arrRep <- subset(gffa, !is.na(arrayID))[,list(
    mem = ofID[!isArrayRep], rep = ofID[isArrayRep]),
    by = c("arrayID", "genome", "chr")]
  arrRep[,`:=`(inPg = rep %in% pg$ofID, ofID = rep)]

  # -- add the arrays to the pangenome with the reps
  pga <- merge(
    subset(pg, ofID %in% arrRep$rep),
    with(arrRep, data.table(ofID = rep, mem = mem)),
    by  = "ofID", allow.cartesian = T)

  # -- rename ofIDs from arrays
  pga[,`:=`(ofID = mem, mem = NULL)]
  with(pga, cat(sprintf(
    "\t... found %s genes, %s OGs and %s/%s placed/unplaced OGs",
    uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))

  # -- add arrays to the pg
  pg <- rbind(pg, pga)

  ##############################################################################
  # -- integrate orthologs
  if(!is.na(gsParam$paths$orthologuesDir)){
    if(verbose)
      cat("\n\tAdding in non-syntenic orthologs ... \n")
    gv <- gffa$genome; names(gv) <- gffa$ofID
    pgs <- merge(
      subset(pg, ofID %in% nsOrtho$ofID),
      with(nsOrtho, data.table(ofID = ofID, mem = orthIDs)),
      by  = "ofID", allow.cartesian = T)
    pgs[,`:=`(ofID = mem, mem = NULL, nonSynOrtho = TRUE)]

    # -- ensure the genomeIDs are right
    pgs[,genome := gv[ofID]]
    pgs <- subset(pgs, !duplicated(pgs))
    pg[,nonSynOrtho := FALSE]
    with(pgs, cat(sprintf(
      "\t... found %s genes, %s OGs and new %s entries\n",
      uniqueN(ofID), uniqueN(og), uniqueN(pgID), uniqueN(og[is.na(ord)]))))

    # -- combine pg with orthologs
    pg <- rbind(pg, pgs)
  }
  setorder(pg, ord, na.last = T)
  pg[,`:=`(pgID = as.numeric(factor(pgID, levels = unique(pgID))),
           og = as.numeric(factor(og, levels = unique(og))))]


  ##############################################################################
  # -- output and write
  if(verbose)
    cat("\tFormating and writing the pangenome ... ")
  pgout <- subset(pg, !is.na(ofID))

  # -- give real names
  iv  <- gffa$id; names(iv) <- gffa$ofID
  pgout[,id := iv[ofID]]

  # -- flag non-syn orthos
  pgout$id[pgout$nonSynOrtho] <- paste0(pgout$id[pgout$nonSynOrtho], "*")

  # -- reshape to wide format
  pgout <- dcast(pgout, pgID + og + chr + ord ~ genome,
                 value.var = "id", fun.aggregate = function(x) list(x))

  # -- order by pg position
  setorder(pgout, ord, na.last = T)

  # -- write text file
  pgf <- file.path(
    gsParam$paths$results,
    sprintf("%s_pangenomeDB.txt.gz", refGenome))
  fwrite(pg, file = pgf, sep = "\t", showProgress = F)
  if(verbose)
    cat(sprintf("Done!\nPangenome written to results/%s_pangenomeDB.txt.gz\n", refGenome))
  return(pgout)
}
