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
#' @param noSecondaryAnchors logical, should "secondary" (see synteny) hits be
#' used as pangenome position anchors? If so, like with polyploids, gene
#' positions may be duplicated
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
                      noSecondaryAnchors = TRUE,
                      refGenome = NULL){


  ##############################################################################
  # -- ad hoc function to pull non-syntenic orthologs
  pull_nonSynOrthologs <- function(gsParam,
                                   gff){

    gen1 <- og2 <- og1 <- orthIDs <- ofID <- id2 <- gen2 <- id1 <- NULL
    idv <- gff$ofID; names(idv) <- with(gff, paste(genome, id))
    ogv <- gff$og; names(ogv) <- gff$ofID

    orths <- rbindlist(lapply(unique(gff$genome), function(i){
      x <- parse_orthologues(
        gsParam = gsParam,
        refGenome = i)
      x[,`:=`(ofID = idv[paste(gen1, id1)],
              orthIDs = idv[paste(gen2, id2)])]
      x[,`:=`(og1 = ogv[ofID], og2 = ogv[orthIDs])]
      x <- subset(x, og1 != og2)
      return(x)
    }))
    return(orths[,c("ofID", "orthIDs", "orthID")])
  }

  pull_arrayRep <- function(gff){
    if(!"arrayID" %in% colnames(gff))
      stop("Cannot find a column named arrayID ... can't add arrayReps\n")

    if(!all(c("start", "ord", "pepLen") %in% colnames(gff)))
      stop("cannot find start, ord and/or pepLen in data.table columns\n")

    # -- split single and multi-member arrays
    n <- isArrayRep <- genome <- ord <- NULL
    gff[,n := .N, by = "arrayID"]
    out <- subset(gff, n == 1)
    out[,isArrayRep := TRUE]
    tmp <- subset(gff, n > 1)
    if(nrow(tmp) > 1){
      setkey(tmp, genome, ord)

      # -- calulate the distance to the median for each gene
      med <- ord <- medbp <- start <- dist2med <- dist2bp <- NULL
      tmp[,`:=`(med = as.numeric(median(ord)),
                medbp = as.numeric(median(start)))]
      tmp[,`:=`(dist2med = abs(med - ord),
                dist2bp = abs(medbp - start))]

      # -- rank and choose representatives
      pepLen <- dist2med <- dist2bp <- rnk <- isArrayRep <- med <- medbp <- NULL
      setorder(tmp, dist2med, dist2bp, -pepLen)
      tmp[,rnk := 1:.N, by = "arrayID"]
      tmp[,`:=`(isArrayRep = rnk == 1, rnk = NULL, dist2med = NULL,
                dist2bp = NULL, med = NULL, medbp = NULL)]

      # -- combine and return
      genome <- ord <- NULL
      out <- rbind(out, tmp)
      out[,n := NULL]
      setkey(out, genome, ord)
    }
    return(out)
  }

  interp_refPos <- function(gsParam,
                            refGenome,
                            nCores,
                            gff,
                            blkSize,
                            noSecondaryAnchors){

    ov <- gff$ord; ogv <- gff$og; names(ogv) <- names(ov) <- gff$ofID

    # -- find the syntenic hit files
    genome1 <- genome2 <- NULL
    writeTo <- gsParam$paths$results
    sp <- subset(gsParam$params$synteny,
                 genome1 == refGenome | genome2 == refGenome)
    synHitsFiles <- with(sp, file.path(
      writeTo, sprintf("%s_%s_synHits.txt.gz", genome1, genome2)))
    synHitsFiles <- synHitsFiles[file.exists(synHitsFiles)]

    # -- for each file, pull the hits
    regID <- blkID <- isOg <- ofID1 <- ofID2 <- isAnchor <- ord1 <- ord2 <-
      tmp1 <- tmp2 <- clus <- ofID <-
    synpos <- rbindlist(mclapply(synHitsFiles, mc.cores = nCores, function(i){
      syh <- fread(i, na.strings = c("", "NA"), showProgress = F)
      syh <- subset(syh, !is.na(regID))
      syh <- subset(syh, isOg & !grepl("self", blkID))
      if(noSecondaryAnchors)
        syh <- subset(syh, !grepl("second", blkID))
      syb <- subset(
        syh, !is.na(blkID) & ofID1 %in% names(ov) & ofID2 %in% names(ov))
      if(nrow(syb) > blkSize){
        if(syb$gen1[1] == refGenome){
          hits <- syb[,c("chr1", "ofID1", "ofID2", "blkID", "isAnchor")]
        }else{
          hits <- syb[,c("chr2", "ofID2", "ofID1", "blkID", "isAnchor")]
          setnames(hits, 1:3, c("chr1", "ofID1", "ofID2"))
        }
        anch <- subset(hits, isAnchor)
        anch[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]
        anch[,`:=`(tmp1 = frank(ord1, ties.method = "dense"),
                   tmp2 = frank(ord2, ties.method = "dense"))]
        anch[,clus := dbscan(frNN(cbind(tmp1, tmp2), eps = blkSize*2),
                             minPts = ceiling(blkSize/2))$cluster,
             by = c("blkID")]
        anch <- subset(anch, clus != 0)
        anch[,`:=`(tmp1 = NULL, tmp2 = NULL, clus = NULL)]
        spl <- split(anch, by = "blkID")
        # -- pull gff for genes without anchor hits to the reference
        interpos <- rbindlist(lapply(spl, function(x){
          chri <- x$chr1[1]
          blki <- x$blkID[1]
          wh <- which(gff$ofID %in% x$ofID2)
          tmp <- subset(gff[min(wh):max(wh),], !ofID %in% x$ofID2)
          if(nrow(tmp) > 1){
            y <- with(tmp, data.table(
              chr1 = chri, ofID1 = NA, ofID2 = ofID, blkID = blki, isAnchor = TRUE,
              ord1 = NA, ord2 = ord))
            x <- rbind(x, y)
            setkey(x, ord2)
            if(any(is.na(x$ord1)) & any(!is.na(x$ord1)))
              suppressWarnings(
                x[,ord1 := interp_linear(x = ord2, y = ord1)])
            out <- with(x, data.table(
              refOrd = ord1, ofID = ofID2, isAnchor = isAnchor,
              blkID = blkID[1], chr = chr1[1]))
            return(out)
          }else{
            return(NULL)
          }
        }))
        return(interpos)
      }else{
        return(NULL)
      }
    }))

    og <- NULL
    synpos[,og := ogv[ofID]]
    return(subset(synpos, complete.cases(synpos)))
  }
  ##############################################################################
  # 1. check the basic parameters
  # -- genomeIDs
  nCores <- gsParam$params$nCores
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
  blkSize <- max(gsParam$params$synteny$blkSize)

  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)
  verbose <- gsParam$params$verbose

  ##############################################################################
  # -- load synteny-constrained orthogroups
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(!file.exists(gffFile))
    stop("can't find the annotated gff-like text file\t\n ... have you run annotate_gff yet?\n")
  if(verbose)
    cat("Building reference pangenome source data ...\n")
  gf <- fread(gffFile, showProgress = F, na.strings = c("NA", ""))
  if(!"og" %in% colnames(gf))
    stop("can't find the og column name in the gff-like text file\n\tHas synteny been run yet?\n")

  if(!refGenome %in% gf$genome)
    stop(sprintf("%s (specified refGenome) not in the gff. Available genomes are: \n\t%s\n",
                 refGenome, paste(unique(gf$genome), collapse = ",")))

  # -- 1. Choose the representatives from the reference genome and build scaff.
  # -- Check if there are any syntenic blocks within the reference genome that
  # hit the same chromosome twice
  gen1 <- gen2 <- chr1 <- chr2 <- blkID <- isSelf <- firstGene1 <-
    firstGene2 <- lastGene2 <- lastGene1 <- NULL
  blks <- fread("results/syntenicBlocks.txt.gz", showProgress = F)
  blks <- subset(blks, !grepl("^NA_", blkID))
  blks[,isSelf := firstGene1 == firstGene2 | firstGene1 == lastGene2 | firstGene2 == lastGene1]
  tmp <- data.table(blks)
  setnames(tmp, gsub("1$","3",colnames(tmp)))
  setnames(tmp, gsub("2$","1",colnames(tmp)))
  setnames(tmp, gsub("3$","2",colnames(tmp)))
  blks <- rbind(blks[,colnames(blks) %in% colnames(tmp), with = F],
                tmp[,colnames(tmp) %in% colnames(blks), with = F], use.names = T)
  blks <- subset(blks, !duplicated(blks))
  refBlks <- subset(blks, gen1 == refGenome & gen1 == gen2)
  nflag <- nrow(subset(refBlks, !isSelf & chr1 == chr2))
  refBlks <- subset(refBlks, chr1 != chr2 | isSelf)
  if(noSecondaryAnchors)
    refBlks <- subset(refBlks, !grepl("second", blkID))
  if(nflag > 0)
    warning(sprintf("\tThe reference genome %s has multiple overlapping blocks \n\ton the same chr. This situation is not yet supported for pangenome \n\tconstruction. Dropping %s blocks from the pangenome\n",
                    refGenome, nflag))

  # -- subset the gff to chrs only in intergenic syntenic blocks
  gen1 <- gen2 <- blkID <- chr1 <- chr2 <- NULL
  genChr <- unlist(lapply(genomeIDs, function(i){
    b1 <- subset(blks, gen1 == i & !grepl("self", blkID) & gen2 != gen1)$chr1
    b2 <- subset(blks, gen2 == i & !grepl("self", blkID) & gen2 != gen1)$chr2
    return(paste(i, unique(c(b1, b2))))
  }))
  refBlks <- subset(
    refBlks, paste(gen1, chr1) %in% genChr & paste(gen2, chr2) %in% genChr)

  # -- recalculate arrayReps
  genome <- chr <- isArrayRep <- og <- arrayID <- NULL
  gchr <- subset(gf, paste(genome, chr) %in% genChr)
  gfrep <- subset(gchr, isArrayRep)
  gfrep[,arrayID := sprintf("%s_%s_%s", genome, chr, og)]
  gfrep <- pull_arrayRep(gfrep)
  gfrep <- subset(gfrep, isArrayRep)
  gref <- subset(gfrep, genome == refGenome)
  ov <- gfrep$ord; names(ov) <- gfrep$ofID

  # -- pull the self hits for the array reps and build a scaffold pg
  pg <- gref[,c("chr","ord","og","ofID")]
  setnames(pg, c("ofID", "ord","chr"), c("repID", "pgOrd", "pgChr"))

  # -- pull other array reps and add to the pangenome
  if(verbose)
    cat(sprintf("\tBuilt the scaffold pangenome with %s entries\n", nrow(pg)))

  ##############################################################################
  # -- 2. Build out a multi-genome pangenome
  # -- Determine the position of all genes against the reference
  if(verbose)
    cat("\tInterpolating syntenic positions of genes against reference ...")
  pgord <- interp_refPos(
    gsParam = gsParam,
    refGenome = refGenome,
    nCores = nCores,
    gff = gfrep,
    blkSize = blkSize,
    noSecondaryAnchors = noSecondaryAnchors)

  # -- choose the anchors and positions for each chr/og combination
  if(verbose)
    cat(" Done!\n\tChecking position accuracy ...")
  genome <- chr <- og <- chr <- refOrd <- med <- d2m <- pgChr <- pgOrd <- repID <- NULL
  uchr <- unique(pg$pgChr)
  pgAnch <- rbindlist(lapply(uchr, function(i){
    pgi <- subset(pg, pgChr == i)
    msk <- subset(gfrep, genome == refGenome & chr == i)
    ordi <- subset(pgord, chr == i & !og %in% pgi$og)
    ordi <- subset(ordi, !ofID %in% msk$ofID)
    setkey(ordi, og, refOrd)
    ordi[,med := median(refOrd), by = "og"]
    ordi[,d2m := abs(refOrd - med)]
    ordi <- subset(ordi, !duplicated(og))
    out <- rbind(pgi, with(ordi, data.table(
      pgChr = i, pgOrd = refOrd, og = og, repID = ofID)))
    out <- subset(out, !duplicated(og))
    return(out)
  }))

  # -- drop unlikely mapped positions
  og <- n <- prop <- keep <- pgChr <- og <- NULL
  for(k in 1:2){
    for(i in unique(pgAnch$pgChr)){
      ogs <- unique(pgAnch$og[pgAnch$pgChr == i])
      tmp <- subset(pgAnch, og %in% ogs)
      tmpn <- tmp[,list(n = uniqueN(og)), by = "pgChr"]
      tmpn[,prop := n / (sum(n))]
      tmpn[,keep := prop >= 0.05 | n > blkSize]
      z <- subset(tmp, pgChr %in% tmpn$pgChr[!tmpn$keep])
      pgAnch <- subset(pgAnch, !paste(pgChr, og) %in% paste(z$pgChr, z$og))
    }
  }

  # -- clean up by dcast/melt/unlist
  if(verbose)
    cat(" Done!\nBuilding full pangenome source data ...\n")
  og <- genome <- pgOrd <- pgID <- ofID <- NULL
  tmp <- dcast(gfrep, og ~ genome, value.var = "ofID", fun.aggregate = list)
  pg <- merge(pgAnch, tmp, by = "og", all = T, allow.cartesian = T)
  setorder(pg, pgOrd, na.last = T)
  pg[,pgID := 1:.N]
  pgm <- melt(
    pg, id.vars = c("pgID", "pgChr", "pgOrd", "og", "repID"),
    measure.vars = genomeIDs, value.name = "ofID", variable.name = "genome")
  pgm <- pgm[,list(ofID = unlist(ofID)),
             by = c("pgID", "pgChr", "pgOrd", "og", "repID", "genome")]


  # -- add in orthogroup reps for genes without placements
  repID <- pgOrd <- pgID <- NULL
  pgmis <- subset(pgm, is.na(repID))
  pgmis[,repID := ofID[1], by = "og"]
  pgm <- rbind(subset(pgm, !is.na(repID)), pgmis)
  setorder(pgm, pgOrd, na.last = T)

  if(verbose){
    cat(sprintf("\tBuilt the pangenome scaffold with %s unique positions ...\n",
                nrow(pg)))
    pgc <- pgm[,list(nGenome = uniqueN(genome)), by = "pgID"]
    cat(sprintf("\tn private entries = %s; all genomes = %s; PAV = %s\n",
                sum(pgc$nGenome == 1),
                sum(pgc$nGenome == length(genomeIDs)),
                sum(pgc$nGenome > 1 & pgc$nGenome< length(genomeIDs))))

    tmp <- subset(pgm, !duplicated(pgID))
    nun <- sum(is.na(tmp$pgChr))
    tab <- table(table(tmp$og))

    cat(sprintf("\tn single-copy entries = %s; 2x = %s; 2+x = %s; unplaced = %s\n",
                as.numeric(tab["1"]), as.numeric(tab["2"]),
                sum(as.numeric(tab[!names(tab) %in% c("1","2")])),
                sum(is.na(tmp$pgChr))))
  }

  # -- add in the array members
  if(verbose)
    cat("Annotating the pangenome ...\n\tAdding collinear array members ...")
  og <- ofID <- isArrayRep <- pgOrd <- NULL
  gfm <- subset(gf, og %in% unique(pgm$og) & !ofID %in% pgm$ofID)
  tmp <- pgm[,c("pgID", "pgChr", "pgOrd", "og", "repID")]
  tmp <- subset(tmp, !duplicated(tmp))
  pgarr <- merge(tmp, gfm[,c("og", "ofID", "genome")], by = "og")
  if(verbose)
    cat(sprintf("found %s genes, %s OGs and %s entires\n",
                uniqueN(pgarr$ofID), uniqueN(pgarr$og), uniqueN(pgarr$pgID)))

  # -- add in non-syntenic orthologs
  isNSortho <- ofID <- orthIDs <- pgOrd <- isArrayRep <- NULL
  pgarr[,isArrayRep := FALSE]
  pgm[,isArrayRep := TRUE]
  pg <- rbind(pgm, pgarr, use.names = T)
  setorder(pg, pgOrd, na.last = T)
  ogdir <- gsParam$paths$orthologuesDir
  pg[,isNSortho := FALSE]
  if(!is.na(ogdir) && dir.exists(ogdir)){
    if(verbose)
      cat("\tAdding non-syntenic orthologs ...")
    nso <- pull_nonSynOrthologs(gsParam = gsParam, gff = gf)
    gv <- gf$genome; names(gv) <- gf$ofID
    nso <- with(subset(nso, ofID %in% pg$repID), data.table(
      repID = ofID, ofID = orthIDs, genome = gv[orthIDs]))
    nso[,genome := gv[ofID]]
    tmp <- pg[,!colnames(pg) %in% c("genome", "ofID"), with = F]
    tmp <- subset(tmp, !duplicated(tmp) & !is.na(pgOrd))
    nso <- merge(tmp, nso, by = "repID")
    nso[,`:=`(isNSortho = TRUE, isArrayRep = FALSE)]
    setcolorder(nso, colnames(pg))
    pg <- rbind(pg, nso)
    if(verbose)
      cat(sprintf("Found %s\n", nrow(nso)))
  }


  # -- add in unplaced bottom drawer orthologs
  og <- repID <- tmp <- ofID <- isArrayRep <- genome <- pgID <- pgChr <-
    pgOrd <- isBottomDrawer <- isNSortho <- NULL
  gfm <- subset(gf, !og %in% unique(pgm$og))

  if(nrow(gfm) > 0){
    if(verbose)
      cat("\tAdding bottom-drawer non-syntenic orthogroups ...")
    gfm[,repID := ofID[isArrayRep][1], by = "og"]
    gfm[,tmp := ofID[1], by = "og"]
    gfm$repID[is.na(gfm$repID)] <- gfm$tmp[is.na(gfm$repID)]
    gfm <- dcast(gfm, og + repID + isArrayRep ~ genome,
                 value.var = "ofID", fun.aggregate = list)
    m <- max(pg$pgID)
    gfm[,pgID := (m + 1):(nrow(gfm) + m)]
    pgm <- melt(
      gfm, id.vars = c("pgID", "og", "repID", "isArrayRep"),
      measure.vars = genomeIDs[genomeIDs %in% colnames(gfm)],
      value.name = "ofID", variable.name = "genome")
    pgm[,`:=`(pgChr = NA, pgOrd = NA)]
    pgm <- pgm[,list(ofID = unlist(ofID)),
               by = c("pgID", "pgChr", "pgOrd", "og", "repID", "genome", "isArrayRep")]
    if(verbose)
      cat(sprintf("found %s entries", uniqueN(pgm$pgID)))
    pgm[,`:=`(isBottomDrawer = TRUE, isNSortho = FALSE)]
    pg <- rbind(pg, pgm[,colnames(pg), with = F])
  }

  ##############################################################################
  # -- output and write
  if(verbose)
    cat("\n\tFormating and writing the pangenome ... ")
  ofID <- genome <- id <- pgID <- pgChr <- pgOrd <- repID <- NULL
  pgout <- subset(pg, !is.na(ofID))
  setkey(pgout, pgID, genome)

  # -- give real names
  iv  <- gf$id; names(iv) <- gf$ofID
  pgout[,id := iv[ofID]]

  # -- flag non-syn orthos
  pgout$id[pgout$isNSortho] <- paste0(pgout$id[pgout$isNSortho], "*")

  # -- reshape to wide format
  pgout <- dcast(pgout, pgID + pgChr + pgOrd + repID ~ genome,
                 value.var = "id", fun.aggregate = function(x) list(x))

  # -- order by pg position
  setorder(pgout, pgOrd, na.last = T)

  # -- write text file
  pgf <- file.path(
    gsParam$paths$results,
    sprintf("%s_pangenomeDB.txt.gz", refGenome))
  fwrite(pg, file = pgf, sep = "\t", showProgress = F)
  if(verbose)
    cat(sprintf("Done!\nPangenome written to results/%s_pangenomeDB.txt.gz\n", refGenome))
  return(pgout)
}

