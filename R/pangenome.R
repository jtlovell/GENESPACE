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
#' @param propAssignThresh numeric of length 1, the minimum proportion of genes
#' in a pangenome entry needed to anchor the physical position to a new location
#' @param maxPlacementsPerRefChr integer of length 1, specifying the maximum
#' number of pan-genome entries that a gene can belong to on a single
#' reference chromosome.
#' @param maxNonSynOrthos2keepPerGenome integer of length 1 specifying the
#' number of non-syntenic orthologs to include in the pan-genome annotation.
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
                      maxPlacementsPerRefChr = 2,
                      propAssignThresh = 0.5,
                      maxNonSynOrthos2keepPerGenome = 5){

  ##############################################################################
  ##############################################################################
  # -- ad hoc function to pull non-syntenic orthologs
  pull_nonSynOrthologs <- function(gsParam,
                                   synParamsDt,
                                   nCores,
                                   gff,
                                   pangenomeDt,
                                   maxNonSynOrthos2keepPerGenome){

   gen1 <- id1 <- gen2 <- id2 <- ofID1 <- pgID1 <- pgID2 <- n <- pgID <-
     ofID <- NULL

    # -- make vector of gene OF IDs
    idv <- gff$ofID
    names(idv) <- with(gff, paste(genome, id))

    pg1 <- with(pangenomeDt, data.table(ofID1 = ofID, pgID1 = pgID))
    pg2 <- with(pangenomeDt, data.table(ofID2 = ofID, pgID2 = pgID))
    pg1 <- subset(pg1, !duplicated(pg1))
    pg2 <- subset(pg2, !duplicated(pg2))
    # -- for each genomeID
    orths <- rbindlist(lapply(genomeIDs, function(i){

      # -- read in all orthologs
      x <- parse_orthologues(
        gsParam = gsParam,
        refGenome = i)

      # -- rename to unique orthofinderIDs
      x[,`:=`(ofID1 = idv[paste(gen1, id1)],
              ofID2 = idv[paste(gen2, id2)],
              genome = gen2)]
      x <- x[,c("ofID1", "ofID2", "genome")]
      x <- subset(x, ofID1 %in% pangenomeDt$ofID[pangenomeDt$isRep])
      x <- merge(pg1, x, by = "ofID1", allow.cartesian = T, all.y = T)
      x <- subset(x, !duplicated(x))
      x <- merge(pg2, x, by = "ofID2", allow.cartesian = T, all.y = T)
      x <- subset(x, !duplicated(x))
      x <- subset(x, pgID1 != pgID2 | is.na(pgID1) | is.na(pgID2))
      print(x)
      x[,n := 1:.N, by = c("ofID1", "genome")]
      print(x)
      x <- subset(x, n <= maxNonSynOrthos2keepPerGenome)
      x <- with(x, data.table(pgID = pgID1, ofID = ofID2, genome = genome))
      x <- subset(x, !duplicated(x))
      # -- subset to only hits that are not in the pangenome
      return(x)
    }))
    orths <- subset(orths, !duplicated(orths))
    pg2fill <- subset(pangenomeDt, pgID %in% orths$pgID)
    pg2fill <- pg2fill[,c("pgChr", "pgOrd", "pgID", "og")]
    pg2fill <- subset(pg2fill, !duplicated(pg2fill))

    pgns <- merge(pg2fill, orths, by = "pgID", allow.cartesian = T)

    ov <- gff$ord; cv <- gff$chr; gv <- gff$genome
    names(ov) <- names(cv) <- names(gv) <- gff$ofID

    pgns[,`:=`(isRep = FALSE, chr = cv[ofID], ord = ov[ofID],
               isDirectSyn = FALSE, isSynOgOnly = FALSE, isArrayRep = FALSE)]

    pangenomeDt <- rbind(pangenomeDt, pgns[,colnames(pangenomeDt), with = F])
    return(pangenomeDt)
  }

  ##############################################################################
  ##############################################################################
  # -- internal function to read in hits against a ref genome
  read_hits4pg <- function(synParamsDt,
                           pgRepOfIDs,
                           refGenome,
                           gff,
                           nCores){

    hitsFile <- genome1 <- genome2 <- regID <- isRep1 <- isRep2 <- ofID1 <-
      ofID2 <- gen1 <- NULL

    # -- check that the hits files are in the synParam obj and keep only those
    # that exist.
    if(!"hitsFile" %in% colnames(synParamsDt))
      stop("need to add hits file path to synParams\n")
    synParamsDt <- subset(synParamsDt, file.exists(hitsFile))
    if(nrow(synParamsDt) < 1)
      stop("no synHits files exist\n")

    # -- get vector of orthogroup IDs
    ov <- gff$og; names(ov) <- gff$ofID

    # -- split into two groups, one with the ref as genome1
    altSyn1 <- subset(synParamsDt, genome1 == refGenome)
    altSyn2 <- subset(synParamsDt, genome2 == refGenome & genome1 != refGenome)
    cols <- c("gen1", "gen2", "ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2",
              "blkID", "regID", "isAnchor", "inBuffer")

    # -- for each group ...
    if(nrow(altSyn1) > 0){
      alt1 <- rbindlist(mclapply(altSyn1$hitsFile, mc.cores = nCores, function(i){

        # -- read in the synHits file
        tmp <- fread(
          i, showProgress = FALSE, na.strings = c("", "NA"),
          select = c(cols, "isRep1", "isRep2"))

        # -- keep only syntenic orthogroups and array reps
        tmp <- subset(tmp, !is.na(regID) & isRep1 & isRep2)
        tmp <- subset(tmp, ov[ofID1] == ov[ofID2])

        # -- keep only hits in the pgRepOfIDs vector
        tmp <- subset(tmp, ofID1 %in% pgRepOfIDs)
        return(tmp[, cols, with = F])
      }))
    }else{
      alt1 <- NULL
    }
    if(nrow(altSyn2) > 0){
      tmpCols <- c("gen2", "gen1", "ofID2", "ofID1", "chr2", "chr1", "ord2",
                   "ord1", "blkID", "regID", "isAnchor", "inBuffer")
      alt2 <- rbindlist(mclapply(altSyn2$hitsFile, mc.cores = nCores, function(i){
        tmp <- fread(i, showProgress = FALSE, na.strings = c("", "NA"),
                     select = c(cols, "isRep1", "isRep2"))
        tmp <- subset(tmp, !is.na(regID) & isRep1 & isRep2)
        tmp <- subset(tmp, ov[ofID1] == ov[ofID2])
        tmp <- subset(tmp, ofID2 %in% pgRepOfIDs)[,tmpCols, with = F]
        setnames(tmp, cols)
        return(tmp)
      }))
    }else{
      alt2 <- NULL
    }

    # -- combine the hits
    hits <- rbind(alt1, alt2)

    # -- change the names to match the pangenome
    setnames(hits, c("ofID1", "chr1", "ord1"),
             c("repOfID", "repChr", "repOrd"))
    setnames(hits, c("ofID2", "chr2", "ord2", "gen2"),
             c("ofID", "chr", "ord", "genome"))
    hits[,gen1 := NULL]
    return(hits)
  }

  ##############################################################################
  ##############################################################################
  # -- internal function to re-cluster syntenic anchors via micro collinearity
  clus_synAnchorPos <- function(ord1, ord2, blkSize, include){

    xc <- x <- yc <- y <- clus <- index <- NULL


    # -- convert to vectors
    ord1 <- as.integer(ord1)
    ord2 <- as.integer(ord2)

    # -- check that orders are specified correctly
    if(all(is.na(ord1)) || all(is.na(ord2)))
      stop("only NA values found ... cannot conduct linear interpolation\n")
    if(length(ord1) != length(ord2))
      stop("order vectors must be of the same length\n")

    # -- make a data table with these data and remove rows with NAs
    ix <- data.table(index = 1:length(ord1))
    dt <- data.table(index = 1:length(ord1), x = ord1, y = ord2)

    if(length(include) != length(ord1) || !any(include))
      stop("problem with include vector, which should be logical of the same length as ord1\n")
    dt <- subset(dt, include)

    # -- set the minimum retained block size and check that its ok
    minJitBlk <- ceiling(blkSize/2)
    if(is.na(minJitBlk) || minJitBlk < 1)
      stop("blkSize must be posive and numeric\n")

    # -- one last pruning and classify so that clus can be used
    dt[,xc := frank(x, ties.method = "dense")]
    dt[,yc := frank(y, ties.method = "dense")]

    # -- if enough observations, do the clustering
    if(nrow(dt) >= minJitBlk){
      dt[,clus := dbscan(
        frNN(cbind(xc, yc), eps = minJitBlk), minPts = minJitBlk)$cluster]
    }else{
      dt[,clus := 0]
    }

    # -- subset to hits in clusters
    dt <- subset(dt, clus > 0)

    # -- get the range of original index's in each cluster
    rng <- dt[,list(index = min(index):max(index)), by = "clus"]

    # -- subset to only unique entries and merge, keeping all in the original
    rng <- subset(rng, !duplicated(index))
    out <- merge(ix, rng, by = "index", all.x = T)
    if(nrow(out) != length(ord1))
      stop("problem matching input and output order string\n")
    return(out$clus)
  }

  ##############################################################################
  ##############################################################################
  # -- internal function to make sure the syntenic anchors are in good order
  check_synAnchorPos <- function(ord1, ord2, blkSize){

    xc <- x <- yc <- y <- clus <- index <- NULL

    # -- convert to vectors
    ord1 <- as.integer(ord1)
    ord2 <- as.integer(ord2)
    excl <- is.na(ord1) | is.na(ord2)

    # -- check that orders are specified correctly
    if(all(is.na(ord1)) || all(is.na(ord2)))
      stop("only NA values found ... cannot conduct linear interpolation\n")
    if(length(ord1) != length(ord2))
      stop("order vectors must be of the same length\n")

    # -- make a data table with these data and remove rows with NAs
    ix <- data.table(index = 1:length(ord1))
    dt <- data.table(index = 1:length(ord1), x = ord1, y = ord2)
    dt <- subset(dt, complete.cases(dt))

    # -- set the minimum retained block size and check that its ok
    minJitBlk <- ceiling(blkSize/2)
    if(is.na(minJitBlk) || minJitBlk < 1)
      stop("blkSize must be posive and numeric\n")

    # -- iteratively remove offending hits that are not in exact synteny
    for(j in 1:minJitBlk){

      # -- rank xy
      dt[,xc := frank(x, ties.method = "dense")]
      dt[,yc := frank(y, ties.method = "dense")]

      # -- cluster
      if(nrow(dt) >= j){
        dt[,clus := dbscan(
          frNN(cbind(xc, yc), eps = j), minPts = j)$cluster]
      }else{
        dt[,clus := 0]
      }

      # -- remove unclustered
      dt <- subset(dt, clus != 0)
    }

    # -- re-rank if any clustered
    if(nrow(dt) >= minJitBlk){
      dt[,xc := frank(x, ties.method = "dense")]
      dt[,yc := frank(y, ties.method = "dense")]

      # -- recluster
      dt[,clus := dbscan(
        frNN(cbind(xc, yc), eps = minJitBlk), minPts = minJitBlk)$cluster]
      dt <- subset(dt, clus != 0)
    }else{
      dt[,clus := NA]
    }

    # -- one last pruning and classify so that clus can be used
    dt <- merge(ix, dt[,c("index", "clus")], by = "index", all.x = T)
    setkey(dt, index)
    if(nrow(dt) != length(ord1))
      stop("problem matching input and output order string\n")
    return(!is.na(dt$clus) & !excl)
  }

  ##############################################################################
  ##############################################################################
  # -- internal function to get reference positions.
  # this function has ad hoc (non-standard) input: anchorHits
  # due to this input, keeping function internal to pangenome to avoid confusion
  # during documentation
  interp_refPos <- function(gff,
                            anchorHits,
                            blkSize,
                            refGenome,
                            nCores,
                            verbose){

    isArrayRep <- index <- isSelf <- ind1 <- ind2 <- ind2 <- interpRefOrd <-
      ord1 <- useAsAnch <- ord2 <- clus <- ord2 <- ofID2 <- interpRefOrd <-
      ord1 <- interpRefOrd <- ofID <- refChr <- NULL

    # -- get indexed position of geneIDs
    gff <- subset(gff, isArrayRep)
    gff[,index := 1:.N]
    indv <- gff$index; names(indv) <- gff$ofID
    ah <- with(anchorHits, data.table(
      gen1 = refGenome, gen2 = genome,
      chr1 = repChr, chr2 = chr,
      ord1 = repOrd, ord2 = ord,
      ind1 = indv[repOfID], ind2 = indv[ofID],
      ofID1 = repOfID, ofID2 = ofID,
      blkID = sprintf("%s_%s", genome, blkID)))

    # -- make block coords with index non self
    ah[,isSelf := any(ind1 == ind2), by = "blkID"]
    blkIndex <- subset(ah, !isSelf)[,list(
      st1 = min(ind1), st2 = min(ind2), en1 = max(ind1), en2 = max(ind2)),
      by = c("blkID", "chr1", "chr2", "gen1", "gen2")]

    # -- for each block, get anchors, non-position anchors
    gtmp1 <- gff[,c("index","ofID", "ord")]
    gtmp2 <- gff[,c("index","ofID", "ord")]
    setnames(gtmp1, c("ind1","ofID1", "ord1"))
    setnames(gtmp2, c("ind2","ofID2", "ord2"))

    splAnchs <- split(ah[,c("blkID", "ofID1", "ofID2", "ord1", "ord2")], by = "blkID")

    # -- for each block, get the interpolated reference order
    splH <- rbindlist(mclapply(1:nrow(blkIndex), mc.cores = nCores, function(i){
      z <- blkIndex[i,]
      x <- data.table(splAnchs[[z$blkID]][,c("ofID1", "ofID2", "ord1", "ord2")])
      y2 <- subset(gtmp2[blkIndex$st2[i]:blkIndex$en2[i],], !ind2 %in% x$ind2)
      y <- rbind(x, y2, fill = T)

      if(all(complete.cases(y[,c("ord1", "ord2")]))){
        y[,interpRefOrd := ord1]
      }else{
        # -- cluster and check the positions
        # -- if not enough anchor regions, just use them all and one cluster
        if(sum(!is.na(y$ord1)) < blkSize){
          y[,`:=`(useAsAnch = !is.na(ord1), clus = 1)]
        }else{
          # -- re-cluster the anchors, make sure they are perfectly collinear
          y[,useAsAnch := check_synAnchorPos(
            ord1 = ord1, ord2 = ord2, blkSize = blkSize)]

          # -- cluster the anchors into groups if multiple regions
          if(any(y$useAsAnch)){
            y[,clus := clus_synAnchorPos(
              ord1 = ord1, ord2 = ord2, blkSize = blkSize, include = useAsAnch)]
          }else{
            y[,clus := NA]
          }
        }

        # -- pull and split clustered hits
        yclus <- subset(y, !is.na(clus) & !is.na(ord1) & !is.na(ord2))
        spl <- split(yclus, by = "clus")
        y2interp <- subset(y, !is.na(ord2))
        setorder(y2interp, ord1, na.last = T)

        # -- get hits within each new micro block and interpolate
        y <- rbindlist(lapply(spl, function(cl){
          tmp <- subset(y2interp, ord2 >= min(cl$ord2) & ord2 <= max(cl$ord2))
          tmp <- subset(tmp, !duplicated(ofID2))
          tmp[,interpRefOrd := interp_linear(refOrd = ord1, toInterpOrd = ord2)]
          return(tmp)
        }))
      }

      intp <- data.table(
        blkID = z$blkID, ofID = y$ofID2,
        refChr = z$chr1, interpRefOrd = y$interpRefOrd)
      return(intp)
    }), fill = T)
    out <- subset(splH, complete.cases(splH))

    # -- make sure we don't have duplicated NAs
    out <- subset(out, !is.na(interpRefOrd))

    # -- merge with gff
    gfm <- merge(gff, out, by = "ofID", all = T, allow.cartesian = T)
    setorder(gfm, interpRefOrd, na.last = T)
    gfY <- subset(gfm, !is.na(interpRefOrd))
    gfN <- subset(gfm, is.na(interpRefOrd))
    gfN <- subset(gfN, !ofID %in% gfY$ofID)
    tab <- table(gfY$ofID)
    if(verbose)
      cat(sprintf(
        "n. genes mapped: 1x = %s, 2+x = %s, 0x = %s\n",
        sum(tab == 1), sum(tab > 1), nrow(gfN)))
    intPos <- rbind(gfY, gfN)
    intPos <- with(subset(intPos, !is.na(refChr) & !is.na(interpRefOrd)),
                   data.table(ofID = ofID, altChr = chr, altOrd = ord,
                              interpRefChr = refChr, interpRefOrd = interpRefOrd))
    return(intPos)
  }

  ##############################################################################
  # -- internal function to build reference pangenome
  build_refPg <- function(refAnchorHits,
                          interpolatedPositionDT,
                          gff,
                          genomeIDs,
                          refGenome,
                          propAssignThresh,
                          maxPlacementsPerRefChr){

    genome <- isArrayRep <- og <- id1 <- id2 <- ofID <- clus <- lev <- med <-
      interpRefChr <- interpRefOrd <- nInChr <- rat <- nOnChr <- tmpi <- nPos <-
      clusThis <- dist2med <- isSameChr <- chr <- n <- pgID <- altID <- NULL

    # -- get full set of positions
    ipInterp <- with(interpolatedPositionDT, data.table(
      ofID = ofID, chr = altChr, ord = altOrd,
      interpRefChr = interpRefChr,  interpRefOrd = interpRefOrd))
    ipGff <- with(subset(gff, genome == refGenome & isArrayRep), data.table(
      ofID = ofID, chr = chr, ord = ord,
      interpRefChr = chr, interpRefOrd = ord))
    ip <- rbind(ipGff, ipInterp)
    ip <- subset(ip, complete.cases(ip))

    # -- get vector of position info
    ov <- gff$ord; cv <- gff$chr; gv <- gff$genome
    names(ov) <- names(cv) <- names(gv) <- gff$ofID

    # -- get vector of orthogroups based on hits (more specific than global)
    tmp <- with(refAnchorHits, data.table(
      id1 = c(ofID, repOfID), id2 = c(repOfID, ofID)))
    tmp <- subset(tmp, !duplicated(tmp))
    tmp[,og := clus_igraph(id1, id2)]
    ogv <- tmp$og; names(ogv) <- tmp$id1
    ip[,og := ogv[ofID]]

    # -- toss chromosomes that are rare (default < 25%)
    ip[,clus := 1]
    setorder(ip, interpRefChr, interpRefOrd, na.last = T)
    tmp <- subset(ip, !is.na(interpRefChr))
    tmp[,nInChr := uniqueN(ofID), by = c("interpRefChr", "og")]
    tmp[,n := uniqueN(ofID), by = "og"]
    tmp[,rat := nInChr / n]
    tmp <- subset(tmp, rat >= propAssignThresh)[,colnames(ip), with = F]

    # -- separate groups that need to be clustered from those that don't
    tmp[,tmpi := round(interpRefOrd)]
    tmp[,nPos := uniqueN(tmpi), by = c("og", "interpRefChr")]

    pg1 <- subset(tmp, nPos == 1)[,colnames(ip), with = F]
    pg2 <- subset(tmp, nPos > 1)[,colnames(ip), with = F]

    # -- cluster multi-position groups
    if(nrow(pg2) > 0){
      pg2[,clusThis := max(diff(interpRefOrd)) > synBuff,
          by = c("og", "interpRefChr")]
      pg1 <- rbind(
        pg1, subset(pg2, !clusThis)[,colnames(ip), with = F])
      pg2 <- subset(pg2, clusThis)[,colnames(ip), with = F]

      # -- cluster if needed
      if(nrow(pg2)){
        pg2[,clus := dbscan(frNN(cbind(interpRefOrd, interpRefOrd),
                                 eps = synBuff), minPts = 1)$cluster,
            by = c("og","interpRefChr")]
        pg1 <- rbind(
          pg1, pg2[,colnames(ip), with = F])
      }
    }

    # -- again drop clusters without enough support
    out <- data.table(pg1)
    out[,nInChr := uniqueN(ofID), by = c("interpRefChr", "og", "clus")]
    out[,n := uniqueN(ofID), by = "og"]
    out[,rat := nInChr / n]
    out <- subset(out, rat >= propAssignThresh)

    # -- re-level the genomes by user specs
    gord <- c(refGenome, genomeIDs[genomeIDs != refGenome])
    out[,genome := gv[ofID]]
    out[,lev := factor(genome, levels = gord)]

    # -- get median position
    out[,med := median(interpRefOrd), by = c("og", "interpRefChr", "clus")]
    out[,dist2med := abs(interpRefOrd - med)]

    # -- assign whether the placement is the true location, then rank
    out[,isSameChr := genome == refGenome & chr == cv[ofID]]
    setorder(out, -isSameChr, -rat, dist2med,  lev)

    # -- subset to only the best n placement for each ref chrs
    scaf <- subset(out, !duplicated(paste(og, interpRefChr, clus)))
    setorder(scaf, -isSameChr, -rat, dist2med,  lev)
    scaf[,nOnChr := 1:.N, by = c("interpRefChr", "og")]
    scaf <- subset(scaf, nOnChr <= maxPlacementsPerRefChr)

    # -- get pangenome IDs
    setkey(scaf, interpRefChr, interpRefOrd)
    scaf[,pgID := 1:.N]
    scaf <- scaf[,c("ofID", "interpRefChr", "interpRefOrd", "pgID", "og")]

    # -- merge with hits
    tmp <- with(refAnchorHits, data.table(
      ofID = c(repOfID, ofID), altID = c(ofID, repOfID)))
    tmp <- subset(tmp, ofID %in% scaf$ofID)
    tmp <- subset(tmp, !duplicated(tmp))
    out <- merge(scaf, tmp, by = "ofID", allow.cartesian = T)
    out[,`:=`(isRep = ofID == altID, ofID = NULL)]
    out <- subset(out, !duplicated(out))
    # -- rename and format
    setnames(out, c("pgChr", "pgOrd", "pgID", "og", "ofID","isRep"))
    out[,`:=`(genome = gv[ofID], chr = cv[ofID], ord = ov[ofID])]
    out[,lev := factor(genome, levels = gord)]
    setkey(out, pgID, lev)
    out[,lev := NULL]
    return(out)
  }

  ##############################################################################
  # -- internal function to add array members
  add_arrayMembers <- function(pangenomeDt, gff){

    arrayID <- arrayMemID <- NULL

    gff <- data.table(gff)
    pangenomeDt <- data.table(pangenomeDt)

    # -- pull non-array representative genes
    gfnonrep <- subset(gff, !isArrayRep)
    gfrep <- subset(gff, arrayID %in% gfnonrep$arrayID & isArrayRep)
    gfr <- merge(
      with(gfrep, data.table(ofID = ofID, arrayID = arrayID)),
      with(gfnonrep, data.table(arrayMemID = ofID, arrayID = arrayID, genome, chr, ord)),
      by = "arrayID", allow.cartesian = T)
    pgr <- subset(pangenomeDt, ofID %in% gfrep$ofID)
    pgr <- pgr[,c("pgChr", "pgOrd", "pgID", "og", "ofID",
                  "isRep", "isDirectSyn", "isSynOgOnly")]
    pgr <- subset(pgr, !duplicated(pgr))
    gfr <- subset(gfr, ofID %in% pgr$ofID)
    gfr[,arrayID := NULL]
    gfr <- subset(gfr, !duplicated(gfr))
    pgr <- merge(
      pgr, gfr,
      all.y = T,
      by = "ofID",
      allow.cartesian = T)
    pgr[,ofID := arrayMemID]
    pgr <- pgr[,colnames(pangenomeDt), with = F]
    pangenomeDt[,isArrayRep := TRUE]
    pgr[,isArrayRep := FALSE]
    pangenomeDt <- rbind(pangenomeDt, pgr)
    return(pangenomeDt)
  }

  ##############################################################################
  # -- internal function to pull missing direct edges
  pull_missingDirectEdges <- function(pangenomeDt,
                                      synParamsDt,
                                      nCores,
                                      gff){

    repOfID <- ofID <- repOfID <- NULL

    ofIDl <- split(pangenomeDt$ofID, pangenomeDt$genome)
    hitsm <- rbindlist(lapply(names(ofIDl), function(i){
      hitsi <- read_hits4pg(
        synParamsDt = synParamsDt,
        pgRepOfIDs = unique(ofIDl[[i]]),
        refGenome = i,
        nCores = nCores,
        gff = gff)
      tmp <- subset(hitsi, repOfID != ofID)[,c("ofID", "repOfID")]
      ptmp1 <- with(subset(pangenomeDt, ofID %in% tmp$repOfID),
                    data.table(repOfID = ofID, pgID1 = pgID))
      ptmp2 <- with(subset(pangenomeDt, ofID %in% tmp$ofID),
                    data.table(ofID = ofID, pgID2 = pgID))
      tmp <- merge(tmp, ptmp1, by = "repOfID", allow.cartesian = T)
      tmp <- merge(tmp, ptmp2, by = "ofID", allow.cartesian = T)
      tmp <- with(tmp, paste(repOfID, ofID))
      hitso <- subset(hitsi, !paste(repOfID, ofID) %in% tmp & repOfID != ofID)
      return(hitso[,c("repOfID", "ofID")])
    }))
    hitsm <- rbind(hitsm, with(hitsm, data.table(repOfID = ofID, ofID = repOfID)))
    hitsm <- subset(hitsm, !duplicated(hitsm))
    ov <- gff$ord; cv <- gff$chr; gv <- gff$genome
    names(ov) <- names(cv) <- names(gv) <- gff$ofID
    hitsm[,`:=`(genome = gv[ofID], chr = cv[ofID], ord = ov[ofID])]

    tmp <- pangenomeDt[,c("pgChr", "pgOrd","pgID", "og", "ofID")]
    tmp[,`:=`(isRep = FALSE, isDirectSynOG = TRUE, repOfID = ofID, ofID = NULL)]
    hmiss <- subset(hitsm, repOfID %in% pangenomeDt$ofID)
    hmiss <- subset(hmiss, !ofID %in% pangenomeDt$ofID)
    hmiss <- hmiss[,c("genome","repOfID", "ofID", "chr", "ord")]
    tmp <- subset(tmp, repOfID %in% hmiss$repOfID)
    hmiss <- subset(hmiss, !duplicated(hmiss))
    tmp <- subset(tmp, !duplicated(tmp))
    hmiss <- merge(tmp, hmiss, by = "repOfID", allow.cartesian = T)
    hmiss <- hmiss[,colnames(pangenomeDt), with = F]
    pangenomeDt <- rbind(pangenomeDt, hmiss)
    return(pangenomeDt)
  }

  ##############################################################################
  # -- internal function to check for missing entries
  chk_missingEntries <- function(pangenomeDt,
                                 refPgHits){

    ofID <- repOfID <- NULL

    tmp <- pangenomeDt[,c("pgChr", "pgOrd","pgID", "og", "ofID")]
    tmp[,`:=`(isRep = FALSE, isDirectSynOG = TRUE, repOfID = ofID, ofID = NULL)]
    hmiss <- subset(refPgHits, repOfID %in% pangenomeDt$ofID)
    hmiss <- subset(hmiss, !ofID %in% pangenomeDt$ofID)
    hmiss <- hmiss[,c("genome","repOfID", "ofID", "chr", "ord")]
    tmp <- subset(tmp, repOfID %in% hmiss$repOfID)
    hmiss <- subset(hmiss, !duplicated(hmiss))
    tmp <- subset(tmp, !duplicated(tmp))
    hmiss <- merge(tmp, hmiss, by = "repOfID", allow.cartesian = T)
    hmiss <- hmiss[,colnames(pangenomeDt), with = F]
    out <- rbind(pangenomeDt, hmiss)
    return(out)
  }

  ##############################################################################
  # -- internal function to pull indirect syntenic edges
  pull_indirectEdges <- function(pangenomeDt,
                                 synParamsDt,
                                 nCores,
                                 gff){

    isArrayRep <- genome <- genome1 <- genome2 <- repOfID <- ofID <- u <-
      repOfID <- ofID <- pgID <- NULL

    pgRef <- data.table(pangenomeDt)
    ofIDl <- subset(gff, isArrayRep & genome != refGenome)
    ofIDl <- split(ofIDl$ofID, ofIDl$genome)
    hitsm <- rbindlist(lapply(names(ofIDl), function(i){
      hitsi <- read_hits4pg(
        synParamsDt = subset(synParamsDt,
                             genome1 != refGenome & genome2 != refGenome),
        pgRepOfIDs = unique(ofIDl[[i]]),
        refGenome = i,
        nCores = nCores,
        gff = gff)
      tmp <- subset(hitsi, repOfID != ofID)[,c("ofID", "repOfID")]
      ptmp1 <- with(subset(pgRef, ofID %in% tmp$repOfID), data.table(repOfID = ofID, pgID1 = pgID))
      ptmp2 <- with(subset(pgRef, ofID %in% tmp$ofID), data.table(ofID = ofID, pgID2 = pgID))
      tmp <- merge(tmp, ptmp1, by = "repOfID", allow.cartesian = T)
      tmp <- merge(tmp, ptmp2, by = "ofID", allow.cartesian = T)
      tmp <- with(tmp, paste(repOfID, ofID))
      hitsi[,u := paste(repOfID, ofID)]
      hitso <- subset(hitsi, !u %in% tmp & !paste(ofID, repOfID) %in% tmp)
      return(hitso[,c("repOfID", "ofID")])
    }))
    hitsOut <- data.table(hitsm)
    hitso <- subset(hitsm, repOfID != ofID)
    # -- pull out hits that are indirect syn OGs
    hitsm <- rbind(hitsm, with(hitsm, data.table(repOfID = ofID, ofID = repOfID)))
    hitsm <- subset(hitsm, !duplicated(hitsm))
    indHits <- subset(hitsm, ofID %in% pgRef$ofID & ofID != repOfID)
    pgInd <- subset(pgRef, ofID %in% indHits$ofID)
    pgInd <- pgInd[,c("pgChr", "pgOrd", "pgID", "ofID", "og")]
    pgInd <- merge(pgInd, indHits, allow.cartesian = T)
    pgInd <- subset(pgInd, !duplicated(pgInd))

    # -- add to the pangenome
    ov <- gff$ord; cv <- gff$chr; gv <- gff$genome
    names(ov) <- names(cv) <- names(gv) <- gff$ofID
    pgInd[,`:=`(genome = gv[ofID], chr = cv[ofID], ord = ov[ofID], isRep = FALSE,
                isDirectSynOG = FALSE, ofID = repOfID)]
    pgInd <- pgInd[,colnames(pgRef), with = F]
    pgInd <- subset(pgInd, !paste(ofID, pgID) %in% paste(pgRef$ofID, pgRef$pgID))
    pgRef <- rbind(pgRef, pgInd)
    return(list(pg = pgRef, missingHits = hitsOut))
  }

  ##############################################################################
  # -- internal function to pull non-reference syntenic orthogroups
  pull_nonRefSynOgs <- function(pangenomeDt,
                                synParamsDt,
                                nCores,
                                gff){

    ofID <- repOfID <- genome <- refg <- og <- genome2 <- ng <- hasInterp <-
      ni <- isRep <- clus <- interpRefChr <- interpRefOrd <- nInChr <- n <-
      rat <- tmpi <- nPos <- clusThis <- lev <- med <- dist2med <- nOnChr <-
      pgID <- altID <- NULL

    # -- invert and copy hits
    u <- unique(pgInd$ofID)
    ov <- gff$ord; cv <- gff$chr; gv <- gff$genome
    names(ov) <- names(cv) <- names(gv) <- gff$ofID
    tmp <- with(subset(hitsMissing, !ofID %in% u & !repOfID %in% u), data.table(
      repOfID = c(repOfID, ofID), ofID = c(ofID, repOfID)))

    levs <- genomeIDs[genomeIDs!=refGenome]
    tmp[,genome := gv[repOfID]]
    anch <- NULL

    # -- hierarchical orthogroup assignment
    for(i in levs){
      anchi <- subset(tmp, genome == i)
      anchi[,refg := i]
      lefto <- subset(tmp, !ofID %in% anch$ofID & !repOfID %in% anch$ofID)
      anchi[,og := as.numeric(factor(repOfID, levels = unique(repOfID)))]
      anch <- rbind(anch, anchi)
    }
    anch[,og := as.numeric(as.factor(paste(refg, og)))]
    anch[,refg := NULL]

    # -- order by positions, genome and interpolated positions to choose reps
    tmp <- data.table(anch)
    levs <- c(refGenome, genomeIDs[genomeIDs!=refGenome])
    tmp[,genome := factor(gv[repOfID], levels = levs)]
    tmp[,genome2 := factor(gv[ofID], levels = levs)]
    tmp[,ng := uniqueN(genome2), by = "repOfID"]
    tmp <- subset(tmp, !duplicated(tmp))
    tmp[,hasInterp := ofID %in% interp$ofID]
    setorder(tmp, -hasInterp, -ng, genome)
    tmp[,ni := 1:.N, by = "repOfID"]
    tmp[,isRep := ni == 1 & repOfID == ofID]
    reps <- tmp$repOfID[tmp$isRep]
    tmp <- subset(tmp, repOfID %in% reps)
    setkey(tmp, og)

    synHits <- data.table(tmp)

    tmp <- with(tmp, rbind(data.table(ofID = repOfID, og = og),
                           data.table(ofID = ofID, og = og)))
    tmp <- subset(tmp, !duplicated(tmp))
    ip <- merge(tmp, interp, by = "ofID", allow.cartesian = T, all.x = T)

    # -- toss chromosomes that are rare (default < 25%)
    ip[,clus := 1]
    setorder(ip, interpRefChr, interpRefOrd, na.last = T)
    tmp <- subset(ip, !is.na(interpRefChr))
    tmp[,nInChr := uniqueN(ofID), by = c("interpRefChr", "og")]
    tmp[,n := uniqueN(ofID), by = "og"]
    tmp[,rat := nInChr / n]
    tmp <- subset(tmp, rat >= propAssignThresh)[,colnames(ip), with = F]

    # -- separate groups that need to be clustered from those that don't
    tmp[,tmpi := round(interpRefOrd)]
    tmp[,nPos := uniqueN(tmpi), by = c("og", "interpRefChr")]

    pg1 <- subset(tmp, nPos == 1)[,colnames(ip), with = F]
    pg2 <- subset(tmp, nPos > 1)[,colnames(ip), with = F]

    # -- cluster multi-position groups
    if(nrow(pg2) > 0){
      pg2[,clusThis := max(diff(interpRefOrd)) > synBuff,
          by = c("og", "interpRefChr")]
      pg1 <- rbind(
        pg1, subset(pg2, !clusThis)[,colnames(ip), with = F])
      pg2 <- subset(pg2, clusThis)[,colnames(ip), with = F]

      # -- cluster if needed
      if(nrow(pg2)){
        pg2[,clus := dbscan(frNN(cbind(interpRefOrd, interpRefOrd),
                                 eps = synBuff), minPts = 1)$cluster,
            by = c("og","interpRefChr")]
        pg1 <- rbind(
          pg1, pg2[,colnames(ip), with = F])
      }
    }

    # -- again drop clusters without enough support
    out <- data.table(pg1)
    out[,nInChr := uniqueN(ofID), by = c("interpRefChr", "og", "clus")]
    out[,n := uniqueN(ofID), by = "og"]
    out[,rat := nInChr / n]
    out <- subset(out, rat >= propAssignThresh)

    # -- re-level the genomes by user specs
    gord <- c(refGenome, genomeIDs[genomeIDs != refGenome])
    out[,genome := gv[ofID]]
    out[,lev := factor(genome, levels = gord)]

    # -- get median position
    out[,med := median(interpRefOrd), by = c("og", "interpRefChr", "clus")]
    out[,dist2med := abs(interpRefOrd - med)]

    setorder(out, -rat, dist2med,  lev)

    # -- subset to only the best n placement for each ref chrs
    scaf <- subset(out, !duplicated(paste(og, interpRefChr, clus)))
    setorder(scaf, -rat, dist2med,  lev)
    scaf[,nOnChr := 1:.N, by = c("interpRefChr", "og")]
    scaf <- subset(scaf, nOnChr <= maxPlacementsPerRefChr)
    scaf <- scaf[,c("ofID", "interpRefChr", "interpRefOrd")]
    scaf <- subset(scaf, !duplicated(scaf))
    # -- get pangenome IDs
    setkey(scaf, interpRefChr, interpRefOrd)
    scaf[,pgID := 1:.N]

    scaf <- subset(scaf, !duplicated(scaf))

    # -- merge with hits
    tmp <- with(synHits, data.table(
      ofID = c(repOfID, ofID), altID = c(ofID, repOfID), og = og))
    tmp <- subset(tmp, ofID %in% scaf$ofID)
    tmp <- subset(tmp, !duplicated(tmp))
    out <- merge(scaf, tmp, by = "ofID", allow.cartesian = T)
    out[,`:=`(isRep = ofID == altID, ofID = NULL)]
    out <- subset(out, !duplicated(out))
    # -- rename and format
    setnames(out, c("pgChr", "pgOrd", "pgID", "ofID", "og", "isRep"))
    out[,`:=`(genome = gv[ofID], chr = cv[ofID], ord = ov[ofID])]
    out[,lev := factor(genome, levels = gord)]
    setkey(out, pgID, lev)
    out[,lev := NULL]

    pgMis <- out[,colnames(pangenomeDt), with = F]
    pgMis <- subset(pgMis, !paste(ofID, pgID) %in% paste(pangenomeDt$ofID, pangenomeDt$pgID))
    pgOut <- rbind(pangenomeDt, pgMis)
    return(pgOut)
  }

  ##############################################################################
  # -- internal function to add missing syntenic orthogroups
  add_missingSynOgs <- function(pangenomeDt,
                                gff){
    ofID <- isArrayRep <- og <- pgID <- synog <- isRep <- NULL

    gfmiss <- subset(gff, !ofID %in% pangenomeDt$ofID & isArrayRep)
    gfin <- subset(gff, og %in% gfmiss$og & ofID %in% pangenomeDt$ofID)
    if(nrow(gfin) > 0){
      synogv <- gff$og; names(synogv) <- gff$ofID
      pgin <- subset(pangenomeDt,
                     pgID %in% subset(pangenomeDt, ofID %in% gfin$ofID)$pgID)
      pgin[,synog := synogv[ofID]]
      pgin <- pgin[,c("pgChr", "pgOrd", "pgID", "og", "synog")]
      pgin <- subset(pgin, !duplicated(pgin))

      gfin <- with(gfmiss, data.table(
        ofID = ofID, chr = chr, ord = ord,
        genome = genome, isRep = FALSE, synog = og))

      pgin <- merge(
        subset(pgin, !duplicated(pgin)),
        subset(gfin, !duplicated(gfin)),
        by = "synog", allow.cartesian = T, all.y = T)
      pgin <- pgin[,colnames(pangenomeDt), with = F]
      return(rbind(pangenomeDt, pgin))
    }else{
      return(pangenomeDt)
    }
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
  if(gsParam$genomes$ploidy[refGenome] > 1)
    cat(sprintf(
      "*NOTE* RefGenome %s has >1x ploidy - this is fine, but will be slower.\n\tDepending on the size of your run, you may run into memory issues.\n",
      refGenome))

  # -- check that synteny params exist
  synp <- gsParam$params$synteny
  if(!is.data.table(synp))
    stop("Must run set_syntenyParams first!\n")
  synp <- data.table(synp)
  blkSize <- max(synp$blkSize)
  synBuff <- max(synp$selfRegionMask)

  # -- check the synhits exist
  hitsFile <- genome1 <- genome2 <- NULL
  synp[, hitsFile := file.path(
    gsParam$paths$results,
    sprintf("%s_%s_synHits.txt.gz", genome1, genome2))]
  synp <- subset(synp, file.exists(hitsFile))
  if(nrow(synp) < 1)
    stop("can't find synHits accompanying synParams ... has synteny been run?\n")

  # -- get the orthofinder directory if needed
  if(is.na(gsParam$paths$orthogroupsDir))
    gsParam <- find_orthofinderResults(gsParam)
  verbose <- gsParam$params$verbose

  # -- check the gff file
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  if(!file.exists(gffFile))
    stop("can't find the annotated gff-like text file\t\n ... have you run annotate_gff yet?\n")
  gf <- fread(gffFile, showProgress = F, na.strings = c("NA", ""))
  if(!"og" %in% colnames(gf))
    stop("can't find the og column name in the gff-like text file\n\tHas synteny been run yet?\n")

  # -- check that the reference is in the gff
  if(!refGenome %in% gf$genome)
    stop(sprintf("%s (specified refGenome) not in the gff. Available genomes are: \n\t%s\n",
                 refGenome, paste(unique(gf$genome), collapse = ",")))

  # -- check that the blocks look right
  blksFile <- file.path(gsParam$paths$results, "syntenicBlocks.txt.gz")
  if(!file.exists(blksFile))
    stop("can't find the syntenic block coordinates file. Has synteny been run yet?\n")
  blks <- fread(blksFile, showProgress = F, na.strings = c("NA", ""))
  blkID <- lastGene1 <- lastGene2 <- firstGene2 <- firstGene1 <- isSelf <- NULL
  blks <- subset(blks, !grepl("^NA_", blkID))
  blks[,isSelf := firstGene1 == firstGene2 | firstGene1 == lastGene2 | firstGene2 == lastGene1]

  # -- get a list of the chromosomes that are ever in synteny on refGenome
  gen2 <- gen1 <- blkID <- NULL
  genChr <- unlist(lapply(genomeIDs, function(i){
    b1 <- subset(blks, gen1 == i & !grepl("self", blkID) & gen2 != gen1)$chr1
    b2 <- subset(blks, gen2 == i & !grepl("self", blkID) & gen2 != gen1)$chr2
    return(paste(i, unique(c(b1, b2))))
  }))

  ##############################################################################
  # 1. use hits to build a reference scaffold
  if(verbose)
    cat(sprintf("Building reference-anchored scaffold against %s\n", refGenome))

  # -- pull all self hits that are array reps
  genome <- isArrayRep <- genome1 <- genome2 <- NULL
  refSynp <- subset(synp, genome1 == genome2 & genome1 == refGenome)
  if(nrow(refSynp) != 1)
    stop("not exactly 1 refgenome entry in synParams. something is wrong\n")

  # -- convert to pangenome format
  pgReps <- subset(gf, genome == refGenome & isArrayRep)$ofID
  if(verbose)
    cat(sprintf("\tn. ref positions = %s\n", length(pgReps)))

  ##############################################################################
  # 2. Interpolate reference positions
  # -- read in all synOg hits against reference
  if(verbose)
    cat(sprintf("\tReading in hits against %s ... ", refGenome))
  hits <- read_hits4pg(
    synParamsDt = synp,
    pgRepOfIDs = pgReps,
    refGenome = refGenome,
    nCores = nCores,
    gff = gf)

  # -- interpolate reference position
  if(verbose)
    cat(sprintf(
      "found %s\n\tInterpolating positions ... ",
      nrow(hits)))

  blkID <- isAnchor <- NULL
  anch <- subset(hits, !is.na(blkID) & isAnchor)
  interp <- interp_refPos(
    gff = gf,
    refGenome = refGenome,
    anchorHits = anch,
    blkSize = blkSize,
    nCores = nCores, verbose = T)


  ##############################################################################
  # 3. Build ref pangenome
  # -- merge with the reference genes to build pan-genome scaffold
  if(verbose)
    cat("\tForming ref.-anchored db ... ")
  pgRef <- build_refPg(
    refAnchorHits = anch,
    interpolatedPositionDT = interp,
    gff = gf,
    genomeIDs = genomeIDs,
    refGenome = refGenome,
    propAssignThresh = propAssignThresh,
    maxPlacementsPerRefChr = maxPlacementsPerRefChr)
  if(verbose)
    cat(sprintf("found %s genes for %s placements\n",
                uniqueN(pgRef$ofID), uniqueN(pgRef$pgID)))

  ##############################################################################
  # 4. Completing pangenome
  # -- go back into the hits, make sure they are all there
  if(verbose)
    cat("Completing the pan-genome annotation ...\n\tAdding non-anchor entries ... ")
  pgChk <- chk_missingEntries(
    pangenomeDt = pgRef,
    refPgHits = hits)
  ofID <- pgID <- NULL
  pgAddedChk <- subset(
    pgChk, !paste(ofID, pgID) %in% paste(pgRef$ofID, pgRef$pgID))
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgAddedChk$ofID), uniqueN(pgAddedChk$pgID)))

  # -- get all missing reference hits against the reference
  if(verbose)
    cat("\tChecking missing direct ref. syn. OGs ... ")
  pgAdd <- pull_missingDirectEdges(
    pangenomeDt = pgChk,
    synParamsDt = synp,
    nCores = nCores,
    gff = gf)
  ofID <- pgID <- NULL
  pgAddedAdd <- subset(
    pgAdd, !paste(ofID, pgID) %in% paste(pgChk$ofID, pgChk$pgID))
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgAddedAdd$ofID), uniqueN(pgAddedAdd$pgID)))

  # get syn OGs that are part of the graph, but not anchored to the reference
  if(verbose)
    cat("\tAdding indirect syn. OGs ... ")
  tmp <- pull_indirectEdges(
    pangenomeDt = pgAdd,
    synParamsDt = synp,
    nCores = nCores,
    gff = gf)
  pgInd <- data.table(tmp$pg)
  hitsMissing <- data.table(tmp$missingHits)
  ofID <- pgID <- NULL
  pgAddedInd <- subset(pgInd, !paste(ofID, pgID) %in% paste(pgAdd$ofID, pgAdd$pgID))
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgAddedInd$ofID), uniqueN(pgAddedInd$pgID)))

  # -- get synOGs that do not have reference genome in the graph
  if(verbose)
    cat("\tAdding syn. OGs without ref. anchor ... ")
  pgNS <- pull_nonRefSynOgs(
    pangenomeDt = pgInd,
    synParamsDt = synp,
    nCores = nCores,
    gff = gf)
  ofID <- pgID <- NULL
  pgAddedNS <- subset(pgNS, !paste(ofID, pgID) %in% paste(pgInd$ofID, pgInd$pgID))
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgAddedNS$ofID), uniqueN(pgAddedNS$pgID)))

  # -- final check and adding by synOgs
  if(verbose)
    cat("\tAdding missing genes by synOG identity ... ")
  pgOut <- add_missingSynOgs(pangenomeDt = pgNS, gff = gf)
  ofID <- pgID <- NULL
  pgAddedSynog <- subset(pgOut, !paste(ofID, pgID) %in% paste(pgNS$ofID, pgNS$pgID))
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgAddedSynog$ofID), uniqueN(pgAddedSynog$pgID)))

  ##############################################################################
  # 5. Add in missing information
  if(verbose)
    cat("Annotating and formatting pan-genome\n\tAdding non-anchor entries ... ")

  # -- classifiers
  ofID <- pgID <- isDirectSyn <- pgOrd <- u <- isSynOgOnly <- NULL
  pgOut[,u := paste(pgID, ofID)]
  genesInPgRef <- with(pgRef, unique(paste(pgID, ofID)))
  genesInPgAddedChk <- with(pgAddedChk, unique(paste(pgID, ofID)))
  genesInPgAddedAdd <- with(pgAddedAdd, unique(paste(pgID, ofID)))
  genesInPgAddedInd <- with(pgAddedInd, unique(paste(pgID, ofID)))
  genesInPgAddedNS <- with(pgAddedNS, unique(paste(pgID, ofID)))
  genesInPgAddedSynog <- with(pgAddedSynog, unique(paste(pgID, ofID)))
  pgOut[,`:=`(
    isDirectSyn = u %in% c(genesInPgRef, genesInPgAddedChk, genesInPgAddedAdd),
    isSynOgOnly = u %in% genesInPgAddedSynog | is.na(pgOrd),
    u = NULL)]

  # -- array representatives
  isArrayRep <- NULL
  pgAll <- add_arrayMembers(pangenomeDt = pgOut, gff = gf)
  pgMem <- subset(pgAll, !isArrayRep)
  if(verbose)
    cat(sprintf("found %s genes and %s placements\n",
                uniqueN(pgMem$ofID), uniqueN(pgMem$pgID)))

  # -- pull non-syntenic ortho
  ogdir <- gsParam$paths$orthologuesDir
  if(!is.na(ogdir) && dir.exists(ogdir)){
    if(verbose)
      cat("\tAdding non-syn. orthologs ... ")
    pgNSog <- pull_nonSynOrthologs(
      gsParam = gsParam,
      synParamsDt = synp,
      nCores = nCores,
      gff = gf,
      pangenomeDt = pgAll,
      maxNonSynOrthos2keepPerGenome = maxNonSynOrthos2keepPerGenome)
    ofID <- pgID <- NULL
    pgAddedNSog <- subset(
      pgNSog, !paste(ofID, pgID) %in% paste(pgAll$ofID, pgAll$pgID))
    if(verbose)
      cat(sprintf("found %s genes and %s placements\n",
                  uniqueN(pgAddedNSog$ofID), uniqueN(pgAddedNSog$pgID)))
    pgAll <- data.table(pgNSog)
    ofID <- pgID <- isNSOrtho <- NULL
    pgAll[,isNSOrtho := paste(ofID, pgID) %in% paste(pgAddedNSog$ofID, pgAddedNSog$pgID)]
  }else{
    pgAll[,isNSOrtho := FALSE]
    if(verbose)
      cat("\tNo orthologue file available. Will ignore\n")
  }

  ##############################################################################
  # -- output and write
  if(verbose)
    cat(sprintf("\tWriting pangenome to results/%s_pangenomeDB.txt.gz\n", refGenome))
  pgChr <- ofID <- pgOrd <- isRep <- isDirectSyn <- isSynOgOnly <- isArrayRep <-
    isNSOrtho <- genome <- chr <- ord <- u <- NULL
  setorder(
    pgAll, pgChr, pgOrd, -isRep, -isDirectSyn, -isSynOgOnly,
    -isArrayRep, isNSOrtho, genome, chr, ord, na.last = T)
  pgAll[,u := paste(pgChr, pgOrd, pgID)]
  pgAll$u[is.na(pgAll$pgOrd)] <- as.numeric(as.factor(pgAll$og[is.na(pgAll$pgOrd)]))
  pgAll[,pgID := as.numeric(factor(u, levels = unique(u)))]
  pgAll[,u := NULL]
  pgout <- data.table(pgAll)

  # -- give real names
  iv  <- gf$id; names(iv) <- gf$ofID
  pgout[,`:=`(id = iv[ofID])]

  # -- write text file
  pgf <- file.path(
    gsParam$paths$results,
    sprintf("%s_pangenomeDB.txt.gz", refGenome))
  fwrite(pgout, file = pgf, sep = "\t", showProgress = F)

  if(verbose)
    cat("\tReturning wide-format with only syntenic array reps\n")
  # -- flag non-syn orthos
  wh <- which(pgout$isNSOrtho)
  pgout$id[wh] <- paste0(pgout$id[wh], "*")

  # -- flag array reps
  wh <- which(!pgout$isArrayRep & !pgout$isNSOrtho)
  pgout$id[wh] <- paste0(pgout$id[wh], "+")

  # -- reshape to wide format
  pgID <- pgChr <- pgOrd <- genome <- NULL
  pgw <- dcast(
    pgout,
    pgID + pgChr + pgOrd ~ genome,
    value.var = "id",
    fun.aggregate = function(x) list(x))

  # -- order by pg position
  setorder(pgw, pgChr, pgOrd, na.last = T)

  if(verbose)
    cat("\tDone!\n")
  return(pgw)
}

