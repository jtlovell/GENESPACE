#' @title Integrate syntenic positions across genomes
#' @description
#' \code{integrate_synteny} Linear interpolation and block coordinate phasing
#' across multiple genomes
#' @name integrate_synteny
#'
#' @param gsParam A list of genespace parameters created by init_genespace.
#'
#' \cr
#' If called, \code{integrate_synteny} returns its own arguments.
#'
#' @details Linear interpolation is conducted by using perfectly collinear
#' syntenic anchors and known fixed positions between any two genomes. The
#' expected syntenic positions of genes between these anchors (without perfect
#' positions) are interpolated from the anchor positions. Syntenic block
#' breakpoints are then calculated for each pair of genomes and also phased
#' against each reference based on the interpolated syntenic location of
#' genes in each block.

#' @title integrate_synteny
#' @description
#' \code{integrate_synteny} integrate_synteny
#' @rdname integrate_synteny
#' @import data.table
#' @export
integrate_synteny <- function(gsParam){

  ##############################################################################
  # -- ad hoc function to add array reps
  add_arrayReps2bed <- function(bed){
    n <- ord <- medOrd <- medDiff <- pepLen <- arrayID <- isRep <-
      isArrayRep <- ofID <- NULL
    bed[,n := .N, by = "arrayID"]
    tmp <- subset(bed, n > 1)
    tmp[,medOrd := as.numeric(median(ord)), by = "arrayID"]
    tmp[,medDiff := as.numeric(abs(medOrd - ord))]
    setorder(tmp, arrayID, medDiff, -pepLen)
    tmp$noAnchor[duplicated(tmp$arrayID)] <- TRUE
    tmp[,isRep := !duplicated(arrayID)]
    bed[,isArrayRep := ofID %in% tmp$ofID[tmp$isRep] | n == 1]
    tmp <- subset(bed, isArrayRep)
    tmp[,ord := frank(ord, ties.method = "dense"), by = "genome"]
    di <- tmp$ord; names(di) <- tmp$ofID
    return(bed)
  }

  ##############################################################################
  # -- ad hoc function to add arrays to bed
  add_array2bed <- function(bed, maxPlaces, synBuff, maxIter = 10){

    # -- we can make arrays for everything except that we want to exclude the
    # arrays that are huge and problematic
    genome <- chr <- id <- arrayID <- nOGPlaces <- tord <- ofID <- NULL
    tmp <- subset(bed, nOGPlaces <= maxPlaces)
    tmp[,tord := as.numeric(ord)]

    # -- set up the iteration
    cnt <- 1
    diffn <- 1
    tmp[,arrayID := sprintf(
      "tmp%s", as.integer(as.factor(paste(genome, chr, id))))]
    while(cnt <= maxIter && diffn > 0){
      # -- for each iteration, calculate clusters by the size of jump between
      # genes larger than the synBuffer
      cnt <- cnt + 1
      initn <- uniqueN(tmp$arrayID)
      genome <- chr <- og <- ord <- jumpLeft <- clus <- n <- NULL
      setkey(tmp, genome, chr, og, tord)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      tmp[,jumpLeft := c(synBuff + 1, diff(tord)), by = c("genome", "chr", "og")]
      tmp[,clus := as.integer(jumpLeft > synBuff), by = c("genome", "chr")]
      tmp[,clus := cumsum(clus), by = c("genome", "chr", "og")]
      tmp[,`:=`(
        arrayID = sprintf("tmp%s",
                          as.integer(as.factor(paste(genome, chr, og, clus)))),
        jumpLeft = NULL, clus = NULL)]

      # -- get the new order of genes based on array ID
      tmp[,n := .N, by = "arrayID"]
      tmp1 <- subset(tmp, n == 1)
      tmp2 <- subset(tmp, n > 1)
      tmp2[,tord := as.numeric(tord)]
      tmp2[,tord := mean(as.numeric(tord)), by = "arrayID"]
      tmp <- rbind(tmp1, tmp2)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      newn <- uniqueN(tmp$arrayID)
      diffn <- initn - newn
    }

    # -- relabel
    ogID <- NULL
    lab <- gsub(" ", "0",
                align_charRight(
                  as.numeric(factor(tmp$arrayID,
                                    levels = unique(tmp$arrayID)))))
    arrv <- sprintf("Arr%s", lab); names(arrv) <- tmp$ofID
    bed[,arrayID := arrv[ofID]]

    tmp <- subset(bed, is.na(arrayID))
    wh <- with(tmp,
               as.numeric(as.factor(paste(genome, chr, og))))
    lab <- gsub(" ", "0", align_charRight(wh))
    wh <- which(is.na(bed$arrayID))
    bed$arrayID[wh] <- sprintf("NoArr%s", lab)
    return(bed)
  }

  ##############################################################################
  # -- ad hoc function to add array info to bed
  add_arrayInfo2bed <- function(gsParam){

    md <- data.table(gsParam$annotBlastMd)
    bedFile <- file.path(gsParam$paths$results, "combBed.txt")
    bed <- read_combBed(bedFile)

    # -- 1.2 find the arrays
    bed <- add_array2bed(
      bed = bed,
      maxPlaces = 4,
      synBuff = gsParam$params$synBuff,
      maxIter = 10)

    # -- 1.3 choose the array representative genes
    bed <- add_arrayReps2bed(bed)
    write_combBed(x = bed, filepath = bedFile)
    return(gsParam)
  }

  ##############################################################################
  # -- ad hoc function to aggregate syntenic positions
  aggregate_synpos <- function(md, bed){

    query <- target <- genome <- isArrayRep <- chr  <- ord <- interpChr <- NULL

    gids <- with(md, unique(c(query, target)))
    interpChrs <- rbindlist(lapply(gids, function(i){

      # -- 2.1 get the input data
      mds <- subset(md, query == i | target == i)
      beds <- subset(bed, genome == i)

      # -- 2.2 read in the interpolated positions (stored in the tmp dir)
      interps <- rbindlist(lapply(which(!is.na(mds$interPosFile)), function(j)
        read_intSynPos(mds$interPosFile[j])))


      # -- 2.3 match the interpolated and bed columns
      beda <- subset(beds, isArrayRep)[,c("genome", "ofID", "chr", "ord", "og")]
      beda[,`:=`(interpGenome = genome, interpChr = chr,
                 interpOrd = ord, isAnchor = TRUE)]

      inta <- subset(interps, !is.na(interpChr) & genome == i)
      inta <- inta[,colnames(beda), with = F]

      # -- 2.4 combine and write to file in the /pangenome directory
      intb <- rbind(beda, inta)
      setkey(intb, ord)

      spFile <- file.path(gsParam$paths$pangenome, sprintf(
        "%sintegratedSynPos.txt", i))
      fwrite(intb, spFile, sep = "\t", quote = FALSE)
      return(intb[,c("ofID", "interpGenome", "interpChr", "interpOrd", "isAnchor")])
    }))
    return(interpChrs)
  }

  ##############################################################################
  # 1. add in array information into bed file, using the new OGs
  lab <- query <- target <- interPosFile <- isAnchor <- lgBlkID <- blkID <-
    interpGenome <- interpChr <- NULL

  gsParam <- add_arrayInfo2bed(gsParam)

  md <- data.table(gsParam$annotBlastMd)
  md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]
  md[,interPosFile := NA]

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))

  ##############################################################################
  # 2. for each pair of genomes, do linear interpolation, save output
  cat("\t##############\n\tLinear interpolation of syntenic positions ... \n")
  gsParam <- interp_synPos(gsParam)
  md <- data.table(gsParam$annotBlastMd)
  md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]

  ##############################################################################
  # 3. Aggregate interpolated positions across all hits to each genome
  cat("\t##############\n\tAggregating positions across genomes ...")
  interpChrs <- aggregate_synpos(bed = bed, md = md)

  # 3. For each genome, calculate block coordinates
  cat("Done!\n\t##############\n\tSplitting syntenic block coordinates by ref. chr. ... ")

  # -- 3.1 organize the interpolated positions for a merge
  # spFile <- file.path(gsParam$paths$pangenome, sprintf(
  #   "%sintegratedSynPos.txt", i))
  interp1 <- with(interpChrs, data.table(
    ofID1 = ofID, interpGenome = interpGenome, interpChr = interpChr))
  interp2 <- with(interpChrs, data.table(
    ofID2 = ofID, interpGenome = interpGenome, interpChr = interpChr))

  # -- 3.2 merge interpolated chrosomes with block IDs
  blkComb <- lapply(1:nrow(md), function(j){
    # -- read in the blasts
    g1 <- md$query[j]
    g2 <- md$target[j]
    hits <- subset(fread(
      md$annotBlastFile[j],
      na.strings = c("", "NA"),
      select = c("ofID1", "chr1", "start1", "end1", "ord1", "genome1",
                 "ofID2", "chr2", "start2", "end2", "ord2", "genome2",
                 "isAnchor", "lgBlkID"),
      showProgress = FALSE),
      isAnchor & !is.na(lgBlkID))

    # -- merge with interp via genome1
    tmp <- merge(interp1, hits, by = "ofID1", allow.cartesian = T)
    tmp[,blkID := sprintf("%sXXXgrpbyXXX%sXXXgrpbyXXX%s",
                          gsub("_", "", lgBlkID), interpGenome, interpChr)]
    tmpb1 <-  calc_blkCoords(tmp, mirror = T)
    tmpb1[,c("blkID", "refGenome", "refChr") := tstrsplit(blkID, "XXXgrpbyXXX")]

    # -- merge with interp via genome2
    tmp <- merge(interp2, hits, by = "ofID2", allow.cartesian = T)
    tmp[,blkID := sprintf("%sXXXgrpbyXXX%sXXXgrpbyXXX%s",
                          gsub("_", "", lgBlkID), interpGenome, interpChr)]
    tmpb2 <-  calc_blkCoords(tmp, mirror = T)
    tmpb2[,c("blkID", "refGenome", "refChr") := tstrsplit(blkID, "XXXgrpbyXXX")]
    hits[,blkID := lgBlkID]
    out <- list(rawBlks = calc_blkCoords(hits),
                phasedBlks = rbind(tmpb1, tmpb2))
    return(out)
  })
  rawBlks <- rbindlist(lapply(blkComb, function(x) x$rawBlks))
  rawBlks <- subset(rawBlks, !duplicated(rawBlks))
  blkComb <- rbindlist(lapply(blkComb, function(x) x$phasedBlks))
  blkComb <- subset(blkComb, !duplicated(blkComb))
  fwrite(
    rawBlks, file = file.path(gsParam$paths$results, "blkCoords.txt"),
    sep = "\t", quote = FALSE, showProgress = FALSE)
  fwrite(
    blkComb, file = file.path(gsParam$paths$riparian, "refPhasedBlkCoords.txt"),
    sep = "\t", quote = FALSE, showProgress = FALSE)
  cat("Done!\n")
  return(gsParam)
}

#' @title Linear interpolation of syntenic positions
#' @description
#' \code{interp_synPos} Find the likely position of a gene based on the syntenic
#' position of its neighbors
#' @rdname integrate_synteny
#' @import data.table
#' @export
interp_synPos <- function(gsParam){

  ##############################################################################
  # -- add hoc function to ensure that syntenic anchors are all collinear
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
  # -- add hoc function to interpolate syntenic hits position
  interp_hitsPos <- function(hits, bed, blkSize = 5, md){

    genome <- isArrayRep <- chr <- start <- end <- sameOg <- isAnchor <-
      blkID <- useAsAnch <- ord1 <- ord2 <- anySelf <- ofID1 <- ofID2 <- chr1 <-
      start1 <- end1 <- chr2 <- start2 <- end2 <- interpRefOrd <- NULL
    # 1. Get the data set up
    # -- 1.1 subset the bed to each genome, key for merge
    bed1 <- subset(bed, genome == hits$genome1[1] & isArrayRep)
    bed2 <- subset(bed, genome == hits$genome2[1] & isArrayRep)
    setkey(bed1, chr, start, end)
    setkey(bed2, chr, start, end)

    # -- 1.2 subset to potential anchor hits
    blkSize2 <- ceiling(blkSize/2)
    anch <- subset(hits, sameOg & isAnchor & !is.na(blkID))

    # -- 1.3 determine fully collinear hits
    anch[,useAsAnch := check_synAnchorPos(
      ord1 = ord1, ord2 = ord2, blkSize = blkSize), by = "blkID"]
    anch <- subset(anch, useAsAnch)
    anch[,anySelf := any(ofID1 == ofID2), by = "blkID"]
    anch <- subset(anch, !anySelf)
    if(nrow(anch) < blkSize){
      cat("no non-self blocks found, skipping interpolation\n")
      return(NULL)
    }else{
      splAnch <- split(anch[,c("ofID1", "ofID2", "blkID")], by = "blkID")

      # -- 1.4 get anchor block coordinates
      anchCoord1 <- anch[,list(chr = chr1[1], start = min(start1), end = max(end1)),
                         by = "blkID"]
      setkey(anchCoord1, chr, start, end)
      anchCoord2 <- anch[,list(chr = chr2[2], start = min(start2), end = max(end2)),
                         by = "blkID"]
      setkey(anchCoord2, chr, start, end)

      # -- 1.5 f overlap join
      bedCoord1 <- split(foverlaps(anchCoord1, bed1), by = "blkID")
      bedCoord2 <- split(foverlaps(anchCoord2, bed2), by = "blkID")
      blkIDs <- intersect(names(bedCoord1), names(bedCoord2))

      # 2. Interpolate by block
      interpd <- rbindlist(lapply(blkIDs, function(i){

        # -- 2.1 merge anchors with all genes in the block
        gen1 <- bedCoord1[[i]]$genome[1]
        gen2 <- bedCoord2[[i]]$genome[1]
        chromo1 <- bedCoord1[[i]]$chr[1]
        chromo2 <- bedCoord2[[i]]$chr[1]
        g1 <- with(bedCoord1[[i]], data.table(ofID1 = ofID, chr1 = chr, ord1 = ord))
        g2 <- with(bedCoord2[[i]], data.table(ofID2 = ofID, chr2 = chr, ord2 = ord))
        g12 <- merge(g1, merge(
          g2, splAnch[[i]], by = "ofID2", all = T), by = "ofID1", all = T)

        # -- 2.2 Do genome 1 (merge, strip leading / trailing nas, interp ofID2)
        tmp1 <- subset(g12, !is.na(ofID1))
        setkey(tmp1, ord1)
        tmp1 <- subset(tmp1, !flag_boundingNAs(ofID2))
        tmp1[,interpRefOrd := interp_linear(x = ord2, y = ord1)]

        # -- 2.3 Do genome 2 (merge, strip leading / trailing nas, interp ofID1)
        tmp2 <- subset(g12, !is.na(ofID2))
        setkey(tmp2, ord2)
        tmp2 <- subset(tmp2, !flag_boundingNAs(ofID1))
        tmp2[,interpRefOrd := interp_linear(x = ord1, y = ord2)]

        o1 <- with(tmp1, data.table(
          ofID = ofID1, interpChr = chromo2,
          interpGenome = gen2, interpOrd = interpRefOrd, isAnchor = !is.na(ofID2)))
        o2 <- with(tmp2, data.table(
          ofID = ofID2, interpChr = chromo1,
          interpGenome = gen1, interpOrd = interpRefOrd, isAnchor = !is.na(ofID1)))
        out <- rbind(o1, o2)
        return(out)
      }))

      interpout <- merge(bed, interpd, by = "ofID", all = T, allow.cartesian = T)
      interpout <- subset(interpout, !duplicated(interpout))
      return(interpout)
    }
  }

  isArrayRep <- lab <- query <- target <- interPosFile <- ofID1 <- ofID2 <-
    genome <- nPlace <- interpOrd <- NULL

  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))
  bedrep <- subset(bed, isArrayRep)

  md <- data.table(gsParam$annotBlastMd)
  md[,lab := align_charLeft(sprintf("%s v. %s: ", query, target))]
  md[,interPosFile := NA]

  for(i in 1:nrow(md)){
    cat("\t...", md$lab[i])
    g1 <- md$query[i]
    g2 <- md$target[i]
    hits <- read_synHits(md$annotBlastFile[i])
    synPos <- interp_hitsPos(
      hits = subset(hits, ofID1 %in% bedrep$ofID & ofID2 %in% bedrep$ofID),
      bed = subset(bedrep, genome %in% c(g1, g2)),
      md = md)
    if(!is.null(synPos)){
      synPos[,nPlace := uniqueN(interpOrd, na.rm = T),
             by = c("ofID", "genome")]
      synPos$nPlace[synPos$nPlace >= 2] <- "2+"
      tab <- synPos[,list(n = .N), by = c("genome", "nPlace")]
      with(tab, cat(sprintf(
        "1x = %s/%s || 2+x = %s/%s || 0x = %s/%s\n",
        n[genome == g1 & nPlace == "1"],  n[genome == g2 & nPlace == "1"],
        n[genome == g1 & nPlace == "2+"],  n[genome == g2 & nPlace == "2+"],
        n[genome == g1 & nPlace == "0"],  n[genome == g2 & nPlace == "0"])))
      spFile <- file.path(gsParam$paths$tmp, sprintf(
        "%s_vs_%s.interpSynPos.txt", g1, g2))
      md$interPosFile[i] <- spFile
      fwrite(synPos, file = spFile, sep = "\t", quote = FALSE)
    }
  }
  gsParam$annotBlastMd <- md
  return(gsParam)
}
