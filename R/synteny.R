#' @title synteny
#' @description
#' \code{synteny} synteny
#' @name synteny
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param synRad integer of length 1, taken as the euclidean distance from
#' synBuff (* sqrt 2).
#' @param dropSmallNonOGBlks logical of length 1, after splitting overlapping
#' non-duplicated blocks, should remaining blocks that are smaller than blkSize
#' and have no sameOG hits be dropped?
#' @param blkSize integer of length 1, minimum block size
#' @param nGaps integer of length 1, number of gaps (-m) to allow
#' @param topn1  integer of length 1, specifying the top n scoring hits in
#' genome 1. Typically taken from ploidy in gsParam
#' @param topn2 integer of length 1, specifying the top n scoring hits in
#' genome 2. Typically taken from ploidy in gsParam
#' @param tmpDir file.path of length 1, pointing to the tmp directory, taken
#' from gsParam
#' @param onlyOgAnchors logical of length 1, should anchors only be in the same
#' OG?
#' @param MCScanX_hCall file.path to MCScanX_h, taken from gsParam
#' @param outDir file.path of length 1, pointing to the directory to write plots
#' @param gridSize integer of length 1,
#' @param appendName character string of length 1, to append to plot file name
#' @param inBufferRadius integer of length 1,
#' @param dropSmallNonOGBlks logical of length 1,
#' @param maxIter integer of length 1, specifying the number of iterations to
#' allow to split overlapping non-duplicated blocks
#' @param hits data.table containing the blast hits, also stored in /synHits
#' \cr
#' If called, \code{synteny} returns its own arguments.
#'
#' @details info here

#' @title synteny
#' @description
#' \code{synteny} synteny
#' @rdname synteny
#' @import data.table
#' @import R.utils
#' @export
synteny <- function(gsParam){
  ##############################################################################
  # 0. setup
  if(!"synRad" %in% colnames(gsParam$annotBlastMd)){
    strwrap(
      "**NOTE** synteny parameters have not been set, using parameters from
      gsParam. If you want to manually adjust parameters for each combination
      of genomes, run `set_synParam()` separately\n", indent = 0, exdent = 8)
    gsParam <- set_syntenyParams(gsParam, overwrite = T)
  }

  blMd <- data.table(gsParam$annotBlastMd)
  ##############################################################################

  query <- target <- lab <- NULL
  blMd[,lab := align_charLeft(sprintf("%s v. %s:", query, target))]
  # -- loop through the metadata
  blMdOut <- rbindlist(lapply(1:nrow(blMd), function(i){
    x <- blMd[i,]
    hits <- fread(x$annotBlastFile, na.strings = c("", "NA"))
    cat(x$lab)
    ############################################################################
    # 1. intragenomic hits
    if(x$query == x$target){
      hits <- find_selfSyn(
        hits = hits, synRad = x$synRad)

      ########################################################################
      # 2. self hits if ploidy > 1
      if(x$ploidy1 > 1){
        regID <- NULL
        hitsMask <- subset(hits, is.na(regID))
        hitsMask <- find_initialAnchors(
          hits = hitsMask,
          nGaps = x$nGaps,
          blkSize = x$blkSize,
          topn1 = x$ploidy2 - 1,
          topn2 = x$ploidy1 - 1,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchors)

        hitsMask <- find_synRegions(
          hits = hitsMask,
          blkSize = x$blkSize,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          nGaps = x$nGaps,
          synRad = x$synRad)

        hitsMask <- find_synBlks(
          hits = hitsMask,
          dropSmallNonOGBlks = x$onlyOgAnchors,
          inBufferRadius = x$inBufferRadius,
          maxIter = 10,
          blkSize = x$blkSize)
        hits <- rbind(subset(hits, !is.na(regID)), hitsMask)
      }
    }else{
      ##########################################################################
      # 3. intergenomic hits
      hits <- find_initialAnchors(
        hits = hits,
        nGaps = x$nGaps,
        blkSize = x$blkSize,
        topn1 = x$ploidy2,
        topn2 = x$ploidy1,
        MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
        tmpDir = gsParam$paths$tmp,
        onlyOgAnchors = x$onlyOgAnchors)

      hits <- find_synRegions(
        hits = hits,
        blkSize = x$blkSize,
        MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
        tmpDir = gsParam$paths$tmp,
        nGaps = x$nGaps,
        synRad = x$synRad)

      hits <- find_synBlks(
        hits = hits,
        dropSmallNonOGBlks = x$onlyOgAnchors,
        inBufferRadius = x$inBufferRadius,
        maxIter = 10,
        blkSize = x$blkSize)
    }

    ##########################################################################
    # 4. secondary hits
    if(x$nSecondaryHits > 0){
      regID <- NULL
      hitsMask <- subset(hits, is.na(regID))
      hitsMask <- find_initialAnchors(
        hits = hitsMask,
        nGaps = x$nGapsSecond,
        blkSize = x$blkSizeSecond,
        topn1 = x$ploidy2 * x$nSecondaryHits,
        topn2 = x$ploidy1 * x$nSecondaryHits,
        MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
        tmpDir = gsParam$paths$tmp,
        onlyOgAnchors = x$onlyOgAnchors)

      hitsMask <- find_synRegions(
        hits = hitsMask,
        blkSize = x$blkSizeSecond,
        MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
        tmpDir = gsParam$paths$tmp,
        nGaps = x$nGapsSecond,
        synRad = x$synRad)

      hitsMask <- find_synBlks(
        hits = hitsMask,
        dropSmallNonOGBlks = x$onlyOgAnchorsSecond,
        inBufferRadius = x$inBufferRadius,
        maxIter = 10,
        blkSize = x$blkSizeSecond)
      hits <- rbind(subset(hits, !is.na(regID)), hitsMask)
    }
    x[,`:=`(
      nRegionHits = sum(!is.na(hits$regID)),
      nRegions =uniqueN(hits$regID, na.rm = T),
      nAnchorHits = sum(hits$isAnchor),
      nBlks = uniqueN(hits$blkID, na.rm = T),
      nSVs = uniqueN(hits$blkID, na.rm = T) -
        uniqueN(paste(hits$chr1, hits$chr2)[!is.na(hits$blkID)]))]

    nRegionHits <- nRegions <- nAnchorHits <- nBlks <- nSVs <- NULL
    with(x, cat(sprintf(
      " %s hits in %s regions || %s anchors (%s blks, %s SVs)\n",
      nRegionHits, nRegions, nAnchorHits, nBlks, nSVs)))

    ggdotplot_blkRegs(
      hits = hits,
      outDir = gsParam$paths$dotplots,
      appendName = "synHits")

    fwrite(
      hits,
      file = x$annotBlastFile,
      quote = F,
      sep = "\t",
      showProgress = FALSE)
    return(x)
  }))
  gsParam$annotBlastMd <- blMdOut
  return(gsParam)
}


#' @title find_selfSyn
#' @description
#' \code{find_selfSyn} find_selfSyn
#' @rdname synteny
#' @import data.table
#' @export
find_selfSyn <- function(hits, synRad){
  ##############################################################################
  # 1. get the minimum distance between two gene as either ancOrd or ord
  ofID1 <- ofID2 <- chr1 <- chr2 <- ord1 <- ord2 <- synRad <- regID <-
    blkID <- isAnchor <- inBuffer <- NULL
  hits[,isAnchor := ofID1 == ofID2]
  tmp <- subset(hits, chr1 == chr2)

  spl <- split(tmp, by = c("chr1", "chr2"))
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  hits[,inBuffer := paste(ofID1, ofID2) %in% ofInBuff]
  hits[,regID := ifelse(inBuffer, paste(chr1, chr2), NA)]
  hits[,blkID := regID]
  return(hits)
}

#' @title split_ovlBlks
#' @description
#' \code{split_ovlBlks} split_ovlBlks
#' @rdname synteny
#' @import data.table
#' @export
split_ovlBlks <- function(hits,
                          dropSmallNonOGBlks,
                          blkSize){
  isAnchor <- ord1 <- ord2 <- blkID <- chr1 <- start1 <- end1 <- i.blkID <- NULL
  tmp <- subset(hits, isAnchor)
  # dupBlks <- find_dupBlks(tmp)
  # u <-  unique(with(dupBlks,c(paste(blk1, blk2), paste(blk2, blk1))))
  bc <- tmp[,list(start1 = min(ord1), end1 = max(ord1),
                  start2 = min(ord2), end2 = max(ord2)),
            by = c("chr1", "chr2", "blkID")]

  # -- for genome 1
  bc1 <- with(bc, data.table(
    blkID = blkID, chr = chr1, start = start1, end = end1,
    key = c("chr","start","end")))
  fo1 <- subset(foverlaps(bc1, bc1), blkID != i.blkID)
  splh <- split(tmp, by = "blkID")
  nOgHits <- NULL
  if(nrow(fo1) > 0){
    splTmp <- rbindlist(lapply(1:nrow(fo1), function(j){
      h1 <- splh[[fo1$blkID[j]]]
      h2 <- splh[[fo1$i.blkID[j]]]
      isDup <- any(h1$ofID1 %in% h2$ofID1 | h1$ofID2 %in% h2$ofID2)
      if(!isDup){
        blkID <- ofID1 <- ofID2 <- sameOg <- n <- rl <- newblkID <- NULL
        hs <- rbind(h1, h2)
        setkey(hs, ord1, ord2)
        hs[,rl := add_rle(blkID, which = "id")]
        if(dropSmallNonOGBlks){
          hs[,`:=`(nOgHits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOgHits > 0)
          hs[,rl := add_rle(blkID, which = "id")]
        }
        hs[,newblkID := paste(blkID, rl)]
        return(hs)
      }
    }))
    if(nrow(splTmp) > 0){
      blkID <- newblkID <- NULL
      otmp <- subset(tmp, !blkID %in% splTmp$blkID)
      splTmp[,blkID := newblkID]
      tmp <- rbind(otmp, splTmp[,colnames(tmp), with = F])
    }
  }

  ord1 <- ord2 <- blkID <- i.blkID <- sameOg <- ofID1 <- ofID2 <-
  bc <- tmp[,list(start1 = min(ord1), end1 = max(ord1),
                  start2 = min(ord2), end2 = max(ord2)),
            by = c("chr1", "chr2", "blkID")]
  bc2 <- with(bc, data.table(
    blkID = blkID, chr = chr2, start = start2, end = end2,
    key = c("chr","start","end")))
  fo2 <- subset(foverlaps(bc2, bc2), blkID != i.blkID)
  splh <- split(tmp, by = "blkID")

  if(nrow(fo2) > 0){
    splTmp <- rbindlist(lapply(1:nrow(fo2), function(j){
      h1 <- splh[[fo2$blkID[j]]]
      h2 <- splh[[fo2$i.blkID[j]]]
      isDup <- any(h1$ofID1 %in% h2$ofID1 | h1$ofID2 %in% h2$ofID2)
      if(!isDup){
        blkID <- ofID1 <- ofID2 <- sameOg <- n <- rl <- newblkID <- NULL
        hs <- rbind(h1, h2)
        setkey(hs, ord2, ord1)
        hs[,rl := add_rle(blkID, which = "id")]
        if(dropSmallNonOGBlks){
          hs[,`:=`(nOgHits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOgHits > 0)
          hs[,rl := add_rle(blkID, which = "id")]
        }
        hs[,newblkID := paste(blkID, rl)]
        return(hs)
      }
    }))
    if(nrow(splTmp) > 0){
      blkID <- newblkID <- NULL
      otmp <- subset(tmp, !blkID %in% splTmp$blkID)
      splTmp[,blkID := newblkID]
      tmp <- rbind(otmp, splTmp[,colnames(tmp), with = F])
    }
  }

  ofID1 <- ofID2 <- NULL
  blku <- tmp$blkID; names(blku) <- with(tmp, paste(ofID1, ofID2))
  hits[,blkID := blku[paste(ofID1, ofID2)]]
  return(hits)
}


#' @title find_initialAnchors
#' @description
#' \code{find_initialAnchors} find_initialAnchors
#' @rdname synteny
#' @import data.table
#' @importFrom stats complete.cases
#' @export
find_initialAnchors <- function(hits,
                                nGaps,
                                blkSize,
                                topn1,
                                topn2,
                                tmpDir,
                                onlyOgAnchors,
                                MCScanX_hCall){

  noAnchor <- isArrayRep1 <- isArrayRep2 <- sameOg <- isAnchor <- NULL
  tmp <- subset(hits, !noAnchor & isArrayRep1 & isArrayRep2)
  if(onlyOgAnchors)
    tmp <- subset(tmp, sameOg)

  # -- 1.2 subset to top n hits / gene
  bitScore <- sr1 <- sr2 <- ord1 <- ord2 <- ofID1 <- ofID2 <- NULL
  tmp[,sr1 := frank(-bitScore, ties.method = "dense"), by = "ofID1"]
  tmp[,sr2 := frank(-bitScore, ties.method = "dense"), by = "ofID2"]
  tmp <- subset(tmp, sr1 <= topn1 & sr2 <= topn2)

  # -- 1.3 rerank order
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- 1.4 pull collinear hits via mcscanx
  tmp <- tmp[,c("ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2", "bitScore")]
  setnames(tmp, "bitScore", "score")
  tmp <- subset(tmp, complete.cases(tmp))
  anchu <- run_mcscanx(
    nGaps = nGaps,
    tmpDir = tmpDir,
    blkSize = blkSize,
    hits = tmp,
    MCScanX_hCall = MCScanX_hCall)
  hits[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]
  return(hits)
}

#' @title find_synRegions
#' @description
#' \code{find_synRegions} find_synRegions
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom stats complete.cases
#' @export
find_synRegions <- function(hits,
                            blkSize,
                            nGaps,
                            synRad,
                            tmpDir,
                            MCScanX_hCall){

  # -- 1 potential anchors are all hits that are not `noAnchor`
  noAnchor <- anySyn <- ord1 <- ord2 <- ofID1 <- ofID2 <- isAnchor <- NULL
  tmp <- subset(hits, !noAnchor)
  tmp[,anySyn := any(isAnchor),by = c("chr1", "chr2")]
  tmp <- subset(tmp, anySyn)
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- 2 pull all hits within synBuff of an anchor
  spl <- split(tmp, by = c("chr1", "chr2"))
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  tmp <- subset(tmp, paste(ofID1, ofID2) %in% ofInBuff)

  # -- 3 run mcscanx
  tmp <- tmp[,c("ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2", "bitScore")]
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
  setnames(tmp, "bitScore", "score")
  tmp <- subset(tmp, complete.cases(tmp))
  anchu <- run_mcscanx(
    nGaps = nGaps,
    tmpDir = tmpDir,
    blkSize = blkSize,
    hits = tmp,
    MCScanX_hCall = MCScanX_hCall)

  # -- add all potential anchors to hits (these are used for the rest)
  hits[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]

  # -- 4 cluster anchors into regions
  # -- rerank anchors
  tmp <- subset(hits, isAnchor)
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- cluster
  tmp[,regID := dbscan(frNN(
    x = cbind(ord1, ord2),
    eps = synRad),
    minPts = blkSize)$cluster,
    by = c("chr1", "chr2")]

  # -- toss blocks without enough hits
  regID <- chr1 <- chr2 <- ofID1 <- ofID2 <- ord1 <- ord2 <- NULL
  tmp <- subset(tmp, regID > 0)

  # -- rename blocks with uniqueID, add to hits
  tmp[,regID := paste(chr1, chr2, regID)]
  regu <- tmp$regID; names(regu) <- with(tmp, paste(ofID1, ofID2))
  hits[,regID := regu[paste(ofID1, ofID2)]]

  # -- 5 get all hits within region coordinates
  u <- with(subset(hits, !is.na(regID)), unique(paste(chr1, chr2)))
  tmp <- subset(hits, paste(chr1, chr2) %in% u)
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- pull all hits within synBuff of an anchor
  spl <- split(tmp, by = c("chr1", "chr2"))
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  tmp <- subset(tmp, paste(ofID1, ofID2) %in% ofInBuff)

  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- get block coordinates
  blks <- subset(tmp, !is.na(regID))[,list(
    start1 = min(ord1), end1 = max(ord1),
    start2 = min(ord2), end2 = max(ord2)),
    by = c("chr1", "chr2", "regID")]
  u <- with(blks, unique(paste(chr1, chr2)))

  # -- split out hits and block coordinates by block ID
  splh <- split(subset(tmp, paste(chr1, chr2) %in% u), by = c("chr1", "chr2"))
  splb <- split(blks, by = c("chr1", "chr2"))
  splh <- splh[names(splb)]

  # -- 6 finalize hits within each region
  tmpReg <- rbindlist(lapply(1:nrow(blks), function(i){
    y <- blks[i,]
    x <- splh[[with(y, paste(chr1, chr2, sep = "."))]]
    x <- subset(x, ord1 >= y$start1 & ord1 <= y$end1 &
                  ord2 >= y$start2 & ord2 <= y$end2)
    x[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    out <- subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")]
    out[,regID := y$regID]
    return(out)
  }))

  regu <- tmpReg$regID; names(regu) <- with(tmpReg, paste(ofID1, ofID2))
  hits[,regID := NULL]
  hits[,regID := regu[paste(ofID1, ofID2)]]
  return(hits)
}

#' @title find_synBlks
#' @description
#' \code{find_synBlks} find_synBlks
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
find_synBlks <- function(hits,
                         blkSize,
                         inBufferRadius,
                         dropSmallNonOGBlks,
                         maxIter){

  noAnchor <- regID <- isArrayRep1 <- isArrayRep2 <- sr1 <- sr2 <-
    bitScore <- NULL
  tmp <- subset(hits, !is.na(regID) & !noAnchor & isArrayRep1 & isArrayRep2)
  tmp[,sr1 := frank(-bitScore, ties.method = "dense"), by = c("regID","ofID1")]
  tmp[,sr2 := frank(-bitScore, ties.method = "dense"), by = c("regID","ofID2")]
  tmp <- subset(tmp, sr1 == 1 & sr2 == 1)

  # -- 3.2 re-rank genes within each region and split by region
  ord1 <- ord2 <- blkID <- regID <- ofID1 <- ofID2 <- NULL
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense")), by = "regID"]
  spl <- split(tmp, by = "regID")

  # -- for each region ...
  tmpBlk <- rbindlist(lapply(spl, function(x){

    # -- 3.3 dbscan prune with slightly larger buffer
    x[,blkID := dbscan(frNN(
      x = cbind(ord1, ord2),
      eps = blkSize),
      minPts = blkSize)$cluster]
    x <- subset(x, blkID > 0)

    # -- 3.4 re-rank order and re-cluster with finest scale blocks
    x[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
    x[,blkID := dbscan(frNN(
      x = cbind(ord1, ord2),
      eps = ceiling((blkSize)/sqrt(2))),
      minPts = blkSize)$cluster]
    x <- subset(x, blkID > 0)

    # -- 3.5 final clustering
    x[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
    x[,blkID := dbscan(frNN(
      x = cbind(ord1, ord2),
      eps = ceiling((blkSize)/sqrt(2))),
      minPts = blkSize)$cluster]
    x <- subset(x, blkID > 0)
    return(x)
  }))

  u <- with(tmpBlk, paste(regID, blkID))
  names(u) <- with(tmpBlk, paste(ofID1, ofID2))
  hits[,blkID := u[paste(ofID1, ofID2)]]
  hits[,isAnchor := !is.na(blkID)]

  ##############################################################################
  # 4. finalize all the data
  # -- 4.1 split non-duplicated overlapping blocks
  isAnchor <- blkID <- NULL
  tmp <- subset(hits, isAnchor)
  iter = 1; nblks = 1; oblks = 0
  while(iter <= maxIter & nblks > oblks){
    iter = iter + 1
    oblks <- uniqueN(tmp$blkID, na.rm = T)
    tmp <- split_ovlBlks(
      tmp,
      dropSmallNonOGBlks = dropSmallNonOGBlks,
      blkSize = blkSize)
    nblks <- uniqueN(tmp$blkID, na.rm = T)
  }

  hits[,blkID := NULL]
  u <- with(tmp, paste(regID, blkID))
  names(u) <- with(tmp, paste(ofID1, ofID2))
  hits[,blkID := u[paste(ofID1, ofID2)]]
  hits[,isAnchor := !is.na(blkID)]

  # -- 4.2 get inBuffer hits (1/2 of the synbuffer)
  anySyn <- inBuffer <- NULL
  tmp <- subset(hits, !is.na(regID))
  tmp[,anySyn := any(isAnchor),by = c("chr1", "chr2")]
  tmp <- subset(tmp, anySyn)
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense")), by = "regID"]
  spl <- split(tmp, by = "regID")
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = inBufferRadius))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  hits[,inBuffer := paste(ofID1, ofID2) %in% ofInBuff]
  return(hits)
}

#' @title ggdotplot_blkRegs
#' @description
#' \code{ggdotplot_blkRegs} ggdotplot_blkRegs
#' @rdname synteny
#' @import data.table
#' @import ggplot2
#' @importFrom grDevices pdf dev.off rgb
#' @export
ggdotplot_blkRegs <- function(hits,
                              outDir,
                              gridSize = 20,
                              appendName = "synHits"){

  blkCols <- sample(
    gs_colors(20),
    uniqueN(hits$blkID, na.rm = T),
    replace = T)
  regCols <- sample(
    gs_colors(20),
    uniqueN(hits$regID, na.rm = T),
    replace = T)

  tp <- data.table(hits)
  un1 <- uniqueN(tp$ofID1)
  un2 <- uniqueN(tp$ofID2)

  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- NULL
  tp[,ngene1 := uniqueN(ofID1), by = "chr1"]
  tp[,ngene2 := uniqueN(ofID2), by = "chr2"]
  tp <- subset(tp, sameOg & ngene1 > gridSize & ngene2 > gridSize)
  tp[,`:=`(rnd1 = round_toInteger(ord1, gridSize),
           rnd2 = round_toInteger(ord2, gridSize))]

  if(un1 > un2){
    ht <- 12
    wd <- ht * (un1/un2)
  }else{
    wd <- 12
    ht <- wd * (un2/un1)
  }

  dpFile <- file.path(
    outDir, sprintf("%s_vs_%s.%s.pdf",
                    hits$genome1[1], hits$genome2[1], appendName))
  pdf(dpFile, height = ht, width = wd)

  isAnchor <- n <- rnd1 <- rnd2 <- chr1 <- chr2 <- ns1 <- ns2 <- blkID <-
    regID <- NULL
  hcBlk <- subset(tp, isAnchor)[,list(
    n = .N),
    by = c("chr1", "chr2", "rnd1", "rnd2", "blkID")]
  hcBlk[,ns2 := sum(n), by = c("chr2")]
  hcBlk[,ns2 := sum(n), by = c("chr2")]
  hcBlk$n[hcBlk$n > 20] <- 20

  p1 <- ggplot(hcBlk, aes(rnd1, rnd2, col = blkID)) +
    geom_point(size = .01) +
    scale_color_manual(values = blkCols, guide = "none") +
    scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd1), by = 1e3))+
    scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd2), by = 1e3))+
    theme(panel.background = element_rect(fill = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = rgb(1,1,1,.2), size = .2, linetype = 2),
          panel.spacing = unit(.1, "mm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, family = "Helvetica", size = 5),
          strip.text.y.left = element_text(angle = 0, family = "Helvetica", size = 5),
          axis.title = element_text( family = "Helvetica", size = 6),
          plot.title = element_text( family = "Helvetica", size = 7)) +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
    labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome1[1], uniqueN(hits$ofID1)),
         y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome2[1], uniqueN(hits$ofID2)),
         title = sprintf("Syntenic block anchor blast hits, colored by blockID, %s-gene window", gridSize))

  hcBlk <- subset(tp, !is.na(regID))[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2", "regID")]
  hcBlk[,ns2 := sum(n), by = c("chr2")]
  hcBlk[,ns2 := sum(n), by = c("chr2")]
  hcBlk$n[hcBlk$n > 10] <- 10
  p2 <- ggplot(hcBlk, aes(rnd1, rnd2, col = regID, alpha = n)) +
    geom_point(pch = ".") +
    scale_color_manual(values = blkCols, guide = "none") +
    scale_alpha_continuous(guide = "none", range = c(.2, 1)) +
    scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd1), by = 1e3))+
    scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd2), by = 1e3))+
    theme(panel.background = element_rect(fill = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = rgb(1,1,1,.2), size = .2, linetype = 2),
          panel.spacing = unit(.1, "mm"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(angle = 90, family = "Helvetica", size = 5),
          strip.text.y.left = element_text(angle = 0, family = "Helvetica", size = 5),
          axis.title = element_text( family = "Helvetica", size = 6),
          plot.title = element_text( family = "Helvetica", size = 7)) +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
    labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome1[1], uniqueN(hits$ofID1)),
         y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome2[1], uniqueN(hits$ofID2)),
         title = sprintf("All blast hits proximate to syntenic regions, colored by regionID, %s-gene window thresholded at 10 hits", gridSize))
  print(p1)
  print(p2)
  dev.off()
}
