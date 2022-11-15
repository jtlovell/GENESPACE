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
#' @param appendName character string of length 1, to append to plot file name
#' @param blkRadius integer of length 1, specifying the search area to cluster
#' hits into blocks
#' @param dropSmallNonOGBlks logical of length 1,
#' @param hits data.table containing the blast hits, also stored in /synHits
#' @param dotsPerIn integer, number of dots to plot per inch.
#' @param plotSize numeric, length in inches of the shortest plot dimension
#' @param minGenes2plot integer, minimum unique orthogroups on a chromosome
#' for it to be plotted
#' @param verbose logical, should updates be printed to the console?
#'
#' \cr
#' If called, \code{synteny} returns its own arguments.
#'
#' @details Details coming soon.

#' @title Flag syntenic hits
#' @description
#' \code{synteny} The main engine to call syntenic blocks and regions
#' @rdname synteny
#' @import data.table
#' @import R.utils
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
synteny <- function(gsParam, verbose = TRUE, overwrite = TRUE){
  ##############################################################################
  # 1. setup
  # -- 1.1 get env vars set up
  query <- target <- lab <- nRegionHits <- nRegions <- nAnchorHits <- nBlks <-
    nSVs <- selfOnly <- queryPloidy <- targetPloidy <- nGlobOgHits <- synHits <-
    nTotalHits <- chunk <- inBuffer <- NULL

  nCores <- gsParam$params$nCores

  # -- 1.2 check that files are all ok
  if(!"synteny" %in% names(gsParam))
    stop("must run set_syntenyParams prior to synteny")

  if(!all(file.exists(gsParam$synteny$blast$synHits)))
    stop("some annotated blast files dont exist, run annotate_blast() first\n")

  # -- 1.1 split the metadata into chunks
  blMd <- data.table(gsParam$synteny$blast)
  blMd[,lab := align_charLeft(sprintf("%s v. %s:", query, target))]
  blMd[,selfOnly := query == target & queryPloidy == 1 & targetPloidy == 1]

  if(!"nGlobOgHits" %in% colnames(blMd))
    blMd[,nGlobOgHits := file.size(synHits)]

  if(!"nTotalHits" %in% colnames(blMd))
    blMd[,nTotalHits := file.size(synHits)]

  setorder(blMd, selfOnly, -nGlobOgHits, -nTotalHits)
  blMd[,chunk := rep(1:.N, each = nCores)[1:.N]]
  synMdSpl <- split(blMd, by = "chunk")

  # -- 1.2 Check and see if all files are there and have good blkID calls
  if(!overwrite){
    tmpFiles <- blMd$synHits
    if(all(file.exists(tmpFiles))){
      chkd <- sapply(tmpFiles, function(x){
        chk <- !is.na(tryCatch(
          read_synHits(x, nrows = 2)$ofID1[1],
          error = function(err) NA))
        if(chk){
          chk <- sum(complete.cases(fread(
            x, select = c("blkID", "isAnchor", "inBuffer"),
            na.strings = c("", "NA")))) >= gsParam$params$blkSize
        }
        return(chk)
      })
      if(all(chkd))
        overwrite <- FALSE
    }else{
      overwrite <- TRUE
    }
  }

  ##############################################################################
  # -- 2. loop through each chunk
  blMdOut <- rbindlist(lapply(1:length(synMdSpl), function(chnki){
    chnk <- data.table(synMdSpl[[chnki]])

    if(nCores > 1)
      cat(sprintf(
        "\t# Chunk %s / %s (%s) ... \n",
        chnki, max(blMd$chunk), format(Sys.time(), "%X")))

    ############################################################################
    # -- loop through each row in each chunk
    outChnk <- rbindlist(mclapply(1:nrow(chnk), mc.cores = nCores, function(i){

      # -- 2.1 read in the metadata and hits
      x <- chnk[i,]
      hits <- read_synHits(x$synHits)

      ############################################################################
      # 1. intragenomic hits
      if(overwrite){
        if(x$query == x$target){
          hits <- find_selfSyn(
            hits = hits, synRad = x$synRad)

          ########################################################################
          # 2. self hits if ploidy > 1
          if(x$queryPloidy > 1){
            hitsMask <- subset(hits, !inBuffer)
            hitsSelf <- subset(hits, inBuffer)
            hitsMask <- synteny_engine(
              hits = hitsMask,
              nGaps = x$nGaps,
              blkRadius = x$blkRadius,
              blkSize = x$blkSize,
              synRad = x$synRad,
              topn1 = x$targetPloidy - 1,
              topn2 = x$queryPloidy - 1,
              MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
              tmpDir = gsParam$paths$tmp,
              onlyOgAnchors = x$onlyOgAnchors)
            hits <- rbind(hitsMask, hitsSelf)
          }
        }else{
          ##########################################################################
          # 3. intergenomic hits
          hits <- synteny_engine(
            hits = hits,
            nGaps = x$nGaps,
            blkSize = x$blkSize,
            blkRadius = x$blkRadius,
            topn1 = x$targetPloidy,
            topn2 = x$queryPloidy,
            synRad = x$synRad,
            MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
            tmpDir = gsParam$paths$tmp,
            onlyOgAnchors = x$onlyOgAnchors)
        }

        ##########################################################################
        # 4. secondary hits
        if(x$nSecondaryHits > 0){
          hitsMask <- subset(hits, !inBuffer)
          hitsPrim <- subset(hits, inBuffer)
          hitsMask <- synteny_engine(
            hits = hitsMask,
            nGaps = x$nGapsSecond,
            blkSize = x$blkSizeSecond,
            blkRadius = x$blkRadiusSecond,
            topn1 = x$targetPloidy,
            topn2 = x$queryPloidy,
            synRad = x$synRad,
            MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
            tmpDir = gsParam$paths$tmp,
            onlyOgAnchors = x$onlyOgAnchorsSecond)
          hits <- rbind(hitsMask, hitsPrim)
        }
      }

      out <- with(hits, data.table(
        lab = chnk$lab[i],
        nAnchorHits = sum(isAnchor),
        nBufferHits = sum(inBuffer),
        nBlks = uniqueN(blkID, na.rm = T),
        nSVs = uniqueN(blkID, na.rm = T) -
          uniqueN(paste(chr1, chr2)[!is.na(blkID)])))
      out <- data.table(x, out)
      out[,`:=`(nTotalHits = nrow(hits),
                nGlobOgHits = sum(hits$sameOg))]
      write_synBlast(x = hits, filepath = x$synHits)
      if(!overwrite)
        ggdotplot(hits = hits, outDir = gsParam$paths$dotplots)
      return(out)
    }))

    with(outChnk, cat(sprintf(
      "\t...%s %s hits (%s anchors) in %s blocks (%s SVs)\n",
      lab, nBufferHits, nAnchorHits, nBlks, nSVs)))

    return(outChnk)
  }))
  gsParam$annotBlastMd <- blMdOut
  return(gsParam)
}


#' @title Flag synteny in self-hits
#' @description
#' \code{find_selfSyn} Self hits where genes are array reps are the anchors.
#' Also pulls all inbuffer hits around the anchors
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
find_selfSyn <- function(hits, synRad){
  ##############################################################################
  # 1. get the minimum distance between two gene as either ancOrd or ord
  ofID1 <- ofID2 <- chr1 <- chr2 <- ord1 <- ord2 <- regID <- genome1 <- n <-
    blkID <- isAnchor <- inBuffer <- lgBlkID <- genome2 <- NULL
  hits[,isAnchor := ofID1 == ofID2]
  tmp <- subset(hits, chr1 == chr2)
  chrs2search <- tmp[,list(n = diff(range(ord1))), by = "chr1"]
  chrNosearch <- subset(chrs2search, n <= synRad)
  chrs2search <- subset(chrs2search, n > synRad)

  spl <- split(subset(tmp, chr1 %in% chrs2search$chr1), by = "chr1")
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  hits[,inBuffer := paste(ofID1, ofID2) %in% ofInBuff |
         (chr1 == chr2 & chr1 %in% chrNosearch$chr1)]
  hits[,blkID := ifelse(inBuffer, sprintf(
    "%s %s selfBlk %s", genome1, genome2, chr1), NA)]
  return(hits)
}

#' @title Split overlapping non-duplicated blocks
#' @description
#' \code{split_ovlBlks} Where syntenic blocks overlap, converts overlapping
#' regions into RLE and splits accordingly.
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
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
        if(nrow(hs) > blkSize*2){
          hs[,rl := add_rle(blkID, which = "id")]
        }else{
          hs[,rl := 1]
        }

        if(dropSmallNonOGBlks){
          hs[,`:=`(nOgHits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOgHits > 0)
          if(nrow(hs) > blkSize*2){
            hs[,rl := add_rle(blkID, which = "id")]
          }else{
            hs[,rl := 1]
          }
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
        if(nrow(hs) > blkSize*2){
          hs[,rl := add_rle(blkID, which = "id")]
        }else{
          hs[,rl := 1]
        }
        if(dropSmallNonOGBlks){
          hs[,`:=`(nOgHits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOgHits > 0)
          if(nrow(hs) > blkSize*2){
            hs[,rl := add_rle(blkID, which = "id")]
          }else{
            hs[,rl := 1]
          }
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

#' @title Find syntenic regions
#' @description
#' \code{synteny_engine} The main engine for synteny discovery. Takes a "hits"
#' data.table and classifies block IDs, whether a gene is an anchor in the block
#' and whether it is within a syntenic buffer of the block anchor hits.
#' @rdname synteny
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
synteny_engine <- function(hits,
                           onlyOgAnchors,
                           blkSize,
                           blkRadius,
                           tmpDir,
                           topn1,
                           topn2,
                           nGaps,
                           synRad,
                           MCScanX_hCall){


  noAnchor <- isArrayRep1 <- isArrayRep2 <- sameOg <- sr1 <- bitScore <- sr2 <-
    ord1 <- ord2 <- isAnchor <- ofID1 <- ofID2 <- chr1 <- chr2 <- blkID <-
    nAnchorHits <- genome1 <- genome2 <- NULL

  ################################################################################
  # 1. Find initial syntenic anchors
  tmp <- subset(hits, !noAnchor & isArrayRep1 & isArrayRep2)
  if(onlyOgAnchors)
    tmp <- subset(tmp, sameOg)

  # -- 1.2 subset to top n hits / gene
  tmp[,sr1 := frank(-round(bitScore, -1), ties.method = "dense"), by = "ofID1"]
  tmp[,sr2 := frank(-round(bitScore, -1), ties.method = "dense"), by = "ofID2"]
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

  ################################################################################
  # 2. pull all hits within synBuff of an anchor
  # -- 2.1 subset to hits that could be in blks
  tmp <- subset(hits, !noAnchor)
  if(onlyOgAnchors)
    tmp <- subset(tmp, sameOg)

  # -- 2.2 subset to hits on chrs with anchors
  anch <- subset(hits, isAnchor)
  anchu <- with(anch, paste(chr1, chr2))
  tmp <- subset(tmp, paste(chr1, chr2) %in% anchu)

  # -- 2.3 re-rank potential anchors
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]

  # -- 2.4 get all hits within the buffer
  spl <- split(tmp, by = c("chr1", "chr2"))
  ofInBuff <- with(rbindlist(lapply(spl, function(x){
    nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
    wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
    return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
  })), paste(ofID1, ofID2))
  tmp <- subset(tmp, paste(ofID1, ofID2) %in% ofInBuff)

  ################################################################################
  # 3. cluster hits into blocks
  # -- 3.1 re-rank and re-run mcscanx
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

  # -- 3.2 re-rank and cluster into regions
  tmp <- subset(hits, paste(ofID1, ofID2) %in% names(anchu))
  tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
  tmp[,blkID := dbscan(frNN(
    x = cbind(ord1, ord2),
    eps = blkRadius),
    minPts = blkSize)$cluster,
    by = c("chr1", "chr2")]

  # -- 3.3 split overlapping blocks
  tmp <- split_ovlBlks(
    tmp,
    dropSmallNonOGBlks = onlyOgAnchors,
    blkSize = blkSize)
  anchu <- tmp$blkID; names(anchu) <- with(tmp, paste(ofID1, ofID2))

  hits[,blkID := anchu[paste(ofID1, ofID2)]]
  hits[,isAnchor := !is.na(blkID)]

  ################################################################################
  # 4. Get all hits in each syntenic block

  # -- 4.1 get hits on chrs with anchors
  chru <- with(subset(hits, isAnchor), unique(paste(chr1, chr2)))
  tmp <- subset(hits, paste(chr1, chr2) %in% chru)

  # -- 4.2 get block coordinates
  blks <- subset(tmp, isAnchor)[,list(
    start1 = min(ord1), end1 = max(ord1),
    start2 = min(ord2), end2 = max(ord2),
    nAnchorHits = .N),
    by = c("chr1", "chr2", "blkID")]
  blks <- subset(blks, complete.cases(blks) & nAnchorHits > blkSize)
  setorder(blks, -nAnchorHits)

  # -- 4.3 split out hits by chr
  tmpOut <- data.table(tmp[0,])
  for(i in 1:nrow(blks)){
    y <- blks[i,]
    x <- subset(tmp, chr1 == y$chr1 & chr2 == y$chr2)
    x[, isAnchor := blkID == y$blkID]

    # -- only hits within coordinates of block
    x <- subset(x, ord1 >= y$start1 & ord1 <= y$end1 &
                  ord2 >= y$start2 & ord2 <= y$end2)
    if(nrow(x) > blkSize){
      # -- only hits within synBuff of anchor
      nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
      wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
      out <- subset(x, 1:nrow(x) %in% wh)
      out[,`:=`(isAnchor = !is.na(blkID),
                blkID = y$blkID,
                inBuffer = TRUE)]
      u <- with(out, paste(ofID1, ofID2))
      tmp <- subset(tmp, !paste(ofID1, ofID2) %in% u)
      tmpOut <- rbind(tmpOut, out)
    }
  }
  tmpOut[,blkID := paste(genome1, genome2, as.numeric(as.factor(paste(chr1, chr2, blkID))))]
  bu <- with(tmpOut, paste(ofID1, ofID2))
  au <- with(subset(tmpOut, isAnchor), paste(ofID1, ofID2))
  bid <- tmpOut$blkID; names(bid) <- with(tmpOut, paste(ofID1, ofID2))
  hits[,`:=`(blkID = bid[paste(ofID1, ofID2)],
             inBuffer = paste(ofID1, ofID2) %in% bu,
             isAnchor = paste(ofID1, ofID2) %in% au)]

  return(hits)
}

#' @title make dotplots of syntenic hits
#' @description
#' \code{ggdotplot} ggplot2 integrated graphics to produce dotplots
#' @rdname synteny
#' @import data.table
#' @import ggplot2
#' @importFrom grDevices pdf dev.off rgb
#' @importFrom dbscan dbscan frNN
#' @export
ggdotplot <- function(hits,
                      outDir = NULL,
                      minGenes2plot = 100,
                      appendName = "synHits",
                      dotsPerIn = 256,
                      plotSize = 12,
                      verbose = is.null(outDir)){
  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- blkID <-
    inBuffer <- rnd2 <- rnd1 <- ns2 <- n <- isArrayRep2 <- isArrayRep1 <-
    noAnchor <- NULL

  tp <- data.table(hits)

  un1 <- uniqueN(tp$ofID1)
  un2 <- uniqueN(tp$ofID2)
  if(un1 > un2){
    ht <- plotSize
    wd <- ht * (un1/un2)
  }else{
    wd <- plotSize
    ht <- wd * (un2/un1)
  }

  x <- max(tp$ord1, na.rm = T)
  y <- max(tp$ord2, na.rm = T)

  ordPerIn <- x / dotsPerIn
  totDots <- wd * dotsPerIn
  xrnd2 <- floor(x / totDots)

  ordPerIn <- y / dotsPerIn
  totDots <- ht * dotsPerIn
  yrnd2 <- floor(y / totDots)

  tp[,`:=`(rnd1 = round_toInteger(ord1, xrnd2),
           rnd2 = round_toInteger(ord2, yrnd2))]

  # pdf(dpFile, height = ht, width = wd*1.07)
  tp[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1]), by = "chr1"]
  tp[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2]), by = "chr2"]

  hc <- subset(tp, sameOg & ngene1 > minGenes2plot & ngene2 > minGenes2plot)
  hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
  setorder(hc, -n)
  hc[,ns2 := sum(n), by = c("chr2")]
  hc$n[hc$n > 20] <- 20

  un1 <- uniqueN(tp$ofID1)
  un2 <- uniqueN(tp$ofID2)

  setorder(hc, n)
  p1 <- ggplot(subset(hc, n > 1), aes(rnd1, rnd2, col = n)) +
    geom_point(pch = ".") +
    scale_color_viridis_c(begin = .1, trans = "log10", guide = "none") +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
    theme_genespace()+
    facet_grid(chr2 ~ chr1, scales = "free",
               space = "free", as.table = F, switch = "both")+
    labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome1[1], uniqueN(hits$ofID1)),
         y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome2[1], uniqueN(hits$ofID2)),
         title = sprintf("Blast hits where query and target are in the same OG, %s/%s-gene x/y window thresholded at 20 hits", xrnd2, yrnd2))

  hcBlk <- subset(tp, sameOg & inBuffer & ngene1 > minGenes2plot & ngene2 > minGenes2plot)
  hcBlk <- hcBlk[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2", "blkID")]
  blkCols <- sample(gs_colors(uniqueN(hcBlk$blkID)))
  p2 <- ggplot(subset(hcBlk, n > 1), aes(x = rnd1, y = rnd2, col = blkID)) +
    geom_point(pch = ".") +
    scale_color_manual(values = blkCols, guide = "none") +
    scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd1), by = 1e3))+
    scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hcBlk$rnd2), by = 1e3))+
    theme_genespace() +
    facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
    labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome1[1], uniqueN(hits$ofID1[hits$isAnchor])),
         y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                     hits$genome2[1], uniqueN(hits$ofID2[hits$isAnchor])),
         title = sprintf("Syntenic anchor blast hits, colored by block ID"))


  if(is.null(outDir)){
    if(verbose)
      cat("writing to the present graphics device")

  }else{
    dpFile <- file.path(outDir,
                        sprintf("%s_vs_%s.rawHits.pdf",
                                tp$genome1[1], tp$genome2[1]))
    if(verbose)
      cat(sprintf("writing to file: %s", dpFile))
    pdf(dpFile, height = ht, width = wd)
    print(p1)
    print(p2)
    dev.off()
  }
}

