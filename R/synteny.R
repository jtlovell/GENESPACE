#' @title synteny
#' @description
#' \code{synteny} synteny
#' @name synteny
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
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
synteny <- function(gsParam, verbose = TRUE){
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
      outMd <- data.table(chnk[i,])
      x <- data.table(outMd)
      rawHits <- read_synHits(x$synHits)

       ############################################################################
      # 1. intragenomic haploid hits
      if(x$query == x$target & x$queryPloidy == 1){
        outHits <- find_selfSyn(
          hits = data.table(rawHits), synRad = x$synRad)
      }

      ############################################################################
      # 2. intragenomic polyploid
      if(x$query == x$target & x$queryPloidy > 1){

        # -- self synteny (as above)
        selfHits <- find_selfSyn(
          hits = data.table(rawHits),
          synRad = x$synRad)

        # -- self synteny + mask (to reduce effects of large tandem dups)
        maskHits <- find_selfSyn(
          hits = data.table(rawHits),
          synRad = gsParam$params$maskBuffer)
        maskHits <- subset(maskHits, !inBuffer)

        # -- run the synteny engine
        maskHits <- synteny_engine(
          hits = data.table(maskHits),
          nGaps = x$nGapsSecond,
          blkRadius = x$blkRadiusSecond,
          blkSize = x$blkSizeSecond,
          synRad = x$synRad,
          topn1 = x$targetPloidy - 1,
          topn2 = x$queryPloidy - 1,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchorsSelf)

        # -- if there is no synteny, return an error
        if(sum(maskHits$isAnchor) < x$blkSize){
          um <- NA
          maskHits <- maskHits[1,]
          maskHits[,note := sprintf("WARNING!! no paralogous syntenic hits in %s vs %s. You may need a better outgroup or set onlyOgAnchorsSelf = FALSE", genome1, genome2)]
        }else{
          um <- with(maskHits, unique(paste(ofID1, ofID2)))
        }
        outHits <- rbind(maskHits, subset(selfHits, !paste(ofID1, ofID2) %in% um), fill = T)
      }

      ##########################################################################
      # 3. intergenomic hits
      if(x$query != x$target){
        outHits <- synteny_engine(
          hits = data.table(rawHits),
          nGaps = x$nGaps,
          blkSize = x$blkSize,
          blkRadius = x$blkRadius,
          topn1 = x$targetPloidy,
          topn2 = x$queryPloidy,
          synRad = x$synRad,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchors)

        # -- if there is no synteny, return an error
        if(sum(outHits$isAnchor) < x$blkSize){
          outHits <- rawHits[1,]
          outHits[,note := sprintf("ERROR!! no primary syntenic hits in %s vs %s. If this isn't expected, look at the raw blast hits to troubleshoot or relax synteny settings", genome1, genome2)]
        }
      }

      ##########################################################################
      # 4. secondary hits
      if(x$nSecondaryHits > 0){
        hitsMask <- subset(outHits, !inBuffer)
        hitsPrim <- subset(outHits, inBuffer)
        hitsMask <- synteny_engine(
          hits = data.table(hitsMask),
          nGaps = x$nGapsSecond,
          blkSize = x$blkSizeSecond,
          blkRadius = x$blkRadiusSecond,
          topn1 = x$targetPloidy,
          topn2 = x$queryPloidy,
          synRad = x$synRad,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchorsSecond)
        outHits <- rbind(hitsMask, hitsPrim)
      }

      out <- with(outHits, data.table(
        lab = chnk$lab[i],
        nAnchorHits = sum(isAnchor),
        nBufferHits = sum(inBuffer),
        nBlks = uniqueN(blkID, na.rm = T),
        nSVs = uniqueN(blkID, na.rm = T) -
          uniqueN(paste(chr1, chr2)[!is.na(blkID)])))

      if("note" %in% colnames(outHits)){
        out[,status := ifelse(any(grepl("WARNING", outHits$note)), "no paralog synteny",
                              ifelse(any(grepl("ERROR", outHits$note)), "no syntenic hits"))]
      }else{
        out[,status := "PASS"]
      }

      if("nAnchorHits" %in% colnames(outMd)){
        outMd[,`:=`(nAnchorHits = NULL, lab = NULL,
                nBufferHits = NULL, nBlks = NULL, nSVs = NULL)]
      }

      out <- data.table(outMd, out)
      out[,`:=`(nTotalHits = nrow(rawHits),
                nGlobOgHits = sum(outHits$sameOg))]
      write_synBlast(x = outHits, filepath = x$synHits)
      return(out)
    }), fill = T)

    if(any(outChnk$status == "no syntenic hits")){
      cat("\n############\nERROR - No synteny found for the following pairs:\n")
      print(subset(outChnk, status == "no syntenic hits")[,c("query", "target")])
      stop("GENESPACE requires synteny to work ... try new parameters or genomes\n")
    }

    if(any(outChnk$status == "no paralog synteny")){
      cat("\n############\nWARNING - No synteny found for paralogs among the following pairs:\n")
      print(subset(outChnk, status == "no paralog synteny")[,c("query", "target")])
      cat(strwrap(
        "If you have polyploids, you really should use an outgroup (see README).
      If this isn't possible, you can relax secondary synteny parameters or set
      onlyOgAnchorsSelf to FALSE\n",  indent = 8, exdent = 8), sep = "\n")
      cat("############\n")
    }

    with(outChnk, cat(sprintf(
      "\t...%s %s hits (%s anchors) in %s blocks (%s SVs)\n",
      lab, nBufferHits, nAnchorHits, nBlks, nSVs)))

    return(outChnk)
  }))

  gsParam$synteny$blast <- blMdOut
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
  ofID1 <- ofID2 <- chr1 <- chr2 <- ord1 <- ord2 <- genome1 <- n <-
    blkID <- isAnchor <- inBuffer <- genome2 <- NULL

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
  return(data.table(hits))
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
  nOghits <- NULL
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
          hs[,`:=`(nOghits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOghits > 0)
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
          hs[,`:=`(nOghits = sum(sameOg),
                   n = min(c(uniqueN(ofID1, na.rm = T),
                             uniqueN(ofID2, na.rm = T)))), by = "rl"]
          hs <- subset(hs, n >= blkSize | nOghits > 0)
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
    nAnchorhits <- genome1 <- genome2 <- start1 <- start2 <- end1 <- end2 <-
    i.blkID <- blkID1 <- blkID2 <- tmp <- NULL

  ################################################################################
  # 1. Find initial syntenic anchors
  tmp <- subset(hits, !noAnchor & isArrayRep1 & isArrayRep2)
  if(nrow(tmp) < blkSize){
    errs <- hits[1,]
    errs[,note := "ERROR no array rep hits"]
    return(errs)
  }else{
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
      hits = data.table(tmp),
      MCScanX_hCall = MCScanX_hCall)
    hits[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]

    if(sum(hits$isAnchor) < blkSize){
      errs <- hits[1,]
      errs[,note := "ERROR no syntenic hits"]
      return(errs)
    }else{

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
        hits = data.table(tmp),
        MCScanX_hCall = MCScanX_hCall)

      # -- 3.2 re-rank and cluster into regions
      tmp <- subset(hits, paste(ofID1, ofID2) %in% names(anchu))

      if(nrow(tmp) < blkSize){
        errs <- hits[1,]
        errs[,note := "ERROR no syntenic hits"]
        return(err)
      }else{
        tmp[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
                  ord2 = frank(ord2, ties.method = "dense"))]
        tmp[,blkID := dbscan(frNN(
          x = cbind(ord1, ord2),
          eps = blkRadius),
          minPts = blkSize)$cluster,
          by = c("chr1", "chr2")]

        tmp[,blkID := paste("blk", as.numeric(as.factor(paste(chr1, chr2, blkID))))]
        # -- 3.3 split overlapping blocks
        tmp <- split_ovlBlks(
          tmp,
          dropSmallNonOGBlks = onlyOgAnchors,
          blkSize = blkSize*2)
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
          start1 = min(start1), end1 = max(end1),
          start2 = min(start2), end2 = max(end2),
          nAnchorhits = .N),
          by = c("chr1", "chr2", "blkID")]
        blks <- subset(blks, complete.cases(blks) & nAnchorhits > blkSize)
        setorder(blks, -nAnchorhits)

        # -- 4.3 split out hits by blk coords
        b1 <- blks[,c("chr1", "start1" ,"end1", "blkID")]
        b2 <- blks[,c("chr2", "start2" ,"end2", "blkID")]
        setkey(b1, chr1, start1, end1)
        setkey(b2, chr2, start2, end2)
        setkey(tmp, chr1, start1, end1)
        t1 <- foverlaps(b1, tmp)
        t1[,`:=`(i.start1 = NULL, i.end1 = NULL, blkID1 = i.blkID, i.blkID = NULL)]
        setkey(t1, chr2, start2, end2)
        t2 <- foverlaps(b2, t1)
        t2[,`:=`(i.start2 = NULL, i.end2 = NULL, blkID2 = i.blkID, i.blkID = NULL)]
        tmpOut <- subset(t2, blkID1 == blkID2)
        tmpOut[,blkID := paste(genome1, genome2, as.numeric(as.factor(paste(chr1, chr2, blkID1))))]
        spl <- split(tmpOut, by = "blkID")
        ofInBuff <- with(rbindlist(lapply(spl, function(x){
          nn <- with(x, frNN(x = data.frame(ord1, ord2), eps = synRad))
          wh <- unique(c(which(x$isAnchor), unlist(nn$id[x$isAnchor])))
          return(subset(x, 1:nrow(x) %in% wh)[,c("ofID1", "ofID2")])
        })), paste(ofID1, ofID2))
        tmpOut <- subset(tmpOut, paste(ofID1, ofID2) %in% ofInBuff)

        bu <- with(tmpOut, paste(ofID1, ofID2))
        au <- with(subset(tmpOut, isAnchor), paste(ofID1, ofID2))
        bid <- tmpOut$blkID; names(bid) <- with(tmpOut, paste(ofID1, ofID2))

        hits[,`:=`(blkID = bid[paste(ofID1, ofID2)],
                   inBuffer = paste(ofID1, ofID2) %in% bu,
                   isAnchor = paste(ofID1, ofID2) %in% au)]

        return(hits)
      }
    }
  }
}
