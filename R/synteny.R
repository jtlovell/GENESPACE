#' @title flag and split syntenic hits
#' @description
#' \code{synteny} from an annotated blast file, assign syntenic (or not) blocks
#' to each hit.
#' @name synteny
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param verbose logical, should updates be printed to the console?
#' @param hits data.table of hits, see read_allBlast
#' @param synRad see init_genespace
#' @param onlyOgAnchors see init_genespace
#' @param blkSize see init_genespace
#' @param blkRadius see init_genespace
#' @param tmpDir see init_genespace
#' @param topn1 integer, the number of best scoring hits per gene in genome 1
#' @param topn2 integer, the number of best scoring hits per gene in genome 1
#' @param nGaps see init_genespace
#' @param onlySameChrs logical, should only hits on chromosomes with the same
#' name be permitted?
#' @param MCScanX_hCall see init_genespace
#'
#' \cr
#' If called, \code{synteny} returns its own arguments.
#'

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
    nTotalHits <- chunk <- inBuffer <- allBlast <- isArrayRep1 <- isArrayRep2 <-
    noAnchor <- ord1 <- ord2 <- isAnchor <- note <- genome1 <- genome2 <-
    ofID1 <- ofID2 <- nSynBlks <- status <- isSyntenic <- NULL

  nCores <- gsParam$params$nCores
  runOFINBLK <- gsParam$params$orthofinderInBlk

  # -- 1.2 check that files are all ok
  if(!"synteny" %in% names(gsParam))
    stop("must run set_syntenyParams prior to synteny")

  if(!all(file.exists(gsParam$synteny$blast$allBlast)))
    stop("some annotated blast files dont exist, run annotate_blast() first\n")

  # -- 1.1 split the metadata into chunks
  blMd <- data.table(gsParam$synteny$blast)
  blMd[,lab := align_charLeft(sprintf("%s v. %s:", query, target))]
  blMd[,selfOnly := query == target & queryPloidy == 1 & targetPloidy == 1]

  if(!"nGlobOgHits" %in% colnames(blMd))
    blMd[,nGlobOgHits := file.size(allBlast)]

  if(!"nTotalHits" %in% colnames(blMd))
    blMd[,nTotalHits := file.size(allBlast)]

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
      rawHits <- read_allBlast(x$allBlast)
      synHits <- subset(rawHits, isArrayRep1 & isArrayRep2 & !noAnchor)

      if(x$query == x$target & x$queryPloidy == 1){
        outHits <- find_selfSyn(
          hits = data.table(rawHits), synRad = x$synRad)
      }

      if(x$query == x$target & x$queryPloidy > 1){

        # -- self synteny (as above)
        selfHits <- find_selfSyn(
          hits = data.table(rawHits),
          synRad = x$synRad)

        # -- self synteny + mask (to reduce effects of large tandem dups)
        selfHits[,inBuffer := flag_hitsInRadius(
          x = ord1, y = ord2, isAnchor = isAnchor, radius = x$synBuff*2),
          by = c("chr1", "chr2")]
        maskHits <- subset(selfHits, !inBuffer)
        # -- run the synteny engine
        maskHits <- synteny_engine(
          hits = data.table(maskHits),
          nGaps = x$nGapsSecond,
          blkRadius = x$blkRadiusSecond,
          blkSize = x$blkSizeSecond,
          synRad = x$synRad,
          topn1 = x$targetPloidy - 1,
          topn2 = x$queryPloidy - 1,
          onlySameChrs = gsParam$params$onlySameChrs,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchorsSelf,
          verbose = FALSE)

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

      if(x$query != x$target){
        ##########################################################################
        # 2.2. primary hits
        outHits <- synteny_engine(
          hits = data.table(synHits),
          nGaps = x$nGaps,
          blkSize = x$blkSize,
          blkRadius = x$blkRadius,
          topn1 = x$targetPloidy,
          topn2 = x$queryPloidy,
          synRad = x$synRad,
          onlySameChrs = gsParam$params$onlySameChrs,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchors,
          verbose = FALSE)
      }

      ##########################################################################
      # 2.3 secondary hits
      if(x$nSecondaryHits > 0){
        hitsMask <- data.table(outHits)
        hitsMask[,inBuffer := flag_hitsInRadius(
          x = ord1, y = ord2, isAnchor = isAnchor, radius = x$synBuff),
          by = c("chr1", "chr2")]
        hitsPrim <- subset(hitsMask, inBuffer)
        hitsMask <- synteny_engine(
          hits = data.table(subset(hitsMask, !inBuffer)),
          nGaps = x$nGapsSecond,
          blkSize = x$blkSizeSecond,
          blkRadius = x$blkRadiusSecond,
          topn1 = x$targetPloidy,
          topn2 = x$queryPloidy,
          synRad = x$synRad,
          onlySameChrs = gsParam$params$onlySameChrs,
          MCScanX_hCall = gsParam$shellCalls$mcscanx_h,
          tmpDir = gsParam$paths$tmp,
          onlyOgAnchors = x$onlyOgAnchorsSecond)
        wh <- which(!is.na(hitsMask$regID))
        hitsMask$regID[wh] <- sprintf("%s_secondary", hitsMask$regID[wh])
        wh <- which(!is.na(hitsMask$blkID))
        hitsMask$blkID[wh] <- sprintf("%s_secondary", hitsMask$blkID[wh])
        outHits <- rbind(hitsMask, hitsPrim)
      }

      ##########################################################################
      # 2.4 check and compile notes
      outHits$isAnchor[is.na(outHits$isAnchor)] <- FALSE
      outHits$isSyntenic[is.na(outHits$isSyntenic)] <- FALSE
      out <- with(outHits, data.table(
        lab = chnk$lab[i],
        nTotalHits = length(chr1),
        nOgHits = sum(sameOG),
        nAnchorHits = sum(isAnchor),
        nSyntenicHits = sum(isSyntenic),
        nSynRegions = uniqueN(regID, na.rm = T),
        nSynBlks = uniqueN(blkID, na.rm = T)))
      nsynchr <- with(subset(outHits, isAnchor), uniqueN(paste(chr1, chr2)))
      out[,nSVs := nSynBlks - nsynchr]
      if("note" %in% colnames(outHits)){
        out[,status := ifelse(
          any(grepl("WARNING", outHits$note)), "no paralog synteny",
          ifelse(any(grepl("ERROR", outHits$note)), "no syntenic hits"))]
      }else{
        out[,status := "PASS"]
      }

      out <- data.table(outMd, out)

      write_allBlast(x = outHits, filepath = x$allBlast)
      write_synHits(x = subset(outHits, isSyntenic), filepath = x$synHits)

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
      "\t...%s %s hits (%s anchors) in %s blocks (%s SVs, %s regions)\n",
      lab, nSyntenicHits, nAnchorHits, nSynBlks, nSVs, nSynRegions)))

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
    blkID <- isAnchor <- isSyntenic <- genome2 <- NULL

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
  hits[,isSyntenic := paste(ofID1, ofID2) %in% ofInBuff |
         (chr1 == chr2 & chr1 %in% chrNosearch$chr1)]
  hits[,blkID := ifelse(isSyntenic, sprintf(
    "%s %s selfBlk %s", genome1, genome2, chr1), NA)]
  return(data.table(hits))
}

#' @title Find syntenic regions
#' @description
#' \code{synteny_engine} The main engine for synteny discovery. Takes a "hits"
#' data.table and classifies block IDs, whether a gene is an anchor in the block
#' and whether it is within a syntenic buffer of the block anchor hits. The
#' steps are as follows:
#' 1) Filter to initial hits - onlyOgAnchors (if TRUE) and topN hits depending
#' on ploidy. Built global collinear hits by piping these hits into MCScanX.
#' 2) Cluster global collinear hits into large regions searching within synRad
#' then pull all (or onlyOG hits) within synRad of global anchor hits.
#' 3) Re-run MCScanX within large regions, re-cluster and re-cull region hits
#' within synRad. Split overlapping regions. These are the "regions" named in
#' "regID" column
#' 4) For each region, pull all hits (regardless of score, OG etc) and re-run
#' MCScanX. Cluster these 'potential' anchors into blocks with blkRadius search
#' radius and blkSize minimum size. The resulting hits are the anchors, flagged
#' 'isAnchor = TRUE.
#' 5) Split anchors of interleaved blocks, then extract hits within the physical
#' bounds of the blocks and within blkRadius distance of an anchor (within-block
#' distance). These hits are flagged 'isSyntenic = TRUE'.
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
                           MCScanX_hCall,
                           onlySameChrs,
                           verbose = FALSE){


  noAnchor <- isArrayRep1 <- isArrayRep2 <- sameOG <- sr1 <- bitScore <- sr2 <-
    ord1 <- ord2 <- isAnchor <- ofID1 <- ofID2 <- chr1 <- chr2 <- blkID <-
    nAnchorhits <- genome1 <- genome2 <- start1 <- start2 <- end1 <- end2 <-
    i.blkID <- blkID1 <- blkID2 <- tmp <- note <- inBuffer <- need2keep <-
    need4anchor <- anyNeed <- hasSameOG1 <- hasSameOG2 <- hasSameOG <- clus <-
    nOG <- rl <- clus1 <- clus2 <- regID <- V1 <- V2 <- regID <- isSyntenic <-
    NULL

  # hits[,regID := NULL]
  raw <- data.table(hits)
  cull <- subset(hits, !noAnchor & isArrayRep1 & isArrayRep2)
  if(onlySameChrs)
    cull <- subset(cull, chr1 == chr2)

  if(verbose)
    cat(sprintf("%s vs %s:\n\tn. total hits = %s\n\tn. rep. hits = %s\n",
                hits$genome1[1], hits$genome2[2], nrow(hits), nrow(cull)))

  ##############################################################################
  ##############################################################################
  # 1. Pull global collinear hits
  # The overall goal here is to find the 'seed' anchors. These are hits that
  # MCScanX thinks are collinear, based on the initial criteria: (a) only OG
  # hits (optional, but default) and (b) the approx top scoring hits based on
  # ploidy (nhits). These initial seed hits are re-ranked by order across the
  # genome, then fed into MCScanX. Those that MCScanX thinks are collinear are
  # flagged
  ##############################################################################

  ##############
  # -- 1.1 [optionally] subset to only OG hits
  h1 <- data.table(cull)
  if(onlyOgAnchors)
    h1 <- subset(h1, sameOG)

  # if no hits are availble, return an note that will be convert to error
  if(nrow(h1) < blkSize){
    errs <- hits[1,]
    errs[,note := "ERROR no array rep hits"]
    return(errs)
  }

  ##############
  # -- 1.2 Subset to top-scoring hits
  # keeping anything within 10 or so, this gives a bit of wiggle in intial
  # anchor
  h1[,sr1 := frank(-round(bitScore, -1), ties.method = "dense"), by = "ofID1"]
  h1[,sr2 := frank(-round(bitScore, -1), ties.method = "dense"), by = "ofID2"]
  h1 <- subset(h1, sr1 <= topn1 & sr2 <= topn2)

  ##############
  # -- 1.3 rerank order
  h1[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
           ord2 = frank(ord2, ties.method = "dense"))]

  ##############
  # -- 1.4 pull collinear hits via mcscanx
  tmp <- h1[,c("ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2", "bitScore")]
  setnames(tmp, "bitScore", "score")
  tmp <- subset(tmp, complete.cases(tmp))
  anchu <- run_mcscanx(
    nGaps = 25,
    tmpDir = tmpDir,
    blkSize = 5,
    hits = data.table(tmp),
    MCScanX_hCall = MCScanX_hCall)
  cull[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]

  if(sum(cull$isAnchor) < blkSize){
    errs <- hits[1,]
    errs[,note := "ERROR no synteny"]
    return(errs)
  }

  if(verbose)
    cat(sprintf("\tn. potenial anchor hits = %s\n\tn. seed anchor hits = %s\n",
                nrow(h1), sum(cull$isAnchor)))

  ##############################################################################
  ##############################################################################
  # 2. Find hits in 'initial' seed regions
  # using the seed anchor hits from step 1, we find large syntenic 'initial
  # regions'. These are determined by re-ranking the position of the seed hits
  # and clustering using the syntenic radius ('synRad') in 2d dbscan. These
  # initial regions are used to physically query any nearby hits as follows: (a)
  # hits must be within the physical bounds, (b) must be within the synteny
  # radius of anchors, except that (c) if the query and target are same
  # orthogroup, can extend outside of the regions (still needs to be within the
  # synBuff of an initial seed anchor)
  ##############################################################################

  ##############
  # -- 2.1 cluster into large regions
  h2 <- subset(cull, isAnchor)
  h2[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
           ord2 = frank(ord2, ties.method = "dense"))]

  h2[,blkID := dbscan(frNN(
    x = cbind(ord1, ord2),
    eps = synRad),
    minPts = blkSize)$cluster,
    by = c("chr1", "chr2")]
  h2 <- subset(h2, blkID != 0)
  h2[,blkID := paste(chr1, chr2, blkID)]
  bu <- h2$blkID; names(bu) <- with(h2, paste(ofID1, ofID2))
  cull[,`:=`(blkID = bu[paste(ofID1, ofID2)],
             isAnchor = paste(ofID1, ofID2) %in% names(bu))]

  ##############
  # -- 2.2 pull hits on same chromosomes as anchors and re-rank
  chru <- with(subset(h2, !is.na(blkID)), unique(paste(chr1, chr2)))
  h2 <- subset(cull, paste(chr1, chr2) %in% chru)
  h2[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
           ord2 = frank(ord2, ties.method = "dense"))]

  ##############
  # -- 2.3 pull hits in buffer to the original anchors
  h2[,inBuffer := flag_hitsInRadius(
    x = ord1, y = ord2,
    isAnchor = isAnchor, radius = synRad),
    by = c("chr1", "chr2")]
  h2 <- subset(h2, inBuffer)

  ##############
  # -- 2.4 flag the block IDs of each hit
  h2[,blkID := flag_hitsInBlk(x = ord1, y = ord2, blkID = blkID)]

  ##############
  # -- 2.5 Let the coordinates expand to any genes in the same orthogroup
  h2[,need2keep := is.na(blkID) & sameOG & inBuffer]
  h2[,need4anchor := !is.na(blkID) & isAnchor & !is.na(blkID)]
  h2[,anyNeed := any(need2keep) & any(need4anchor), by = c("chr1", "chr2")]
  spl <- split(subset(h2, anyNeed), by = c("chr1", "chr2"))
  nnout <- rbindlist(lapply(spl, function(inchr){

    hasAnch <- subset(inchr, need4anchor)
    out <- subset(inchr, need2keep)
    fn <- frNN(
      x = hasAnch[,c("ord1", "ord2")],
      query = out[,c("ord1", "ord2")],
      eps = synRad)

    index <- sapply(fn$id, function(x) x[1])
    out <- out[!is.na(index)]
    index <- index[!is.na(index)]

    out[,blkID := hasAnch$blkID[index]]
    out <- subset(out, !is.na(blkID))
    if(nrow(out) > 0)
      return(out)
  }))
  u <- with(nnout, paste(ofID1, ofID2))
  h2 <- rbind(subset(h2, !paste(ofID1, ofID2) %in% u), nnout)

  ##############
  # -- 2.6 ensure we have synteny
  h3 <- subset(h2, !is.na(blkID))
  h2 <- NULL

  if(nrow(h3) < blkSize){
    errs <- hits[1,]
    errs[,note := "ERROR no synteny"]
    return(errs)
  }

  if(verbose)
    cat(sprintf("\tn. inital anchors hits / blocks = %s / %s\n",
                nrow(h3), uniqueN(h3$blkID)))
  ##############################################################################
  ##############################################################################
  # 3. Find Syntenic anchors
  # Here, we change our view of synteny from that of the genome wide comparison
  # (sections 1-2) to purely within the large inital regions (built in step 2).
  # The main logic here is that when considered in isolation, collinearity
  # within each syntenic region can be inferred with the same parameters, even
  # if they have very different evolutionary histories (this subgenome
  # dominance). There are four steps here: (a) potential anchors are defined
  # as hits with the top scoring hit for both the query and target gene OR in
  # the same orthogroup within each large region, (b) gene rank order position
  # is re-calculated within each region, (c) mcscanx is re-run within each
  # region, (d) hits within a fixed radius of anchors from (c) are pulled and
  # iteratively clustered into blocks using a radius from blkRadius to blkSize.
  # These anchors are cleaned up into finalized blocks in the next section.
  ##############################################################################

  ##############
  # -- 3.1 just to the top scoring or OG hits in each block
  h3[,sr1 := frank(-bitScore, ties.method = "first"),
     by = c("blkID", "ofID1")]
  h3[,sr2 := frank(-bitScore, ties.method = "first"),
     by = c("blkID", "ofID2")]
  h3 <- subset(h3, (sr1 == 1 & sr2 == 1) | sameOG)

  ##############
  # -- 3.2 drop hits that are top scoring but are not OG and either the
  # query or target has an OG hit in the block
  h3[,hasSameOG1 := any(sameOG), by = c("blkID", "ofID1")]
  h3[,hasSameOG2 := any(sameOG), by = c("blkID", "ofID2")]
  h3[,hasSameOG := hasSameOG1 | hasSameOG2]
  h3 <- subset(h3, !hasSameOG | (hasSameOG & sameOG))

  ##############
  # -- 3.3 re-rank gene order within blocks
  h3[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
           ord2 = frank(ord2, ties.method = "dense")),
     by = "blkID"]

  ##############
  # -- 3.4 mcscanx within block to find good proximate hits
  spl <- split(h3, by = "blkID")
  outu <- rbindlist(lapply(spl, function(h3spl){
    tmp <- h3spl[,c("ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2", "bitScore")]
    setnames(tmp, "bitScore", "score")
    tmp <- subset(tmp, complete.cases(tmp))
    anchu <- run_mcscanx(
      nGaps = nGaps,
      tmpDir = tmpDir,
      blkSize = blkSize,
      hits = data.table(tmp),
      MCScanX_hCall = MCScanX_hCall)
    tmp[,isAnchor := paste(ofID1, ofID2) %in% names(anchu)]
    return(subset(tmp, isAnchor)[,c("ofID1", "ofID2")])
  }))
  anchu <- with(outu, paste(ofID1, ofID2))
  h3[,isAnchor := paste(ofID1, ofID2) %in% anchu]

  ##############
  # -- 3.5 pull inbuffer hits
  h3[,inBuffer := flag_hitsInRadius(
    x = ord1, y = ord2,
    isAnchor = isAnchor, radius = blkRadius),
    by = c("chr1", "chr2", "blkID")]
  h3 <- subset(h3, inBuffer)

  ##############
  # -- 3.6 iterative pruning to drop non-collinear hits and form new anchors
  steps <- seq(from = blkRadius, to = blkSize)
  steps <- floor(steps)
  steps <- steps[!duplicated(steps)]
  for(i in steps){
    h3[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
              ord2 = frank(ord2, ties.method = "dense")), by = "blkID"]
    h3[,clus := dbscan(frNN(
      x = cbind(ord1, ord2),
      eps = i),
      minPts = blkSize)$cluster,
      by = c("chr1", "chr2", "blkID")]
    h3 <- subset(h3, clus > 0)
    if(onlyOgAnchors){
      h3[,nOG := sum(sameOG), by = c("chr1", "chr2", "blkID", "clus")]
      h3 <- subset(h3, nOG > 0)
    }
  }

  ##############
  # -- 3.7 drop no-OG blocks
  setkey(h3, blkID, ord1)
  h3[,rl := add_rle(clus, which = "id"), by = "blkID"]
  setkey(h3, blkID, ord2)
  h3[,rl := add_rle(rl, which = "id"), by = "blkID"]
  if(onlyOgAnchors){
    h3[,nOG := sum(sameOG), by = c("chr1", "chr2", "blkID", "rl")]
    h3 <- subset(h3, nOG > 0)
  }

  ##############
  # -- 3.8 final clustering of anchors within regions
  h3[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense")), by = "blkID"]
  h3[,clus := dbscan(frNN(
    x = cbind(ord1, ord2),
    eps = i),
    minPts = 1)$cluster,
    by = c("chr1", "chr2", "blkID")]
  h3 <- subset(h3, clus > 0)
  if(onlyOgAnchors){
    h3[,nOG := sum(sameOG), by = c("chr1", "chr2", "blkID", "clus")]
    h3 <- subset(h3, nOG > 0)
  }

  ##############
  # -- 3.9 final splitting within regions by RLEs
  setkey(h3, blkID, ord1)
  h3[,rl := add_rle(clus, which = "id"), by = "blkID"]
  setkey(h3, blkID, ord2)
  h3[,rl := add_rle(rl, which = "id"), by = "blkID"]
  if(onlyOgAnchors){
    h3[,nOG := sum(sameOG), by = c("chr1", "chr2", "blkID", "rl")]
    h3 <- subset(h3, nOG > 0)
  }

  h3[,blkID := paste(chr1, chr2, blkID, clus, rl)]
  setkey(h3, ord1)
  h3[,blkID := as.integer(factor(blkID, levels = unique(blkID)))]
  bv <- h3$blkID; names(bv) <- with(h3, paste(ofID1, ofID2))

  cull[,blkID := bv[paste(ofID1, ofID2)]]

  ##############################################################################
  ##############################################################################
  # 4. Final clean up and checking of blocks/regions
  #
  ##############################################################################

  ##############
  # -- 4.1 re-rank gene order for syntenic anchors
  h4 <- subset(cull, isAnchor)
  h4 <- subset(h4, !is.na(blkID))
  h4[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
           ord2 = frank(ord2, ties.method = "dense"))]

  ##############
  # -- 4.2 check for consecutive gaps inside blocks > blkRadius
  setkey(h4, ord1)
  h4[,clus1 := cumsum(c(1, as.numeric(diff(ord1) > blkRadius))),
     by = c("chr1", "blkID")]
  setkey(h4, ord2)
  h4[,clus2 := cumsum(c(1, as.numeric(diff(ord2) > blkRadius))),
     by = c("chr2", "blkID")]
  h4[,blkID := paste(blkID, clus1, clus2)]

  ##############
  # -- 4.3 add new block IDs back to the original data
  h4[,blkID := as.integer(factor(blkID, levels = unique(blkID)))]
  h4[,blkID := sprintf("%s_vs_%s: %s", genome1, genome2, blkID)]
  bv <- h4$blkID; names(bv) <- with(h4, paste(ofID1, ofID2))

  cull[,blkID := bv[paste(ofID1, ofID2)]]
  cull[,isAnchor := !is.na(blkID)]

  ##############
  # -- 4.4 initial cluster blocks into regions
  h4 <- subset(cull, isAnchor)
  h4[,`:=`(ord1 = frank(ord1, ties.method = "dense"),
            ord2 = frank(ord2, ties.method = "dense"))]
  h4[,clus := dbscan(frNN(
    x = cbind(ord1, ord2),
    eps = blkRadius),
    minPts = blkSize)$cluster,
    by = c("chr1", "chr2")]
  h4[,regID := paste(chr1, chr2, clus)]
  h4[,regID := as.integer(factor(regID, levels = unique(regID)))]

  ##############
  # -- 4.5 merge regions that have overlapping blocks
  regMerge <- h4[,CJ(unique(blkID), unique(blkID)), by = "regID"]
  regMerge[,clus := clus_igraph(V1, V2)[V1]]
  regMerge <- subset(regMerge, !duplicated(V1))
  rv <- regMerge$V2; names(rv) <- regMerge$V1
  h4[,regID := rv[blkID]]
  h4 <- subset(h4, !is.na(blkID))

  ##############
  # -- 4.6 add regions back to the data
  setkey(h4, ord1)
  h4 <- subset(h4, !is.na(regID))
  h4[,regID := as.integer(factor(regID, levels = unique(regID)))]
  h4[,regID := sprintf("%s_vs_%s: %s", genome1, genome2, regID)]
  rv <- h4$regID; names(rv) <- with(h4, paste(ofID1, ofID2))
  cull[,regID := rv[paste(ofID1, ofID2)]]

  ################################################################################
  # 5. Flag syntenic hits

  ##############
  # -- 5.1 find global hits in buffer
  h5 <- data.table(cull)
  h5[,inBuffer := flag_hitsInRadius(
    x = ord1, y = ord2,
    isAnchor = isAnchor, radius = blkRadius),
    by = c("chr1", "chr2")]

  ##############
  # -- 5.2 find hits in physical limits of blocks
  h5 <- subset(h5, inBuffer)
  h5[,regID := flag_hitsInBlk(x = ord1, y = ord2, blkID = regID)]
  h5[,blkID := flag_hitsInBlk(x = ord1, y = ord2, blkID = blkID)]

  ##############
  # -- 5.3 assign block IDs to OG hits just outside of block coords
  h5[,need2keep := is.na(regID) & sameOG & inBuffer]
  h5[,need4anchor := !is.na(regID) & isAnchor & !is.na(blkID)]
  h5[,anyNeed := any(need2keep) & any(need4anchor), by = c("chr1", "chr2")]
  if(any(h5$anyNeed)){
    spl <- split(subset(h5, anyNeed), by = c("chr1", "chr2"))
    nnout <- rbindlist(lapply(spl, function(inchr){

      hasAnch <- subset(inchr, need4anchor)
      out <- subset(inchr, need2keep)
      fn <- frNN(
        x = hasAnch[,c("ord1", "ord2")],
        query = out[,c("ord1", "ord2")],
        eps = synRad)

      index <- sapply(fn$id, function(x) x[1])
      out <- out[!is.na(index)]
      index <- index[!is.na(index)]

      out[,blkID := hasAnch$blkID[index]]
      out[,regID := hasAnch$regID[index]]
      out <- subset(out, !is.na(blkID) | !is.na(regID))
      if(nrow(out) > 0)
        return(out)
    }))
    u <- with(nnout, paste(ofID1, ofID2))
    h5 <- rbind(subset(h5, !paste(ofID1, ofID2) %in% u), nnout)
  }

  rv <- h5$regID; bv <- h5$blkID
  names(rv) <- names(bv) <- with(h5, paste(ofID1, ofID2))

  cull[,`:=`(regID = rv[paste(ofID1, ofID2)],
             blkID = bv[paste(ofID1, ofID2)])]
  cull[,isSyntenic := !is.na(regID) | !is.na(blkID)]

  if(verbose)
    cat(sprintf("\tn. final anchors = %s \n\tn. final syntenic hits = %s\n",
                sum(cull$isAnchor), sum((cull$isSyntenic))))
  if(verbose)
    cat(sprintf("\tn. final blocks = %s \n\tn. final regions = %s\n",
                uniqueN(cull$blkID, na.rm = T),
                uniqueN(cull$regID, na.rm = T)))
  u <- with(cull, paste(ofID1, ofID2))
  out <- subset(raw, !paste(ofID1, ofID2) %in% u)
  outSyn <- cull[,colnames(out), with = F]
  out <- rbind(outSyn, out)
  return(out)
}

