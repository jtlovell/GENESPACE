#' @title Genespace plotting routines
#' @description
#' \code{plot_genespace} Genespace plotting routines
#' @name plot_genespace
#'
#' @param synParam file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param gsAnnot named character vector of file.paths to gff-like annotation
#' files. The names in this vector must match genome1/genome2 columns in the
#' syntenyParams data.table.
#'
#' @details ...
#'
#' @note \code{plot_genespace} is a generic name for the functions documented.
#' \cr
#' If called, \code{plot_genespace} returns its own arguments.
#'
#' @title Riparian plot
#' @description
#' \code{plot_riparian} Make a riparian plot
#' @rdname plot_genespace
#' @import data.table
#' @export
plot_riparian <- function(gsParam,
                          gsAnnot,
                          synParam,
                          hitsFile = NULL,
                          hits = NULL,
                          blks = NULL,
                          braidColor = "blue",
                          inversionColor = "blue",
                          braidAlpha = .8,
                          braidBorderLwd = NULL,
                          blksFile = NULL,
                          genomeIDs = NULL,
                          refGenome = NULL,
                          reorderChrs = TRUE,
                          minGenes = 5,
                          gapProp = .005,
                          useOrder = TRUE,
                          axisTitleCex = 1,
                          chrLabCex = .5,
                          chrRectWidth = NULL,
                          chrRectBuffer = 1.5,
                          refAtTop = F,
                          chrBorder = "black",
                          chrFill = "white",
                          genomesWithChrIDs = NULL,
                          colByChrs = c("#6B2701", "#ED9004", "#F9C70E", "#EAE075", "#BAE0DB", "#8BEDF9", "#74B8FC", "#4871F9", "#040DC9","#0E004C","#5E09A3","#C054F9","#E6BDFC"),
                          labelChrBiggerThan = ifelse(useOrder, 500, 5e6),
                          xlab = sprintf("Chromosomes and %s", ifelse(useOrder, "gene rank order", "physical position")),
                          ylab = "genomes",
                          chrNameFun = function(x) gsub("^chr|^scaffold","",tolower(x)),
                          title = "Riparian plot"){
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]

  if(is.null(refGenome))
    refGenome <- genomeIDs[1]

  genomeIDs <- c(refGenome, genomeIDs[genomeIDs != refGenome])

  if(!is.null(genomesWithChrIDs))
    genomesWithChrIDs <- genomesWithChrIDs[genomesWithChrIDs %in% genomeIDs]
  if(is.null(genomesWithChrIDs))
    genomesWithChrIDs <- genomeIDs



  if(is.null(blksFile) & is.null(blks)){
    blksFile <- file.path(gsParam$results, "blockCoordinates_geneOrder.txt.gz")
    if(!file.exists(blksFile))
      stop("no blks data.table or blksFile file.path specied. Also could not find blockCoordinates file in /results")
  }
  if(is.null(blks))
    blks <- fread(blksFile)
  tmp <- blks[,c("chr2","chr1","blkID","ordStart2","ordEnd2","ordStart1", "ordEnd1",
                 "bpStart2","bpEnd2","bpStart1","bpEnd1","nHits","orient",
                 "firstGene2", "firstGene1", "lastGene2", "lastGene1",
                 "genome2","genome1")]
  setnames(tmp, names(blks))
  blks <- rbind(blks, tmp)
  blks <- subset(blks, !duplicated(blks))

  chrl <- match_synChrs(
    refGenome = refGenome,
    gsParam = gsParam,
    gsAnnot = gsAnnot,
    synParam = synParam,
    minGenes = minGenes,
    reorderChrs = reorderChrs)

  ucl <- with(blks, data.table(genome = c(genome1, genome2), chr = c(chr1, chr2)))
  ucl <- with(subset(ucl, !duplicated(ucl)), split(chr, genome))
  for(i in names(chrl))
    chrl[[i]] <- chrl[[i]][chrl[[i]] %in% ucl[[i]]]

  genePos <- calc_linearGenePos(
    gsAnnot = gsAnnot,
    gsParam = gsParam,
    chrList = chrl,
    gapProp = gapProp)

  if(useOrder){
    xv <- genePos$xord
    chrPos <- genePos[,list(chrStart = min(xord),
                            chrEnd = max(xord)),
                      by = c("genome","y","chr")]
  }else{
    xv <- genePos$xbp
    chrPos <- genePos[,list(chrStart = min(xbp),
                            chrEnd = max(xbp)),
                      by = c("genome","y","chr")]
  }
  setkey(chrPos, y, chrStart)
  if(refAtTop){
    yv <- -genePos$y
    names(xv) <- names(yv) <- genePos$ofID
    chrPos[,y := -y]
  }else{
    yv <- genePos$y
    names(xv) <- names(yv) <- genePos$ofID
  }


  # make the plot
  with(chrPos,
       plot(
         1, 1, type = "n", axes = F, bty = "n",
         xlim = c(min(chrStart) - (max(chrEnd)/10),
                  max(chrEnd)),
         ylim = c(min(y), max(y)),
         ylab = ylab,
         xlab = xlab,
         main = title))
  if (is.null(chrRectWidth))
    chrRectWidth <- strheight("chr", cex = chrLabCex)*chrRectBuffer

  setkey(blks, genome1, genome2, chr1, chr2, orient)

  if(length(colByChrs) > 0){
    cpBlks <- make_blksByRefChr(
      refGenome = refGenome,
      gsAnnot = gsAnnot,
      blks = blks,
      genePos = genePos,
      gsParam = gsParam,
      chrList = chrl,
      synParam = synParam,
      genomeIDs = genomeIDs)
    if(length(colByChrs) == length(chrl[[1]])){
      cols <- colByChrs
    }else{
      cols <- colorRampPalette(colByChrs)(length(chrl[[1]]))
    }
    names(cols) <- chrl[[1]]

    chrBlks <- with(cpBlks, data.table(
      g1 = genome1, g2 = genome2, blkID = blkID,
      c1 = chr1, c2 = chr2, y1 = yv[firstGene1], y2 = yv[firstGene2],
      f1 = xv[firstGene1], l1 = xv[lastGene1],
      f2 = ifelse(orient == "+", xv[firstGene2], xv[lastGene2]),
      l2 = ifelse(orient == "+", xv[lastGene2], xv[firstGene2]),
      col = add_alpha(cols[refChr], braidAlpha), u = 1:length(genome1)))

    chrPolygonList <- mclapply(split(chrBlks, by = "u"), mc.cores = gsParam$nCores, function(x)
      with(x, calc_curvePolygon(
        start1 = f1, start2 = f2,
        end1 = l1, end2 = l2,
        y1 = y1 + chrRectWidth/4, y2 = y2 - chrRectWidth/4)))
    nu <- sapply(1:length(chrPolygonList), function(i){
      if(!is.null(braidBorderLwd)){
        polygon(chrPolygonList[[i]],
                border = ifelse(is.null(braidBorderLwd), NA, chrBlks$col[i]),
                lwd = ifelse(is.null(braidBorderLwd), NA, braidBorderLwd),
                col = chrBlks$col[i])
      }else{
        polygon(chrPolygonList[[i]],
                border = NA,
                col = chrBlks$col[i])
      }
    })
  }else{
    tpBlks <- subset(with(blks, data.table(
      g1 = genome1, g2 = genome2, blkID = blkID,
      c1 = chr1, c2 = chr2, y1 = yv[firstGene1], y2 = yv[firstGene2],
      f1 = xv[firstGene1], l1 = xv[lastGene1],
      f2 = ifelse(orient == "+", xv[firstGene2], xv[lastGene2]),
      l2 = ifelse(orient == "+", xv[lastGene2], xv[firstGene2]),
      col = ifelse(orient == "+",
                   add_alpha(braidColor, braidAlpha),
                   add_alpha(inversionColor, braidAlpha)))),
      y2 - y1 == 1)
    polygonList <- lapply(split(tpBlks, by = "blkID"), function(x)
      with(x, calc_curvePolygon(
        start1 = f1, start2 = f2,
        end1 = l1, end2 = l2,
        y1 = y1 + chrRectWidth/4, y2 = y2 - chrRectWidth/4)))
    nu <- sapply(1:length(polygonList), function(i){
      if(!is.null(braidBorderLwd)){
        polygon(polygonList[[i]],
                border = ifelse(is.null(braidBorderLwd), NA, tpBlks$col[i]),
                lwd = ifelse(is.null(braidBorderLwd), NA, braidBorderLwd),
                col = tpBlks$col[i])
      }else{
        polygon(polygonList[[i]],
                border = NA,
                col = tpBlks$col[i])
      }
    })
  }

  # chrPos <<- chrPos
  for (i in 1:nrow(chrPos)) {
    polygon(with(chrPos[i,], round_rect(
      xleft = chrStart, xright = chrEnd,
      ybottom = y - chrRectWidth/2,
      ytop = y + chrRectWidth/2)),
      border = chrBorder,
      col = chrFill,
      lwd = .5)
    if(with(chrPos[i,], chrEnd - chrStart) > labelChrBiggerThan){
      with(chrPos[i,],
           text((chrStart + chrEnd)/2, y, labels = chrNameFun(chr), cex = chrLabCex))
    }
  }

  cp <- chrPos[,list(x = min(chrStart)), by = c("genome","y")]
  with(cp, text(x = x, y = y, label= genome, adj = c(1.2,.5), cex = axisTitleCex))
}

#' @title make_blksByRefChr
#' @description
#' \code{make_blksByRefChr} make_blksByRefChr
#' @rdname plot_genespace
#' @import data.table
#' @export
make_blksByRefChr <- function(refGenome,
                              blks,
                              chrList,
                              gsAnnot,
                              genePos,
                              gsParam,
                              synParam,
                              genomeIDs){
  gffAll <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), blastDir = gsParam$blast)
  ov <- gffAll$ord; names(ov) <- gffAll$ofID
  cv <- gffAll$chr; names(cv) <- gffAll$ofID
  sv <- gffAll$start; names(sv) <- gffAll$ofID
  ev <- gffAll$end; names(ev) <- gffAll$ofID
  gv <- gffAll$genome; names(gv) <- gffAll$ofID

  sp <- subset(synParam, runBlast)
  sp[,sogf := file.path(gsParam$synteny, sprintf("%s_vs_%s_synog.txt.gz", query, target))]
  sp[,`:=`(c1 = match(query, genomeIDs), c2 = match(target, genomeIDs))]
  sp[,invert := c2 < c1]

  ug <- with(rbindlist(lapply(names(chrList), function(i)
    data.table(genome = i, chr = chrList[[i]]))),
    paste(genome, chr))

  # pull ref chr of all genes
  splg <- split(gffAll, by = c("genome","chr"))
  bc <- subset(blks, genome1 == refGenome & genome2 != refGenome &
                      paste(genome1, chr1) %in% ug & paste(genome2, chr2) %in% ug)
  refcg <- rbindlist(lapply(1:nrow(bc), function(i){
    x <- bc[i,]
    y2 <- splg[[sprintf("%s.%s", x$genome2, x$chr2)]]
    wh2s <- which(y2$ofID == x$firstGene2)
    wh2e <- which(y2$ofID == x$lastGene2)
    y2 <- y2[wh2s:wh2e,]
    return(data.table(refChr = x$chr1, ofID = y2$ofID))
  }))
  refcg <- subset(refcg, !duplicated(refcg))

  bl <- data.table(blks)
  bl[,`:=`(c1 = match(genome1, genomeIDs), c2 = match(genome2, genomeIDs))]
  bl <- subset(bl, c2 - c1 == 1)
  blRef <- subset(bl, genome1 == refGenome)
  blRef[,refChr := chr1]

  blAlt <- subset(bl, genome1 != refGenome)
  if(nrow(blAlt) == 0){
    blkGenes <- data.table(bl)
    blkGenes[,refChr := chr1]
    return(blkGenes)
  }else{
    blkHits <- rbindlist(lapply(1:nrow(blAlt), function(i){
      x <- blAlt[i,]
      y1 <- splg[[sprintf("%s.%s", x$genome1, x$chr1)]]
      wh1s <- which(y1$ofID == x$firstGene1)
      wh1e <- which(y1$ofID == x$lastGene1)
      y2 <- splg[[sprintf("%s.%s", x$genome2, x$chr2)]]
      wh2s <- which(y2$ofID == x$firstGene2)
      wh2e <- which(y2$ofID == x$lastGene2)
      y <- rbind(y1[wh1s:wh1e,], y2[wh2s:wh2e,])
      y[,`:=`(blkID = x$blkID, genome1 = x$genome1, genome2 = x$genome2)]
      return(y)
    }))
    blkGenes <- merge(blkHits, refcg, by = "ofID", allow.cartesian = T)

    blkGenes <- subset(blkGenes, complete.cases(blkGenes))
    blkGenes <- subset(blkGenes, !duplicated(paste(genome1, genome2, blkID, refChr)))
    refw <- blkGenes[,list(refChr = list(unique(refChr))),
                     by = c("genome1", "genome2", "blkID")]
    refout <- subset(refw, sapply(refChr, length) == 1)
    out <- merge(blks, refout, by = c("genome1", "genome2","blkID"))

    refspl <- subset(refw, sapply(refChr, length) > 1)
    setnames(refspl, "refChr", "refChrs")
    tospl <- merge(refspl, blkHits, by = c("genome1", "genome2", "blkID"))
    tospl <- merge(tospl, refcg, by = "ofID", allow.cartesian = T)
    tospl <- subset(tospl, !is.na(ofID))
    ispl <- tospl[,list(ofID = unique(ofID)), by = c("genome1","genome2","blkID","refChr")]
    ispl[,`:=`(ord = ov[ofID], start = sv[ofID], end = ev[ofID], genome = gv[ofID])]
    setkey(ispl, ord)
    spl <- split(ispl, by = c("blkID","refChr"))
    ospl <- rbindlist(lapply(spl, function(x){
      if(uniqueN(x$genome) == 2){
        y <- subset(blks, blkID == x$blkID[1] & genome1 == x$genome1[1] & genome2 == x$genome2[1])
        x1 <- subset(x, genome == genome1)
        x2 <- subset(x, genome == genome2)
        # setkey(x1, ord)
        # setkey(x2, ord)
        o <- with(y, data.table(
          chr1 = chr1, chr2 = chr2, blkID = sprintf("%s_%s", blkID, x1$refChr[1]),
          ordStart1 = min(x1$ord), ordEnd1 = max(x1$ord),
          ordStart2 = min(x2$ord), ordEnd2 = max(x2$ord),
          bpStart1 = min(x1$start), bpEnd1 = max(x1$end),
          bpStart2 = min(x2$start), bpEnd2 = max(x2$end),
          nHits = nrow(x1) + nrow(x2), orient = orient,
          firstGene1 = x1$ofID[1], firstGene2 = x2$ofID[1],
          lastGene1 = x1$ofID[nrow(x1)], lastGene2 = x2$ofID[nrow(x2)],
          genome1 = genome1, genome2 = genome2, refChr = x$refChr[1]))
        return(o)
      }
    }))
    blRef[,`:=`(c1 = NULL, c2 = NULL)]
    out[,refChr := sapply(refChr, function(x) x[1])]
    blRef[,refChr := sapply(refChr, function(x) x[1])]
    return(rbind(blRef, out, ospl))
  }
}


#' @title Plot hits
#' @description
#' \code{plot_hits} Make a pairwise hits dotplot
#' @rdname plot_genespace
#' @import data.table
#' @export
plot_hits <- function(gsParam,
                      gsAnnot,
                      synParam,
                      hitsFile = NULL,
                      hits = NULL,
                      reorderChrs = TRUE,
                      gapProp = .005,
                      useOrder = TRUE,
                      minGenes = 5,
                      axisTitleCex = .5,
                      chrLabCex = .25,
                      ptCex = .15,
                      ptBgCex = ptCex*3.5,
                      nHeatmapBins = 1000,
                      darkChrFill = "grey60",
                      lightChrFill = "grey85",
                      emptyChrFill = "grey97",
                      heatmapCols = c("green","blue4","black"),
                      blkCols = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                  "#66A61E", "#E6AB02", "#A6761D", "#666666",
                                  "darkred","darkblue"),
                      returnSourceData = FALSE){
  if(is.null(hitsFile) & is.null(hits))
    stop("must provide either a file path or a hits data.table\n")

  if(!is.null(hitsFile))
    hits <- fread(hitsFile)
  if(!any(c("chr1","chr2","genome1","genome2") %in% names(hits))){
    gff <- add_ofID2gff(read_gff(gsAnnot$gff), blastDir = gsParam$blast)
    cv <- gff$chr; gv <- gff$genome; names(cv) <- names(gv) <- gff$ofID
    ov <- gff$ord; names(ov) <- gff$ofID
    hits[,`:=`(genome1 = gv[ofID1], genome2 = gv[ofID2],
               chr1 = cv[ofID1], chr2 = cv[ofID2],
               ord1 = ov[ofID1], ord2 = ov[ofID2])]
    hits <- subset(hits, complete.cases(hits))
  }
  refGenome <- hits$genome1[1]
  genomeIDs <- c(refGenome, hits$genome2[1])
  if(uniqueN(genomeIDs) == 1){
    h1 <- factor(hits$chr1, levels = unique(hits$chr1[order(hits$ord1)]))
    h2 <- factor(hits$chr12, levels = unique(hits$chr2[order(hits$ord2)]))
    chrl <- list(names(table(h1)[table(h1) >= minGenes]))
    names(chrl) <- genomeIDs[1]
  }else{
    chrl <- match_synChrs(
      refGenome = refGenome,
      gsParam = gsParam,
      gsAnnot = gsAnnot,
      synParam = synParam,
      genomeIDs = genomeIDs,
      minGenes = minGenes,
      reorderChrs = reorderChrs)
  }

  xPos <- calc_linearGenePos(
    gsAnnot = gsAnnot,
    gsParam = gsParam,
    chrList = chrl,
    gapProp = gapProp)
  if(useOrder)
    xPosv <- xPos$xord
  if(!useOrder)
    xPosv <- xPos$xbp
  names(xPosv) <- xPos$ofID
  hits[,`:=`(x = xPosv[ofID1], y = xPosv[ofID2])]
  g1 <- hits$genome1[1]; g2 <- hits$genome2[1]
  hits <- subset(hits, complete.cases(hits))

  # get chr bound/mean coordinates
  chrBnd <- calc_chrBounds(
    hits = hits,
    genome1 = g1,
    genome2 = g2)
  chrBnd[[1]][,mini := min(chrStart1) - 1]
  chrBnd[[2]][,mini := min(chrStart2) - 1]
  chrBnd[[1]][,`:=`(chrStart1 = chrStart1 - mini,
                    chrEnd1 = chrEnd1 - mini,
                    chrMid1 = chrMid1 - mini)]
  chrBnd[[2]][,`:=`(chrStart2 = chrStart2 - mini,
                    chrEnd2 = chrEnd2 - mini,
                    chrMid2 = chrMid2 - mini)]
  # get chr combination background colors
  chrBgCol <- color_chrByHits(
    chrList = chrl,
    hits = hits,
    genome1 = g1,
    genome2 = g2,
    darkChrFill = darkChrFill,
    lightChrFill = lightChrFill,
    emptyChrFill= emptyChrFill)
  hits[,`:=`(x = x - min(x) + 1,
             y = y - min(y) + 1)]
  make_dotplot(
    chrBounds = chrBnd,
    chrBg = chrBgCol,
    hits = hits,
    genome1 = g1, genome2 = g2,
    axisTitleCex = axisTitleCex,
    chrLabCex = chrLabCex,
    ptCex = ptCex,
    useOrder = useOrder,
    ptBgCex = ptBgCex,
    blkCols = blkCols,
    gapProp = gapProp,
    heatmapCols = heatmapCols,
    nbinsx = nHeatmapBins)
  if(returnSourceData)
    return(hits)
}

#' @title calculate chromosome bounds from hits
#' @description
#' \code{calc_chrBounds} calculate chromosome bounds from hits
#' @rdname plot_genespace
#' @import data.table
#' @export
calc_chrBounds <- function(hits, genome1, genome2){
  chrMd1 <- hits[,list(
    chrStart1 = min(x), chrEnd1 = max(x), chrMid1 = (max(x) + min(x))/2),
    by = "chr1"]
  chrMd2 <- hits[,list(
    chrStart2 = min(y), chrEnd2 = max(y), chrMid2 = (max(y) + min(y))/2),
    by = "chr2"]
  out <- list(chrBound1 = chrMd1, chrBound2 = chrMd2)
  names(out) <- c(genome1, genome2)
  return(out)
}

#' @title get colors for chromosome backgrounds
#' @description
#' \code{color_chrByHits}  get colors for chromosome backgrounds
#' @rdname plot_genespace
#' @import data.table
#' @export
color_chrByHits <- function(chrList,
                            hits,
                            genome1,
                            genome2,
                            darkChrFill,
                            lightChrFill,
                            emptyChrFill){
  tmp <- data.table(CJ(chr1 = chrList[[genome1]],
                       chr2 = chrList[[genome2]]))
  hitsCnt <- hits[,list(n = uniqueN(c(ofID1, ofID2))), by = c("chr1", "chr2")]
  hitsCnt <- merge(hitsCnt, tmp, all = T, by = c("chr1", "chr2"))
  hitsCnt <- subset(hitsCnt, !duplicated(paste(chr1, chr2)))
  hitsCnt$n[is.na(hitsCnt$n)] <- 0

  # add colors to background
  if(min(hitsCnt$n) == 0){
    col20 <- rev(colorRampPalette(colors = c(darkChrFill,lightChrFill))(21))
    col20 <- c(col20, emptyChrFill)
    names(col20) <- as.character(0:20)
  }else{
    col20 <- rev(colorRampPalette(colors = c(darkChrFill,lightChrFill))(21))
    names(col20) <- as.character(1:20)
  }
  if(uniqueN(hits$chr1) > uniqueN(hits$chr2)){
    hitsCnt[,mmax := max(n), by = "chr2"]
  }else{
    hitsCnt[,mmax := max(n), by = "chr1"]
  }
  hitsCnt[,propMax := n / mmax]
  hitsCnt[,propclass := ceiling(propMax * 20)]
  hitsCnt[,col := col20[as.character(propclass)]]
  hitsCnt$col[hitsCnt$n == 0] <- emptyChrFill
  return(hitsCnt)
}

#' @title make_dotplot
#' @description
#' \code{make_dotplot}  make_dotplot
#' @rdname plot_genespace
#' @import data.table
#' @export
make_dotplot <- function(chrBounds,
                         chrBg,
                         hits,
                         genome1,
                         genome2,
                         gapProp,
                         axisTitleCex,
                         chrLabCex,
                         ptCex = ptCex,
                         useOrder,
                         ptBgCex,
                         blkCols,
                         nbinsx,
                         heatmapCols){

  xoffset <- min(hits$x) - (diff(range(hits$x))/20)
  yoffset <- min(hits$y) - (diff(range(hits$y))/20)
  plot(
    NA, NA,
    xlim = c(xoffset, max(hits$x)),
    ylim = c(yoffset, max(hits$y)),
    type = "n", axes = F,
    xlab = "", ylab = "")

  if(useOrder){
    title(
      xlab = paste(genome1,"chromosomes (gene rank order)"),
      ylab = paste(genome2,"chromosomes (gene rank order)"),
      line = 0, cex.lab = axisTitleCex)
  }else{
    title(
      xlab = paste(genome1,"chromosomes (physical gene position)"),
      ylab = paste(genome2,"chromosomes (physical gene position)"),
      line = 0, cex.lab = axisTitleCex)
  }

  cb <- merge(chrBg, chrBounds[[1]], by = "chr1", all.x = T)
  cb <- merge(cb, chrBounds[[2]], by = "chr2", all.x = T)
  with(cb, rect(
    xleft = chrStart1, xright = chrEnd1,
    ybottom = chrStart2, ytop = chrEnd2,
    col = col, border = NA))

  with(chrBounds[[1]],
       text(
         x = chrMid1, y = 0,
         labels = chr1,
         cex = chrLabCex, adj = c(1.05,.5), srt = 90))
  with(chrBounds[[2]],
       text(
         x = 0, cex = chrLabCex, labels = chr2,
         y = chrMid2,  adj = c(1.05,.5)))

  if(!"blkID" %in% colnames(hits)){
    colScale <- colorRampPalette(heatmapCols)
    binSize <- max(hits$x)/nbinsx
    hits[,`:=`(rndx = ceiling(x/binSize)*binSize, rndy = ceiling(y/binSize)*binSize)]
    cnts <- hits[,list(n = uniqueN(c(ofID1, ofID2))),
                 by = c("rndx", "rndy")]
    levs <- 1:max(cnts$n)
    colPal <- colScale(length(levs))
    cnts[,col := colPal[n]]
    setkey(cnts, n)
    with(cnts, points(
      x = rndx, y = rndy,
      pch = 19, cex = ptBgCex, col = "white", lwd = 0))
    with(cnts, points(
      x = rndx, y = rndy,
      pch = 19, cex = ptCex*2, col = col, lwd = .1))
  }else{
    pal <-  rep(blkCols, uniqueN(hits$blkID))
    cols <- pal[1:uniqueN(hits$blkID)]
    names(cols) <- unique(hits$blkID)
    hits[,col := cols[blkID]]
    hits[,n := uniqueN(c(ofID1, ofID2)), by = "blkID"]
    setorder(hits, n)
    with(hits, points(
      x = x, y = y,
      pch = 19, cex = ptBgCex, col = "white", lwd = 0))
    with(hits, points(
      x = x, y = y,
      pch = 19, cex = ptCex, col = col, lwd = .1))
  }
}

#' @title calc_linearGenePos
#' @description
#' \code{calc_linearGenePos} calc_linearGenePos
#' @rdname plot_utils
#' @export
calc_linearGenePos <- function(chrList,
                               gsParam,
                               gsAnnot,
                               gapProp){
  chrdt <- rbindlist(lapply(names(chrList), function(i)
    data.table(genome = i, chr = chrList[[i]], chrOrd = 1:length(chrList[[i]]))))
  gff <- merge(
    chrdt,
    add_ofID2gff(read_gff(gsAnnot$gff[names(chrList)]), gsParam$blast),
    by = c("genome", "chr"))
  gff[,y := as.numeric(factor(genome, levels = names(chrList)))]
  gff <- subset(gff, complete.cases(gff))
  setkey(gff, y, chrOrd, ord)
  gff[,`:=`(ord = 1:.N, bp = (start + end)/2), by = c("genome","chr")]
  gffc <- gff[,list(maxBp = max(end), maxOrd = max(ord)), by = c("genome","chrOrd")]
  gffc[,`:=`(totBp = sum(as.numeric(maxBp)), totOrd = sum(as.numeric(maxOrd))), by = c("genome")]
  gffc[,`:=`(longestBp = max(totBp), longestOrd = max(totOrd))]
  gffc[,`:=`(pLongestBp = totBp/longestBp, pLongestOrd = totOrd/longestOrd)]
  minGapBp <- with(gffc, min(ceiling(totBp * gapProp)))
  minGapOrd <- with(gffc, min(ceiling(totOrd * gapProp)))
  gffc[,`:=`(gapSizeBp = ceiling(minGapBp/pLongestBp),
             gapSizeOrd = ceiling(minGapOrd/pLongestOrd),
             pLongestBp = NULL, pLongestOrd = NULL)]
  gffc[,`:=`(totGapBp = gapSizeBp*.N,
             totGapOrd = gapSizeOrd*.N), by = "genome"]
  gffc[,`:=`(chrStartBp = floor(-(totBp+totGapBp)/2),
             chrStartOrd = floor(-(totOrd+totGapOrd)/2))]
  gffc[,`:=`(cumBp = c(1, cumsum(maxBp[1:(.N-1)]+gapSizeBp[-1])),
             cumOrd = c(1, cumsum(maxOrd[1:(.N-1)]+gapSizeOrd[-1]))), by = "genome"]
  gffc[,`:=`(xStartBp = cumBp + chrStartBp,
             xStartOrd = cumOrd + chrStartOrd), by = "genome"]
  gffc <- gffc[,c("genome","chrOrd","xStartBp", "xStartOrd")]
  gffo <- merge(
    gff[,c("genome","ofID","id","y","chrOrd","chr","bp","ord")],
    gffc,  by = c("genome","chrOrd"))
  gffo[,`:=`(xbp = xStartBp+bp, xord = xStartOrd + ord,
             ord = NULL, bp = NULL, xStartBp = NULL, xStartOrd = NULL)]
  return(gffo)
}


#' @title match_synChrs
#' @description
#' \code{match_synChrs} match_synChrs
#' @rdname plot_utils
#' @export
match_synChrs <- function(refGenome,
                          synParam,
                          minGenes = 5,
                          gsParam,
                          gsAnnot,
                          reorderChrs = TRUE,
                          genomeIDs = NULL){

  # get genomeIDs ironed out
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]
  if(!refGenome %in% genomeIDs)
    refGenome <- genomeIDs[1]
  genomeIDs <- c(refGenome, genomeIDs[genomeIDs != refGenome])

  # read gff and parse to just chrs that are big enough
  gff <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), blastDir = gsParam$blast)
  gff[,nGenesChr := .N, by = c("genome", "chr")]
  gff <- subset(gff, nGenesChr >= minGenes)
  gu <- with(gff, unique(paste(genome, chr)))
  if(!reorderChrs){
    gff[,chrn := as.numeric(gsub("[^0-9]", "", chr))]
    setorder(gff, chrn, chr, na.last = T)
    spl <- split(gff, by = "genome")
    chrList <- sapply(spl, USE.NAMES = T, simplify = F, function(x) unique(x$chr))
  }else{

    refChrs <- subset(gff, genome == refGenome)
    refChrs[,chrn := as.numeric(gsub("[^0-9]", "", chr))]
    setorder(refChrs, chrn, chr, na.last = T)
    tmp <- unique(refChrs$chr); refChrOrd <- 1:length(tmp); names(refChrOrd) <- tmp

    ov <- gff$ord; names(ov) <- gff$ofID
    cv <- gff$chr; names(cv) <- gff$ofID
    sv <- gff$start; names(sv) <- gff$ofID
    ev <- gff$end; names(ev) <- gff$ofID

    sp <- subset(synParam, runBlast)
    sp[,sogf := file.path(gsParam$synteny, sprintf("%s_vs_%s_synog.txt.gz", query, target))]
    sp[,`:=`(c1 = match(query, genomeIDs), c2 = match(target, genomeIDs))]
    sp[,invert := c2 < c1]

    # metadata for hits against the reference
    spRef <- subset(sp, (c1 == 1 | c2 == 1) & c1 != c2)

    allChrGenes <- rbindlist(lapply(1:nrow(spRef), function(i){
      if(!spRef$invert[i]){
        d <- fread(spRef$sogf[i])
        d[,genome2 := spRef$genome2[i]]
      }else{
        d <- fread(spRef$sogf[i], select = c(2,1,3), col.names = c("ofID1", "ofID2", "blkID"))
        d[,genome2 := spRef$genome1[i]]
      }
      return(d)
    }))
    allChrGenes[,`:=`(chr1 = cv[ofID1], chr2 = cv[ofID2])]
    allChrGenes[,chr1Ord := refChrOrd[chr1]]
    allChrGenes <- subset(allChrGenes, chr1 %in% names(refChrOrd) &
                            paste(genome2, chr2) %in% gu)
    allChrGenes[,chrn := as.numeric(gsub("[^0-9]", "", chr2))]
    allChrGenes[,chrord := as.numeric(as.factor(chr2))]
    altChrOrd <- allChrGenes[,list(medOrd = median(chr1Ord)),
                             by = c("genome2","chr2","chrn","chrord")]
    setkey(altChrOrd, genome2, medOrd, chrn, chrord)
    out <- rbind(data.table(genome = refGenome, chr = names(refChrOrd)),
                 with(altChrOrd, data.table(genome = genome2, chr = chr2)))
    chrList <- with(out, split(chr, genome))
  }
  return(chrList[genomeIDs])
}

#' @title calc_chrOffset
#' @description
#' \code{calc_chrOffset} calc_chrOffset
#' @rdname plot_utils
#' @export
calc_chrOffset <- function(gsAnnot,
                           chrList,
                           useOrder,
                           gapProp){
  genomeIDs <- names(chrl)
  gff <- read_gff(gsAnnot$gff[genomeIDs])
  gff <- rbindlist(lapply(split(gff, by = "genome"), function(x){
    x[,chrn := as.numeric(factor(chr, levels = chrList[[x$genome[1]]]))]
    x <- subset(x, !is.na(chrn))
    setkey(x, chrn, start, end)
    return(x)
  }))
  if(useOrder){
    gff <- gff[,list(start = min(ord), end = max(ord)),
               by = c("genome", "chr")]
  }else{
    gff <- gff[,list(start = min(start), end = max(end)),
               by = c("genome", "chr")]
  }
  gff[,`:=`(eo = end, end = end - start)]
  gff[,totLen := sum(end), by = "genome"]
  gff[,maxLen := max(totLen)]
  gff[,propOfLongest := totLen/maxLen]
  minGap <- with(gff, min(ceiling(totLen * gapProp)))
  gff[,gapSize := ceiling(minGap/propOfLongest)]
  gff[,totGap := gapSize*.N, by = "genome"]
  gff[,chrStart := floor(-(totLen+totGap)/2)]
  gff[,offset := chrStart + c(0, (cumsum(gapSize) + cumsum(end))[1:(.N - 1)]),
      by = "genome"]
  gff[,y := as.numeric(factor(genome, levels = names(chrList)))]
  gff[,end := eo]
  return(gff[,c("genome", "chr", "y","start","end","offset")])
}

#' @title calc_linearBlkCoords
#' @description
#' \code{calc_linearBlkCoords} calc_linearBlkCoords
#' @rdname plot_utils
#' @export
calc_linearBlkCoords <- function(chrOffset,
                                 gsParam,
                                 gsAnnot,
                                 chrList,
                                 blks,
                                 useOrder){
  blk <- data.table(blks)
  cvo <- chrOffset$offset; names(cvo) <- with(chrOffset, paste(genome, chr))
  so <-
    yo <- chrOffset$y; names(yo) <- with(chrOffset, genome); yo <- yo[!duplicated(names(yo))]
  if(useOrder){
    blk[,`:=`(start1 = bpStart1 + cvo[paste(genome1, chr1)],
              end1 = bpEnd1 + cvo[paste(genome1, chr1)],
              start2 = bpStart2 + cvo[paste(genome2, chr2)],
              end2 = bpEnd1 + cvo[paste(genome2, chr2)])]
  }else{
    blk[,`:=`(start1 = bpStart1 + cvo[paste(genome1, chr1)],
              end1 = bpEnd1 + cvo[paste(genome1, chr1)],
              start2 = bpStart2 + cvo[paste(genome2, chr2)],
              end2 = bpEnd1 + cvo[paste(genome2, chr2)])]
  }
  blk[,`:=`(y1 = yo[genome1], y2 = yo[genome2])]
  blk <- blk[,c("genome1","genome2","y1","y2","chr1","chr2",
                "start1","end1","start2","end2","orient")]
  blk <- subset(blk, complete.cases(blk))

  if(any(w$orient == "-")){
    tmp <- subset(blk, orient == "-")
    blk <- subset(blk, orient == "+")
    tmp[,`:=`(tm = end2, end2 = start2)]
    tmp[,`:=`(start2 = tm, tm = NULL)]
    blk <- rbind(blk, tmp)
  }
  setorder(blk, y1, y2, start1, end1, start2, end2)
  return(w)
}


#' @title ...
#' @description
#' \code{round_rect} ...
#' @rdname utils
#' @importFrom graphics par
#' @importFrom grDevices dev.size
#' @export
round_rect <- function(xleft, ybottom, xright, ytop){

  if (ytop <= ybottom)
    stop("ytop must be > ybottom")
  if (xleft >= xright)
    stop("xleft must be < xright")

  # measure graphics device
  asp <- diff(par("usr")[3:4]) / diff(par("usr")[1:2])
  dev <- dev.size()[1] / dev.size()[2]

  # make a curve and split into left and right
  radius <- (ytop - ybottom) / 2
  centerY <- ytop - radius
  centerX <- mean(c(xleft, xright))
  theta <- seq(0, 2 * pi, length = 200)
  circX <- cos(theta)
  circY <- sin(theta)
  leftC <- which(circX <= 0)
  rightC <- which(circX >= 0)

  xR <- circX[rightC]
  yR <- circY[rightC]
  ordYR <- rev(order(yR))
  xR <- xR[ordYR]
  yR <- yR[ordYR]

  xL <- circX[leftC]
  yL <- circY[leftC]
  ordYL <- order(yL)
  xL <- xL[ordYL]
  yL <- yL[ordYL]

  # project onto graphics device and scale
  xRightS <- xright - (radius / asp / dev)
  xLeftS <- xleft + (radius / asp / dev)
  if (centerX < xLeftS)
    xLeftS <- centerX
  if (centerX > xRightS)
    xRightS <- centerX
  xLS <- scale_between(xL, xleft, xLeftS)
  xRS <- scale_between(xR, xRightS, xright)
  yLS <- scale_between(yL, ybottom, ytop)
  yRS <- scale_between(yR, ybottom, ytop)
  return(data.table(x = c(xRS,xLS), y = c(yRS, yLS)))
}


#' @title ...
#' @description
#' \code{calc_curvePolygon} ...
#' @rdname utils
#' @export
calc_curvePolygon <- function(start1,
                              end1 = NULL,
                              start2,
                              end2 = NULL,
                              y1,
                              y2){
  scaledCurve <- cosine_points()
  if (!is.null(end1) | !is.null(end2)) {
    tp <- rbind(
      start1 = data.table(x = start1, y = y1),
      poly1 = data.table(
        x = scale_between(
          x = scaledCurve[,1],
          min = start1,
          max = start2),
        y = scale_between(
          x = scaledCurve[,2],
          min = y1,
          max = y2)),
      start2 = data.table(
        x = start2,
        y = y2),
      end2 = data.table(
        x = end2,
        y = y2),
      poly2 = data.table(
        x = scale_between(
          x = scaledCurve[,1],
          min = end2,
          max = end1),
        y = scale_between(
          x = scaledCurve[,2],
          min = y2,
          max = y1)),
      end1 = data.table(
        x = end1,
        y = y1))
  }else{
    tp <- data.table(
      x = scale_between(
        x = scaledCurve[,1],
        min = start1,
        max = start2),
      y = scale_between(
        x = scaledCurve[,2],
        min = y1,
        max = y2))
  }
  return(tp)
}

#' @title ...
#' @description
#' \code{scale_between} ...
#' @rdname utils
#' @export
scale_between <- function(x, min, max, scale1toMean = TRUE){
  if(length(unique(x)) > 1){
    return((x - min(x)) / (max(x) - min(x)) * (max - min) + min)
  }else{
    if(scale1toMean){
      return(mean(c(min, max)))
    }else{
      return(max)
    }
  }
}

#' @title ...
#' @description
#' \code{cosine_points} ...
#' @rdname utils
#' @export
cosine_points <- function(){
  npts = 1e4 # initial number of points
  keepat = round(npts / 20) # grid to keep always
  grid <- seq(from = 0, to = pi, length.out = npts) # grid
  x <- (1 - cos(grid)) / max((1 - cos(grid))) # scaled cosine
  y <- grid / max(grid) # scaled grid
  # calculate slope for each point
  x1 <- x[-1];  y1 <- y[-1]
  x2 <- x[-length(x)];  y2 <- y[-length(y)]
  s <-  (y1 - y2) / (x1 - x2)
  # choose points that capture changes in slope
  ds <- cumsum(abs(diff(s)))*5
  wh <- c(1,which(!duplicated(round(ds))), length(x))
  wh2 <- c(wh, seq(from = 0, to = length(x), by = round(keepat)))
  wh <- c(wh, wh2)[!duplicated(c(wh, wh2))]
  wh <- wh[order(wh)]
  return(cbind(x[wh], y[wh]))
}
