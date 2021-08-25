#' @title Genespace plotting routines
#' @description
#' \code{plot_pav} Genespace plotting routines
#' @name plot_pav
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param genomeIDs character vector specifying the genomeIDs to plot
#' @param maxComb2plot integer, maximum number of combinations to plot
#' @param excudePrivate logical, should private PAV be ignored?
#' @param gff data.table of gff-like annotations
#' @param genomeLabCex character expansion genome labels
#' @param ptCex character expansion of points
#' @param leftBuffer area to left of plot to be left empty
#' @param barWidth width of bars
#' @param nCols number of columns to plot
#' @param nLabCex character expansion of number label
#' @param ylabCex character expansion of y axis label
#' @param xlabCex character expansion of x axis label
#' @param cols vector of colors
#' @details ...
#'
#' @note \code{plot_pav} is a generic name for the functions documented.
#' \cr
#' If called, \code{plot_pav} returns its own arguments.
#'
#' @title plot_pav plot
#' @description
#' \code{plot_pav} plot_pav
#' @rdname plot_pav
#' @import data.table
#' @export
plot_pav <- function(gsParam,
                     genomeIDs = NULL,
                     maxComb2plot = 50,
                     excudePrivate = F){

  ##############################################################################
  # 0. Initial setup / checking
  ##############################################################################
  genome <- NULL
  if(is.null(genomeIDs)){
    genomeIDs <- gsParam$genomes$genomeIDs
    genomeIDs <- genomeIDs[!genomeIDs %in% gsParam$genomes$outgroup]
  }else{
    tmp <- gsParam$genomes$genomeIDs
    tmp <- tmp[!tmp %in% gsParam$genomes$outgroup]
    genomeIDs <- genomeIDs[genomeIDs %in% tmp]
  }

  ##############################################################################
  # 1. Read gff
  ##############################################################################
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- fread(gffFile)
  gff <- subset(gff, genome %in% genomeIDs)
  if(!"synOg" %in% colnames(gff))
    stop("Not seeing syntenic OGs in gff - has synteny been run yet?\n")

  ##############################################################################
  # 2. plot the overall orthogroup type counts
  ##############################################################################
  ogcnt <- plot_ogCnts(
    gff = gff,
    genomeIDs = genomeIDs,
    ylabCex = .75,
    xlabCex = .75)

  ##############################################################################
  # 3. plot unique pav
  ##############################################################################
  l = length(genomeIDs)
  ptc <- ifelse(l > 20, .5, ifelse(l > 15, .7, ifelse(l > 10, 1, ifelse(l > 5, 1.5, 1.8))))
  pavcnt <- plot_pavCnts(
    gff = gff,
    ptCex = ptc,
    genomeIDs = genomeIDs,
    maxComb2plot = maxComb2plot,
    excudePrivate = excudePrivate)

}

#' @title plot_pavCnts plot
#' @description
#' \code{plot_pavCnts} plot_pavCnts
#' @rdname plot_pav
#' @import data.table
#' @importFrom graphics title rect text segments points
#' @export
plot_pavCnts <- function(gff,
                         genomeIDs,
                         excudePrivate = F,
                         genomeLabCex = .75,
                         ptCex = 1.5,
                         leftBuffer = .25,
                         barWidth = .8,
                         nCols = 1,
                         nLabCex = .5,
                         maxComb2plot = 20){
  ng <- n <- genome <- x <- scl <- type <- NULL
  g <- data.table(gff)
  gi <- unique(g$genome)
  if(is.null(genomeIDs))
    genomeIDs <- gi
  gid <- data.table(genome = genomeIDs)
  gid[,`:=`(ord = 1:.N,
            y = scale_between(1:.N, min = -.1, max = -.9))]
  if(excudePrivate){
    g[,ng := uniqueN(genome), by = "synOg"]
    g <- subset(g, !ng %in% c(1, length(gi)))
  }
  d <- dcast(
    g, synOg ~ genome, value.var = "ofID",
    fun.aggregate = function(x) ifelse(length(x) == 0, 0, 1))
  dc <- d[,list(n = .N), by =gi]
  setorder(dc, -n)
  dc <- dc[1:maxComb2plot, ]
  sclSt <- ceiling(nrow(dc) * leftBuffer)
  dc[,x := sclSt:(.N + sclSt - 1)]
  mi <- min(dc$n) / max(dc$n)
  dc[,scl := scale_between(n, min = mi, max = 1)]
  nk <- round(max(dc$n)/3, -2)
  scl10k <- scale_between(c(nk, dc$n), min = mi, max = 1)[1]
  r <- range(c(gid$y, dc$scl))
  r[1] <- r[1] - .1; r[2] <- r[2] + .1
  plot(
    NA, NA,
    xlim = c(1-barWidth, max(dc$x)+barWidth),
    ylim = r,
    xlab = "",
    bty = "n", ylab = "", axes = F)
  text(x = min(dc$x)-barWidth, y = 0,
       label = "n syn. ogs", srt = 90, cex = genomeLabCex, adj = c(-.125,.25))
  segments(
    x0 = sclSt - (barWidth/2), x1 = max(dc$x) + (barWidth/2), y0 = scl10k, y1 = scl10k,
    lty = 3, col = "black")
  text(
    y = scl10k, x = max(dc$x),
    sprintf("%s genes",nk),
    adj = c(1.5, -.25))
  for(i in 1:nrow(dc)){
    with(dc[i,], rect(
      xleft = x-(barWidth/2), xright = x+(barWidth/2), ybottom = 0, ytop = scl,
      col = "lightgrey"))
    with(dc[i,], text(
      x = x-(barWidth/2), y = scl, labels = n, cex = nLabCex, adj = c(0,-.25)))
  }
  for(i in 1:nrow(gid)){
    with(gid[i, ], text(
      x = sclSt-(barWidth/2)-.1, y = y, labels = genome,
      adj = c(1,.5), cex = genomeLabCex))
    points(
      dc$x, rep(gid$y[i], nrow(dc)),
      pch = ifelse(dc[[gid$genome[i]]] == 1, 16, 1),
      cex = ptCex)
  }
  dco <- data.table(dc, type = "syntenic")

  g <- data.table(gff)
  gi <- unique(g$genome)
  if(is.null(genomeIDs))
    genomeIDs <- gi
  gid <- data.table(genome = genomeIDs)
  gid[,`:=`(ord = 1:.N,
            y = scale_between(1:.N, min = -.1, max = -.9))]
  if(excudePrivate){
    g[,ng := uniqueN(genome), by = "og"]
    g <- subset(g, !ng %in% c(1, length(gi)))
  }
  d <- dcast(
    g, og ~ genome, value.var = "ofID",
    fun.aggregate = function(x) ifelse(length(x) == 0, 0, 1))
  dc <- d[,list(n = .N), by =gi]
  setorder(dc, -n)
  dc <- dc[1:maxComb2plot, ]
  sclSt <- ceiling(nrow(dc) * leftBuffer)
  dc[,x := sclSt:(.N + sclSt - 1)]
  mi <- min(dc$n) / max(dc$n)
  dc[,scl := scale_between(n, min = mi, max = 1)]
  nk <- round(max(dc$n)/3, -2)
  scl10k <- scale_between(c(nk, dc$n), min = mi, max = 1)[1]
  r <- range(c(gid$y, dc$scl))
  r[1] <- r[1] - .1; r[2] <- r[2] + .1
  plot(
    NA, NA,
    xlim = c(1-barWidth, max(dc$x)+barWidth),
    ylim = r,
    bty = "n", ylab = "", xlab = "", axes = F)
  text(x = min(dc$x)-barWidth, y = 0,
       label = "n glob. ogs", srt = 90, cex = genomeLabCex, adj = c(-.125,.25))
  segments(
    x0 = sclSt - (barWidth/2), x1 = max(dc$x) + (barWidth/2), y0 = scl10k, y1 = scl10k,
    lty = 3, col = "black")
  text(
    y = scl10k, x = max(dc$x),
    sprintf("%s genes",nk),
    adj = c(1.5, -.25))
  for(i in 1:nrow(dc)){
    with(dc[i,], rect(
      xleft = x-(barWidth/2), xright = x+(barWidth/2), ybottom = 0, ytop = scl,
      col = "lightgrey"))
    with(dc[i,], text(
      x = x-(barWidth/2), y = scl, labels = n, cex = nLabCex, adj = c(0,-.25)))
  }
  for(i in 1:nrow(gid)){
    with(gid[i, ], text(
      x = sclSt-(barWidth/2)-.1, y = y, labels = genome,
      adj = c(1,.5), cex = genomeLabCex))
    points(
      dc$x, rep(gid$y[i], nrow(dc)),
      pch = ifelse(dc[[gid$genome[i]]] == 1, 16, 1),
      cex = ptCex)
  }
  dc[,type := "global"]
  # par(mfrow = pm)
  # par(mar = pmar)
  return(rbind(dco, dc))
}

#' @title plot_ogCnts plot
#' @description
#' \code{plot_ogCnts} plot_ogCnts
#' @rdname plot_pav
#' @import data.table
#' @importFrom graphics title barplot
#' @export
plot_ogCnts <- function(gff,
                        genomeIDs,
                        ylabCex = .75,
                        xlabCex = .75,
                        cols = c("lightgoldenrod1", "#1FA187FF", "#365C8DFF", "#443A83FF")){
  genome <- ngenome <- og <- synOg <- NULL
  g <- subset(gff, genome %in% genomeIDs)
  g[,genome := factor(genome, levels = genomeIDs)]
  g[,ngenome := uniqueN(genome), by = "og"]
  gc <- subset(g, !duplicated(paste(og, genome)))
  gc <- gc[,list(ng = .N), by = c("genome", "ngenome")]
  dc <- dcast(gc, ngenome ~ genome, value.var = "ng")
  colr <- colorRampPalette(cols)(nrow(dc))
  setorder(dc, -ngenome)
  bp <- barplot(
    data.matrix(dc[,-1])/1000,
    col = colr,
    xlab = "n global orthogroups (x1000)",
    sub = "yellow = complete, darkblue = private",
    ylab = "", axes = T, axisnames = F,
    main = "",
    horiz = T,
    las = 2,
    cex.axis = ylabCex,
    cex.names = xlabCex)
  text(
    0, bp, labels = colnames(dc)[-1],
    srt = 0, adj = c(1.1, .5), xpd = TRUE, cex = ylabCex)

  dcs <- data.table(dc)

  g[,ngenome := uniqueN(genome), by = "synOg"]
  gc <- subset(g, !duplicated(paste(synOg, genome)))
  gc <- gc[,list(ng = .N), by = c("genome", "ngenome")]
  dc <- dcast(gc, ngenome ~ genome, value.var = "ng")
  colr <- colorRampPalette(cols)(nrow(dc))
  setorder(dc, -ngenome)
  barplot(
    data.matrix(dc[,-1])/1000,
    col = colr,
    xlab = "n syntenic orthogroups (x1000)",
    sub = "yellow = complete, darkblue = private",
    ylab = "", axes = T, axisnames = F,
    main = "",
    horiz = T,
    las = 2,
    cex.axis = ylabCex,
    cex.names = xlabCex)
  text(
    0, bp, labels = colnames(dc)[-1],
    srt = 0, adj = c(1.1,0.5), xpd = TRUE, cex = ylabCex)
  return(dc)
}

