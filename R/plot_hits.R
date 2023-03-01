#' @title plot hits as a xy dotplot
#' @description
#' \code{plot_hits} routines to make visually appealing dotplots
#' @name plot_hits
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param type character string of "all", "raw", or "syntenic" specifying which
#' type of dotplot to generate.
#' @param verbose logical, should updates be printed to the console?
#' @param hits data.table containg hits. See read_allBlast.
#' @param outDir file.path where pdf should be written
#' @param minGenes2plot integer specifying the minimum number of genes that can
#' be plotted
#' @param appendName character with text to append to the file name
#' @param dotsPerIn integer specifying how fine-scaled the heatmap is
#' @param quantileThresh integer specifying the top quantile to be thresholded
#' @param plotSize numeric smalled dimension of the plot
#' @param minScore numeric the minimum score to be permitted in the first plot
#' @param maxFacets integer the maximum number of facets to plot (doesn't plot
#' for genomes with lots of small chrs/scaffolds)
#' @param colorByBlks logical, should blocks be colored?
#' @param alpha numeric [0-1], specifying the transparency of the points
#' @param useOrder logical, should gene order or bp position be used?
#' @param minScore numeric, the minimum scoring hit to plot
#' @param minGenes2plot integer, the minimum number of hits to plot a chromosome
#' combination.
#' @param print2currentDevice logical, should the plot be printed to a pdf or
#' to the current device?
#'
#' \cr
#' If called, \code{plot_hits} returns its own arguments.
#'
#' @details Dotplots here aggregate across proximate positions to reduce file
#' size, especially for very large genomes. XY positions are always in gene-rank
#' order positions. Graphics are built with ggplot2.

#' @title Plot syntenic hits
#' @description
#' \code{plot_hits} The pipeline to plot syntenic hits in parallel
#' @rdname plot_hits
#' @import data.table
#' @import R.utils
#' @importFrom dbscan dbscan frNN
#' @importFrom parallel mclapply
#' @export
plot_hits <- function(gsParam,
                      verbose = TRUE,
                      type,
                      dotsPerIn = 256,
                      quantileThresh = .5,
                      plotSize = 12,
                      minScore = 50){
  ##############################################################################
  # 1. setup
  # -- 1.1 get env vars set up
  query <- target <- lab <- nRegionHits <- nRegions <- nAnchorHits <- nBlks <-
    nSVs <- selfOnly <- queryPloidy <- targetPloidy <- nGlobOgHits <- synHits <-
    nTotalHits <- chunk <- inBuffer <- NULL

  if(!"synteny" %in% names(gsParam))
    gsParam <- set_syntenyParams(gsParam)

  nCores <- gsParam$params$nCores

  # -- 1.2 check that files are all ok
  if(!"synteny" %in% names(gsParam))
    stop("must run set_syntenyParams prior to synteny")

  if(!all(file.exists(gsParam$synteny$blast$allBlast)) && type %in% c("all", "raw"))
    stop("some annotated blast files dont exist, run annotate_blast() first\n")
  if(!all(file.exists(gsParam$synteny$blast$synHits)) && type %in% c("all", "syntenic"))
    stop("some syntenic blast files dont exist, run synteny() first\n")
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
  blMdOut <- lapply(1:length(synMdSpl), function(chnki){

    chnk <- data.table(synMdSpl[[chnki]])

    ############################################################################
    # -- loop through each row in each chunk
    outChnk <- mclapply(1:nrow(chnk), mc.cores = nCores, function(i){

      # -- 2.1 read in the metadata and hits
      outMd <- data.table(chnk[i,])
      x <- data.table(outMd)
      rawHits <- read_allBlast(x$allBlast)

      l1 <- mean(table(rawHits$chr1[rawHits$sameOG]))/10
      l2 <- mean(table(rawHits$chr2[rawHits$sameOG]))/10

      dps <- gsParam$params$dotplots
      if(dps == "check"){
        ggdotplot(
          hits = data.table(rawHits),
          outDir = gsParam$paths$dotplots,
          minGenes2plot = min(c(l1, l2)),
          maxFacets = 10000,
          type = type,
          dotsPerIn = dotsPerIn,
          quantileThresh = quantileThresh,
          plotSize = plotSize,
          minScore = minScore)
      }else{
        if(dps == "always"){
          ggdotplot(
            hits = data.table(rawHits),
            outDir = gsParam$paths$dotplots,
            minGenes2plot = min(c(l1, l2)),
            maxFacets = Inf,
            type = type,
            dotsPerIn = dotsPerIn,
            quantileThresh = quantileThresh,
            plotSize = plotSize,
            minScore = minScore)
        }
      }
    })
  })
  return(gsParam)
}



#' @title make dotplots of syntenic hits
#' @description
#' \code{ggdotplot} ggplot2 integrated graphics to produce dotplots
#' @rdname plot_hits
#' @import data.table
#' @import ggplot2
#' @importFrom grDevices pdf dev.off rgb
#' @importFrom dbscan dbscan frNN
#' @export
ggdotplot <- function(hits,
                      type,
                      outDir = NULL,
                      minGenes2plot = 100,
                      appendName = "synHits",
                      dotsPerIn = 256,
                      quantileThresh = .5,
                      plotSize = 12,
                      minScore = 50,
                      maxFacets = 10000,
                      verbose = is.null(outDir),
                      print2currentDevice = FALSE){
  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- blkID <-
    inBuffer <- rnd2 <- rnd1 <- n <- isArrayRep2 <- isArrayRep1 <- chr1 <-
    noAnchor <- bitScore <- quantile <- chr2 <- sameOG <- isAnchor <- NULL

  ##############################################################################
  # 1. Get the plot size figured out
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
  xrnd2 <- floor(x / totDots)+1

  ordPerIn <- y / dotsPerIn
  totDots <- ht * dotsPerIn
  yrnd2 <- floor(y / totDots)+1


  tp[,`:=`(rnd1 = round_toInteger(ord1, xrnd2),
           rnd2 = round_toInteger(ord2, yrnd2))]
  tp <- subset(tp, complete.cases(tp[,c("rnd1", "rnd2", "chr1", "chr2")]))

  ##############################################################################
  # 2. Make the plot with all hits, regardless of og

  # -- 2.1 subset the hits to those with high enough score
  if(type %in% c("all", "raw")){
    hc <- subset(tp, bitScore > minScore)
    ng1 <- as.integer(uniqueN(hc$ofID1))
    ng2 <- as.integer(uniqueN(hc$ofID2))

    # -- 2.2 get axis labels
    xlab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome1[1], ng1)
    ylab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome2[1], ng2)

    # -- 2.3 subset to chrs with enough genes on them
    hc[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1], na.rm = T), by = "chr1"]
    hc[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2], na.rm = T), by = "chr2"]
    hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)

    # -- 2.4 count n hits in each aggregated position
    hc <- hc[,c("chr1", "chr2", "rnd1", "rnd2")]
    hc <- subset(hc, complete.cases(hc))
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    setorder(hc, -n)
    hc <- subset(hc, !is.na(n))

    # -- 2.5 threshold n to not highlight super strong regions
    qthresh <- quantile(hc$n, quantileThresh)
    if(qthresh > 20)
      qthresh <- 20
    if(qthresh < 5)
      qthresh <- 5
    hc$n[hc$n > qthresh] <- qthresh

    # -- 2.6 get plot title
    titlab <- sprintf(
      "All blast hits with score > %s, %s/%s-gene x/y windows (heatmap range: 2-%s+ hits/window)",
      minScore, xrnd2, yrnd2, round(qthresh))

    # -- 2.7 make the plot
    setorder(hc, n)
    hc <- subset(hc, n > 1)
    nfacets <- nrow(with(hc, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(tp$chr1[order(tp$rnd1)])
      chrOrd2 <- unique(tp$chr2[order(tp$rnd2)])
      hc[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
               chr2 = factor(chr2, levels = chrOrd2))]
      p0 <- ggplot(hc, aes(rnd1, rnd2, col = n)) +
        geom_point(pch = ".") +
        scale_color_viridis_c(begin = .1, trans = "log10", guide = "none") +
        scale_x_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
        scale_y_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
        theme_genespace()+
        facet_grid(chr2 ~ chr1, scales = "free",
                   space = "free", as.table = F, switch = "both")+
        labs(x = xlab, y = ylab, title = titlab)
    }else{
      p0 <- NULL
    }

    ##############################################################################
    # 3. Make the plot with just OG hits

    # -- 2.1 subset the hits to those with high enough score
    hc <- subset(tp, sameOG)
    ng1 <- as.integer(uniqueN(hc$ofID1))
    ng2 <- as.integer(uniqueN(hc$ofID2))

    # -- 2.2 get axis labels
    xlab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome1[1], ng1)
    ylab <- sprintf(
      "%s: gene rank order position (%s genes w/ blast hits), grids every 1000 genes",
      hits$genome2[1], ng2)

    # -- 2.3 subset to chrs with enough genes on them
    hc[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1]), by = "chr1"]
    hc[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2]), by = "chr2"]
    hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)

    # -- 2.4 count n hits in each aggregated position
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    setorder(hc, -n)
    hc <- subset(hc, !is.na(n))

    # -- 2.5 threshold n to not highlight super strong regions
    qthresh <- quantile(hc$n, quantileThresh)
    if(qthresh > 20)
      qthresh <- 20
    if(qthresh < 5)
      qthresh <- 5
    hc$n[hc$n > qthresh] <- qthresh

    # -- 2.6 get plot title
    titlab <- sprintf(
      "Blast hits where query and target are in the same orthogroup, %s/%s-gene x/y windows (heatmap range: 1-%s+ hits/window)",
      xrnd2, yrnd2, round(qthresh))

    # -- 2.7 make the plot
    setorder(hc, n)
    nfacets <- nrow(with(hc, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(tp$chr1[order(tp$rnd1)])
      chrOrd2 <- unique(tp$chr2[order(tp$rnd2)])
      hc[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
               chr2 = factor(chr2, levels = chrOrd2))]
      p1 <- ggplot(hc, aes(rnd1, rnd2, col = n)) +
        geom_point(pch = ".") +
        scale_color_viridis_c(begin = .1, trans = "log10", guide = "none") +
        scale_x_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
        scale_y_continuous(expand = c(0,0),
                           breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
        theme_genespace()+
        facet_grid(chr2 ~ chr1, scales = "free",
                   space = "free", as.table = F, switch = "both")+
        labs(x = xlab, y = ylab, title = titlab)
    }else{
      p1 <- NULL
    }
  }else{
    p1 <- p0 <- NULL
  }

  if(type %in% c("all", "syntenic")){
    ##############################################################################
    # 4. Make the plot with just anchors
    hcBlk <- subset(tp, isAnchor)
    hcBlk[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1]), by = "chr1"]
    hcBlk[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2]), by = "chr2"]
    hcBlk <- subset(hcBlk, ngene1 > minGenes2plot & ngene2 > minGenes2plot)
    hcBlk <- hcBlk[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2", "blkID")]
    blkCols <- sample(gs_colors(uniqueN(hcBlk$blkID)))

    ng1 <- as.integer(uniqueN(hcBlk$ofID1))
    ng2 <- as.integer(uniqueN(hcBlk$ofID2))

    nfacets <- nrow(with(hcBlk, expand.grid(unique(chr1), unique(chr2))))
    if(nfacets < maxFacets){
      chrOrd1 <- unique(hcBlk$chr1[order(hcBlk$rnd1)])
      chrOrd2 <- unique(hcBlk$chr2[order(hcBlk$rnd2)])
      hcBlk[,`:=`(chr1 = factor(chr1, levels = chrOrd1),
                  chr2 = factor(chr2, levels = chrOrd2))]
      hcBlk <- subset(hcBlk, !is.na(rnd1) & !is.na(rnd2))
      if(nrow(hcBlk) < 1){
        warning(sprintf("no syntenic hits found for %s vs. %s",
                        hits$genome1[1], hits$genome2[1]))
        p2 <- NULL
      }else{
        p2 <- ggplot(hcBlk, aes(x = rnd1, y = rnd2, col = blkID)) +
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
      }
    }else{
      p2 <- NULL
    }
  }else{
    p2 <- NULL
  }


  if(is.null(outDir)){
    if(verbose)
      cat("writing to the present graphics device")
    if(!is.null(p0))
      print(p0)
    if(!is.null(p1))
      print(p1)
    if(!is.null(p2))
      print(p2)
  }else{
    dpFile <- file.path(outDir,
                        sprintf("%s_vs_%s.%sHits.pdf",
                                tp$genome1[1], tp$genome2[1], type))
    if(verbose)
      cat(sprintf("writing to file: %s", dpFile))
    if(!print2currentDevice)
      pdf(dpFile, height = ht, width = wd)
    if(!is.null(p0))
      print(p0)
    if(!is.null(p1))
      print(p1)
    if(!is.null(p2))
      print(p2)
    if(!print2currentDevice)
      dev.off()
  }
}


#' @title simple dotplots from a hits data.table
#' @description
#' \code{gghits} ggplot2 integrated graphics to produce dotplots
#' @rdname plot_hits
#' @import data.table
#' @import ggplot2
#' @export
gghits <- function(hits,
                   colorByBlks = TRUE,
                   alpha = ifelse(colorByBlks, 1, .25),
                   useOrder = TRUE,
                   minScore = 0,
                   minGenes2plot = 0){
  ofID1 <- ofID2 <- sameOg <- ngene1 <- ngene2 <- ord1 <- ord2 <- blkID <-
    inBuffer <- rnd2 <- rnd1 <- n <- isArrayRep2 <- isArrayRep1 <- chr1 <-
    noAnchor <- bitScore <- quantile <- chr2 <- sameOG <- isAnchor <-
    start1 <- start2 <- x <- y <- NULL

  tp <- data.table(hits)

  if(colorByBlks){
    hc <- subset(tp, !is.na(blkID) & isAnchor)
  }else{
    hc <- subset(tp, bitScore > minScore)
  }

  ng1 <- as.integer(uniqueN(hc$ofID1))
  ng2 <- as.integer(uniqueN(hc$ofID2))

  hc[,ngene1 := uniqueN(ofID1[!noAnchor & isArrayRep1], na.rm = T), by = "chr1"]
  hc[,ngene2 := uniqueN(ofID2[!noAnchor & isArrayRep2], na.rm = T), by = "chr2"]
  hc <- subset(hc, ngene1 > minGenes2plot & ngene2 > minGenes2plot)

  if(useOrder){
    hc[,`:=`(x = ord1, y = ord2)]
    xlab <- "query gene rank order position"
    ylab <- "target gene rank order position"
  }else{
    hc[,`:=`(x = start1/1e6, y = start2/1e6)]
    xlab <- "query physical (Mb) gene position"
    ylab <- "target physical (Mb) gene position"
  }

  if(!colorByBlks){
    p <- ggplot(hc, aes(x = x, y = y)) +
      geom_point(pch = ".", alpha = alpha) +
      scale_x_continuous(expand = c(0,0), breaks = pretty(hc$x, n = 10))+
      scale_y_continuous(expand = c(0,0), breaks = pretty(hc$y, n = 10))+
      facet_grid(genome2 + chr2 ~ genome1 + chr1, scales = "free",
                 space = "free", as.table = F)+
      labs(x = xlab, y = ylab)+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(
              color = rgb(1, 1, 1, .2), size = .2, linetype = 2),
            panel.spacing = unit(.1, "mm"),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(family = "Helvetica", size = 5),
            axis.title = element_text(family = "Helvetica", size = 6),
            plot.title = element_text(family = "Helvetica", size = 7))
  }else{
    blkCols <- sample(gs_colors(uniqueN(hc$blkID)))
    p <- ggplot(hc, aes(x = x, y = y, col = blkID)) +
      geom_point(pch = ".", alpha = alpha) +
      scale_color_manual(values = blkCols, guide = "none") +
      scale_x_continuous(expand = c(0,0), breaks = pretty(hc$x, n = 10))+
      scale_y_continuous(expand = c(0,0), breaks = pretty(hc$y, n = 10))+
      facet_grid(genome2 + chr2 ~ genome1 + chr1, scales = "free", space = "free",
                 as.table = F)+
      labs(x = xlab, y = ylab)+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(
              color = rgb(1, 1, 1, .2), size = .2, linetype = 2),
            panel.spacing = unit(.1, "mm"),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            axis.text = element_text(family = "Helvetica", size = 5),
            axis.title = element_text(family = "Helvetica", size = 6),
            plot.title = element_text(family = "Helvetica", size = 7))
  }
  print(p)
}
