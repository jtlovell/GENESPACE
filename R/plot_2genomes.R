#' @title Comparative genome structure plots
#' @description
#' \code{plot_2genomes} Single function to find synteny between two genomes
#' and plot the density of genes and repeats along with syntenic links
#'
#' @param genomeIDs character vector specifying the genome IDs
#' @param faFiles character vector coercible to file paths specifying the
#' locations of the fasta assembly files
#' @param wd character string coercible to file path where the results should
#' be stored
#' @param geneGffFiles character vector coercible to file paths specifying the
#' locations of the gene gff3 annotation files
#' @param repeatGffFiles character vector coercible to file paths specifying the
#' locations of the repeatmasker (or similar) gff3 annotation files
#' @param minChrSize integer specifying the minimum chromosome size to use
#' @param verbose logical specifying whether updates should be printed to the
#' console
#' @param kmers character specifying kmers that should be classified, defualt is
#' to not plot kmer density
#' @param nCores integer specifying the number of parallel processes to run
#' @param kmerMisMatch integer specifying the number of mismatches allowed when
#' search for kmer matches
#' @param slidingwindowSize integer specifying the sliding window size
#' @param slidingwindowStep integer specifying the step between windows
#' @param plotGapSize numeric (0-1) specifying the size of gaps between the
#' largest genomes chromosomes as a fraction of the total genome size
#' @param repeatClassColumnName character specifying which column should be used
#' to grep for repeat classes
#' @param repeatGrep1 character specifying the string to grep for the first
#' repeat density
#' @param repeatGrep2 character specifying the string to grep for the second
#' repeat density
#' @param cdsGrep character specifying the string to grep for CDS in the gene
#' gff3 annotation file
#' @param transcriptGrep character specifying the string to grep for transcript
#' in the gene gff3 annotation file
#' @param plotCols character vector specifying the plot colors
#' @param plotTheme ggplot2 theme to add to the plot
#' @param overwrite logical specifying whether results should be overwritted
#' @param pdfFile character string coercible to a file path where the plot
#' should be written
#' @param returnSourceData logical, should source data be returned
#' @param forceCleanWindows logical, should clean windows output be overwritten?
#' @param ... additional arguments passed to clean_windows
#'
#' @details Coming soon
#'
#' @return A plot, written to file.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom Biostrings readDNAStringSet
#' @export
plot_2genomes <- function(genomeIDs,
                          faFiles,
                          wd,
                          geneGffFiles,
                          repeatGffFiles,
                          minChrSize = 1e6,
                          verbose = TRUE,
                          kmers = NULL,
                          nCores = 1,
                          kmerMisMatch = 0,
                          slidingwindowSize = 1e6,
                          slidingwindowStep = 1e5,
                          plotGapSize = .1,
                          repeatClassColumnName = "class",
                          repeatGrep1 = "Gypsy",
                          repeatGrep2 = "Copia",
                          cdsGrep = "CDS",
                          transcriptGrep = "mRNA",
                          plotCols = NULL,
                          plotTheme = NULL,
                          overwrite = FALSE,
                          pdfFile = NULL,
                          returnSourceData = FALSE,
                          forceCleanWindows = FALSE,
                          ...){

  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to slide genome, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to slide genome, install BiocGenerics from bioconductor\n")
  if(!requireNamespace("rtracklayer", quietly = TRUE))
    stop("to slide genome, install rtracklayer from bioconductor\n")
  if(!requireNamespace("gridExtra", quietly = TRUE))
    stop("to slide genome, install gridExtra from CRAN\n")
  requireNamespace("gridExtra", quietly = TRUE)
  type <- genome <- isLarger <- start <- len <- leftGap <- end <- y <-
    genome1 <- genome2 <- start1 <- chr1 <- start2 <- chr2 <- end1 <- end2 <-
    x <- chr <- propWind <- id <- index <- NULL

  if(is.null(plotCols))
    plotCols <- c(
      "#CC6828", "#F4A460", "#FFFFFF", "#0F4F8B",
      "#4C86C6", "#AED8E6", "#CC2027")
  if(is.null(plotTheme))
    plotTheme <- theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      panel.spacing = unit(0, "cm"),
      axis.title.y = element_blank(),
      plot.title = element_blank())
  if(!is.null(pdfFile)){
    if(!dir.exists(dirname(pdfFile)))
      stop("path to pdf file is not valid")
  }
  ##############################################################################
  # -- 1. get synteny map
  blkFile <- file.path(wd, "output", sprintf(
    "%s_vs_%s.blockCoords.txt", genomeIDs[1], genomeIDs[2]))
  if(verbose)
    cat("Building synteny map ...\n")
  if(file.exists(blkFile) && !overwrite && !forceCleanWindows){
    if(verbose)
      cat("\tblocks file exists and !overwrite, so not re-running clean_windows\n")
  }else{
    winds <- clean_windows(
      faFiles = faFiles,
      wd = wd,
      nCores = nCores,
      genomeIDs = genomeIDs,
      verbose = FALSE,
      overwrite = forceCleanWindows,
      ...)
    if(verbose)
      cat("\tDone!\n")
  }

  if(verbose)
    cat("Reading in assembly fasta files ... ")


  dnass1 <- readDNAStringSet(faFiles[1])
  dnass2 <- readDNAStringSet(faFiles[2])
  seqInfo1 <- pull_seqInfo(dnass1)
  seqInfo2 <- pull_seqInfo(dnass2)

  ##############################################################################
  # -- 2. if kmers are specified, find the positions of these
  if(!is.null(kmers)){
    if(verbose)
      cat("Done!\nFinding kmers ... ")
    if(kmerMisMatch > 0 || length(kmers) < 100){
      kmers1 <- find_fewKmers(
        dnass = dnass1,
        kmers = kmers,
        nCores = nCores,
        max.mismatch = kmerMisMatch)
      if(!is.null(kmers1)){
        kmers1 <- with(as.data.frame(GenomicRanges::reduce(kmers1)), data.table(
          chr = seqnames, start = start, end = end))
      }else{
        kmers1 <- data.table(chr = names(dnass1), start = 1, end = 2)
      }
      kmers2 <- find_fewKmers(
        dnass = dnass2,
        kmers = kmers,
        nCores = nCores,
        max.mismatch = kmerMisMatch)
      if(!is.null(kmers2)){
        kmers2 <- with(as.data.frame(GenomicRanges::reduce(kmers2)), data.table(
          chr = seqnames, start = start, end = end))
      }else{
        kmers2 <- data.table(chr = names(dnass2), start = 1, end = 2)
      }
    }else{
      kmers1 <- find_manyKmers(
        dnass = dnass1,
        kmers = kmers,
        nCores = nCores)

      if(!is.null(kmers1)){
        kmers1 <- with(as.data.frame(GenomicRanges::reduce(kmers1)), data.table(
          chr = seqnames, start = start, end = end))
      }else{
        kmers1 <- data.table(chr = names(dnass1), start = 1, end = 2)
      }
      kmers2 <- find_manyKmers(
        dnass = dnass2,
        kmers = kmers,
        nCores = nCores)

      if(!is.null(kmers2)){
        kmers2 <- with(as.data.frame(GenomicRanges::reduce(kmers2)), data.table(
          chr = seqnames, start = start, end = end))
      }else{
        kmers2 <- data.table(chr = names(dnass2), start = 1, end = 2)
      }
    }
  }else{
    kmers1 <- data.table(chr = names(dnass1), start = 1, end = 2)
    kmers2 <- data.table(chr = names(dnass2), start = 1, end = 2)
  }


  if(verbose)
    cat("\tDone!\nClassifying the genomes ... ")

  ##############################################################################
  # -- 3. Classify the genome
  genes1 <- rtracklayer::readGFF(geneGffFiles[1])
  genes2 <- rtracklayer::readGFF(geneGffFiles[2])
  repeats1 <- rtracklayer::readGFF(repeatGffFiles[1])
  repeats2 <- rtracklayer::readGFF(repeatGffFiles[2])

  nogrp <- sprintf("%s|%s",repeatGrep1, repeatGrep2)

  rep11 <- as.data.frame(subset(repeats1, grepl(repeatGrep1, repeats1[[repeatClassColumnName]])))
  rep12 <- as.data.frame(subset(repeats1, grepl(repeatGrep2, repeats1[[repeatClassColumnName]])))
  rep13 <- as.data.frame(subset(repeats1, !grepl(nogrp, repeats1[[repeatClassColumnName]])))
  cds1 <- as.data.frame(subset(genes1, type == cdsGrep))
  gene1 <- as.data.frame(subset(genes1, type == transcriptGrep))

  rep21 <- subset(repeats2, grepl(repeatGrep1, repeats2[[repeatClassColumnName]]))
  rep22 <- subset(repeats2, grepl(repeatGrep2, repeats2[[repeatClassColumnName]]))
  rep23 <- subset(repeats2, !grepl(nogrp, repeats2[[repeatClassColumnName]]))
  cds2 <- as.data.frame(subset(genes2, type == cdsGrep))
  gene2 <- as.data.frame(subset(genes2, type == transcriptGrep))

  ##############################################################################
  # - hierarchical classification


  ##############################################################################
  # 4. Sliding windows
  if(is.null(kmers)){
    beds1 <- list(
      cds = with(cds1, data.table(chr = seqid, start = start, end = end)),
      rep1 = with(rep11, data.table(chr = seqid, start = start, end = end)),
      rep2 = with(rep12, data.table(chr = seqid, start = start, end = end)),
      otherRepeat = with(rep13, data.table(chr = seqid, start = start, end = end)),
      introns = with(gene1, data.table(chr = seqid, start = start, end = end)))

    beds2 <- list(
      cds = with(cds2, data.table(chr = seqid, start = start, end = end)),
      rep1 = with(rep21, data.table(chr = seqid, start = start, end = end)),
      rep2 = with(rep22, data.table(chr = seqid, start = start, end = end)),
      otherRepeat = with(rep23, data.table(chr = seqid, start = start, end = end)),
      introns = with(gene2, data.table(chr = seqid, start = start, end = end)))
  }else{
    beds1 <- list(
      kmers = kmers1,
      cds = with(cds1, data.table(chr = seqid, start = start, end = end)),
      rep1 = with(rep11, data.table(chr = seqid, start = start, end = end)),
      rep2 = with(rep12, data.table(chr = seqid, start = start, end = end)),
      otherRepeat = with(rep13, data.table(chr = seqid, start = start, end = end)),
      introns = with(gene1, data.table(chr = seqid, start = start, end = end)))

    beds2 <- list(
      kmers = kmers2,
      cds = with(cds2, data.table(chr = seqid, start = start, end = end)),
      rep1 = with(rep21, data.table(chr = seqid, start = start, end = end)),
      rep2 = with(rep22, data.table(chr = seqid, start = start, end = end)),
      otherRepeat = with(rep23, data.table(chr = seqid, start = start, end = end)),
      introns = with(gene2, data.table(chr = seqid, start = start, end = end)))
  }

  classes1 <- classify_genome(
    dnaSS = dnass1, listOfBeds = beds1, verbose = T)
  classes2 <- classify_genome(
    dnaSS = dnass2, listOfBeds = beds2, verbose = T)

  if(is.null(kmers)){
    sw1 <- slide_genome(
      seqInfo = seqInfo1,
      listOfGrs = classes1[c(1,5,6,2,3,4)],
      windowSize = slidingwindowSize,
      stepSize = slidingwindowStep)
    sw2 <- slide_genome(
      seqInfo = seqInfo2,
      listOfGrs = classes2[c(1,5,6,2,3,4)],
      windowSize = slidingwindowSize,
      stepSize = slidingwindowStep)
  }else{
    sw1 <- slide_genome(
      seqInfo = seqInfo1,
      listOfGrs = classes1[c(2,6,7,3,4,5,1)],
      windowSize = slidingwindowSize,
      stepSize = slidingwindowStep)
    sw2 <- slide_genome(
      seqInfo = seqInfo2,
      listOfGrs = classes2[c(2,6,7,3,4,5,1)],
      windowSize = slidingwindowSize,
      stepSize = slidingwindowStep)
  }

  sw1[,genome := genomeIDs[1]]
  sw2[,genome := genomeIDs[2]]
  swtp <- rbind(sw1, sw2)

  ##############################################################################
  # 5. get linear coordinates in order
  dnass1 <- dnass1[width(dnass1) >= minChrSize]
  dnass2 <- dnass2[width(dnass2) >= minChrSize]
  chrs1 <- names(dnass1)
  chrs2 <- names(dnass2)
  clens1 <- width(dnass1); names(clens1) <- chrs1
  clens2 <- width(dnass2); names(clens2) <- chrs2

  blks <- fread(blkFile)

  len1 <- sum(clens1)
  len2 <- sum(clens2)

  plotGapSize <- (max(c(len1, len2)) * plotGapSize) / (length(c(chrs1, chrs1)) / 2)
  size1 <- ((length(clens1) - 1) * plotGapSize) + len1
  size2 <- ((length(clens2) - 1) * plotGapSize) + len2
  if(size1 > size2){
    gapSize1 <- plotGapSize
    gapSize2 <- (size1 - len2) / (length(clens2) + 1)
    smaller <- "genome2"
    gpSize <- gapSize2
    mpSize <- gapSize1
  }else{
    gapSize2 <- plotGapSize
    gapSize1 <- (size2 - len1) / (length(clens1) + 1)
    smaller <- "genome1"
    gpSize <- gapSize1
    mpSize <- gapSize2
  }

  pmd <- rbind(
    data.table(genome = genomeIDs[1], chr = chrs1, len = clens1,
               leftGap = gapSize1),
    data.table(genome = genomeIDs[2], chr = chrs2, len = clens2,
               leftGap = gapSize2))
  pmd[,isLarger := ifelse(size1 > size2, genomeIDs[1], genomeIDs[2]) == genome]
  pmd[,start := c(0, (cumsum(len + leftGap))[-.N]), by = "genome"]
  pmd$start[!pmd$isLarger] <- pmd$start[!pmd$isLarger] + gpSize
  pmd[,end := start + len]
  pmd[,`:=`(y = match(genome, genomeIDs))]
  pmd[,`:=`(y1 = y + 0.04, y2 = y - 0.04)]


  lin1 <- pmd$start[pmd$genome == genomeIDs[1]]
  lin2 <- pmd$start[pmd$genome == genomeIDs[2]]
  names(lin1) <- pmd$chr[pmd$genome == genomeIDs[1]]
  names(lin2) <- pmd$chr[pmd$genome == genomeIDs[2]]
  lin <- pmd$start; names(lin) <- with(pmd, paste(genome, chr))

  ################################################################################
  # -- Parse the synteny map
  dat <- blks[,c("chr1", "chr2", "start1", "start2", "end1", "end2", "orient")]
  dat[,genome1 := genomeIDs[1]]
  dat[,genome2 := genomeIDs[2]]
  dat[,`:=`(start1 = start1 + lin[paste(genome1, chr1)],
            start2 = start2 + lin[paste(genome2, chr2)],
            end1 = end1 + lin[paste(genome1, chr1)],
            end2 = end2 + lin[paste(genome2, chr2)])]
  dat <- subset(dat, complete.cases(dat))

  braidPolygons <- rbindlist(lapply(1:nrow(dat), function(i){
    x <- with(dat[i, ], calc_curvePolygon(
      start1 = start1, end1 = end1, start2 = start2,
      end2 = end2, y1 = 1.05, y2 = 1.95))
    x[,`:=`(chr = dat$chr1[i], index = i, type = dat$orient[i])]
    return(x)
  }))

  chrPolygons <- rbindlist(lapply(1:nrow(pmd), function(i){
    z <- pmd[i,]
    out <- data.table(z[,c("genome","chr")], with(z, round_rect(
      xleft = start, xright = end, ybottom = y2,
      ytop = y1, yrange = range(c(pmd$y1, pmd$y2)),
      xrange = range(c(pmd$start, pmd$end)),
      plotWidth = 18,
      plotHeight = 1)))
    return(out)
  }))

  swtp[,x := (start + end - 1)/2 + lin[paste(genome, chr)]]
  swtp <- subset(swtp, complete.cases(swtp))

  # 6. get data parse to make the plot
  dictext <- paste(c(cdsGrep, transcriptGrep, "unannotated", repeatGrep1, repeatGrep2, "otherRepeat"),
                   collapse = ", ")
  if(!is.null(kmers))
    dictext <- sprintf("%s, kmers", dictext)
  xlab1 <- sprintf(
    "%s chromosomes (%s total Mb, %sMb-overlapping %sMb windows)\nColors (top to bottom): %s",
    genomeIDs[1],
    round(sum(clens1)/1e6, 1),
    round((slidingwindowSize - slidingwindowStep)/1e6, 2),
    round(slidingwindowSize/1e6, 2),
    dictext)
  xlab2 <- sprintf(
    "%s chromosomes (%s total Mb, %sMb-overlapping %sMb windows)\nColors (top to bottom): %s",
    genomeIDs[2],
    round(sum(clens2)/1e6, 1),
    round((slidingwindowSize - slidingwindowStep)/1e6, 2),
    round(slidingwindowSize/1e6, 2),
    dictext)
  swtp[,genome := factor(genome, levels = rev(genomeIDs))]
  p1 <- ggplot(swtp, aes(x = x, y = propWind, fill = id))+
    scale_fill_manual(values = plotCols, guide = "none") +
    scale_x_continuous(expand = c(0,0), name = xlab1, position = "top",
                       limits = c(0, max(pmd$end) + mpSize))+
    scale_y_continuous(expand = c(0,0))+
    plotTheme +
    theme(plot.margin = unit(c(.1,0,0,.2), "lines"),
          axis.title.x = element_text(family = "Helvetica", size = 6, vjust = -3))
  for(i in chrs1){
    p1 <- p1 + geom_area(data = subset(swtp, chr == i & genome == genomeIDs[1]))
  }

  p2 <- ggplot(swtp, aes(x = x, y = propWind, fill = id))+
    scale_fill_manual(values = plotCols, guide = "none") +
    scale_x_continuous(expand = c(0,0), name = xlab2, limits = c(0, max(pmd$end) + mpSize))+
    scale_y_continuous(expand = c(0,0))+
    plotTheme+
    theme(plot.margin = unit(c(-.1,0,.1,.2), "lines"),
          axis.title.x = element_text(family = "Helvetica", size = 6, vjust = 3))
  for(i in chrs2){
    p2 <- p2 + geom_area(data = subset(swtp, chr == i & genome == genomeIDs[2]))
  }

  p3 <- ggplot()+
    geom_polygon(
      data = braidPolygons,
      aes(x = x, y = y, group = index, fill = type),
      alpha = 1,  lwd = .01, col = NA)+
    geom_polygon(
      data = chrPolygons,
      aes(x = x, y = y, group = genome),
      alpha = 1,  lwd = .01, fill = "black", col = "white")+
    geom_text(
      data = pmd,
      aes(x = (start + end)/2, y = y, label = gsub("Gm0|Gm","",chr)),
      alpha = 1, col = "white", size = 2)+
    scale_x_continuous(expand = c(0,0), limits = c(0, max(pmd$end) + mpSize))+
    scale_y_reverse(expand = c(0,0))+
    scale_fill_manual(values = c("darkred", "grey"), guide = "none")+
    plotTheme+
    theme(axis.title.x = element_blank(),
          plot.margin = unit(c(0,0,0,.2), "lines"))


  if(is.null(pdfFile))
    pdfFile <- file.path(wd, sprintf("%s_vs_%s_swPlot.pdf", genomeIDs[1], genomeIDs[2]))
  pdf(pdfFile, height = 3.5, width = 9)
  gridExtra::grid.arrange(p1, p3, p2, nrow = 3)
  dev.off()

  if(returnSourceData)
    return(list(braidPolygons = braidPolygons, chrPolygons = chrPolygons,
                plotMetadata = pmd, slidingWindows = swtp))
}






