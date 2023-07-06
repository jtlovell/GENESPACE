#' @title Syntenic windowed whole-genome alignments
#' @description
#' \code{clean_windows} Method to use genomic sequences instead of peptides to
#' estimate synteny across genomes. This is completely distinct from the main
#' GENESPACE routines, relying on minimap2 alignments of windowed sequences.
#' Windowing increases speed relative to standard WGA-->SYRI methods, but comes
#' at a cost of precision. Like the main GENESPACE functions, this is to provide
#' a coarse view of synteny and subsequent local alignments should be used to
#' test for basepair-level differences.
#' @name clean_windows
#'
#' @param faFiles character vector coercible to a file path pointing to the
#' genome assembly fasta files to use
#' @param genomeIDs character vector with the names for the genomes stored
#' in fasta files
#' @param ploidy integer vector corresponding to the ploidies of the genomes
#' stored in the fasta files. Currently polyploids are not allowed. Setting to
#' > 1 will result in an error.
#' @param stripFaNames function applied to basename(faFiles) to get a shorter
#' ID for each genome. Ignored if genomeIDs != NULL.
#' @param windowSize integer, specifying the width of windows. Smaller windows
#' are marginally faster and offer more resolution in unique sequences. Larger
#' windows can be slower and are more precise in repetitive genomes, but also
#' may miss more SVs.
#' @param syntenicBlkSize integer, specifying the number of windows needed to
#' form a syntenic block. By default, the smallest block that can be detected is
#' 40kb, so 40 1kb blocks. If not specified, this will scale so that blocks
#' are ~ 40kb.
#' @param syntenicHitRadius integer, specifying the number of hits away from an
#' initial anchor a syntenic anchor hit can exist.
#' @param syntenicBlkBpRadius integer, the maximum basepair distance between two
#' adjacent anchor hits for which those hits can be in the same syntenic block.
#' @param minChrSize integer, the minimum size in basepairs for a sequence to be
#' considered.
#' @param nCores integer, number of parallel processes to run.
#' @param asmPreset character string matching one of the minimap -x options
#' @param minimap2call character string coercible to a file path point to the
#' minimap2 program call. If in the path, can leave black or set as "minimap2"
#' @param keepSecondary logical, specifying whether secondary hits should be
#' retained.
#' @param mm2mode character string either "fast" or "default". If fast,
#' parameters minimap2 -k 25 -w 20 are added.
#' @param quantileThresh numeric [0-1] specifying the stringency of culling for
#' initial syntenic anchor hits. Lower numbers are more inclusive.
#' @param overwrite logical, if results exist, should they be overwritten.
#' @param verbose logical, should updates be printed to the console?
#' @param outDir file.path for internal functions, not meant for direct use.
#' @param wd file.path for internal functions, not meant for direct use.
#' @param onlySameChrs logical, for internal functions, not meant for direct use.
#' @param faFile1 character, for internal functions, not meant for direct use.
#' @param faFile2 character, for internal functions, not meant for direct use.
#' @param genomeID1 character, for internal functions, not meant for direct use.
#' @param genomeID2 character, for internal functions, not meant for direct use.
#' @param mm2Dir file.path for internal functions, not meant for direct use.
#' @param inFile character, for internal functions, not meant for direct use.
#' @param outFile character, for internal functions, not meant for direct use.
#' @param paf data.table, for internal functions, not meant for direct use.
#' @param x character, for internal functions, not meant for direct use.
#' @param stepSize numeric, for internal functions, not meant for direct use.
#'
#' \cr
#' If called, \code{clean_windows} returns its own arguments.
#'
#' @examples
#' \dontrun{
#' st <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
#' hg38url <- file.path(
#'   st,
#'   "000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz")
#' t2turl <- file.path(
#'   st,
#'   "009/914/755/GCA_009914755.3_T2T-CHM13v1.1/GCA_009914755.3_T2T-CHM13v1.1_genomic.fna.gz")
#' wd <- file.path("~/Downloads/test_genespace")
#' hg38file <- file.path(wd, "hg38.fa.gz")
#' t2tfile <- file.path(wd, "t2t.fa.gz")
#' hg38file1 <- file.path(wd, "hg38Chr1.fa.gz")
#' t2tfile1 <- file.path(wd, "t2tChr1.fa.gz")
#'
#'
#' if(!dir.exists(wd))
#'   dir.create(wd)
#' options(timeout=1000)
#' download.file(url = hg38url, destfile = hg38file)
#' download.file(url = t2turl, destfile = t2tfile)
#'
#'
#' tmp <- Biostrings::readDNAStringSet(hg38file)[1]
#' names(tmp) <- "chr1"
#' Biostrings::writeXStringSet(tmp, filepath = hg38file1, compress = T)
#'
#' tmp <- Biostrings::readDNAStringSet(t2tfile)[1]
#' names(tmp) <- "chr1"
#' Biostrings::writeXStringSet(tmp, filepath = t2tfile1, compress = T)
#'
#' test <- clean_windows(
#'   wd = wd,
#'   nCores = 12,
#'   faFiles = c(hg38file1, t2tfile1))
#' }
#'
#'
#' @title clean_windows
#' @description
#' \code{clean_windows} clean_windows
#' @rdname clean_windows
#' @import data.table
#' @export
clean_windows <- function(faFiles,
                          wd,
                          onlySameChrs = FALSE,
                          genomeIDs = NULL,
                          ploidy = 1,
                          stripFaNames = function(x) gsub("\\.fa(?:\\.gz)?$", "", basename(x)),
                          windowSize = 1e3,
                          syntenicBlkSize = ceiling((200e3 / windowSize) / 5),
                          syntenicHitRadius = 5,
                          syntenicBlkBpRadius = 250e3,
                          minChrSize = windowSize * syntenicBlkSize,
                          nCores = 1,
                          asmPreset = "asm5",
                          minimap2call = "minimap2",
                          keepSecondary = FALSE,
                          mm2mode = "fast",
                          quantileThresh = 0.25,
                          overwrite = FALSE,
                          verbose = TRUE){

  genomeID <- genome1 <- genome2 <- NULL
  # -- 0. Check arguments
  cat("##############", strwrap(
    "**NOTE** `clean_windows()` is still in development - feel free to use it
    (and report errors), but know that many features and error catches are not
    currently fully functional."), "##############\n", sep = "\n")

  if(any(ploidy > 1))
    stop("only haploid genomes are currently permitted\n")

  if(windowSize < 1e3)
    warning(sprintf("window size is set to %s ... using a small window can reduce precision\n", windowSize))

  wd <- path.expand(wd)
  mm2Dir <- file.path(wd, "minimap2")
  outDir <- file.path(wd, "output")

  if(!dir.exists(mm2Dir))
    dir.create(mm2Dir)
  if(!dir.exists(outDir))
    dir.create(outDir)

  if(!all(file.exists(faFiles)))
    stop("some faFiles do not exist\n")

  # -- 1. Make genome path and metadata information
  md <- data.table(rawFile = faFiles)
  if(!is.null(genomeIDs)){
    if(length(genomeIDs) == length(faFiles)){
      md[,genomeID := genomeIDs]
    }else{
      warning("genomeIDs are not NULL and do not match the length of faFiles; parsing using stripFaNames function\n")
      md[,genomeID := stripFaNames(faFiles)]
    }
  }else{
    md[,genomeID := stripFaNames(faFiles)]
  }

  md[,`:=`(fastaFile = file.path(mm2Dir, sprintf("%s.fa", genomeID)),
           windowFile = file.path(mm2Dir, sprintf("%s.window.fa", genomeID)))]
  gv <- md$rawFile; names(gv) <- md$genomeID
  # -- 2. Make the mapping paths and metadata
  cmd <- as.data.table(t(combn(md$genomeID, 2)))
  setnames(cmd, c("genome1", "genome2"))

  cmd[,`:=`(
    faFile1 = file.path(mm2Dir, sprintf("%s.fa", genome1)),
    faFile2 = file.path(mm2Dir, sprintf("%s.fa", genome2)))]

  # -- 3. Run for each combination

  for(i in 1:nrow(cmd)){
    if(verbose)
      cat(sprintf("%s vs %s ... ", cmd$genome1[i], cmd$genome2[i]))
    test <- minimap_4synteny(
      faFile1 = gv[cmd$genome1[i]],
      faFile2 =  gv[cmd$genome2[i]],
      genomeID1 =  cmd$genome1[i],
      genomeID2 = cmd$genome2[i],
      verbose = FALSE,
      windowSize = windowSize,
      syntenicBlkSize = syntenicBlkSize,
      syntenicHitRadius = syntenicHitRadius,
      syntenicBlkBpRadius = syntenicBlkBpRadius,
      minChrSize = minChrSize,
      nCores = nCores,
      outDir = outDir,
      mm2Dir = mm2Dir,
      asmPreset = asmPreset,
      minimap2call = minimap2call,
      keepSecondary = keepSecondary,
      mm2mode = mm2mode,
      quantileThresh = quantileThresh,
      overwrite = overwrite)
    if(verbose)
      cat("Done!\n")
  }
}

#' @title minimap_4synteny
#' @description
#' \code{minimap_4synteny} minimap_4synteny
#' @rdname clean_windows
#' @import data.table
#' @import ggplot2
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stats quantile
#' @export
minimap_4synteny <- function(
    faFile1,
    faFile2,
    genomeID1,
    genomeID2,
    windowSize = 1e3,
    minChrSize = 1e6,
    mm2Dir,
    outDir,
    asmPreset = "asm5",
    minimap2call = "minimap2",
    keepSecondary = FALSE,
    nCores = 12,
    mm2mode = "fast",
    overwrite = FALSE,
    verbose = TRUE,
    quantileThresh = 0.25,
    syntenicBlkSize = ceiling((100e3 / windowSize) / 5),
    syntenicHitRadius = 5,
    syntenicBlkBpRadius = 100e3){

  blkID <- isAnchor <- pos1 <- pos2 <- chr1 <- chr2 <- rnd1 <- rnd2 <- NULL
  ##############################################################################
  # ad hoc function to read and combine the reciprocal paf hits
  combine_pafs <- function(pafFile1, pafFile2, quantileThresh){

    noAnchor <- pmap <- mapq <- windowName <- mapScr <- NULL

    getThresholds <- function(paf, quantileThresh = 0.25){
      mapq <- pmap <- NULL
      tmp1 <- subset(paf, mapq == max(mapq) & pmatch > .5)
      tmp2 <- subset(paf, mapq == max(mapq) & pmap > .5)
      tmp3 <- subset(paf, pmatch > .5 & pmap > .5)
      out <- c(pmap = as.numeric(quantile(tmp1$pmap, quantileThresh)),
               pmatch = as.numeric(quantile(tmp2$pmatch, quantileThresh)),
               mapq = as.numeric(quantile(tmp3$mapq, quantileThresh)))
      return(out)
    }

    parse_pafCols <- function(p){
      nmatches <- windowLength <- pmap <- windowEnd <- windowStart <-
        windowName <- qtmp <- istart <- iend <- queryChr <- windowOrd <-
        queryStart <- queryEnd <- refStart <- refEnd <- mapScr <- mapq <-
        windowRank <- refChr <- refPos <- ord1 <- rank1 <- index <- NULL
      p[,pmatch := nmatches/windowLength]
      p[,pmap := ((windowEnd - windowStart) + 1)/windowLength]
      p[,c("queryChr", "qtmp") := tstrsplit(windowName, ":")]
      p[,c("istart", "iend") := tstrsplit(qtmp, "_")]
      p[, `:=`(istart = as.numeric(istart),
               iend = as.numeric(iend))]
      setkey(p, queryChr, istart, iend)
      p[,windowOrd := as.integer(factor(windowName, levels = unique(windowName)))]
      p[,`:=`(queryStart = istart + windowStart,
              queryEnd = istart + windowEnd,
              qtmp = NULL, istart = NULL, iend = NULL)]
      p[,`:=`(queryPos = (queryStart + queryEnd)/2,
              refPos = (refStart + refEnd)/2,
              mapScr = pmatch * pmap)]
      setorder(p, windowName, -mapScr, -mapq)
      p[,windowRank := 1:.N, by = "windowName"]
      setkey(p, refChr, refPos)

      tmp <- with(p, data.table(
        windowName = windowName, rank1 = windowRank, ord1 = windowOrd,
        chr1 = queryChr, start1 = queryStart, end1 = queryEnd, pos1 = queryPos,
        chr2 = refChr, start2 = refStart, end2 = refEnd, pos2 = refPos,
        mapScr = mapScr, mapq = mapq, pmatch = pmatch, pmap = pmap))
      setkey(tmp, ord1, rank1)
      tmp[,index := sprintf("i%s", 1:.N)]
      return(tmp)
    }

    # genome2 as the reference
    p <- parse_pafCols(read_paf(pafFile2))
    thr <- getThresholds(paf = p, quantileThresh)
    p[,noAnchor := pmap < thr["pmap"] | pmatch < thr["pmatch"] | mapq < thr["mapq"]]
    setorder(p, windowName, -mapScr)
    p1 <- with(subset(p, !duplicated(windowName)), data.table(
      chr1 = chr1, start1 = start1, end1 = end1, pos1 = pos1,
      chr2 = chr2, start2 = start2, end2 = end2, pos2 = pos2,
      index = sprintf("%s_%s", 1, index), noAnchor = noAnchor, mapScr = mapScr, mapq = mapq))

    # genome1 as the reference
    p <- parse_pafCols(read_paf(pafFile1))
    thr <- getThresholds(paf = p, quantileThresh)
    p[,noAnchor := pmap < thr["pmap"] | pmatch < thr["pmatch"] | mapq < thr["mapq"]]
    setorder(p, windowName, -mapScr)
    p2 <- with(subset(p, !duplicated(windowName)), data.table(
      chr1 = chr2, start1 = start2, end1 = end2, pos1 = pos2,
      chr2 = chr1, start2 = start1, end2 = end1, pos2 = pos1,
      index = sprintf("%s_%s", 2, index), noAnchor = noAnchor, mapScr = mapScr, mapq = mapq))
    setnames(p2, colnames(p2))

    p <- rbind(p1, p2)
  }

  ##############################################################################
  # ad hoc function to generate the mm2 command
  run_mm2 <- function(minimap2call, mm2comms, faFile1, faFile2, outFile, overwrite, verbose){
    of <- file.path(basename(dirname(outFile)), basename(outFile))
    of1 <- file.path(basename(dirname(faFile1)), basename(faFile1))
    of2 <- file.path(basename(dirname(faFile2)), basename(faFile2))
    if(file.exists(outFile) && !overwrite){
      if(verbose)
        cat(sprintf(
          "\t%s: file exists and !overwrite, so not re-running\n", of))
    }else{
      coms <- sprintf("%s %s %s -o %s", mm2comms, faFile1, faFile2, outFile)
      comst <- sprintf("%s %s %s %s -o %s", minimap2call, mm2comms, of1, of2, of)
      if(verbose)
        cat(sprintf("\tRunning: %s ... ", comst))
      mm2o <- system2(minimap2call, coms, stdout = TRUE, stderr = TRUE)
      if(verbose)
        cat("Done!\n")
    }
  }

  ##############################################################################
  # 1. parameter checking
  # -- 1.1 minimap2 call
  if(verbose)
    cat("Checking environment ...\n\tminimap2 install ... ")
  stepSize <- windowSize
  asmPreset <- match.arg(asmPreset, choices = c("asm5", "asm10", "asm20"))
  mm2comms <- sprintf("-x %s", asmPreset)
  mm2mode <- match.arg(mm2mode, choices = c("default", "fast"))
  if(mm2mode == "fast")
    mm2comms <- sprintf("%s -k 25 -w 20", mm2comms)
  if(!keepSecondary)
    mm2comms <- sprintf("%s --secondary=no", mm2comms)
  mm2comms <- sprintf("%s -t %s", mm2comms, nCores)
  mmpass <- check_minimap2install(minimap2call)
  if(!mmpass)
    stop("cant find a valid minimap2 install at ", minimap2call, "\n")
  if(verbose)
    cat(sprintf("PASS\n\tminimap2 params: %s\n", mm2comms))

  # -- 1.2 other parameters
  if(!file.exists(faFile1))
    stop("fafile1: ", faFile1, " does not exist\n")
  if(!file.exists(faFile2))
    stop("fafile2: ", faFile2, " does not exist\n")
  genomeID1 <- as.character(genomeID1)
  if(is.na(genomeID1) || is.null(genomeID1))
    stop("genomeID1 must be a unique character string\n")
  if(is.na(genomeID2) || is.null(genomeID2))
    stop("genomeID2 must be a unique character string\n")
  if(genomeID1 == genomeID2)
    stop("genomeID1 must not be the same as genomeID2\n")

  if(!dir.exists(mm2Dir))
    dir.create(mm2Dir)

  if(!dir.exists(outDir))
    dir.create(outDir)

  # -- 1.3 make file.paths
  outRef1 <- file.path(
    mm2Dir, sprintf("%s_minChrSize%s.fa", genomeID1, minChrSize))
  outWin1 <- file.path(
    mm2Dir, sprintf("%s_%swind%sstep.fa", genomeID1, windowSize, stepSize))
  outPaf1 <- file.path(
    mm2Dir, sprintf("%sRef_vs_%sWindow.paf", genomeID1, genomeID2))

  outRef2 <- file.path(
    mm2Dir, sprintf("%s_minChrSize%s.fa", genomeID2, minChrSize))
  outWin2 <- file.path(
    mm2Dir, sprintf("%s_%swind%sstep.fa", genomeID2, windowSize, stepSize))
  outPaf2 <- file.path(
    mm2Dir, sprintf("%sRef_vs_%sWindow.paf", genomeID2, genomeID1))

  outSyn <- file.path(
    outDir, sprintf("%s_vs_%s.synhits.txt", genomeID1, genomeID2))
  outBlk <- file.path(
    outDir, sprintf("%s_vs_%s.blockCoords.txt", genomeID1, genomeID2))
  outDotplot <- file.path(
    outDir, sprintf("%s_vs_%s.dotplot.pdf", genomeID1, genomeID2))

  ##############################################################################
  # 2. [optionally] prepare / parse the genomes
  if(verbose)
    cat(sprintf("Parsing fasta files to sequences with > %s bases ...\n",
                minChrSize))
  tmp <- subset_minChrSize(
    inFile = faFile1, outFile = outRef1, minChrSize = minChrSize,
    overwrite = overwrite, verbose = verbose)
  tmp <- subset_minChrSize(
    inFile = faFile2, outFile = outRef2, minChrSize = minChrSize,
    overwrite = overwrite, verbose = verbose)

  ##############################################################################
  # 3. [optionally] window the genomes
  if(verbose)
    cat(sprintf("Windowing fasta files, step = %s size = %s ...\n",
                windowSize, stepSize))
  tmp <- window_fasta(
    inFile = outRef1, outFile = outWin1,
    windowSize = windowSize, stepSize = stepSize,
    overwrite = overwrite,  verbose = verbose)
  tmp <- window_fasta(
    inFile = outRef2, outFile = outWin2,
    windowSize = windowSize, stepSize = stepSize,
    overwrite = overwrite,  verbose = verbose)

  ##############################################################################
  # 3. [optionally] build the pafs
  if(verbose)
    cat("Running minimap2 and build pafs ...\n")
  tmp <- run_mm2(
    minimap2call, mm2comms = mm2comms,
    faFile1 = outRef1, faFile2 = outWin2, outFile = outPaf1,
    overwrite = overwrite, verbose = verbose)
  tmp <- run_mm2(
    minimap2call, mm2comms = mm2comms,
    faFile1 = outRef2, faFile2 = outWin1, outFile = outPaf2,
    overwrite = overwrite, verbose = verbose)

  ##############################################################################
  # 4. parse the pafs to synteny
  if(verbose)
    cat("Building synteny maps ... ")
  pafComb <- combine_pafs(
    pafFile1 = outPaf1,
    pafFile2 = outPaf2,
    quantileThresh = quantileThresh)

  ps <- flag_pafSynteny(
    paf = pafComb,
    syntenicBlkSize = syntenicBlkSize,
    syntenicHitRadius = syntenicHitRadius,
    windowSize = windowSize,
    syntenicBlkBpRadius = syntenicBlkBpRadius)
  ps[,`:=`(genome1 = genomeID1, genome2 = genomeID2)]

  fwrite(ps, file = outSyn, sep = "\t")
  blks <- calc_pafBlkCoords(ps)
  fwrite(blks, file = outBlk, sep = "\t")

  ##############################################################################
  # 5. Make dotplot

  # 5.1 subset to more-or-less unique hits
  tp <- subset(ps, !is.na(blkID) & isAnchor)
  tp[,`:=`(rnd1 = round_toInteger(pos1, windowSize*2),
           rnd2 = round_toInteger(pos2, windowSize*2))]
  tp <- subset(tp, !duplicated(paste(chr1, chr2, rnd1, rnd2)))

  # 5.2 add in the chromomome ends into the hits
  tmp <- readDNAStringSet(outRef1)
  clens1 <- width(tmp); names(clens1) <- names(tmp)
  tmp <- readDNAStringSet(outRef2)
  clens2 <- width(tmp); names(clens2) <- names(tmp)

  tp0 <- subset(subset(tp, chr1 == chr2), !duplicated(chr1))
  tp0[,`:=`(pos1 = 0, pos2 = 0)]
  tp1 <- data.table(tp0)
  tp1[,`:=`(pos1 = clens1[chr1], pos2 = clens2[chr2])]
  tp <- rbind(tp, tp0, tp1)

  blkCols <- sample(gs_colors(uniqueN(ps$blkID)))
  p1 <- ggplot(tp, aes(x = pos1/1e6, y = pos2/1e6, col = as.factor(blkID)))+
    geom_point(pch = ".")+
    scale_color_manual(values = blkCols, guide = "none")+
    facet_grid(chr1 ~ chr2, scales = "free", space = "free", as.table = F, switch = "both")+
    theme_genespace() +
    scale_x_continuous(expand = c(0,0), breaks = seq(from = 10, to = max(tp$pos1/1e6), by = 10))+
    scale_y_continuous(expand = c(0,0), breaks = seq(from = 10, to = max(tp$pos2)/1e6, by = 10))+
    labs(x = sprintf("%s physical position (grids every 10Mb)", genomeID1),
         y = sprintf("%s physical position (grids every 10Mb)", genomeID2))

  pdf(outDotplot, height = 16, width = 16)
  print(p1)
  dev.off()

  if(verbose)
    cat("Done!\n")
  return(c(synHits = outSyn,
           blkCoords = outBlk,
           dotplots = outDotplot))
}


#' @title check_minimap2install
#' @description
#' \code{check_minimap2install} check_minimap2install
#' @rdname clean_windows
#' @import data.table
#' @export
check_minimap2install <- function(minimap2call){
  path <- path.expand(minimap2call)
  wh <- Sys.which(as.character(path))
  isThere <- basename(wh) == "minimap2"
  if(isThere){
    sys <- system2(minimap2call, "-h", stdout = TRUE)[1]
    isThere <- grepl("Usage: minimap2", sys)
  }
  return(isThere)
}

#' @title subset_minChrSize
#' @description
#' \code{subset_minChrSize} subset_minChrSize
#' @rdname clean_windows
#' @import data.table
#' @importFrom Biostrings readDNAStringSet writeXStringSet width
#' @export
subset_minChrSize <- function(inFile,
                              outFile,
                              minChrSize,
                              overwrite,
                              verbose){
  of <- file.path(basename(dirname(outFile)), basename(outFile))
  if(file.exists(outFile) && !overwrite){
    if(verbose)
      cat(sprintf(
        "\t%s: file exists and !overwrite, so not re-running\n",
        of))
  }else{
    ss <- readDNAStringSet(inFile)
    ss <- ss[width(ss) >= minChrSize]
    writeXStringSet(ss, filepath = outFile)
    if(verbose)
      cat(sprintf(
        "\t%s: %s sequences, %s Mbp\n",
        of, length(ss), round(sum(width(ss)/1e6, 1))))
  }
}

#' @title flag_pafSynteny
#' @description
#' \code{flag_pafSynteny} flag_pafSynteny
#' @rdname clean_windows
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
flag_pafSynteny <- function(paf,
                            syntenicBlkSize,
                            syntenicHitRadius,
                            syntenicBlkBpRadius,
                            windowSize){

  ord1 <- pos1 <- ord2 <- pos2 <- noAnchor <- chr2 <- mapScr <- mapq <- blkID <-
    index <- chr1 <- inBuffer <- isAnchor <- rnd1 <- rnd2 <- m1 <- m2 <- u <- NULL

  # -- subset to potential anchors
  paf <- data.table(paf)
  rnd2n <- (syntenicHitRadius / 2) * windowSize
  paf[,rnd1 := round_toInteger(pos1, rnd2n), by = "chr1"]
  paf[,rnd2 := round_toInteger(pos2, rnd2n), by = "chr2"]
  paf[,ord1 := frank(rnd1, ties.method = "dense"), by = "chr1"]
  paf[,ord2 := frank(rnd2, ties.method = "dense"), by = "chr2"]

  paf <- data.table(paf)

  anch <- subset(paf, !noAnchor)

  setorder(anch, chr2, ord2, -mapScr, -mapq)
  anch <- subset(anch, !duplicated(paste(chr2, ord2)))
  setorder(anch, chr1, ord1, -mapScr, -mapq)
  anch <- subset(anch, !duplicated(paste(chr1, ord1)))

  # -- cull to collinear hits using iterative dbscan
  seqs <- round(seq(from = 2, to = syntenicBlkSize, length.out = 5))

  for(i in seqs){
    anch <- subset(anch, !is.na(ord1) & !is.na(ord2))
    anch[,ord1 := frank(ord1, ties.method = "dense"), by = "chr1"]
    anch[,ord2 := frank(ord2, ties.method = "dense"), by = "chr2"]
    anch[,m1 := max(ord1), by = c("chr1", "chr2")]
    anch[,m2 := max(ord2), by = c("chr1", "chr2")]
    anch <- subset(anch, m1 >= syntenicBlkSize & m2 >= syntenicBlkSize)
    if(max(c(anch$ord1), anch$ord2) >= syntenicBlkSize){
      anch[,blkID := dbscan(frNN(
        x = cbind(ord1, ord2),
        eps = syntenicBlkSize),
        minPts = i)$cluster,
        by = c("chr1", "chr2")]
      anch <- subset(anch, blkID != 0)
    }else{
      print(anch)
      stop()
    }
  }

  ancv <- anch$blkID
  ancu <- with(anch, paste(chr1, chr2, rnd1, rnd2))
  names(ancv) <- names(ancu) <- anch$index
  paf[,u := paste(chr1, chr2, rnd1, rnd2)]
  paf[,`:=`(isAnchor = u %in% ancu,
            blkID = ancu[u])]
  chrs <- with(anch, unique(paste(chr1, chr2)))
  buff <- subset(paf, paste(chr1, chr2) %in% chrs)

  # -- re-rank and pull any hits within a buffer of the anchors
  buffu <- subset(buff, !duplicated(u))
  buffu[,inBuffer := flag_hitsInRadius(
    x = pos1, y = pos2,
    isAnchor = isAnchor,
    radius = syntenicBlkBpRadius),
    by = c("chr1", "chr2")]
  ancu <- unique(subset(buffu, inBuffer)$u)
  paf[,inBuffer := u %in% ancu]

  anch <- subset(paf, inBuffer)
  anch <- subset(anch, !duplicated(u))
  anch[,blkID := dbscan(frNN(
    x = cbind(pos1, pos2),
    eps = syntenicBlkBpRadius),
    minPts = 0)$cluster,
    by = c("chr1", "chr2")]
  anch[,blkID := sprintf("%s_%s_%s",chr1, chr2, blkID)]
  ancv <- anch$blkID; names(ancv) <- anch$u
  paf[,`:=`(isAnchor = u %in% names(ancv),
            blkID = ancv[u])]

  return(paf)
}

#' @title calc_pafBlkCoords
#' @description
#' \code{calc_pafBlkCoords} calc_pafBlkCoords
#' @rdname clean_windows
#' @import data.table
#' @export
calc_pafBlkCoords <- function(paf){

  ord1 <- start1 <- end1 <- start2 <- end2 <-
    index <- ord2 <- max2 <- min2 <- NULL

  setDTthreads(1)

  # -- get the columns and complete observations for these
  hcols <- c("blkID", "start1", "start2", "end1", "end2", "ord1", "ord2",
             "chr1", "chr2", "genome1", "genome2")
  bhits <- subset(paf, complete.cases(paf[,hcols, with = F]))

  # -- get the genome1 coordinates
  setkey(bhits, ord1)
  blks <- bhits[,list(
    start1 = min(start1), end1 = max(end1),
    min2 = min(start2), max2 = max(end2),
    nHits = uniqueN(index),
    orient = ifelse(length(ord1) <= 1, "+",
                    ifelse(cor(jitter(ord1),
                               jitter(ord2)) > 0,"+", "-"))),
    by = c("blkID", "genome1","genome2", "chr1", "chr2")]

  # -- fix the coordinates for inverted blocks
  orient <- NULL
  bgfor <- subset(blks, orient == "+")
  bgrev <- subset(blks, orient == "-")

  maxBp2 <- minBp2 <- maxOrd2 <- minOrd2 <- maxGene2 <- minGene2 <- NULL
  bgrev[,`:=`(start2 = max2, end2 = min2)]
  bgfor[,`:=`(start2 = min2, end2 = max2)]
  blks <- rbind(bgfor, bgrev)
  return(blks)
}

#' @title read_paf
#' @description
#' \code{read_paf} read_paf
#' @rdname clean_windows
#' @import data.table
#' @export
read_paf <- function(x){

  windowName <- mapq <- nmatches <- blkLen <- NULL

  pafNames <- c(
    "windowName", "windowLength", "windowStart", "windowEnd", "windowStrand",
    "refChr", "refLength", "refStart", "refEnd", "nmatches", "blkLen", "mapq")

  # -- read in and give a header to the paf file
  p <- fread(
    x, select = 1:12, col.names = pafNames,
    header = F, showProgress = F, sep = "\t", fill = T)
  setorder(p, windowName, -mapq, -nmatches, -blkLen)
  return(p)
}

#' @title window_fasta
#' @description
#' \code{window_fasta} window_fasta
#' @rdname clean_windows
#' @import data.table
#' @importFrom Biostrings readDNAStringSet width extractAt writeXStringSet
#' @export
window_fasta <- function(inFile,
                         outFile,
                         windowSize,
                         stepSize,
                         verbose,
                         overwrite){

  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to find kmers, install IRanges from bioconductor\n")

  start <- chrlen <- end <- NULL

  of <- file.path(basename(dirname(outFile)), basename(outFile))
  if(file.exists(outFile) && !overwrite){
    if(verbose)
      cat(sprintf(
        "\t%s: file exists and !overwrite, so not re-running\n",
        of))
  }else{
    if(!file.exists(inFile))
      stop("can't find fasta file:", inFile)
    ss <- readDNAStringSet(inFile)

    # -- build the windows start/end coordinate data.table
    windows <- data.table(
      chr = names(ss),
      start = 1,
      chrlen = width(ss))

    windows <- windows[,list(
      start = seq(from = min(start),
                  to = max(chrlen),
                  by = stepSize)),
      by = c("chr", "chrlen")]

    windows[,end := (start + windowSize) - 1]
    windows <- subset(windows, end <= chrlen)
    windows[,chrlen := NULL]

    # -- loop through chromosomes, pulling sequences from the windows
    spl <- split(windows, by = "chr")
    windowsir <- lapply(names(spl), function(i){
      tmp <- extractAt(
        ss[[i]],
        IRanges::IRanges(start = spl[[i]]$start, end = spl[[i]]$end))
      names(tmp) <- sprintf("%s:%s_%s", i, spl[[i]]$start, spl[[i]]$end)
      return(tmp)
    })

    # -- condense into a single dnastringset
    names(windowsir) <- NULL
    windowsir <- do.call(c, windowsir)
    writeXStringSet(windowsir, filepath = outFile)

    if(verbose)
      cat(sprintf(
        "\t%s: %s sequences, %s Mbp\n",
        of, length(windowsir), round(sum(width(windowsir)/1e6, 1))))
  }
}
