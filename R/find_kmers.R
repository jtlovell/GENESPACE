#' @title Generic internal functions used by genespace
#' @description
#' \code{find_kmers} Convenience functions for gscTools, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name find_kmers
#'
#' @param dnass DNAStringSet
#' @param kmers DNAStringSet with kmers to find
#' @param nCores integer, then number of cores to use
#' @param minRunLength numeric, the minimum size for a kmer run to be counted
#' @param maxDistBtwTelo numeric, maximum distance between two telomere kmers
#' for the telomere run to be considered a single run
#' @param minTeloSize numeric, the minimum number of bases for a telomere run
#' to be counted
#' @param minTeloDens numeric [0-1], the minimum density of telomere kmers in
#' the run
#' @param maxDist2end numeric, the maximum distance to a chromosome end for a
#' telomere to be called terminal (and not interstitial)
#' @param ... additional arguments passed to vmatchpattern
#' \cr
#' If called, \code{find_kmers} returns its own arguments.
#'
#'
#' @title Fast method to map many kmers
#' @description
#' \code{find_manyKmers} Finds exact DNA-DNA matches and aggregates across all
#' kmers in the kmer DNAstring
#' @rdname find_kmers
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings PDict matchPDict
#' @export
find_manyKmers <- function(dnass,
                           kmers,
                           nCores = 1){

  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to find kmers, install IRanges from bioconductor\n")
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to find kmers, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to find kmers, install BiocGenerics from bioconductor\n")

  pdictKmer <- PDict(kmers)
  grlist <- mclapply(names(dnass), mc.cores = nCores, function(i){
    tmp <- GenomicRanges::GRanges(matchPDict(
      pdict = pdictKmer,
      subject = dnass[[i]]))
    tmp <- GenomicRanges::reduce(
      tmp,
      drop.empty.ranges = T)
    if(length(tmp) > 1){
      out <- GenomicRanges::GRanges(
        seqnames = i,
        IRanges::IRanges(
          start = BiocGenerics::start(tmp),
          end = BiocGenerics::end(tmp)))
      return(out)
    }
  })
  if(length(grlist) == 0){
    stop("could not find any matches to the kmers\n")
  }else{
    suppressWarnings(outgr <- BiocGenerics::do.call(
      c,
      grlist[sapply(grlist, length) > 0]))
    return(GenomicRanges::reduce(outgr))
  }
}

#' @title Fast method to map many kmers
#' @description
#' \code{find_manyKmers} Finds exact DNA-DNA matches and aggregates across all
#' kmers in the kmer DNAstring
#' @rdname find_kmers
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings PDict matchPDict vmatchPattern
#' @export
find_fewKmers <- function(dnass,
                          kmers,
                          nCores = 1,
                          ...){

  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to find kmers, install BiocGenerics from bioconductor\n")

  si <- pull_seqInfo(dnass)

  kmerposList <- mclapply(kmers, mc.cores = nCores, function(i)
    vmatchPattern(pattern = as.character(i), subject = dnass, ...))

  # -- remove any chromosomes without any Ns
  kmerposList <- kmerposList[sapply(kmerposList, length) > 0]

  # -- combine grs with Ns
  kmerpos <- BiocGenerics::do.call(c, lapply(kmerposList, function(kmeri){
    kmeri <- kmeri[sapply(kmeri, length) > 0]
    kmeri <- BiocGenerics::do.call(c, lapply(names(kmeri), function(i)
      GenomicRanges::GRanges(i, kmeri[[i]], seqinfo = si)))
    return(kmeri)
  }))
  return(kmerpos)
}

#' @title Find contig gap positions as a run of N's
#' @description
#' \code{find_runsOfNs} Exact string matching of Ns to a DNAStringSet to find
#' the position of gaps between contigs
#' @rdname find_kmers
#' @import data.table
#' @importFrom Biostrings vmatchPattern
#' @export
find_runsOfNs <- function(dnass, minRunLength){

  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to find kmers, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to find kmers, install BiocGenerics from bioconductor\n")

  # -- get sequence lengths
  si <- pull_seqInfo(dnass)

  # -- find all occurances of N
  nposList <- vmatchPattern(pattern = "N", dnass)

  # -- remove any chromosomes without any Ns
  nposList <- nposList[sapply(nposList, length) > 0]

  # -- combine grs with Ns
  npos <- BiocGenerics::do.call(c, lapply(names(nposList), function(i)
    GenomicRanges::GRanges(i, nposList[[i]], seqinfo = si)))
  if(!is.null(npos)){
    npos <- GenomicRanges::reduce(
      npos,
      min.gapwidth = 2,
      ignore.strand = TRUE,
      drop.empty.ranges = TRUE)
  }

  return(npos)
}

#' @title find_telomeres
#' @description
#' \code{find_telomeres} find_telomeres
#' @rdname find_kmers
#' @import data.table
#' @importFrom Biostrings vmatchPattern seqinfo
#' @export
find_telomeres <- function(dnass,
                           kmers,
                           maxDistBtwTelo = 20,
                           minTeloSize = 200,
                           minTeloDens = 0.75,
                           maxDist2end = 1e4,
                           nCores = 1){

  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to find kmers, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to find kmers, install BiocGenerics from bioconductor\n")
  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to find kmers, install IRanges from bioconductor\n")
  if(!requireNamespace("GenomeInfoDb", quietly = TRUE))
    stop("to find kmers, install GenomeInfoDb from bioconductor\n")


  si <- seqinfo(dnass)
  kmerpos <- find_fewKmers(
    dnass = dnass,
    kmers = kmers,
    nCores = nCores)

  if(is.null(kmerpos))
    return(NULL)

  kmers <- as.character(kmers)
  if(uniqueN(nchar(kmers)) != 1)
    warning("variable kmer lengths ... using longest for density calculations\n")
  kmerout <- GenomicRanges::reduce(
    kmerpos,
    min.gapwidth = maxDistBtwTelo,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE,
    with.revmap = TRUE)

  dens <- sapply(kmerout@elementMetadata$revmap, length) * max(nchar(kmers))
  wid <- BiocGenerics::width(kmerout)
  dens <- dens / wid
  kmersub <- subset(kmerout, dens >= minTeloDens & wid >= minTeloSize)
  isLeft <- BiocGenerics::start(kmersub) <= maxDist2end
  rightBound <- (GenomeInfoDb::seqlengths(si) - maxDist2end)[as.character(GenomeInfoDb::seqnames(kmersub))]
  isRight <- BiocGenerics::end(kmersub) >= rightBound
  teloclass <- ifelse(isLeft, "left", ifelse(isRight, "right", "inter"))
  kmerout <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(kmersub),
    IRanges::IRanges(start = BiocGenerics::start(kmersub),
                     end = BiocGenerics::end(kmersub)),
    position = teloclass,
    seqinfo = si)
  return(kmerout)
}




