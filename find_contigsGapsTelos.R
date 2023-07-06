#' @title get coordinates aggregated across kmers
#'
#' @description
#' \code{find_contigsGapsTelos} find all positions with an exact match to any kmers
#' and join them
#'
#' @param dnass DNAStringSet containing the assembly fasta
#' @param teloKmers character vector with the kmers in telomers
#' @param minContigGapSize numeric, minimum size of a run of N's to be called a
#' gap
#' @param maxDistBtwTelo numeric, minimum distance between neighboring telomere
#'  kmers
#' @param minTeloSize numeric, minimum size for telomere kmers cluster
#' @param minTeloDens numeric (0-1), minimum density for telomere kmers
#' @param minChrSize integer, minimum size of a scaffold to be included
#' @param maxDist2end integer, maximum distance to the chromosome end for a
#' telomere to be called at that arm.
#' @param verbose logical, should updates be printed to the console?
#'
#' @details finds kmers and reduces overlapping intervals
#'
#' @return a granges object containing any sequences masked by the kmers
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings PDict matchPDict N50 width reverseComplement DNAStringSet
#' @export
find_contigsGapsTelos <- function(dnass,
                                  teloKmers,
                                  minContigGapSize = 100,
                                  maxDistBtwTelo = 20,
                                  minTeloSize = 200,
                                  minTeloDens = 0.75,
                                  minChrSize = 0,
                                  maxDist2end = 1e4,
                                  verbose = TRUE){

  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to join ranges, install IRanges from bioconductor\n")
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to join ranges, install GenomicRanges from bioconductor\n")

  ##############################################################################
  # 1. prep the data.
  # -- 1.1 read in the fasta as a DNAStringSet
  pos <- NULL
  ss <- dnass

  # -- 1.2 subset to chromosomes bigger than minChrSize
  nchr <- length(ss)
  ss <- ss[width(ss) > minChrSize]
  if(verbose && minChrSize > 0 && (nchr - length(ss)) > 0)
    cat(sprintf(
      "Dropping %s of %s chromosomes that are smaller than %sMb\n",
      nchr - length(ss), nchr, round(minChrSize/1e6, 2)))

  # # -- 1.3 re-name chromosomes if necessary
  # ss <- rename_chrs(ss)

  # -- 1.4 get chromosome lengths in the order they are in the fasta
  chrv <- width(ss); names(chrv) <- names(ss)

  # -- 1.5 print initial scaffold stats
  if(verbose)
    cat(strwrap(sprintf(
      "\tN. scaffolds = %s (N50 = %sMb); total length = %sMb\n",
      length(chrv),
      round(N50(chrv / 1e6), 2),
      round(sum(chrv / 1e6, na.rm = T), 2)),
      indent = 8, exdent = 16), sep = "\n")

  ##############################################################################
  # 2. Get the gaps and report contig-level stats
  gapGr <- find_runsOfNs(
    dnass = ss,
    minRunLength = minContigGapSize)
  if(is.null(gapGr)){
    cat("\tThis genome is gapless, scaffold = contig stats\n")
    contigGr <- NULL
  }else{
    gapGr <- gapGr[width(gapGr) >= minContigGapSize]
    contigGr <- find_gaps(gapGr)

    cat(strwrap(sprintf(
      "\tN. contigs = %s (N50 = %sMb); total length = %sMb\n",
      length(contigGr),
      round(N50(width(contigGr) / 1e6), 2),
      round(sum(width(contigGr) / 1e6, na.rm = T), 2)),
      indent = 8, exdent = 16), sep = "\n")
  }

  ##############################################################################
  # 3. Get the telomere kmer positions and cluster
  teloKmers <- c(
    reverseComplement(DNAStringSet(teloKmers)),
    teloKmers)
  teloGr <- find_telomeres(
    dnass = ss,
    kmers = teloKmers,
    maxDistBtwTelo = maxDistBtwTelo,
    minTeloSize = minTeloSize,
    minTeloDens = minTeloDens,
    maxDist2end = maxDist2end)
  if(is.null(teloGr)){
    cat("\tCould not find any telomeres\n")
  }else{
    teloMd <- data.table(
      chr = as.character(GenomeInfoDb::seqnames(teloGr)),
      pos = toupper(substr(teloGr@elementMetadata$position, 1, 1)))
    teloMd <- teloMd[,list(pos = paste(pos, collapse = "")), by = "chr"]
    cat(strwrap(sprintf(
      "\tN. telomeres = %s (%skb) on chrs: %s\n",
      length(teloGr),
      round(N50(width(teloGr) / 1e3), 2),
      paste(sprintf("%s%s", teloMd$chr, teloMd$pos), collapse = ", ")),
      indent = 8, exdent = 16), sep = "\n")
  }

  ##############################################################################
  # 4. Return the results
  return(list(gaps = gapGr, contigs = contigGr, telomeres = teloGr))
}
