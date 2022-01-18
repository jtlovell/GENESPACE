#' @title Sliding window and kmer exploration of genomes
#'
#' @description
#' \code{analyze_genome} Sliding window and kmer exploration of genomes
#'
#' @name analyze_genome
#'
#' @param path2fasta file.path to a fast file
#' @param kmerSize integer of length 1, specifying the size of the kmer to use
#' @param kmers character vector of kmers to search for
#' @param nCores integer of length 1, specifying the number of parallel
#' processes to run
#' @param revComp logical, should the kmers be reverse completemented too?
#' @param minSize2batch integer of length 1, specifying how many kmers should
#' be run in a single batch
#' @param verbose logical, should updates be printed to the console?
#' @param bedObj data.table or data.frame with chr, start, end columns
#' @param step integer of length 1, specifying the distance between adjacent
#' windows in bps.
#' @param size integer of length 1, specifying sliding window size in bps.
#' @details ...
#'
#' @return ...
#'
#' @examples
#' \dontrun{
#' teloKmer <- c("CCCGAAA", "CCCTAAA")
#' teloLoc <- find_kmer(
#'   kmers = teloKmer,
#'   path2fasta = fallaxFa)
#' swtelo <- slide_genome(
#'   path2fasta = fallaxFa,
#'   windowSize = 5e3,
#'   stepSize = 500,
#'   bedObj = teloLoc)
#' ggplot(subset(swtelo$sw, grepl("LG", chr)), aes(x = Mbp, y = prop))+
#'   geom_line()+facet_wrap(~chr)
#' }
#'
#' @note \code{analyze_genome} is a generic name for the functions documented.
#' \cr
#' If called, \code{analyze_genome} returns its own arguments.
#'
#' @title kmerize_fa
#' @description
#' \code{kmerize_fa} kmerize_fa
#' @rdname analyze_genome
#' @import data.table
#' @importFrom Biostrings DNAString readDNAStringSet
#' @export
kmerize_fa <- function(path2fasta, kmerSize){
  # -- read fasta into memory
  dnass <- readDNAStringSet(path2fasta)
  seqList <- lapply(1:length(dnass), function(i) as.character(DNAString(dnass[[i]])))
  ssl <- unlist(lapply(seqList, function(x){
    ind <- data.table(start = 1:(nchar(x)-(kmerSize-1)),
                      end = kmerSize:nchar(x))
    seqs <- sapply(1:nrow(ind), function(i)
      substr(x, ind$start[i], ind$end[i]))
    return(seqs[!duplicated(seqs)])
  }))
  ssl <- ssl[!duplicated(ssl)]
  return(ssl)
}

#' @title find_kmer
#' @description
#' \code{find_kmer} find_kmer
#' @rdname analyze_genome
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings DNAString reverseComplement PDict matchPDict readDNAStringSet
#' @export
find_kmer <- function(path2fasta,
                      kmers,
                      nCores = 1,
                      revComp = TRUE,
                      minSize2batch = 50,
                      verbose = TRUE){

  find_fewKmers <- function(dnass, kmers, nCores, revComp, verbose){
    # -- convert kmer strings into more accessible format
    kmerList <- lapply(1:length(kmers), function(i) DNAString(kmers[i]))
    # -- optionally add reverse complement seq
    if(revComp){
      kmerListr <- lapply(1:length(kmers), function(i) reverseComplement(DNAString(kmers[i])))
      kmerList <- c(kmerList, kmerListr)
    }
    # -- only keep unique kmers
    ndup <- which(!duplicated(sapply(kmerList, as.character)))
    kmerList <- kmerList[ndup]
    # -- for each kmer ...
    out <- rbindlist(mclapply(kmerList, mc.cores = nCores, function(x){
      if(verbose)
        cat(sprintf("\tKmer %s: ", as.character(x)))
      # -- get the index of positions of exact hits
      m <- vmatchPattern(pattern = x, subject = dnass)
      # -- write as a matrix
      o <- rbindlist(lapply(names(m), function(i)
        if(nrow(as.data.table(m[[i]])) > 0)
          data.table(chr = i, start = start(m[[i]]), end = end(m[[i]]))))
      o[,kmer := as.character(x)]
      if(verbose)
        cat(sprintf("found %s hits on %s sequences\n",
                    nrow(o), uniqueN(o$chr)))
      return(o)
    }))
    return(out)
  }

  find_manyKmers <- function(dnass, kmers, nCores, revComp, verbose){
    # -- convert kmer strings into more accessible format
    kmerList <- lapply(1:length(kmers), function(i) DNAString(kmers[i]))
    # -- optionally add reverse complement seq
    if(revComp){
      kmerListr <- lapply(1:length(kmers), function(i)
        reverseComplement(DNAString(kmers[i])))
      kmerList <- c(kmerList, kmerListr)
    }
    # -- only keep unique kmers
    ndup <- which(!duplicated(sapply(kmerList, as.character)))
    kmerSS <- DNAStringSet(kmerList[ndup])
    # -- make the preprocessed pattern dictionary
    pdict0 <- PDict(kmerSS)
    pmn <- as.character(kmerSS)
    outa <- rbindlist(mclapply(names(dnass), mc.cores = nCores, function(j){
      if(verbose)
        cat(sprintf("\tSequence %s: ", j))
      ss1 <- dnass[[j]]
      m <- matchPDict(pdict0, ss1)
      whm <- which(sapply(m, length) > 0)
      if(length(whm) > 0){
        out <- rbindlist(lapply(whm,function(i)
          data.table(
            chr = j,
            kmer = pmn[[i]],
            start = start(m[[i]]),
            end = end(m[[i]]))))
        if(verbose)
          cat(sprintf("found %s hits of %s unique kmers\n",
                      nrow(out), uniqueN(out$kmer)))
        return(out)
      }else{
        if(verbose)
          cat("found no exact matches\n")
      }
    }))
    return(outa)
  }

  if(verbose)
    cat("Reading fasta into memory ... \n")
  dnass <- readDNAStringSet(path2fasta)
  if(verbose)
    cat(sprintf("Looking for exact matches of %s kmers on %s sequences ... \n",
                length(kmers), length(dnass)))
  if(length(kmers) > minSize2batch){
    if(verbose)
      cat(sprintf("Since n. kmers > %s, running by chromosome ...\n",
                  minSize2batch))
    out <- find_manyKmers(
      dnass = dnass,
      kmers = kmers,
      nCores = nCores,
      revComp = revComp,
      verbose = verbose)
  }else{
    if(verbose)
      cat(sprintf("Since n. kmers <= %s, running by kmer ...\n",
                  minSize2batch))
    out <- find_fewKmers(
      dnass = dnass,
      kmers = kmers,
      nCores = nCores,
      revComp = revComp,
      verbose = verbose)
  }
  return(out)
}

#' @title slide_genome
#' @description
#' \code{slide_genome} slide_genome
#' @rdname analyze_genome
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom Biostrings width readDNAStringSet
#' @export
slide_genome <- function(path2fasta,
                         bedObj,
                         step,
                         size,
                         nCores = 1){

  # -- calculate chr sizes
  dnass <- readDNAStringSet(path2fasta)
  chrSize <- width(dnass); names(chrSize) <- names(dnass)
  chrSize <- chrSize[chrSize >= size]

  # -- convert bed to 1-index
  reg1 <- with(bedObj, data.table(chr = chr, start = start - 1, end = end - 1))
  setkey(reg1, chr, start, end)

  # -- make genome into windows, by chr
  out <- rbindlist(mclapply(names(chrSize), mc.cores = nCores, function(i){

    # -- make a data.table with the window coordinates
    mxst <- chrSize[i] - size
    winds <- data.table(chr = i, start = seq(from = 1, to = mxst, by = step))
    winds[,end := start + size -1]

    # -- calculate window width
    winds[,width := end - start]

    # -- set the window key to match the chr region key
    setkey(winds, chr, start, end)
    reg <- subset(reg1, chr == i)

    if(nrow(reg) > 0 && nrow(winds) > 0){

      # -- calculate overlapping intervals
      fo <- foverlaps(reg, winds)

      # -- drop non-overlapping intervals
      fo <- subset(fo, complete.cases(fo))

      # -- summarize the number of bp in intervals
      fo <- fo[,list(nbp = sum(i.end - i.start)),
               by = c("chr", "start", "end", "width")]

      # -- calculate the proportion in intervals and return
      fo[,prop := nbp / width]
      return(fo)
    }
  }))
  setkey(out, chr, start, end)
  return(out)
}

