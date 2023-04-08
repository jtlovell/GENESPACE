#' @title Generic internal functions used by genespace
#' @description
#' \code{gviz_utils} Convenience functions for gscTools, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name gviz_utils
#'
#' @param bed data.table or data.frame with at least chr, start, end columns
#' @param dnass DNAStringSet
#' @param seqInfo seqInfo ranges extracted from dnass
#' @param gr granges object
#' @param existingGr existing granges to mask with
#' @param newGr new granges to mask
#' @param gapGr granges with gaps
#' \cr
#' If called, \code{gviz_utils} returns its own arguments.
#'
#'
#' @title join_ranges
#' @description
#' \code{join_ranges} join_ranges
#' @rdname gviz_utils
#' @export
join_ranges <- function(bed, seqInfo){

  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to join ranges, install IRanges from bioconductor\n")
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to join ranges, install GenomicRanges from bioconductor\n")

  gr <- with(bed, GenomicRanges::GRanges(
    chr,
    IRanges::IRanges(start = start, end = end),
    seqinfo = seqInfo))
  return(GenomicRanges::reduce(gr))
}

#' @title pull_seqInfo
#' @description
#' \code{pull_seqInfo} pull_seqInfo
#' @rdname gviz_utils
#' @export
pull_seqInfo <- function(dnass){

  if(!requireNamespace("GenomeInfoDb", quietly = TRUE))
    stop("to pull seqInfo, install GenomeInfoDb from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to pull seqInfo, install BiocGenerics from bioconductor\n")

  ss <- dnass
  fai <- data.table(
    chr = names(ss),
    start = 1,
    end = BiocGenerics::width(ss))
  seqi <- GenomeInfoDb::Seqinfo(fai$chr, seqlengths = fai$end)
  return(seqi)
}

#' @title find_gaps
#' @description
#' \code{find_gaps} find_gaps
#' @rdname gviz_utils
#' @export
find_gaps <- function(gr){
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to find gaps, install GenomicRanges from bioconductor\n")
  strand <- NULL
  gapr <- GenomicRanges::gaps(GenomicRanges::reduce(gr))
  return(GenomicRanges::reduce(subset(gapr, strand == "*")))
}

#' @title subset_bedInGap
#' @description
#' \code{subset_bedInGap} subset_bedInGap
#' @rdname gviz_utils
#' @export
subset_bedInGap <- function(gr, gapGr){
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to subset bed within gaps, install GenomicRanges from bioconductor\n")
  return(GenomicRanges::intersect(
    x = GenomicRanges::reduce(gr),
    y = GenomicRanges::reduce(gapGr)))
}

#' @title get_gappedRanges
#' @description
#' \code{get_gappedRanges} get_gappedRanges
#' @rdname gviz_utils
#' @export
get_gappedRanges <- function(existingGr, newGr){
  ingap <- subset_bedInGap(
    gr = newGr,
    gapGr = find_gaps(existingGr))
  return(ingap)
}

#' @title convert_si2gr
#' @description
#' \code{convert_si2gr} convert_si2gr
#' @rdname gviz_utils
#' @export
convert_si2gr <- function(seqInfo){
  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to convert2granges, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("IRanges", quietly = TRUE))
    stop("to convert2granges, install IRanges from bioconductor\n")
  if(!requireNamespace("GenomeInfoDb", quietly = TRUE))
    stop("to convert2granges, install GenomeInfoDb from bioconductor\n")
  gr <- GenomicRanges::GRanges(
    names(seqInfo),
    IRanges::IRanges(start = 1, end = GenomeInfoDb::seqlengths(seqInfo)),
    seqinfo = seqInfo)
  return(GenomicRanges::reduce(gr))
}

#' @title rename chromosomes
#' @description
#' \code{rename_chrs} Parse common formats for chromosomes, specifically
#' splitting whitespace separated fields and finding the field like to give
#' the chromosome ID
#' @rdname gviz_utils
#' @import data.table
#' @export
rename_chrs <- function(dnass){

  ss <- dnass
  ns <- as.data.table(tstrsplit(trimws(names(ss)), " "))
  hasna <- sapply(ns, function(x) any(is.na(x)))
  if(all(hasna)){
    return(ss)
  }

  ns <- ns[,!hasna, with = F]
  if(ncol(ns) == 1){
    # -- if one column (no whitespace) return names
    chrns <- ns[[1]]
  }else{
    # -- if the first column is all unique and the first entries are "chrx" or numeric, use it
    u1 <- ns[[1]]
    hasChr <- gsub("chr|chromosome|scaffold|chr0", "", tolower(u1)[[1]]) == "1"
    if(!any(duplicated(u1)) && hasChr){
      chrns <- ns[[1]]
    }else{
      # -- if there is a column with "chr" in it and it is all unique use that
      wh <- grepl("chr", tolower(unlist(ns[1,])))
      us <- apply(ns, 2, function(x) all(!duplicated(x)))
      if(sum(wh & us) == 1){
        chrns <- as.character(ns[[wh & us]])
      }else{
        # -- if there is a column with "chr", followed by one totally unique, use that
        if(sum(wh) == 1){
          index <- which(wh) + 1
          if(all(!duplicated(ns[[index]]))){
            chrns <- ns[[index]]
          }else{
            # -- otherwise just return the names
            chrns <- names(ss)
          }
        }else{
          chrns <- names(ss)
        }
      }
    }
  }
  names(ss) <- chrns
  return(ss)
}

