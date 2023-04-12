#' @title classify the genome based on bed coordinates
#'
#' @description
#' \code{classify_genome} Hierarchically classify the genome by coordinates in
#' the list of bed files. Sequences in the first bed file are masked in the
#' second and so on.
#'
#' @param dnaSS DNAstringset object representing the genome assembly
#' @param listOfBeds list of bed-like data.frames/data.tables with at least
#' the columns chr, start and end.
#' @param verbose logical, should updates be printed to the console?
#'
#' @details intersects, gaps and reduces overlapping intervals
#'
#' @return a list of granges matchin the listOfBeds
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @export
classify_genome <- function(dnaSS,
                            listOfBeds,
                            verbose){

  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to classify genome, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("GenomeInfoDb", quietly = TRUE))
    stop("to classify genome, install GenomeInfoDb from bioconductor\n")

  if(length(listOfBeds) < 1)
    stop("listOfBeds must contain at least one element\n")

  if(!is.list(listOfBeds))
    stop("listOfBeds must be a list of data.frames\n")

  if(!all(sapply(listOfBeds, is.data.frame)))
    stop("all elements in listOfBeds must be a data.frame\n")

  if(!all(sapply(listOfBeds, function(x)
    all(c("start", "end", "chr") %in% colnames(x)))))
    stop("all elements in listOfBeds must have column names: chr, start, end\n")

  if(is.null(names(listOfBeds)))
    stop("listOfBeds must be a named list\n")

  if(sum(duplicated(names(listOfBeds))) > 0)
    stop("listOfBeds had duplicated names\n")

  ##############################################################################
  # 1. Read the sequence lengths
  if(verbose)
    cat(sprintf(
      " (%s Mb in %s chrs)\n",
      round(sum(GenomeInfoDb::seqlengths(dnaSS))/1e6, 2),
      length(dnaSS)))

  ##############################################################################
  # 2. Loop through the bedfiles
  if(verbose)
    cat(sprintf("Classifying the genome by %s bed files ... \n", length(listOfBeds)))
  maskThis <- NULL
  outList <- list()
  for(i in names(listOfBeds)){
    inGr <- join_ranges(listOfBeds[[i]],
                        seqInfo = GenomeInfoDb::seqinfo(dnaSS))
    if(is.null(maskThis)){
      maskThis <- inGr
      outGr <- inGr
    }else{
      maskThis <- GenomicRanges::reduce(maskThis)
      outGr <- get_gappedRanges(
        existingGr = maskThis, newGr = inGr)
      tmp <- GenomicRanges::union(outGr, maskThis)
      maskThis <- GenomicRanges::reduce(tmp)
    }
    outList[[i]] <- outGr
    if(verbose)
      cat(sprintf(
        "\t%sMb %s sequence (%sMb masked)\n",
        sum(round(sum(width(outGr))/1e6, 2)),
        i,
        sum(round(sum(width(maskThis))/1e6, 2))))
  }

  ##############################################################################
  # 3. Classify un-annotated sequence
  unk <- find_gaps(maskThis)
  if(verbose)
    cat(sprintf("\t%sMb un-annotated sequence\n",
                sum(round(sum(width(unk))/1e6, 2))))
  outList[["unknown"]] <- unk
  return(outList)
}
