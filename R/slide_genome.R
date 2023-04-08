#' @title sliding window on a list of granges
#'
#' @description
#' \code{slide_genome} Hierarchically classify the genome by coordinates in
#' the list of granges objects, returning the proportion of sequence each
#' sliding window intersects
#'
#' @param seqInfo seqinfo extracted from the genome assembly DNAstringset
#' @param listOfGrs list of granges containing the regions to intersect
#' @param windowSize numeric, giving the size in basepairs of the windows
#' @param stepSize numeric, giving the step in basepairs between the windows
#'
#' @details intersects, gaps and reduces overlapping intervals within windows
#'
#' @return a data.table containing the amount of sequence covered in each window
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import data.table
#' @importFrom parallel mclapply
#' @export
slide_genome <- function(seqInfo, listOfGrs, windowSize, stepSize){

  if(!requireNamespace("GenomicRanges", quietly = TRUE))
    stop("to slide genome, install GenomicRanges from bioconductor\n")
  if(!requireNamespace("BiocGenerics", quietly = TRUE))
    stop("to slide genome, install BiocGenerics from bioconductor\n")

  seqnames <- chr <- start <- end <- i.end <- i.start <- id <- nBpWindow <-
    propWind <- nBp <- id <- cumProp <- NULL

  # -- get the windows from the chromosome sizes
  siGr <- convert_si2gr(seqInfo)
  sw <- GenomicRanges::slidingWindows(
    siGr,
    width = windowSize,
    step = stepSize)
  sws <- BiocGenerics::do.call(c, sw)

  # -- convert to data.table
  swDt <- data.table(as.data.frame(sws))
  swDt[,`:=`(chr = seqnames, seqnames = NULL,
             strand = NULL, width = NULL)]
  setkey(swDt, chr, start, end)

  # -- count bases by each granges
  outDt <- rbindlist(lapply(names(listOfGrs), function(i){
    xdt <- with(as.data.frame(listOfGrs[[i]]), data.table(
      chr = seqnames,
      start = start,
      end = end,
      key = c("chr", "start", "end")))
    fo <- foverlaps(
      xdt, swDt,
      type = "any", mult = "all")

    fo$i.end[fo$i.end > fo$end] <- fo$end[fo$i.end > fo$end]
    fo$i.start[fo$i.start < fo$start] <- fo$start[fo$i.start < fo$start]

    xout <- fo[,list(
      nBp = sum((i.end + 1) - i.start)),
      by = c("chr", "start", "end")]

    xout <- merge(
      xout, swDt,
      by = c("chr", "start", "end"), all.y = T)
    xout$nBp[is.na(xout$nBp)] <- 0
    xout[,id := i]
    return(xout)
  }))

  outDt[,nBpWindow := (1 + end) - start]
  outDt[,propWind := nBp / nBpWindow]

  outDt[,id := factor(id, levels = names(listOfGrs))]
  setkey(outDt, chr, start, end, id)

  outDt[,cumProp := cumsum(propWind),
        by = c("chr", "start", "end")]
  outDt$cumProp[outDt$cumProp > 1] <- 1

  return(outDt)
}
