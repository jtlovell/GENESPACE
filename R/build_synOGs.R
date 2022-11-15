#' @title construct syntenic orthogroups from pairwise blast
#' @description
#' \code{build_synOGs} integrates many pairwise results from synteny into
#' vectors of syntenic orthogroups. Also can re-run orthofinder within blocks
#' and re-calculate orthogroups from those results.
#'
#' @param gsParam A list of genespace parameters created by init_genespace.

#' @details info here

#' @import data.table
#' @export
build_synOGs <- function(gsParam){

  pull_synOgs <- function(gsParam, onlyInBuffer = TRUE){
    blkID <- sameOg <- sameInblkOg <- inBuffer <- blkID <- NULL
    ##############################################################################
    # -- 1. Get metadata together
    md <- data.table(gsParam$synteny$blast)
    # -- for each line in metadata
    hitsInOgs <- rbindlist(lapply(1:nrow(md), function(i){
      allBlast <- fread(
        md$synHits[i], showProgress = F, na.strings = c("NA", ""))

      # -- subset to only hits in regions
      out <- subset(allBlast, !is.na(blkID) & sameOg)
      out[,sameInblkOg := TRUE]
      if(onlyInBuffer)
        out <- subset(out, inBuffer)
      out <- subset(out, (sameOg | sameInblkOg))
      out <- out[,c("ofID1", "ofID2")]
      return(out)
    }))
    return(hitsInOgs)
  }

  # -- read in the combined bed file
  bed <- read_combBed(gsParam$synteny$combBed)

  ofID <- ofID1 <- ofID2 <- synOG <- bedFile <- NULL
  soh <- pull_synOgs(gsParam = gsParam)
  ic <- with(soh, clus_igraph(id1 = ofID1, id2 = ofID2))
  bed[,synOG := ic[ofID]]
  whna <- which(is.na(bed$synOG))
  mx <- max(bed$synOG, na.rm = T)
  bed$synOG[whna] <- paste((mx + 1): (mx + length(whna)))

  write_combBed(x = bed, filepath = gsParam$synteny$combBed)
  return(gsParam)
}



