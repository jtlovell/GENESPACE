#' @title match_syntenicChrs
#'
#' @description
#' \code{match_syntenicChrs} match_syntenicChrs
#'
#' @param map A data.table with at least the 10 necessary blast columns
#' @param gff A gff-like data.table
#' @param genome1.chrs names for chromosomes in the first element of the list
#' @param genomeIDs Character vector giving the genome IDs.
#' @param smallestchrsize2keep What is the largest chr that should be dropped?
#' @param ... Not currently in use
#' @details None yet

#' @return List of different genomeID levels of genes
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
match_synChrs <- function(map,
                          gff,
                          genomeIDs,
                          genome1.chrs = NULL,
                          smallestchrsize2keep = NULL){

  g1 <- genomeIDs[1]
  if(is.null(genome1.chrs)){
    g1c <- unique(c(map$chr1[map$genome1 == g1],
                    map$chr2[map$genome2 == g1]))
    tab <- table(subset(gff, genome == g1)$chr)
    if(!is.null(smallestchrsize2keep)){
      genome1.chrs <- names(tab)[tab > smallestchrsize2keep]
    }else{
      genome1.chrs <- names(tab)
    }

  }

  m1 <- subset(map,
               genome1 == g1 &
                 genome2 %in% genomeIDs[-1])
  m2 <- subset(map,
               genome2 == g1 &
                 genome1 %in% genomeIDs[-1])
  m <- data.table(genome1 = c(m1$genome1, m2$genome2),
                  genome2 = c(m1$genome2, m2$genome1),
                  chr1 = c(m1$chr1, m2$chr2),
                  chr2 = c(m1$chr2, m2$chr1),
                  id1 = c(m1$id1, m2$id2),
                  id2 = c(m1$id2, m2$id1))
  m <- m[!duplicated(m),]
  map.chrs <- m[,list(n.map = length(unique(id2))),
                by = list(genome1, genome2, chr1, chr2)]
  genome1.chrs <- genome1.chrs[genome1.chrs %in% map.chrs$chr1]

  map.chrs[,ns := -n.map]
  setkey(map.chrs, ns)
  map.chrs <- map.chrs[!duplicated(map.chrs[,c("genome2","chr2")]),]
  map.chrs$ns <- NULL

  if(!is.null(smallestchrsize2keep)){
    map.chrs <- subset(map.chrs, n.map >= smallestchrsize2keep)
  }

  map.chrs[,chr1 := factor(chr1, levels = genome1.chrs)]
  map.chrs[,genome2 := factor(genome2, levels = genomeIDs[-1])]
  setkey(map.chrs, chr1)
  l1 <- list(genome1.chrs)
  names(l1) <- g1
  return(c(l1, split(map.chrs$chr2, map.chrs$genome2)))
}
