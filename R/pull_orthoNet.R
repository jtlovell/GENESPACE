#' @title Pull out gene ids of orthogroups
#'
#' @description
#' \code{pull_orthoNet} Pull out gene ids of orthogroups
#'
#' @param blast A data.table with at least the 10 necessary blast columns
#' @param drop.tandems Logical, should genes in tandem arrays be condensed
#' into a single representatitve hit?
#' @param genomeIDs Character vector giving the genome IDs.
#' @param verbose Logical, should updates be printed?
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
pull_orthoNet <- function(blast,
                          genomeIDs,
                          drop.tandems = F,
                          verbose = T){

  blast <- blast[blast$genome1 %in% genomeIDs &
                   blast$genome2 %in% genomeIDs,]

  if (!"array.id" %in% colnames(blast)) {
    drop.tandems <- FALSE
  }

  if (drop.tandems) {
    if (verbose)
      cat("Condensing tandem arrays into single hit with highest score ... ")
    blast$is.ta <- !is.na(blast$array.id)
    blta <- blast[blast$is.ta,]
    blno <- blast[!blast$is.ta,]
    blta$score1 <- blta$score * (-1)
    setkey(blta, genome1, genome2, score1)
    blta.out <- blta[, head(.SD[1]), list(genome1, genome2, array.id)]
    blta.out <- blta.out[,colnames(blast), with = F]
    blast <- rbind(blno, blta.out)
    if (verbose)
      cat("Done!\n")
  }

  if (verbose)
    cat("Building database of unique hits across genomes ... ")
  bl <- with(blast,
             data.table(genome = c(genome1, genome2),
                        id = c(id1, id2),
                        og = c(og1, og2),
                        score = c(score, score)))
  bl <- bl[!duplicated(bl[,-4,with = F]),]
  bl$genome <- factor(bl$genome, levels = genomeIDs)
  if (verbose)
    cat("Done!\n")

  if (verbose)
    cat("Ranking hits by score ... ")
  bl[,rank := frank(score*(-1), ties.method = "random"),
     by = list(genome, og)]
  if (verbose)
    cat("Done!\n")

  if (verbose)
    cat("Counting hits by genome and orthogroup ... ")
  bll <- bl[,list(count = .N),
            by = list(genome, og)]

  out <- dcast(bll, og ~ genome, value.var = "count")
  out[is.na(out)] <- 0

  out$comb <- apply(out[,-1,with = F], 1, function(x) paste(x, collapse = "_"))
  if (verbose)
    cat("Done!\n")

  if (verbose)
    cat("Building database of orthogroups, split by number of hits ... ")
  tab <- table(out$comb)
  ns <- names(tab)[order(-tab)]
  out$comb <- factor(out$comb, levels = ns)
  setkey(out, comb)
  spl <- split(out, by = "comb")
  names(spl) <- ns

  outl <- lapply(spl, function(x){
    tmp <- bl[bl$og %in% x$og,]
    tmpw <- dcast(tmp, og ~ genome + rank, value.var = "id")
    return(tmpw)
  })
  if(verbose)
    cat("Done!\n")

  return(outl)
}
