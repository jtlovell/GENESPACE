#' @title Run dbscan
#'
#' @description
#' \code{cull_blastByDBS} Cull blast hits by 2d density of gene order ranks
#'
#' @param blast blast results data.table
#' @param n.mappingWithinRadius numeric, number of hits required to be in the radius
#' @param eps.radius numeric, size of the radius
#' @param run.it logical, should dbscan be run, or just order the genes?
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return culled blast data.table
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
cull_blastByDBS <- function(blast,
                            n.mappingWithinRadius = c(5,5,5),
                            eps.radius = c(100,50,25),
                            verbose = T,
                            run.it = T){
  #######################################################
  #######################################################
  run_dbs <- function(y,
                      eps.radius,
                      mappings){
    nn <- frNN(data.frame(y[, c("rank1", "rank2"), with = F]),
               eps = eps.radius)
    dbs <- dbscan(nn,
                  minPts = mappings)
    y$cluster <- dbs$cluster
    return(y)
  }
  #######################################################
  #######################################################

  #######################################################
  blast.cull <- blast
  if (verbose)
    cat("Culling", nrow(blast.cull),
        "BLAST hits by 2d Density\n")

  if (length(n.mappingWithinRadius) != length(eps.radius)) {
    warning("eps.radius and n.mappingWithinRadius not of same length")
    if (length(n.mappingWithinRadius) > length(eps.radius)) {
      n.mappingWithinRadius <- n.mappingWithinRadius[1:length(eps.radius)]
    }else{
      eps.radius <- eps.radius[1:length(n.mappingWithinRadius)]
    }
  }
  #######################################################

  #######################################################
  map <- data.table(blast.cull)
  map$unique <- with(map, paste(genome1, genome2))
  spl <- split(map, "unique")
  out <- rbindlist(lapply(spl, function(x){
    if (verbose)
      cat(paste0("\t",
                 x$genome1[1]), "-->", x$genome2[1],
          paste0("(initial hits = ", nrow(x),") "))

    if (run.it){
      for (i in 1:length(eps.radius)) {
        x$rank1 <- frank(x, "chr1", "start1",
                         ties.method = "dense")
        x$rank2 <- frank(x, "chr2", "start2",
                         ties.method = "dense")
        x <- run_dbs(y = x,
                     eps.radius = eps.radius[i],
                     mappings = n.mappingWithinRadius[i])
        x <- x[x$cluster != 0,]
        if (nrow(x) < min(n.mappingWithinRadius)) {
          break
        }
      }
    } else {
      x$rank1 <- frank(x, "chr1", "start1",
                       ties.method = "dense")
      x$rank2 <- frank(x, "chr2", "start2",
                       ties.method = "dense")
    }
    if (verbose)
      cat("culled hits =", nrow(x), "\n")
    return(x)
  }))
  #######################################################

  #######################################################
  if (verbose)
    cat("\tDone!\n")
  return(out)
}
