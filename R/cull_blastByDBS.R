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
#' @export
cull_blastByDBS <- function(blast,
                            n.mappingWithinRadius = c(5,5,5),
                            eps.radius = c(100,50,25),
                            verbose = T,
                            run.it = T,
                            ...){
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

  map$unique <- with(map, paste(genome1, genome2, chr1, chr2))
  map$unique1 <- with(map, paste(genome1, genome2))

  map[,rank1 := frank(start1,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr1)]
  map[,rank2 := frank(start2,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr2)]

  spl <- split(map, "unique1")
  out <- rbindlist(lapply(spl, function(x){
    if (verbose)
      cat(paste0("\t",
                 x$genome1[1]), "-->", x$genome2[1],
          paste0("(initial hits = ", nrow(x),") "))

    if (run.it) {
      tmp <- split(x, "unique")
      xo = rbindlist(lapply(tmp, function(z){
        z <- run_dbs(y = z,
                     eps.radius = eps.radius,
                     mappings = n.mappingWithinRadius)
        z <- z[z$cluster != 0,]
        return(z)
      }))
    } else {
      xo <- rbindlist(x)
    }
    if (verbose)
      cat("culled hits =", nrow(xo), "\n")
    return(xo)
  }))
  out[,rank1 := frank(start1,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr1)]
  out[,rank2 := frank(start2,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr2)]
  #######################################################

  #######################################################
  if (verbose)
    cat("\tDone!\n")
  return(out)
}
