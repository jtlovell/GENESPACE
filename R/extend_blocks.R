#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param map The map data.frame or data.table
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param blast the blast dataset to screen for syntenic hits
#' @param n.iter Number of iterations to run
#' @param plotit Logical, should plots be made? Will not work with
#' n.core > 1.
#' @param rank.buffer The buffer, in gene rank order.
#' @param n.cores Number of parallel processes to run, when possible
#' @param radius Numeric, length of 1, specifiying the
#' radius to search for hits. Passed to clean_blocks
#' @param n.mappings Numeric, length of 1, specifiying the
#' minimum number of hits in a radius. Passed to clean_blocks
#' @param verbose Logical, should updates be printed
#' @param ... Not currently in use
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom parallel mclapply
#' @export
extend_blocks <- function(map,
                          gff,
                          blast,
                          n.iter = 3,
                          rank.buffer = 250,
                          verbose = TRUE,
                          plotit = FALSE,
                          n.cores = 1,
                          radius = 10,
                          n.mappings = 5){
  actually.plotit <- FALSE
  if (verbose)
    cat("Culling by pairwise genome comparison ...\n")
  for (i in 1:n.iter) {
    if (verbose)
      cat("Iteration", i, ":\n")

    if (plotit & i == n.iter) {
      actually.plotit <- TRUE
    }
    toclean <- cull_syntenicBlast(gff = gff,
                                  map = map,
                                  blast = blast,
                                  plotit = actually.plotit,
                                  rank.buffer = rank.buffer,
                                  verbose = verbose,
                                  n.cores = n.cores)
    toclean <- toclean
    toclean$rank1 <- NULL
    toclean$rank2 <- NULL
    out <- clean_blocks(map = toclean,
                        radius = radius,
                        n.mappings = n.mappings,
                        n.cores = n.cores)
    map <- out$map
  }

  if (verbose)
    cat("Done!\n")

  return(out)
}
