#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{extend_blocks} Subset blast hits to syntenic regions and
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
#' @param clean.it Should the cleaning step be run?
#' @param verbose Logical, should updates be printed
#' @param ... Additional arguments passed to cull_syntenicBlast
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
extend_blocks <- function(map,
                          gff,
                          blast,
                          n.iter = 1,
                          rank.buffer = 250,
                          verbose = TRUE,
                          plotit = FALSE,
                          clean.it = TRUE,
                          n.cores = 1,
                          radius = 10,
                          n.mappings = 5,
                          ...){

  actually.plotit <- FALSE
  if (!clean.it) {
    n.iter <- 1
  }

  if (verbose)
    cat("Culling by pairwise genome comparison ...\n")
  for (i in 1:n.iter) {
    if (verbose & n.iter > 1)
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
                                  n.cores = n.cores,
                                  ...)
    toclean$rank1 <- NULL
    toclean$rank2 <- NULL
    if (clean.it) {
      out <- clean_blocks(map = toclean,
                          radius = radius,
                          n.mappings = n.mappings,
                          n.cores = n.cores)
      map <- out$map
    } else {
      setkey(blast, id1, id2)
      syn.ids <- toclean[ ,c("id1", "id2")]
      setkey(blast, id1, id2)
      setkey(syn.ids, id1, id2)
      out <- merge(syn.ids, blast)
    }
  }

  return(out)
}
