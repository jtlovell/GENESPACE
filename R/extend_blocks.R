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
#' @param rank.buffer The buffer, in gene rank order.
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
                          genomeIDs,
                          n.iter = 1,
                          rank.buffer = 250,
                          verbose = TRUE,
                          clean.it = FALSE,
                          radius = 10,
                          n.mappings = 5){

  if (verbose)
    cat("Culling by pairwise genome comparison ...\n")
  for (i in 1:n.iter) {
    if (verbose & n.iter > 1)
      cat("Iteration", i, ":\n")

    toclean <- cull_syntenicBlast(gff = gff,
                                  map = map,
                                  blast = blast,
                                  rank.buffer = rank.buffer,
                                  verbose = verbose)
    toclean <- toclean[!duplicated(toclean[,c("id1","id2")]),]
    out <- make_blocks(toclean, rename.blocks = F)
    if (clean.it) {
      out <- clean_blocks(map = out$map,
                          radius = radius,
                          genomeIDs = genomeIDs,
                          n.mappings = n.mappings)
    }
    map <- data.table(out$map)

    wh1 <- grep("1$", colnames(map))
    wh2 <- grep("2$", colnames(map))
    who <- which(!grepl("1$|2$", colnames(map)))

    n1 <- colnames(map)[wh1]
    n2 <- colnames(map)[wh2]
    no <- colnames(map)[who]
    m1 <- map[,c(n1,n2,no), with = F]
    m2 <- map[,c(n2,n1,no), with = F]
    setnames(m2, c(n1,n2,no))
    mt <- rbind(m1, m2)

    map <- make_blocks(mt[!duplicated(mt[,c("id1","id2")]),], clean.columns = F)
  }
  ids.syn <- paste(map$id1, map$id2)
  ns <- subset(blast, !(paste(id1, id2) %in% ids.syn |
                         paste(id2, id1) %in% ids.syn))

  return(list(map = out$map,
              block = out$block,
              non.syn.blast = ns))
}
