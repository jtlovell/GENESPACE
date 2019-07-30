#' @title Clean form_syntenicBlocks output
#'
#' @description
#' \code{clean_blocks} Clusters hits and drops low-confidence blocks.
#'
#' @param map the map data.table or data.frame
#' @param rerank logical, should the ranks be re-calculated prior to cleaning?
#' @param radius numeric, what should the radius of 2d density clustering be?
#' @param n.mappings numeric, how many mappings are required for a cluster?
#' @param clean.by.unique.genes Logical, should blocks with few unique genes
#' be culled?
#' @param clean.by.og Logical, should blocks with few unique orthogroups
#' be culled?
#' @param min.unique.genes numeric, legnth 1, if clean.by.unique.genes,
#' this is the number of unique genes needed for a block to be kept.
#' @param min.unique.og numeric, legnth 1, if clean.by.og,
#' this is the number of unique orthogroups needed for a block to be kept.
#' @param clean.columns logical, should extrac columns be dropped when
#' blocks are generated? Passed to make_blocks. Can speed things up.
#' @param n.cores The number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
clean_blocks <- function(map,
                         genomeIDs,
                         rerank = TRUE,
                         radius = 100,
                         n.mappings = 10,
                         verbose = TRUE,
                         by.blk = F,
                         final.n.mapping = 2,
                         final.radius = 2,
                         ...){

  #######################################################
  if ((length(radius) != length(n.mappings)))
    stop("radius and n.mappings must be of same length\n")
  #######################################################

  #######################################################
  # -- Prep the map object
  map <- subset(map, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
  map$unique <- with(map, paste(genome1, genome2, chr1, chr2))
  map$unique.genome <- with(map, paste(genome1, genome2))

  setkey(map, chr1, chr2, start1, start2)
  #######################################################

  #######################################################
  # -- Iteratively (or not) run the cleaning
  for (i in 1:length(n.mappings)) {
    n.map <- n.mappings[i]
    rad <- radius[i]
    if(verbose){
      if(length(n.mappings) == 1){
        cat("Cleaning mappings to", n.map, "hits within",
            rad, "gene-rank radius\n")
      }else{
        cat(paste0("Step",i,":"),
            "blocks must have", n.map, "hits within",
            rad, "gene-rank radius\n")
      }
    }

    cleaned <- clean_it(map = map,
                        genomeIDs = genomeIDs,
                        rerank = rerank,
                        radius = rad,
                        n.mappings = n.map,
                        verbose = verbose)
    map <- cleaned$map

    if (verbose)
      cat("Cleaned n blocks / mappings =",
          nrow(cleaned$block),
          "/",
          nrow(cleaned$map),"\n")
  }
  #######################################################

  if (by.blk) {
    if(length(final.n.mapping) != length(final.radius)){
      final.n.mapping <- rep(final.n.mapping[1], length(final.radius))
    }

    if(verbose)
      cat("Cleaning within blocks ...\n")
    for(i in 1:length(final.radius)){
      if(verbose)
        cat("\tIteration",i,final.n.mapping[i],"hits in",final.radius[i],"hit radius ...")
      cleaned <- clean_it(map = cleaned$map,
                           genomeIDs = genomeIDs,
                           rerank = T,
                           radius = final.radius[i],
                           n.mappings = final.n.mapping[i],
                           verbose = F,
                           by.blk = T)
      if(verbose)
        cat(" kept", nrow(cleaned$map),"hits in", nrow(cleaned$block),"blocks\n")
    }
    if(verbose)
      cat("\tDone!\n")
  }
  cleaned$map[,block.id := paste0("blk_",as.numeric(as.factor(paste(genome1, genome2, chr1, chr2, block.id))))]
  #######################################################
  return(make_blocks(map = cleaned$map,
                     rename.blocks = F))
}

