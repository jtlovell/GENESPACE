#' @title Clean form_syntenicBlocks output
#'
#' @description
#' \code{clean_blocks} Clusters hits and drops low-confidence blocks.
#'
#' @param map the map data.table or data.frame
#' @param blk the block data.table or data.frame
#' @param rerank logical, should the ranks be re-calculated prior to cleaning?
#' @param radius numeric, what should the radius of 2d density clustering be?
#' @param n.mappings numeric, how many mappings are required for a cluster?
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
#' @import parallel
#' @import dbscan
#' @export
clean_blocks <- function(blk,
                         map,
                         rerank = TRUE,
                         radius = NULL,
                         n.mappings = 5,
                         n.cores = 1,
                         verbose = T){

  blk <- data.table(blk)
  map <- data.table(map)

  if (is.null(radius))
    radius <- .1 + ceiling(sqrt((n.mappings^2) * 2))

  if (verbose)
    cat("Cleaning / merging via dbscan\n\tInitial n. blocks / mappings =",
        nrow(blk),
        "/",
        nrow(map),"\n\t")

  map$unique = with(map, paste(genome1, genome2, chr1, chr2))

  setkey(map, chr1, chr2, start1, start2)

  if (rerank) {
    map[, rank1 := frank(start1,
                         ties.method = "random"),
        by = list(genome1, genome2, chr1)]
    map[, rank2 := frank(start2,
                         ties.method = "random"),
        by = list(genome1, genome2, chr2)]
  }

  spl.map <- split(map, "unique")

  merged_map <- rbindlist(mclapply(spl.map, mc.cores = n.cores, mc.preschedule = F, function(tmp){
    x <- run_dbs(y = tmp[, c("rank1","rank2"), with = F],
                 eps.radius = radius,
                 mappings = n.mappings)
    tmp$block.id <- x$cluster
    return(tmp)
  }))
  merged_map <- merged_map[merged_map$block.id != 0, ]

  merged_map$block.id <- with(merged_map,
                              as.numeric(as.factor(paste(unique, block.id))))
  merged_blk <- make_blocks(merged_map,
                            rename.blocks = F,
                            rerank = T)

  if (verbose)
    cat("Cleaned n blocks / mappings =",
        nrow(merged_blk$block),
        "/",
        nrow(merged_blk$map),"\n\tDone!\n")

  return(merged_blk)
}
