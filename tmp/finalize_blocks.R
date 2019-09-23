#' @title Clean form_syntenicBlocks output
#'
#' @description
#' \code{finalize_blocks} Clusters hits and drops low-confidence blocks.
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
finalize_blocks <- function(map,
                            gff,
                            genomeIDs,
                            dir.list,
                            max.blockSize2merge,
                            extend.rank.buffer,
                            min.blockSize,
                            MCScanX.m.param){
  comp<- complete_graph(map = map,
                        gff = gff,
                        verbose = TRUE)
  exted <- extend_mcsBlks(map = map,
                          dir.list = dir.list,
                          gff = gff,
                          genomeIDs = genomeIDs,
                          radius = extend.rank.buffer,
                          min.blockSize = min.blockSize,
                          m.param = MCScanX.m.param)

  pexted <- proc_MCScanBlocks(
    map = exted$map,
    genomeIDs = genomeIDs,
    max.blockSize2merge = max.blockSize2merge,
    gff = gff)

  return(pexted)
}
