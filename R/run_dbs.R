#' @title Clusters hits by dbscan
#'
#' @description
#' \code{run_dbs} Clusters hits by dbscan
#'
#' @param y the map data.table or data.frame, or any object with two columns,
#' 'rank1' and 'rank2'
#' @param eps.radius numeric, what should the radius of 2d density clustering be?
#' @param mappings numeric, how many mappings are required for a cluster?
#' @param n.cores The number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A data.frame / data.table with an additional cluster column.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
run_dbs <- function(y,
                    eps.radius,
                    mappings){

  nn <- frNN(data.frame(x = y$rank1, y = y$rank2),
             eps = eps.radius)
  dbs <- dbscan(nn,
                minPts = mappings)
  y$cluster <- dbs$clusters
  return(y)
}
