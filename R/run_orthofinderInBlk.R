#' @title run_orthofinderInBlk
#'
#' @description
#' \code{run_orthofinderInBlk} run_orthofinderInBlk
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#'
#' @details xxx
#'
#' @return xxx
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame clusters
#' @export
run_orthofinderInBlk <- function(gsParam){
  ##############################################################################
  # -- 1. Get metadata together

  # -- 1.1 synHit files (0-0, 1-1, etc.)

  # -- 1.2 peptide files

  # -- 1.3 tmp directories

  # -- 1.4 write peptides

  # -- 1.5 prep files (-op)

  ##############################################################################
  # -- 2. write the blast files to file
  # -- 2.1 subset to just hits of interest for 0-1, 0-0, etc.
  inreg <- subset(hits, !is.na(regionID) & inBuffer)

  # -- 2.2 get new ids from the tmp directories

  # -- 2.3 make dictionary

  # -- 2.4 re-name hits

  # -- 2.5 split hits by region
  spl <- split(inreg)

  # -- 2.6 write orthofinder-formatted blast to the correct directories

  ##############################################################################
  # -- 3. get within block orthofinder results
  # -- 3.1 call orthofinder for each region

  # -- 3.2 read in and parse HOGs for each region

  ##############################################################################
  # -- 4. combine orthogroups
  # -- 4.1 make pairwise graph of HOGs

  # -- 4.2 parse via igraph to single cluster vector

  # -- 4.3 add to inblkOG column in combBed and replace og with this

  ##############################################################################
  # -- 5. re-run synteny using the new inBlkOG column as og
}
