#' @title Simple orthogroup-constrained tandem array inference.
#'
#' @description
#' \code{find_arrayClusters} Simple orthogroup-constrained tandem array inference
#' @param map the blast dataset to screen for syntenic hits
#' @param rerank Logical, should ranks be remade before each step?
#' @param clean.radius Passed on to clean_blocks
#' @param clean.mappings Passed on to clean_blocks
#' @param merge.buffer Passed on to merge_blocks
#' @param min.hits1 The minimum amount of hits in either genome
#' to be counted as an array
#' @param min.hits2 The minimum amount of hits in the over-represented
#' genome to be counted as an array
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
