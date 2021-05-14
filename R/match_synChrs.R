#' @title Match syntenic chrs
#'
#' @description
#' \code{match_synChrs} Match syntenic chrs
#'
#' @param gsAnnot list of length 2, containing the genespace annotation paths
#' 'gff' and 'peptide' -- file path character vector with the locations of
#' the peptide and gff-like annotation files. Each element is named by the
#' associated genomeID in gsParam
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param blkFile file.path pointing to the blk coordinate text file.
#' @param refGenome character string matching one of the genomeIDs.
#' @param minGenes2plot numeric of length 1 specifying the minimum number of
#' genes on a chromosome to be included in the list.
#' @param refChrs if specified, the order of chromosomes in the ref
#'
#' @details Build
#'
#' @return A
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @import R.utils
#' @import data.table
#' @export

