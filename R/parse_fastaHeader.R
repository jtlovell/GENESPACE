#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{make_inputFileMatrix} Utility function to build metadata
#'
#' @param map The map object (data.frame or data.table)
#' @param blk The block object (data.frame)
#' @param buffer Numeric, the overlapping distance between two blocks.
#' 0 indicates that blocks that overlap by >=0 should be merged.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Primarily used in the run_MCScanX pipeline.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
parse_fastaHeader = function(fasta.dir, is.peptide = T,
                             pattern = "fa", verbose = T){

  files = list.files(fasta.dir,
                      pattern = pattern,
                      full.names = T)

  if(verbose)
    cat("Renaming fasta headers ...\n")
  ss = lapply(files, function(i){
    if(verbose)
      cat("...",i,"\n\t")
    if(is.peptide){
      x = Biostrings::readAAStringSet(i)
    }else{
      x = Biostrings::readDNAStringSet(i)
    }
    if(verbose)
      cat("original names (e.g.):", names(x)[1])
    names(x)<-sapply(gsub(".*locus=","",names(x)),
                     function(y) strsplit(y," ")[[1]][1])
    if(verbose)
      cat("\n\tparsed names (e.g.):", names(x)[1],"\n")
    Biostrings::writeXStringSet(x, filepath = i)
  })
}

