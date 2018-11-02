#' @title Utility to parse the header line of fasta file
#'
#' @description
#' \code{parse_fastaHeader} Returns a fasta file with the headers
#' specifying only the locus name and no other info.
#'
#' @param fasta.dir The path to the directory holding the fasta files to parse.
#' @param is.peptide Logical, should the sequence be treated as a peptide or DNA
#' @param pattern Character, the pattern to grep in the files.
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
                             pattern = "fa", verbose = T, ...){

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

