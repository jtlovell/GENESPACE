#' @title Rename fasta header
#'
#' @description
#' \code{parse_fastaHeader} Rename fasta header
#'
#' @param fasta.dir file path, where fasta files are stored
#' @param is.peptide logical, are the sequences amino acids
#' @param pattern character, string in file names to find
#' @param gff data.table, with parsed gff-formatted annotations
#' @param only.gffGenes logical, should only genes in the gff be read? If true,
#' files MUST be named genomeID.fa, in order to match the gff genome identifier
#' with the gasta file.
#' @param parse_fastaHeader.FUN function, with rules to parse fasta header
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Reads, then overwrites the blast file
#'
#' @return Nothing, just reads, then overwrites the blast file
#'
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
#' @export
parse_fastaHeader <- function(fasta.dir,
                              is.peptide = T,
                              pattern = "fa",
                              gff = NULL,
                              only.gffGenes = TRUE,
                              parse_fastaHeader.FUN,
                              verbose = T){

  files <- list.files(fasta.dir,
                      pattern = pattern,
                      full.names = T)

  if (verbose)
    cat("Renaming fasta headers ...\n")
  ss <- lapply(files, function(i){
    if (verbose)
      cat("...", i, "\n\t")
    if (is.peptide) {
      x <- readAAStringSet(i)
    }else{
      x <- readDNAStringSet(i)
    }
    if (verbose)
      cat("original names (e.g.):",
          names(x)[1])
    names(x) <- sapply(names(x), parse_fastaHeader.FUN)
    if (verbose)
      cat("\n\tparsed names (e.g.):",
          names(x)[1],"\n")
    if(!is.null(gff) & only.gffGenes){
      g <- subset(gff, id %in% names(x) &
                    genome == gsub(".fa","",basename(i), fixed = T))
      x <- x[g$id]
    }
    writeXStringSet(x, filepath = i)
  })
}
