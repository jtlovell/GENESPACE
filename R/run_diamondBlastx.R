#' @title Run diamond blastx
#'
#' @description
#' \code{run_diamondBlastx} A simple wrapper to run diamond blastx
#'
#' @param db.file name of diamond database file
#' @param pep.fa peptide fasta file path
#' @param fa.file assembly fasta file path
#' @param blast.file output blast8-formatted txt file
#' @param max.target.seqs numeric, the maximum number of blast hits
#' to return
#' @param min.score numeric, the minimum blast score to retain
#' @param diamond.blastx.param Additional parameters to pass to
#' diamond blastx
#' @param ... Not currently in use
#' @details  ...

#' @return A data.table with the same format as bed.dt, but with an
#' additional column with the file path to the fasta file.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
run_diamondBlastx <- function(db.file,
                              pep.fa,
                              fa.file,
                              blast.file,
                              min.score = 20,
                              top = 10,
                              diamond.blastx.param = "--quiet",
                              ...){
  system(paste("diamond makedb --quiet",
               "--in", pep.fa,
               "-d", db.file))
  system(paste("diamond blastx",
               "--top", top,
               "--min-score", min.score,
               diamond.blastx.param,
               "-d", db.file,
               "-q", fa.file,
               "-o", blast.file))
}
