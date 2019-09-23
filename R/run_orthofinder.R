#' @title Run orthofinder from R
#'
#' @description
#' \code{run_orthofinder} Utility function to run orthofinder
#'
#' @param peptide.dir file.path, to the subdirectory containing
#' the parsed peptide files
#' @param tmp.dir file.path, to the subdirectory
#' where the temporary results should be written
#' @param output.dir file.path, to the subdirectory
#' where the results should be written
#' @param blast.threads Integer, specifying the number of parallel threads to use for
#' the blasts algorithm
#' @param og.threads Integer, specifying the number of parallel threads to use for
#' the orthofinder algorithm
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#'
#' @details ...
#'
#' @return Nothing
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
run_orthofinder <- function(peptide.dir,
                            tmp.dir,
                            output.dir,
                            blast.threads = 16,
                            og.threads = 1,
                            verbose = T){
  if (verbose)
    cat("Copying peptide fasta files to",
        tmp.dir, "\n")
  if (dir.exists(tmp.dir))
    unlink(tmp.dir, recursive = T)
  dir.create(tmp.dir)
  pep.files <- list.files(peptide.dir,
                          full.names = T)
  nu <- file.copy(pep.files,
                  tmp.dir)
  if (verbose)
    cat("Running blasts within Orthofinder\n")
  com <- paste("orthofinder",
               "-f", tmp.dir,
               "-t", blast.threads,
               "-a", og.threads,
               "-S diamond -og")
  system(com)

  if (verbose)
    cat("Moving blasts results to",
        output.dir,
        "\n")
  blast.loc <- dirname(list.files(tmp.dir,
                                  pattern = "SequenceIDs",
                                  recursive = T,
                                  full.names = T)[1])
  ortho.loc <- dirname(list.files(tmp.dir,
                                  pattern = "Orthogroups.txt",
                                  recursive = T,
                                  full.names = T)[1])
  blast.files <- list.files(blast.loc,
                            pattern = "Blast*",
                            full.names = T)
  fa.files <- list.files(blast.loc,
                         pattern = "Species*",
                         full.names = T)
  fa.files <- fa.files[grep(".fa$", fa.files)]
  dmnd.files <- list.files(blast.loc,
                           pattern = "diamondDBSpecies*",
                           full.names = T)
  og.files <- file.path(ortho.loc,
                        "Orthogroups.txt")
  sp.id.files <- file.path(blast.loc, "SpeciesIDs.txt")
  seq.id.files <- file.path(blast.loc, "SequenceIDs.txt")
  files <- c(blast.files,
             fa.files,
             dmnd.files,
             og.files,
             sp.id.files,
             seq.id.files)
  if (dir.exists(output.dir))
    unlink(output.dir, recursive = T)
  dir.create(output.dir)
  nu <- file.copy(files,
                  output.dir)
  if (verbose)
    cat("\tDecompressing blast results\n")
  system(paste("gunzip -f",
               file.path(output.dir,
                         "*.gz")))
}
