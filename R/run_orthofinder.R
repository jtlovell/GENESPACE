#' @title Run orthofinder from R
#'
#' @description
#' \code{run_orthofinder} Utility function to run orthofinder
#'
#' @param peptide.dir file.path, to the subdirectory containing
#' the parsed peptide files
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
                            output.dir,
                            blast.threads = 16,
                            og.threads = 1,
                            overwrite.output.dir = FALSE,
                            verbose = T){

  # -- Check the output directory
  check_gsDir(dir2check = output.dir,
              overwrite.output.dir = overwrite.output.dir)

  # -- Check the peptide directory
  if (!dir.exists(peptide.dir))
    stop("Cannot find",peptide.dir,"\n")

  # -- Make tmp directory
  tmp.dir <- file.path(getwd(),"gs.of.tmp")
  if (verbose)
    cat("Copying peptide fasta files to",
        tmp.dir, "\n")
  if (dir.exists(tmp.dir))
    unlink(tmp.dir,
           recursive = T)
  dir.create(tmp.dir)

  on.exit(
    expr = unlink(tmp.dir,
                  recursive = T))

  # -- Move peptide files
  pep.files <- list.files(peptide.dir,
                          full.names = T)
  nu <- file.copy(pep.files,
                  tmp.dir)

  # -- Run orthofinder
  if (verbose)
    cat("Running blasts within Orthofinder\n")
  com <- paste("orthofinder",
               "-f", tmp.dir,
               "-t", blast.threads,
               "-a", og.threads,
               "-S diamond -og")
  system(com)

  # -- Find location and copy blast files
  if (verbose)
    cat("Moving blasts results to",
        output.dir,
        "\n")
  files <- find_ofFiles(of.dir = tmp.dir)
  nu <- file.copy(from = files,
                  to = output.dir)

  # -- Decompress blast results
  if (verbose)
    cat("\tDecompressing blast results\n")
  system(paste("gunzip -f",
               file.path(output.dir,
                         "*.gz")))
  if (verbose)
    cat("\tDone!\n")
}
