#' @title Run the orthofinder program
#'
#' @description
#' \code{run_orthofinder} A simple wrapper to run orthofinder from R.
#'
#' @param peptide.dir The path to the directory containing the peptide fasta sequence files
#' @param blast.dir The path to the directory where the blast results should be stored
#' @param tmp.dir The path to the directory where temporary files will be stored
#' then deleted
#' @param blast.threads The number of threads for blast run
#' @param og.threads The number of threads used for orthogroup construction
#' @param og.silent Should orthofinder stout be suppressed?
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details To run successfully, the peptide directory must be populated with a set of
#' primary gene model peptide sequences. R must also be called from an environment
#' with orthofinder and diamond program paths specified. Diamond paramters (e.g. number of
#' cores to use, etc.) must be specified within orthofinders parameter file.

#' @return Nothing, writes results to the blast.dir directory
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
run_orthofinder <- function(peptide.dir = NULL,
                            tmp.dir,
                            blast.dir,
                            blast.threads = 16,
                            og.threads = 1,
                            verbose = T,
                            og.silent = F,
                            ...){

  if (!is.null(peptide.dir)) {

    if (verbose)
      cat("Copying peptide fasta files to",
          tmp.dir,
          "\n")

    if (dir.exists(tmp.dir))
      unlink(tmp.dir, recursive = T)
    if (!dir.exists(tmp.dir))
      dir.create(tmp.dir)

    pep.files <- list.files(peptide.dir,
                            full.names = T)
    nu <- file.copy(pep.files,
                    tmp.dir)

    if (verbose)
      cat("Running blasts within Orthofinder\n")
    if (og.silent) {
      system(paste("orthofinder -f", tmp.dir,
                   "-t", blast.threads,
                   "-a", og.threads,
                   "-S diamond -og 1>/dev/null 2>&1"))
    }else{
      system(paste("orthofinder -f", tmp.dir,
                   "-t", blast.threads,
                   "-a", og.threads,
                   "-S diamond -og"))
    }
  }else{
    if (og.silent) {
      system(paste("orthofinder -b", tmp.dir,
                 "-t", og.threads,
                 "-a", og.threads,
                 "-S diamond -og 1>/dev/null 2>&1"))
    }else{
      system(paste("orthofinder -b", tmp.dir,
                   "-t", og.threads,
                   "-a", og.threads,
                   "-S diamond -og"))
    }
  }


  ################   ################   ################
  ################   ################   ################


  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Moving blasts results to",
        blast.dir,
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

  sp.id.files <- file.path(blast.loc,"SpeciesIDs.txt")
  seq.id.files <- file.path(blast.loc,"SequenceIDs.txt")

  files <- c(blast.files,
            fa.files,
            dmnd.files,
            og.files,
            sp.id.files,
            seq.id.files)

  if (dir.exists(blast.dir))
    unlink(blast.dir, recursive = T)
  if (!dir.exists(blast.dir))
    dir.create(blast.dir)

  nu <- file.copy(files,
                  blast.dir)

}
