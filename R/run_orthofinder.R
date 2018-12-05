#' @title Run the orthofinder program
#'
#' @description
#' \code{run_orthofinder} A simple wrapper to run orthofinder from R.
#'
#' @param peptide.dir The path to the directory containing the peptide fasta sequence files
#' @param blast.dir The path to the directory where the blast results should be stored
#' @param tmp.dir The path to the directory where temporary files will be stored then deleted
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
  ################   ################   ################
  ################   ################   ################
  if(file.exists(blast.dir)){
    system(paste("rm -r",
                 blast.dir))
  }
  system(paste("mkdir",
               blast.dir))

  if(!is.null(peptide.dir)){
    if(verbose)
      cat("Copying peptide fasta files to",
          tmp.dir,
          "\n")
    if(file.exists(tmp.dir)){
      system(paste("rm -r",
                   tmp.dir))
    }
    system(paste("mkdir",
                 tmp.dir))
    system(paste("cp",
                 file.path(peptide.dir, "*"),
                 tmp.dir))
    if(verbose)
      cat("Running blasts within Orthofinder\n")
    if(og.silent){
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
    system(paste("orthofinder -b", tmp.dir,
                 "-a", og.threads,
                 "-S diamond -og"))
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

  system(paste("cp",
               file.path(blast.loc,
                         "SequenceIDs.txt"),
               blast.dir))
  system(paste("cp",
               file.path(blast.loc,
                         "SpeciesIDs.txt"),
               blast.dir))
  system(paste("cp",
               file.path(blast.loc,
                         "Blast*"),
               blast.dir))
  system(paste("cp",
               file.path(ortho.loc,
                         "Orthogroups.txt"),
               blast.dir))
  system(paste("cp",
               file.path(blast.loc,
                         "Species*.fa"),
               blast.dir))
  system(paste("cp",
               file.path(blast.loc,
                         "diamondDBSpecies*.dmnd"),
               blast.dir))
}
