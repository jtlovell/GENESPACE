#' @title Make sure everything is setup
#'
#' @description
#' \code{check_environment} Utility function to make sure environment is
#' correct and files exist
#'
#' @param directory File path, where the subdirectory folders (e.g. raw_assemblies)
#' are found.
#' @param genomeIDs Character, the vector of genome identifiers to consider
#' for analysis
#' @param check.genomes Logical, should the presence of genomes be checked?
#' @param check.pkgs Logical, should the package installs be checked?
#' @param clean Logical, should all directories and files be deleted, except
#' the raw_annotations and raw_assemblies? Use with caution.
#' @param peptide.only Logical, should only peptides be converted/checked?
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return List of directory names and paths.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
check_env <- function(directory,
                      genomeIDs,
                      clean = FALSE,
                      peptide.only = FALSE,
                      check.pkgs = TRUE,
                      check.genomes = TRUE){

  cat("Making R variables and necessary directories... ")

  dirs <- file.path(directory,
                    c("tmp", "results", "blast", "mcscanx",
                      "cull.blast", "genome", "block",
                      "cull.score", "syn.blast"))
  names(dirs) <- c("tmp", "results", "blast", "mcscanx",
                   "cull.blast", "genome", "block",
                   "cull.score.blast", "syn.blast")

  dir.out <- list()
  for (j in names(dirs)) {

    i <- dirs[[j]]
    if (dir.exists(i) & clean & j != "genome")
      unlink(i, recursive = TRUE)

    if (!dir.exists(i))
      dir.create(i)

    dir.out[[j]] <- i
  }

  genome.dir <- file.path(directory,"genome")

  gff.dir <- file.path(genome.dir, "gff")
  peptide.dir <- file.path(genome.dir, "peptide")
  if(!peptide.only){
    cds.dir <- file.path(genome.dir, "cds")
    assembly.dir <- file.path(genome.dir, "assembly")
    dir.out[["cds"]] <- cds.dir
    dir.out[["assembly"]] <- assembly.dir
  }
  dir.out[["gff"]] <- gff.dir
  dir.out[["peptide"]] <- peptide.dir

  cat("Done!\n")

  if (check.genomes) {
    cat("Checking genomes ... ")
    if (!file.exists(genome.dir))
      stop("Genome directory does not exist, build first\n")

    if(peptide.only){
      ftypes <- file.path(genome.dir,
                          c("peptide",
                            "gff"))
    }else{
      ftypes <- file.path(genome.dir,
                          c("peptide",
                            "cds",
                            "gff",
                            "assembly"))
    }

    for (i in ftypes) {
      if (!file.exists(i)) {
        stop("Must build directory:",i,"first\n")
      }
      sap <- sapply(genomeIDs, function(x) any(grepl(x, dir(i))))
      if (any(!sap)) {
        stop("Directory", i,
             "does not contain files for",
             genomeIDs[!sap], "\n")
      }
    }
    cat("Pass!\n")
  }


  if (check.pkgs) {
    programs <- c("orthofinder")

    cat("Checking for program dependencies ... ")
    fi <- sapply(programs, function(x) Sys.which(x) != "")
    if (all(fi)) {
      cat("Pass!\n")
    } else {
      cat("Fail!\nThe following programs need to be installed and added to the path:\n",
          paste(programs[!fi], collapse = "\n"))
    }
  }

  return(dir.out)
}
