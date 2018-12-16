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
#' @param clean Logical, should the existing directories be cleaned out?
#' @param check.genomes Logical, should the presence of genomes be checked?
#' @param check.pkgs Logical, should the package installs be checked?
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return List of directory names and paths.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
check_environment <- function(directory,
                              genomeIDs,
                              clean = FALSE,
                              check.pkgs = TRUE,
                              check.genomes = TRUE){

  cat("Making R variables and necessary directories... ")

  dirs = file.path(directory,
                   c("tmp","results","blast","mcscanx",
                     "cull.blast","genome","block"))
  names(dirs)<-c("tmp","results","blast","mcscanx",
                 "cull.blast","genome","block")

  dir.out = list()
  for (j in names(dirs)) {

    i = dirs[[j]]
    if (dir.exists(i) & clean & j != "genome")
      unlink(i, recursive = TRUE)

    if (!dir.exists(i))
      dir.create(i)

    dir.out[[j]] <- i
  }

  genome.dir <- file.path(directory,"genome")

  gff.dir <- file.path(genome.dir, "gff")
  peptide.dir <- file.path(genome.dir, "peptide")
  cds.dir <- file.path(genome.dir, "cds")
  transcript.dir <- file.path(genome.dir, "transcript")
  assembly.dir <- file.path(genome.dir, "assembly")

  dir.out[["gff"]] <- gff.dir
  dir.out[["peptide"]] <- peptide.dir
  dir.out[["cds"]] <- cds.dir
  dir.out[["transcript"]] <- transcript.dir
  dir.out[["assembly"]] <- assembly.dir
  dir.out[["gff"]] <- gff.dir

  cat("Done!\n")

  if (check.genomes) {
    cat("Checking genomes ... ")
    if(!file.exists(genome.dir))
      stop("Genome directory does not exist, build first\n")

    ftypes <- file.path(genome.dir,
                        c("transcript",
                          "peptide",
                          "cds",
                          "gff",
                          "assembly"))

    for(i in ftypes){
      if(!file.exists(i)){
        stop("Must build directory:",i,"first\n")
      }
      sap = sapply(genomeIDs, function(x) any(grepl(x, dir(i))))
      if(any(!sap)){
        stop("Directory", i,
             "does not contain files for",
             genomeIDs[!sap], "\n")
      }
    }
    cat("Pass!\n")
  }


  if (check.pkgs) {
    programs <- c("bedtools",
                  "MCScanX",
                  "Diamond",
                  "orthofinder",
                  "exonerate",
                  "samtools")

    packages <- c("data.table",
                  "Biostrings",
                  "dbscan")

    cat("Checking for R Package dependencies ... ")
    suppressPackageStartupMessages(
      fi <- sapply(packages, require, quietly = T, character.only = T))
    if(all(fi)){
      cat("Pass!\n")
    }else{
      cat("Fail!\nThe following packages need to be installed:\n",
          paste(packages[!fi], collapse = "\n"))
    }

    cat("Checking for program dependencies ... ")
    fi <- sapply(programs, function(x) Sys.which(x) != "")
    if(all(fi)){
      cat("Pass!\n")
    }else{
      cat("Fail!\nThe following programs need to be installed and added to the path:\n",
          paste(programs[!fi], collapse = "\n"))
    }
  }

  return(dir.out)
}
