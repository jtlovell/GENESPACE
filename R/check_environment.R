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
#' @param abbrevs Character vector of length equal to genomeIDs. Must be unique
#' 2-character words to identify each genome
#' @param ploidy The ploidy of each genome
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
check_environment <- function(directory,
                              genomeIDs,
                              abbrevs,
                              ploidy,
                              clean = T){

  cat("Making R variables and necessary directories... ")

  tmp.dir <<- file.path(directory, "tmp")
  if(file.exists(tmp.dir))
    system(paste("rm -r", tmp.dir))
  system(paste("mkdir", tmp.dir))

  results.dir <<- file.path(directory, "results")
  if(file.exists(results.dir) & clean)
    system(paste("rm -r", results.dir))
  system(paste("mkdir",results.dir))

  blast.dir <<- file.path(directory, "blast")
  if(file.exists(blast.dir) & clean)
    system(paste("rm -r", blast.dir))
  system(paste("mkdir", blast.dir))

  blast.dir <<- file.path(directory, "block")
  if(file.exists(blast.dir) & clean)
    system(paste("rm -r", blast.dir))
  system(paste("mkdir", blast.dir))

  mcscan.dir <<- file.path(directory, "mcscanx")
  if(file.exists(mcscan.dir) & clean)
    system(paste("rm -r", mcscan.dir))
  system(paste("mkdir", mcscan.dir))

  genome.dir <<- file.path(directory,"genome")

  gff.dir <<- file.path(genome.dir, "gff")
  peptide.dir <<- file.path(genome.dir, "peptide")
  cds.dir <<- file.path(genome.dir, "cds")
  transcript.dir <<- file.path(genome.dir, "transcript")
  assembly.dir <<- file.path(genome.dir, "assembly")


  abbrevs <<- abbrevs
  ploidy <<- ploidy

  cat("Done!\n")

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

  programs <- c("bedtools",
                "MCScanX",
                "Diamond",
                "orthofinder",
                "exonerate",
                "samtools")

  packages <- c("data.table",
                "Biostrings",
                "dbscan")


  cat("Checking ploidies ... ")
  if(length(ploidy) == length(genomeIDs)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nPloidy specification is not of the same length as genome IDs\n")
  }

  cat("Checking abbreviations ... ")
  if(length(abbrevs) == length(genomeIDs)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nAbbreviation specification is not of the same length as genome IDs\n")
  }

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
  names(ploidy) <- genomeIDs
  names(abbrevs) <- genomeIDs
  ploidy <<- ploidy
  abbrevs <<- abbrevs
}
