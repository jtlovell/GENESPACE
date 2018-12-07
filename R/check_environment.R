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
                              clean = F){

  cat("Making R variables and necessary directories... ")

  tmp.dir <- file.path(directory, "tmp")
  if(file.exists(tmp.dir))
    system(paste("rm -r", tmp.dir))
  if (clean)
    system(paste("mkdir", tmp.dir))

  results.dir <- file.path(directory, "results")
  if(file.exists(results.dir) & clean)
    system(paste("rm -r", results.dir))
  if (clean)
    system(paste("mkdir",results.dir))

  blast.dir <- file.path(directory, "blast")
  if(file.exists(blast.dir) & clean)
    system(paste("rm -r", blast.dir))
  if (clean)
    system(paste("mkdir", blast.dir))

  mcscan.dir <- file.path(directory, "mcscanx")
  if(file.exists(mcscan.dir) & clean)
    system(paste("rm -r", mcscan.dir))
  if (clean)
    system(paste("mkdir", mcscan.dir))

  cull.blast.dir <- file.path(directory, "cull.blast")
  if(file.exists(cull.blast.dir) & clean)
    system(paste("rm -r", cull.blast.dir))
  if (clean)
    system(paste("mkdir", cull.blast.dir))

  genome.dir <- file.path(directory,"genome")

  gff.dir <- file.path(genome.dir, "gff")
  peptide.dir <- file.path(genome.dir, "peptide")
  cds.dir <- file.path(genome.dir, "cds")
  transcript.dir <- file.path(genome.dir, "transcript")
  assembly.dir <- file.path(genome.dir, "assembly")

  dirs = list(gff = gff.dir,
              peptide = peptide.dir,
              cds = cds.dir,
              transcript = transcript.dir,
              assembly = assembly.dir,
              mcscan = mcscan.dir,
              genome = genome.dir,
              blast = blast.dir,
              results = results.dir,
              tmp = tmp.dir,
              cull.blast = cull.blast.dir)

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

  return(dirs)
}
