#' @title Make sure everything is setup
#'
#' @description
#' \code{check_environment} Utility function to make sure environment is
#' correct and files exist
#'
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
check_environment = function(){

  obj = c("input.dir", "results.dir", "blast.dir", "results.dir", "genomeIDs",  "ploidy",  "abbrevs")
  fi = sapply(obj, exists)

  cat("Checking R objects ... ")
  if(all(fi)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nThe following required objects do not exist\n",paste(obj[!fi], collapse = "\n"))
  }

  programs = c("bedtools","MCScanX","Diamond","orthofinder","exonerate")
  packages = c("data.table", "Biostrings", "dbscan")
  directories = file.path(input.dir,c("peptide","cds","assembly","gff"))

  input.files = c(file.path(directories[1], paste0(genomeIDs,".fa")),
                  file.path(directories[2], paste0(genomeIDs,".fa")),
                  file.path(directories[3], paste0(genomeIDs,".fa")),
                  file.path(directories[4], paste0(genomeIDs,".gff3")))

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

  cat("Checking input.files ... ")
  fi = file.exists(input.files)
  if(all(fi)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nThe following files were expected, but do not exist:\n", paste(input.files[!fi], collapse = "\n"))
  }

  cat("Checking for R Package dependencies ... ")
  suppressPackageStartupMessages(fi <- sapply(packages, require ,quietly = T, character.only = T))
  if(all(fi)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nThe following packages need to be installed:\n", paste(packages[!fi], collapse = "\n"))
  }

  cat("Checking for program dependencies ... ")
  fi = sapply(programs, function(x) Sys.which(x) !="")
  if(all(fi)){
    cat("Pass!\n")
  }else{
    cat("Fail!\nThe following programs need to be installed and added to the path:\n", paste(programs[!fi], collapse = "\n"))
  }
  names(ploidy)<-genomeIDs
  names(abbrevs)<-genomeIDs
  ploidy<<-ploidy
  abbrevs<<-abbrevs
}
