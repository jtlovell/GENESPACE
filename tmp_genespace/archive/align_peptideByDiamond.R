#' @title Conduct reciprocal protein alignments between a pair of genomes
#'
#' @description
#' \code{align_peptideByDiamond} Uses the Diamond Blast program to generate
#' databases, then reciprocally align.
#'
#' @param path_to_diamond The location of the Diamond program executable
#' @param id1 ID for genome 1
#' @param id2 ID for genome 2
#' @param pep1 Peptide file for genome 1
#' @param pep2 Peptide file for genome 2
#' @param blast1 blast output filename for genome 1
#' @param blast2 blast output filename for genome 2
#' @param sensitive.mode Diamond sensitivity mode.
#' @param topPerc Genes this percent from the maximum hit are retained
#' @param minScore The minimum mapping score to be considered
#' @param ortherDiamondOpts Other diamond options to pass the command
#' @param nthreads Number of parallel threads to run blast with
#' @param verbose Logical, should status updates be printed?
#' @param ... Not currently in use
#' @details See pipe_Diamond2MCScanX for more information.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
align_peptideByDiamond<- function(path_to_diamond = NULL,
                                  pep1, pep2,
                                  id1,id2,
                                  blast1, blast2,
                                  sensitive.mode = "--more-sensitive",
                                  topPerc = 50,
                                  minScore = 50,
                                  ortherDiamondOpts = "--quiet",
                                  nthreads= 1, verbose = T){
  db1 = gsub(".fa$","",pep1)
  db2 = gsub(".fa$","",pep2)
  if(is.null(pep1))
    stop("must provide at least one peptide fasta path\n")
  if(is.null(path_to_diamond))
    cat("using diamond in path\n")

  if(verbose)
    cat("1. making diamond commands ...")

  com1 = paste(file.path(path_to_diamond,"diamond"), "makedb --in", pep1,
               "-d", db1, "--quiet",
               "-p", nthreads)
  com2 = paste(file.path(path_to_diamond,"diamond"), "makedb --in", pep2,
               "-d", db2, "--quiet",
               "-p", nthreads)

  com3 = paste(file.path(path_to_diamond,"diamond"), "blastp",
               "--query",pep2,
               "--db",db1,
               "--top",topPerc,
               "--min-score",minScore,
               ifelse(is.na(sensitive.mode),"",sensitive.mode),
               ifelse(is.na(ortherDiamondOpts),"",ortherDiamondOpts),
               "-p",nthreads,
               "--out", blast1)
  com4 = paste(file.path(path_to_diamond,"diamond"), "blastp",
               "--query",pep1,
               "--db",db2,
               "--top",topPerc,
               "--min-score",minScore,
               ifelse(is.na(sensitive.mode),"",sensitive.mode),
               ifelse(is.na(ortherDiamondOpts),"",ortherDiamondOpts),
               "-p",nthreads,
               "--out", blast2)

  if(verbose)
    cat("\n2. making databases\t",id1,"... ")
  system(com1)

  if(verbose)
    cat("Done!", id2,"... ")
  system(com2)
  if(verbose)
    cat("Done!\n3. Running blasts\t",id1,"... ")
  system(com3)
  if(verbose)
    cat("Done!", id2,"... ")
  system(com4)
  if(verbose)
    cat("Done!\n### Diamond Blast Complete! ###")
  system(paste0("rm ", file.path(peptide.dir,"*.dmnd")))
  return(c(com1,com2,com3,com4))
}
