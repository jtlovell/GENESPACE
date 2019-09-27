#' @title re-make orthofinder input
#'
#' @description
#' \code{remake_ofInput} re-make orthofinder input, culling
#' by top hits for each genome and optionally by score.
#'
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param verbose logical, should updates be printed to the console?
#' @param ploidy named integer vector, of length equal to genomeIDs, and named
#' with each genome ID.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param max.dup integer, n top hits to retain for each gene.
#' @param n.cores integer length 1, the number of parallel processes
#' @param min.score numeric, the minimum blast bit score to retain
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param ... Additional arguments passed on to import_gff
#' @details ...
#' @return The function does not return anything to the R console.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
remake_ofInput <- function(dir.list,
                           genomeIDs,
                           ploidy = NULL,
                           gff,
                           max.dup = 2,
                           min.score = 50,
                           verbose = T,
                           n.cores = 1,
                           overwrite.output.dir = FALSE,
                           ...){

  peptide.dir <- dir.list$peptide
  blast.dir <- dir.list$blast
  cull.blast.dir <- dir.list$cull.blast
  cull.score.blast.dir <- dir.list$cull.score.blast
  tmp.dir <- file.path(getwd(),"gs.remake.tmp")

  # -- Make tmp directory
  if (dir.exists(tmp.dir))
    unlink(tmp.dir,
           recursive = T)
  dir.create(tmp.dir)
  on.exit(
    expr = unlink(tmp.dir,
                  recursive = T))

  # -- Check the output directories
  check_gsDir(dir2check = cull.blast.dir,
              overwrite.output.dir = overwrite.output.dir)
  check_gsDir(dir2check = cull.score.blast.dir,
              overwrite.output.dir = overwrite.output.dir)

  # -- Make the new databaase
  map.db <- prep_ofDB(
    tmp.dir = tmp.dir,
    orig.dir = blast.dir,
    peptide.dir = peptide.dir,
    output.dir = cull.blast.dir,
    n.cores = n.cores,
    genomeIDs = genomeIDs,
    min.score = min.score,
    copy2dir = cull.score.blast.dir,
    verbose = verbose)
  if (verbose)
    cat("\tDone!\n")

  #######################################################

  #######################################################
  if (verbose)
    cat("Culling blast by score ... \n")
  map.db$score.filename <- file.path(cull.score.blast.dir,
                                     basename(map.db$new.filename))
  bl <- lapply(1:nrow(map.db), function(i){
    s1 <- map.db$genome1[i]
    s2 <- map.db$genome2[i]
    maxn <- (ploidy[s2]*max.dup)/2
    if (verbose)
      cat(paste0("\t",s1),"-->",s2,"... ")
    blast.file.in = map.db$new.filename[i]
    blast.file.out = map.db$score.filename[i]
    cull_blastByScore(blast.file.in = blast.file.in,
                      blast.file.out = blast.file.out,
                      maxn = maxn,
                      verbose = T)
    if (verbose)
      cat("Done!\n")
  })

  if (verbose)
    cat("\tDone!\n")
  return(list(files = map.db, gff = gff))
}
