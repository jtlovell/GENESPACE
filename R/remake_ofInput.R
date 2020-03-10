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
                           verbose = TRUE,
                           n.cores = 1,
                           overwrite.output.dir = FALSE,
                           rbh.only = F,
                           ...){

  peptide.dir <- dir.list$peptide
  blast.dir <- dir.list$blast
  cull.blast.dir <- dir.list$cull.blast
  cull.score.blast.dir <- dir.list$cull.score.blast

  tmp.dir <- file.path(getwd(), "gs.remake.tmp")
  if (dir.exists(tmp.dir))
    unlink(tmp.dir,
           recursive = T)
  dir.create(tmp.dir)
  on.exit(
    expr = unlink(tmp.dir,
                  recursive = T))

  if (is.null(ploidy))
    ploidy <- rep(2, length(genomeIDs))
  if (length(ploidy) == 1)
    ploidy <- rep(ploidy, length(genomeIDs))

  if (length(min.score) > 1)
    min.score <- min.score[1]

  ################################################
  ################################################
  ################################################
  if (!is.data.table(gff))
    stop("gff must be a data table containing annotation information\n")

  if (!is.character(genomeIDs) |
      length(genomeIDs) == 1)
    stop("genomeIDs must be a character vector of length > 1\n")

  if (!is.list(dir.list))
    stop("dir.list must be a list of directories\n")

  if (!all(dir.exists(c(
    blast.dir,
    cull.blast.dir,
    cull.score.blast.dir,
    tmp.dir,
    peptide.dir))))
    stop("it looks like there is something wrong with directories,
         have you called check_env?")

  if (!all(is.numeric(max.dup),
           length(max.dup) == 1,
           max.dup > 0))
    stop("max.dup must be a numeric vector > 0 of length 1\n")


  if (!all(is.numeric(ploidy),
           length(ploidy) == length(genomeIDs)))
    stop("ploidy must be a numeric vector with the same length as genomeIDs\n")

  if (min(ploidy) < 1) {
    warning("some ploidies are less than haploid, these will be set to haploid (ploidy = 1)")
    ploidy[ploidy < 1 | is.na(ploidy)] <- 1
  }

  if (!is.logical(verbose))
    verbose <- FALSE

  ################################################
  ################################################
  ################################################



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
  map.db$score.filename <- file.path(
    cull.score.blast.dir,
    basename(map.db$new.filename))

  bl <- lapply(1:nrow(map.db), function(i){
    s1 <- map.db$genome1[i]
    s2 <- map.db$genome2[i]
    maxn <- (ploidy[s2]*max.dup)/2
    if (verbose)
      cat(paste0("\t",s1),"-->",s2,"... ")
    blast.file.in = map.db$new.filename[i]
    blast.file.out = map.db$score.filename[i]
    cull_blastByScore(
      blast.file.in = blast.file.in,
      blast.file.out = blast.file.out,
      maxn = maxn,
      rbh.only = rbh.only,
      verbose = T)
    if (verbose)
      cat("Done!\n")
  })

  if (verbose)
    cat("\tDone!\n")
}
