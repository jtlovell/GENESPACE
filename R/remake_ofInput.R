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
                           gff = NULL,
                           max.dup = 2,
                           min.score = 50,
                           verbose = T,
                           n.cores = 1,
                           ...){

  if (dir.exists(dir.list$cull.blast))
    unlink(dir.list$cull.blast, recursive = T)
  dir.create(dir.list$cull.blast)

  if (dir.exists(dir.list$cull.score.blast))
    unlink(dir.list$cull.score.blast, recursive = T)
  dir.create(dir.list$cull.score.blast)

  if (dir.exists(dir.list$tmp))
    unlink(dir.list$tmp, recursive = T)
  dir.create(dir.list$tmp)

  if (verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")
  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = dir.list$cull.blast,
               peptide.dir = dir.list$peptide,
               n.cores = n.cores,
               genomeIDs = genomeIDs)
  fs <- list.files(dir.list$cull.blast, full.names = T)
  cpd <- file.copy(fs, dir.list$cull.score.blast)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (is.null(gff)) {
    if (verbose)
      cat("Importing gff annotations as a data.table ... \n")
    gff <- import_gff(gff.dir = dir.list$gff,
                      genomeIDs = genomeIDs,
                      ...)
    if (verbose)
      cat("\tDone!\n")
  }

  #######################################################

  #######################################################
  if (verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = dir.list$blast,
                             genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = dir.list$cull.blast,
                             genomeIDs = genomeIDs)
  id.db <- merge(old.ids, new.ids, by = "genome")
  setnames(id.db,2:3,c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       of.dir = dir.list$blast,
                       cull.blast.dir = dir.list$cull.blast)
  #
  old.genes <- read_geneIDs(of.dir = dir.list$blast,
                            gff = gff)
  new.genes <- read_geneIDs(of.dir = dir.list$cull.blast,
                            gff = gff)
  genes <- merge(old.genes[,c("genome","id","gene.num")],
                 new.genes[,c("genome","id","gene.num")],
                 by = c("genome","id"))
  g1 <- with(genes,
             data.table(V1 = gene.num.x,
                        new1 = gene.num.y,
                        key = "V1"))
  g2 <- with(genes,
             data.table(V2 = gene.num.x,
                        new2 = gene.num.y,
                        key = "V2"))
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")

  d <- lapply(1:nrow(map.db), function(i){
    if (verbose)
      cat(paste0("\t",map.db$genome1[i]),"-->",map.db$genome2[i],"... ")
    bl <- readRename_blastGenes(gene.dict1 = g1,
                                gene.dict2 = g2,
                                min.score = min.score,
                                blast.file.in = map.db$filename[i],
                                blast.file.out = map.db$new.filename[i])
    if (verbose)
      cat("Done!\n")
    return(bl)
  })
  if (verbose)
    cat("\tDone!\n")

  #######################################################

  #######################################################
  if (verbose)
    cat("Culling blast by score ... \n")
  map.db$score.filename <- file.path(dir.list$cull.score.blast,
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
