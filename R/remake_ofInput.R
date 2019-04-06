#' @title remake_ofInput
#'
#' @description
#' \code{remake_ofInput} remake_ofInput
#'
#' @param dir.list dir.list
#' @param genomeIDs genomeIDs
#' @param ploidy named numeric vector corresponding to the genomeIDs
#' @param cull.blastByScore should blasts be culled to the top hits?
#' @param max.dup maximum ploidy / 2 multiplier
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom utils combn
#' @export
remake_ofInput <- function(dir.list,
                           genomeIDs,
                           ploidy,
                           cull.blastByScore = T,
                           max.dup = 2,
                           verbose = T){
  # -- Step 1. Reformat peptides etc.
  if (dir.exists(dir.list$cull.score.blast))
    unlink(dir.list$cull.score.blast,
           recursive = T)
  dir.create(dir.list$cull.score.blast)

  if (verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")
  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = dir.list$cull.score.blast,
               peptide.dir = dir.list$peptide,
               genomeIDs = genomeIDs)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing gff annotations as a data.table ... \n")
  gff <- import_gff(gff.dir = dir.list$gff,
                    genomeIDs = genomeIDs)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = dir.list$blast,
                             genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = dir.list$cull.score.blast,
                             genomeIDs = genomeIDs)
  id.db <- merge(old.ids, new.ids,
                 by = "genome")
  setnames(id.db, 2:3, c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       blast.dir = dir.list$blast,
                       cull.blast.dir = dir.list$cull.score.blast)
  #
  old.genes <- read_geneIDs(of.dir = dir.list$blast, gff = gff)
  new.genes <- read_geneIDs(of.dir = dir.list$cull.score.blast, gff = gff)
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
  if(verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")

  d <- lapply(1:nrow(map.db), function(i){
    if (verbose)
      cat(paste0("\t", map.db$genome1[i]),
          "-->", map.db$genome2[i], "... ")
    bl <- readRename_blastGenes(gene.dict1 = g1,
                                gene.dict2 = g2,
                                blast.file.in = map.db$filename[i],
                                blast.file.out = map.db$new.filename[i])
    if (verbose)
      cat("Done!\n")
    return(bl)
  })
  #######################################################

  #######################################################
  if (verbose)
    cat("Culling blast by score ... \n")

  map.db$score.filename <- file.path(dir.list$cull.score.blast,
                                     basename(map.db$new.filename))
  bl <- lapply(1:nrow(map.db), function(i){
    s1 <- map.db$genome1[i]
    s2 <- map.db$genome2[i]
    maxn <- (ploidy[s1] * max.dup) / 2
    if(verbose)
      cat(paste0("\t", s1), "-->", s2, "... ")
    blast.file.in <- map.db$new.filename[i]
    blast.file.out <- map.db$score.filename[i]
    cull_blastByScore(blast.file.in = blast.file.in,
                      blast.file.out = blast.file.out,
                      maxn = maxn,
                      verbose = T)
    if (verbose)
      cat("Done!\n")
  })
  if (verbose)
    cat("\tDone!\n")
  return(map.db)
}
