#' @title Import orthofinder results
#'
#' @description
#' \code{import_ofResults} Read in orthogroups, species mappings and gene IDs
#' produced from the orthofinder algorithm
#'
#' @param gff gff annotation data.table
#' @param blast.dir path to orthofinder output directory
#' @param genomeIDs genome identifiers
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return culled blast data.table
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
import_ofResults <- function(gff,
                             blast.dir,
                             verbose = T,
                             genomeIDs,
                             ...){
  gff <- data.table(gff)
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing orthofinder results:\n\t")
  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("Decompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(blast.dir,
                           "*.gz")))
  }
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading Species IDs\n\t")
  si <- read.delim(file.path(blast.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  rownames(si) <- si$genome
  #######################################################

  #######################################################
  sm <- cbind(expand.grid(si$genome.num,
                          si$genome.num,
                          stringsAsFactors = F),
              expand.grid(si$genome,
                          si$genome,
                          stringsAsFactors = F))
  colnames(sm) <- c("n1","n2","genome1","genome2")
  for (i in rev(genomeIDs))
    sm$ref[sm$genome1 == i | sm$genome2 == i] <- i
  sm$alt <- ifelse(sm$genome1 == sm$ref, sm$genome2, sm$genome1)
  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = genomeIDs))
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading gene IDs\n\t")
  sequence.index <- fread(file.path(blast.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading orthogroup networks\n\t")
  og <- readLines(file.path(blast.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])

  if (verbose)
    cat("Building orthogroup data.table\n\t")
  og2 <- readLines(file.path(blast.dir,
                             "Orthogroups.txt"))
  og2 <- lapply(og2, function(x) strsplit(x, " ")[[1]])
  og.name <- sapply(og2, function(x) x[1])
  og.length <- sapply(og2, length) - 1
  og.ids <- sapply(og2, function(x) x[-1])
  og2 <- data.table(block.id = NA,
                    og = rep(og.name, og.length),
                    id = unlist(og.ids),
                    stringsAsFactors = F)
  #######################################################

  #######################################################
  if (verbose)
    cat("Compiling metadata\n")

  gffi <- gff[, c("id","genome")]
  setkey(gffi, "id")
  setkey(og2, "id")
  ogo <- merge(gffi, og2)
  setkey(ogo, "block.id", "genome", "og", "id")
  ogo[, og.n.genes := length(unique(id)), by = list(block.id, og)]
  ogo[, og.n.genomes := length(unique(genome)), by = list(block.id, og)]

  ogo.meta = ogo[!duplicated(ogo[, -c(1:2), with = F]),-c(1:2), with = F]
  #######################################################

  #######################################################
  if (verbose)
    cat("\tDone!\n")
  return(list(orthogroups = og,
              species.mappings = sm,
              species.index = si,
              orthogroup.metadata = ogo.meta,
              orthogroup.data.table = ogo,
              gene.index = sequence.index))
}
