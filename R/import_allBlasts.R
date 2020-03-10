#' @title Read all blast files into memory
#'
#' @description
#' \code{import_allBlasts} Read all blast files in an
#' orthofinder output directory.
#'
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param orthofinder.dir file.path, to the subdirectory containing
#' blast results
#' @param add.orthogroups logical, should the orthogroups be considered at all?
#' @param add.gff logical, should annotation data be added?
#' @param keep.geneNum logical, should the orthofinder gene ID be retained?
#' @param blast.files vector of file paths, optional, if specific
#' files are desired.
#' @param ... Not currently in use
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' all.blast <- read_allBlasts(
#' gff = gff,
#' genomeIDs = genomeIDs,
#' orthofinder.dir = dir.locs$cull.score.blast)
#' }
#' @import data.table
#' @export
import_allBlasts <- function(gff,
                           genomeIDs,
                           orthofinder.dir,
                           add.orthogroups = F,
                           keep.geneNum = F,
                           blast.files = NULL){

  # Functions call
  read_ofBlast <- function(gids,
                           blast.files){
    all.blast <- rbindlist(lapply(blast.files, function(x)
      fread(x,
            col.names = c("gn1", "gn2",
                          "perc.iden", "align.length",
                          "n.mismatch", "n.gapOpen",
                          "q.start", "q.end",
                          "s.start", "s.end",
                          "eval", "score"),
            key = "gn2")))
    if ("gene.num" %in% colnames(gids))
      setnames(gids, "gene.num","gn")
    gids1 <- data.table(gids)
    setnames(gids1, paste0(colnames(gids),"1"))
    gids2 <- data.table(gids)
    setnames(gids2, paste0(colnames(gids),"2"))
    all.blast <- merge(gids1,
                       merge(gids2,
                             all.blast,
                             by = "gn2"),
                       by = "gn1")
    return(all.blast)
  }

  of.dir <- orthofinder.dir

  # -- check the params
  ################################################
  ################################################
  ################################################
  if (!is.data.table(gff))
    stop("gff must be a data table containing annotation information\n")

  if (!is.character(genomeIDs) | length(genomeIDs) == 1)
    stop("genomeIDs must be a character vector of length > 1\n")

  if (!dir.exists(orthofinder.dir))
    stop("can't find orthofinder.dir ... \n")

  if (!is.logical(add.orthogroups) |
     !is.logical(keep.geneNum))
    stop("both keep.geneNum and add.orthogroups need to be logical argumnents\n")

  if (!"OrthoFinder" %in% list.files(of.dir) & add.orthogroups) {
    warning("It doesn't look like the orthofinder program has been run in ",
            of.dir,"\n\tSetting add.orthogroups to FALSE")
    add.orthogroups <- FALSE
  }

  if (!is.null(blast.files))
    if (any(!file.exists(blast.files)))
      stop("Cant find some blast files \n")

  spids <- gsub(" |.fa$","", unlist(
    fread(file.path(of.dir, "SpeciesIDs.txt"),
          select = 2, header = F, sep = ":")))
  if (!all(genomeIDs %in% spids))
    stop("Can't find all genomeIDs in",
         file.path(of.dir,
                   "SpeciesIDs.txt"),"\n")
  ################################################
  ################################################
  ################################################

  if ("gn" %in% colnames(gff))
    gff[,gn := NULL]

  ################################################
  # -- read in the genome IDs
  sids <- read_speciesIDs(
    of.dir = of.dir,
    genomeIDs = genomeIDs)

  ################################################
  # -- read in the gene IDs
  gids <- read_geneIDs(
    of.dir = of.dir,
    species.num.id = sids,
    gff = gff)

  ################################################
  # -- optionally find and read in all blast files
  if (is.null(blast.files))
    blast.files = list.files(
      path = of.dir,
      pattern = "^Blast",
      full.names = T)

  all.blast <- read_ofBlast(
    gids = gids,
    blast.files = list.files(
      path = of.dir,
      pattern = "^Blast",
      full.names = T))

  if (!keep.geneNum) {
    if ("gene.num1" %in% colnames(all.blast))
      all.blast[, gene.num1 := NULL]
    if ("gene.num2" %in% colnames(all.blast))
      all.blast[, gene.num2 := NULL]
    if ("genome.num1" %in% colnames(all.blast))
      all.blast[, genome.num1 := NULL]
    if ("genome.num2" %in% colnames(all.blast))
      all.blast[, genome.num2 := NULL]
  }


  ################################################
  # -- optionally add in orthogroup IDs
  if (add.orthogroups) {
    orthogroup.dt <- read_ogs(
      of.dir = of.dir,
      gff = gff)
    all.blast <- merge(
      all.blast,
      orthogroup.dt,
      by = c("genome1", "genome2", "id1", "id2"))
  }else{
    all.blast[,og := NA]
  }

  return(all.blast)
}
