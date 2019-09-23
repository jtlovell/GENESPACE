#' @title pipe_mcscanx
#'
#' @description
#' \code{pipe_mcscanx} pipe_mcscanx
#'
#' @param blast data.table, containing the merged gff and blast results.
#' Unlike the 'map' object, which really just needs merged gff coordinates,
#' this must have all the blast8 columns. See details.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param MCScanX.s.param numeric of length 1, that specifies the 's'
#' (block size) parameter for MCScanX.
#' @param MCScanX.m.param numeric of length 1, that specifies the 'm'
#' (n. gaps) parameter for MCScanX.
#' @param MCScanX.path file.path, specifying the location of the
#' MCScanX program. This directory must contain the executable
#' '/MCScanX'.
#' @param silent.mcs logical, should MCScanX progress be reported?
#'
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
pipe_mcscanx <- function(blast,
                         genomeIDs,
                         dir.list,
                         gff,
                         MCScanX.s.param,
                         MCScanX.m.param,
                         MCScanX.path,
                         silent.mcs = T,
                         verbose = T){

  blast <- subset(blast,
                genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  blast[, genome1 := factor(genome1, levels = genomeIDs)]
  blast[, genome2 := factor(genome2, levels = genomeIDs)]
  setkey(blast, genome1, genome2)
  blast[, genome1 := as.character(genome1)]
  blast[, genome2 := as.character(genome2)]

  #######################################################
  if (verbose)
    cat("Preparing the gff and blast objects ... ")
  MCScanX.tool <- file.path(MCScanX.path, "MCScanX")
  mcscan.dir <- dir.list$mcscanx
  mcscan.param <- paste("-a -s", MCScanX.s.param,
                        "-m", MCScanX.m.param)
  #######################################################

  gff.tmp <- data.table(gff)
  gff.tmp$genome <- paste0(gff.tmp$genome,"xxxx")
  gff.tmp$id <- paste0(gff.tmp$id,"xxxx")
  gff <- data.table(rbind(gff, gff.tmp))

  genomeIDs <- c(genomeIDs, paste0(genomeIDs, "xxxx"))

  bl.dif <- subset(blast, genome1 != genome2)
  bl.same <- subset(blast, genome1 == genome2)
  bl.same[, id2 := paste0(id2, "xxxx")]
  bl.same[, genome2 := paste0(genome2, "xxxx")]
  blast <- rbind(bl.dif, bl.same)
  #######################################################
  blast <- blast[!duplicated(blast[,c("id1","id2")]),]
  blast[, unique := paste(genome1, genome2)]

  if ("block.id" %in% colnames(blast))
    blast[, block.id := NULL]

  if (verbose)
    cat("Done!\nCulling pairwise blast hits to collinear blocks ... \n")
  spl <- split(blast, by = "unique")
  #######################################################
  out <- rbindlist(lapply(spl, function(x){
    genomes <- c(x$genome1[1], x$genome2[1])

    if (verbose)
      cat(paste0("\t", gsub("xxxx", "", genomes[1])),
          "-->", gsub("xxxx", "", genomes[2]),
          paste0("(initial hits = ", nrow(x), ")"))

    gff.x <- subset(gff, genome %in% genomes)

    tmp <- run_mcs(blast = x,
                   gff = gff.x,
                   genomeIDs = genomeIDs,
                   MCScanX.tool = MCScanX.tool,
                   mcscan.dir = mcscan.dir,
                   mcscan.param = mcscan.param,
                   silent.mcs = silent.mcs)
    if (nrow(tmp) == 0) {
      if (verbose)
        cat(" culled hits = 0!\n")
      return(tmp)
    }else{

      tmp[, block.id := paste0(unique, block.id)]

      if (verbose)
        cat(" culled hits =", nrow(tmp), "\n")

      return(tmp)
    }
  }), fill = T)
  #######################################################
  out[, block.id := as.numeric(as.factor(block.id))]
  out[, genome2 := gsub("xxxx", "", genome2)]
  out[, id2 := gsub("xxxx", "",  out$id2)]
  if (verbose)
    cat("\tDone!\nCompiling block coordinates ... ")
  alo <- make_blocks(out,
                     clean.columns = F,
                     rename.blocks = T)
  if (verbose)
    cat("Done!\n")
  #######################################################
  return(alo)
}
