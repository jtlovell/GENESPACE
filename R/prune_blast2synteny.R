#' @title prune_blast2synteny
#'
#' @description
#' \code{prune_blast2synteny} prune_blast2synteny
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
prune_blast2synteny <- function(blast,
                                gff,
                                genomeIDs,
                                dir.list,
                                MCScanX.path,
                                MCScanX.s.param = 10,
                                MCScanX.m.param = 100,
                                dbs.radius = 100,
                                dbs.hits = 20,
                                overwrite.output.dir,
                                verbose = T){
  if(verbose)
    cat("Pruning inter-genomic blast hits with MCScanX ... ")
  mcs1 <- run_MCScanX(
    blast = blast,
    gff = gff,
    mcscan.dir = dir.list$mcscanx,
    overwrite.output.dir = overwrite.output.dir,
    genomeIDs = genomeIDs,
    MCScanX.path = MCScanX.path,
    MCScanX.s.param = MCScanX.s.param,
    MCScanX.m.param = MCScanX.m.param,
    verbose = F)

  if(verbose)
    cat("Done!\nCompleting unconnected inter-genomic sub-graphs ... ")
  comp1 <- complete_graph(
    map = mcs1,
    gff = gff,
    ignore.self = T,
    verbose = F)

  comp1 <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = comp1,
    intergenome.only = T)

  if(verbose)
    cat("Done!\nInitial fixed-radius pruning of inter-genomic hits ... ")
  cln1 <- clean_it(
    map = comp1,
    genomeIDs = genomeIDs,
    rerank = T,
    radius = dbs.radius,
    n.mappings = dbs.hits,
    verbose = F)$map

  if(verbose)
    cat("Done!\nCompleting unconnected sub-graphs ... ")
  comp.final <- complete_graph(
    map = cln1,
    gff = gff,
    ignore.self = T,
    verbose = F)

  comp.final <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = comp.final,
    intergenome.only = F)

  if(verbose)
    cat("Done!\nFinal fixed-radius pruning via dbscan ... \n")
  clean.final <- clean_it(
    map = comp.final,
    genomeIDs = genomeIDs,
    rerank = T,
    radius = dbs.radius,
    n.mappings = dbs.hits,
    verbose = verbose)$map

  if(verbose)
    cat("\tDone!\n")
  return(mirror_map(clean.final))
}






