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
  sblast <- subset(blast, genome1 == genome2 & id1 == id2)
  mcs1 <- rbind(mcs1,sblast, fill = T)
  mcs1 <- mcs1[!duplicated(mcs1),]
  if(verbose)
    cat("Done!\nCompleting unconnected inter-genomic sub-graphs ... ")
  comp1 <- complete_graph(
    map = mcs1,
    gff = gff,
    ignore.self = T,
    verbose = F)
  comp1 <- rbind(comp1,sblast, fill = T)
  comp1 <- comp1[!duplicated(comp1),]

  comp1 <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = comp1,
    intergenome.only = T)
  comp1 <- rbind(comp1,sblast, fill = T)
  comp1 <- comp1[!duplicated(comp1),]

  if(verbose)
    cat("Done!\nInitial fixed-radius pruning of inter-genomic hits ... ")
  cln1 <- clean_it(
    map = comp1,
    genomeIDs = genomeIDs,
    rerank = T,
    radius = dbs.radius,
    n.mappings = dbs.hits,
    verbose = F)$map
  cln1 <- rbind(cln1, sblast, fill = T)
  cln1 <- cln1[!duplicated(cln1),]

  if(verbose)
    cat("Done!\nCompleting unconnected sub-graphs ... ")
  comp.final <- complete_graph(
    map = cln1,
    gff = gff,
    ignore.self = T,
    verbose = F)
  comp.final <- rbind(comp.final,sblast, fill = T)
  comp.final <- comp.final[!duplicated(comp.final),]

  comp.final <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = comp.final,
    intergenome.only = F)
  comp.final <- rbind(comp.final,sblast, fill = T)
  comp.final <- comp.final[!duplicated(comp.final),]

  if(verbose)
    cat("Done!\nFinal fixed-radius pruning via dbscan ... \n")
  clean.final <- clean_it(
    map = comp.final,
    genomeIDs = genomeIDs,
    rerank = T,
    radius = dbs.radius,
    n.mappings = dbs.hits,
    verbose = verbose)$map
  clean.final <- rbind(clean.final,sblast, fill = T)
  clean.final <- clean.final[!duplicated(clean.final),]

  if(verbose)
    cat("Done!\nCompleting unconnected sub-graphs ... ")
  comp.out <- complete_graph(
    map = clean.final,
    gff = gff,
    ignore.self = T,
    verbose = F)
  comp.out <- rbind(comp.out,sblast, fill = T)
  comp.out <- comp.out[!duplicated(comp.out),]

  if(verbose)
    cat("\tDone!\n")
  return(comp.out)
}






