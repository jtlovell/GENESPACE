#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param map The map data.frame or data.table
#' @param dir.list The directory list produced by check_environment
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param gene.index orthofinder geneID data.table or data.frame
#' giving a dictionary between the 'id' column in the gff object
#' and the 'gene.num' numeric geneIDs in the orthofinder-formatted
#' blast files.
#' @param species.mappings The 'species.mappings' data.table or data.frame
#' from form_syntenicBlocks, giving a dictionary between pairwise
#' blast genome IDs from orthofinder, the blast file locations and
#' the genomeIDs.
#' @param orthogroups The orthogroups list.
#' @param plotit Logical, should plots be made? Will not work with
#' n.core > 1.
#' @param rank.buffer The buffer, in gene rank order.
#' @param n.cores Number of parallel processes to run, when possible
#' @param min.block.size Numeric of length 1, specifying the minimum
#' block size in the cleaning procedure at the end.
#' @param verbose Logical, should updates be printed
#' @param ... Not currently in use
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom parallel mclapply
#' @export
find_syntenicOrthogs <- function(map,
                                 blast,
                                 dir.list,
                                 gff,
                                 gene.index,
                                 genomeIDs,
                                 species.mappings,
                                 n.cores = 1,
                                 orthogroups,
                                 plotit = n.cores == 1,
                                 rank.buffer = 250,
                                 verbose = T,
                                 min.block.size = 5,
                                 ...){

  #######################################################
  if (verbose)
    cat("Cleaning out tmp and culled blast directories ... ")

  unlink(dir.list$tmp, recursive = T)
  dir.create(dir.list$tmp)

  unlink(dir.list$cull.blast, recursive = T)
  dir.create(dir.list$cull.blast)

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  write_ofData(blast = blast,
               genomeIDs = genomeIDs,
               of.dir = dir.list$tmp,
               peptide.dir = dir.list$peptide,
               verbose = verbose)
  #######################################################

  #######################################################
  if (verbose)
    cat("Re-running orthofinder on culled blast hits ...\n")

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = dir.list$tmp,
    blast.dir = dir.list$cull.blast,
    og.threads = 6,
    og.silent = F,
    verbose = T)

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  of.blast <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = dir.list$cull.blast,
    verbose = T)
  #######################################################

  #######################################################
  all.blast <- import_ofBlast(
    species.mappings = of.blast$species.mappings,
    genomeIDs = genomeIDs,
    orthogroups = of.blast$orthogroups,
    gff = gff,
    gene.index = of.blast$gene.index,
    verbose = T)

  all.blast$unique = with(all.blast, paste0(genome1, "_", genome2))

  #######################################################
  return(list(blast = all.blast,
              of.results = of.blast))
}
