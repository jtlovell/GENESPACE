#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param dir.list The directory list produced by check_environment
#' @param blast data.table with at least the 10 necessary columns for blast format.
#' @param genomeIDs character vector giving genome IDs to consider.
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param n.cores Number of parallel processes to run, when possible
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
find_syntenicOrthogs <- function(blast,
                                 dir.list,
                                 gff,
                                 genomeIDs,
                                 n.cores = 1,
                                 verbose = T,
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
    og.threads = n.cores,
    og.silent = verbose,
    verbose = verbose)

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  of.blast <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = dir.list$cull.blast,
    verbose = verbose)
  #######################################################

  #######################################################
  all.blast <- import_ofBlast(
    species.mappings = of.blast$species.mappings,
    genomeIDs = genomeIDs,
    orthogroups = of.blast$orthogroups,
    only.orthologs = TRUE,
    gff = gff,
    gene.index = of.blast$gene.index,
    verbose = verbose)

  all.blast$unique = with(all.blast, paste0(genome1, "_", genome2))

  #######################################################
  return(list(blast = all.blast,
              of.results = of.blast))
}
