#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_syntenicOrthogs2} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param blk The block data.frame or data.table
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
#' @param species.index The 'species.mappings' data.table or data.frame
#' from form_syntenicBlocks, giving a dictionary between pairwise
#' blast genome IDs from orthofinder, the blast file locations and
#' the genomeIDs.
#' @param n.cores Number of parallel processes to run, when possible
#' @param min.block.size The minimum block size to retain.
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
#' @export
find_syntenicOrthogs2 <- function(map,
                                  dir.list,
                                  gff,
                                  gene.index,
                                  species.mappings,
                                  n.cores = 1,
                                  plotit = F,
                                  rerank.gff = T,
                                  rank.buffer = 100,
                                  verbose = T){
  #######################################################

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
  if (verbose)
    cat("Copying orthofiner results to tmp directory ... ")

  files <- list.files(dir.list$blast,
                      full.names = T)
  files <- files[-grep("Blast|Orthogroups.txt", files)]
  nu <- file.copy(files,
                  dir.list$tmp)
  #######################################################

  #######################################################
  if (verbose)
    cat("Done!\nParsing blast results and moving to tmp directory\n")

  rewrite_cullBlast(species.mappings = species.mappings,
                    gff = gff,
                    rerank.gff = rerank.gff,
                    gene.index = gene.index,
                    map = map,
                    rank.buffer = rank.buffer,
                    n.cores = n.cores,
                    tmp.dir = dir.list$tmp, verbose = T,
                    plotit = F)
  #######################################################

  #######################################################
  if (verbose)
    cat("Done!\nRe-running orthofinder on culled blast hits ...\n")

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
  if("gene.num" %in% colnames(gff))
    gff$gene.num <- NULL
  if("rank"  %in% colnames(gff))
    gff$rank <- NULL

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

  return(list(blast = all.blast, of.info = of.blast))
}
