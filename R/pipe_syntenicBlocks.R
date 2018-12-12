#' @title pipe_syntenicBlocks
#'
#' @description
#' \code{pipe_syntenicBlocks} pipe_syntenicBlocks
#'
#' @param genomeIDs Character vector of genome IDs to consider
#' @param dir.list Directory list, produced by `check_environment`.
#' @param peptide.dir Directory containing peptide fasta files
#' @param gff.dir Directory containing peptide fasta files
#' @param tmp.dir Directory containing peptide fasta files
#' @param blast.dir Directory containing peptide fasta files
#' @param block.dir Directory containing peptide fasta files
#' @param cull.blast.dir Directory containing peptide fasta files
#' @param mcscan.dir Directory containing peptide fasta files
#' @param results.dir Directory containing peptide fasta files
#' @param min.block.size Numeric, minimum required block size
#' @param plotit Logical, should plots be generated?
#' @param min.hit.density Numeric, the proportion of hits possible
#' in a block for that block to be retained.
#' @param min.unique.hits The number of unique hits required in each genome
#' for a block to be retained.
#' @param max.hit.density The maximum hit density in a block
#' @param initial.mergeBuffer Numeric, the buffer for an initial scan
#' to merge overlapping blocks. Should be set to either 0 or 1, unless
#' one is ok with sacraficing speed for accuracy, in which case,
#' it can be set anywhere up to min.block.size - 1.
#' @param final.mergeBuffer Numeric, the buffer to search for blocks to merge
#' in the last merge step. Should not exceed min.block.size -1.
#' @param max.propSmallChr2Merge Numeric, the maximum size of two blocks to be
#' merged is set based on this parameter. If set to 1, entire chromosomes
#' may be merged, if set to 0, no merging ever happens. Values between 0 and
#' 1, result in a maximum block size, equal to the proportion of the number of
#' blast hits on the smallest chromosome.
#' @param n.mergeIter Numeric, the number of iterations that the merge should
#' go through.
#' @param return.initial Logical, should the initial orthofinder results
#' be returned?
#' @param return.mcscanRaw Logical, should the unmerged results from MCScanX
#' be returned?
#' @param return.overlapMerged Logical, should the overlap-merged
#' orthofinder results be returned?
#' @param verbose Logical, should updates be printed?
#' @param ... Not in use yet.
#' @details More here
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
pipe_syntenicBlocks <- function(genomeIDs,
                               dir.list = NULL,
                               peptide.dir = NULL,
                               gff.dir = NULL,
                               tmp.dir = NULL,
                               blast.dir = NULL,
                               block.dir = NULL,
                               cull.blast.dir = NULL,
                               mcscan.dir = NULL,
                               results.dir = NULL,
                               verbose = TRUE,
                               min.block.size = 5,
                               min.unique.hits = min.block.size,
                               plotit = T,
                               min.hit.ratio = .2,
                               initial.mergeBuffer = 1,
                               n.cores = 1,
                               final.mergeBuffer = NULL,
                               return.initial = TRUE,
                               return.mcscanRaw = FALSE,
                               return.overlapMerged = FALSE,
                               max.size2merge = 100,
                               ...){

  if (is.null(final.mergeBuffer)) {
    final.mergeBuffer <- (sqrt((min.block.size^2) + (min.block.size^2)) - .1)
  }

  #######################################################
  all.dirs <- !any(is.null(peptide.dir),
                  is.null(gff.dir),
                  is.null(tmp.dir),
                  is.null(blast.dir),
                  is.null(cull.blast.dir),
                  is.null(mcscan.dir),
                  is.null(block.dir))
  if (is.null(dir.list) & !all.dirs) {
    stop("Either a directory list or each directory parameter must be specified\n")
  }

  if (!is.null(dir.list)) {
    peptide.dir <- dir.list$peptide
    gff.dir <- dir.list$gff
    tmp.dir <- dir.list$tmp
    blast.dir <- dir.list$blast
    cull.blast.dir <- dir.list$cull.blast
    block.dir <- dir.list$block
    results.dir <- dir.list$results
    mcscan.dir <- dir.list$mcscan
  }
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 1: Parsing orthofinder ...\n#\tProgress:\n")

  mcsp <- paste("-a -s", min.block.size,
                "-m", min.block.size*20,
                "-w 2 -e 1")
  init.results <- process_orthofinder(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    mcscanx.input.dir = mcscan.dir,
    n.mappingWithinRadius = c(min.block.size,
                              min.block.size,
                              min.block.size/2),
    MCScanX.param = mcsp)
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 2: Building collinear blocks with MCScanX...\n")

  mcsp <- paste("-a -s",min.block.size,
                "-m", min.block.size,
                "-w 2")
  synteny.results <- make_MCSBlocks(
    blast = init.results$blast,
    genomeIDs = genomeIDs,
    mcscanx.input.dir = mcscan.dir,
    MCScanX.params = mcsp,
    min.hit.ratio = min.hit.ratio,
    min.unique.hits = min.unique.hits)
  synteny.results.out = synteny.results

  if (plotit)
    plot_blocksAndMapping(
      map = data.frame(synteny.results$map),
      blk = data.frame(synteny.results$block),
      ref.id = genomeIDs[1],
      altGenome2plot = genomeIDs[2])

  blk <- synteny.results$block
  map <- synteny.results$map

  synteny.results <- make_blocks(map)

  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 3: Merging overlapping blocks...\n")

  max.hits <- min(c(table(synteny.results$map$chr1),
                   table(synteny.results$map$chr2)))

  merged.overlaps <- merge_blocks(
    blk = synteny.results$block,
    map = synteny.results$map,
    buffer = initial.mergeBuffer,
    n.cores =  n.cores,
    max.size2merge = 1e6)
  if (plotit)
    plot_blocksAndMapping(
      map = data.frame(merged.overlaps$map),
      blk = data.frame(merged.overlaps$block),
      ref.id = genomeIDs[1],
      altGenome2plot = genomeIDs[2])
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 4: Merging adjacent blocks...\n")

  merged.close <- merge_blocks(
    blk = data.table(merged.overlaps$block),
    map = data.table(merged.overlaps$map),
    buffer = final.mergeBuffer,
    n.cores = n.cores,
    max.size2merge = max.size2merge)

  if (plotit)
    plot_blocksAndMapping(
      map = data.frame(merged.overlaps$map),
      blk = data.frame(merged.overlaps$block),
      ref.id = genomeIDs[1],
      altGenome2plot = genomeIDs[2])

  #######################################################
  if (verbose)
    cat("##########\n# - Part 5: Re-running orthofinder constrained within known syntenic blocks ...\n")
  rerun.close <- rerun_orthofinderInBlk(
    genomeIDs = genomeIDs,
    blk = merged.close$block,
    map = merged.close$map,
    tmp.dir = tmp.dir,
    gff.dir = gff.dir,
    blast.dir = blast.dir,
    block.dir = block.dir,
    init.results = init.results,
    cull.blast.dir = cull.blast.dir)

  if (plotit)
    plot_blocksAndMapping(
      map = data.frame(rerun.close$map),
      blk = data.frame(rerun.close$block),
      ref.id = genomeIDs[1],
      altGenome2plot = genomeIDs[2])

  out <- list(merged.results = merged.close,
             rerun.results = rerun.close)

  if (return.initial) {
    out[["initial.results"]] <- init.results
  }
  if (return.mcscanRaw) {
    out[["mcscanx.results"]] <- synteny.results.out
  }
  if (return.overlapMerged) {
    out[["overlap.results"]] <- merged.overlaps
  }

  return(out)
}
