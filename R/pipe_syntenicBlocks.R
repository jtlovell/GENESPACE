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
#' @param verbose Logical, should updates be printed?
#' @param min.block.size Numeric, minimum required block size
#' @param plotit Logical, should plots be generated?
#' @param min.hit.density Numeric, the proportion of hits possible
#' in a block for that block to be retained.
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
pipe_syntenicBlocks = function(genomeIDs,
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
                               plotit = T,
                               min.hit.density = .2,
                               ...){

  #######################################################
  all.dirs = !any(is.null(peptide.dir),
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
    peptide.dir = dir.list$peptide
    gff.dir = dir.list$gff
    tmp.dir = dir.list$tmp
    blast.dir = dir.list$blast
    cull.blast.dir = dir.list$cull.blast
    block.dir = dir.list$block
    results.dir = dir.list$results
    mcscan.dir = dir.list$mcscan
  }
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 1: Parsing orthofinder ...\n#\tProgress:\n")

  mcsp <- paste("-a -s",min.block.size,
                "-m", min.block.size*20,
                "-w 2 -e 1")
  init.results <- process_orthofinder(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    mcscanx.input.dir = mcscan.dir,
    n.mappingWithinRadius = c(min.block.size,min.block.size,min.block.size/2),
    MCScanX.param = mcsp)
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 2: Building collinear blocks with MCScanX...\n")

  mcsp <- paste("-a -s",min.block.size,
                "-m", min.block.size,
                "-w 2")
  synteny.results = make_MCSBlocks(
    blast = init.results$blast,
    genomeIDs = genomeIDs,
    mcscanx.input.dir = mcscan.dir,
    MCScanX.params = mcsp)
  if(plotit)
    plot_blocksAndMapping(map = data.frame(synteny.results$map),
                          blk = data.frame(synteny.results$block),
                          ref.id = genomeIDs[1], altGenome2plot = genomeIDs[2])

  blk = synteny.results$block
  map = synteny.results$map
  blk$length = apply(blk[,c("rankend1","rankstart1","rankend2","rankstart2")],1,
                     function(x) max(c(x[1]-x[2],x[3]-x[4])))

  good.blocks = blk$block.id[with(blk, (n.mapping/length)>min.hit.density)]
  nmap = map[map$block.id %in% good.blocks,]
  synteny.results = make_blocks(nmap)
  if(verbose)
    cat("\tRetaining", nrow(synteny.results$map),"hits in",
        nrow(synteny.results$block),"blocks where >",
        min.hit.density*100,"% of possible hits are found in block\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 3: Merging overlapping blocks...\n")

  max.hits = max(c(table(synteny.results$map$chr1),
                   table(synteny.results$map$chr2)))

  merged.overlaps = make_mergedBlocks(
    blk = synteny.results$block,
    map = synteny.results$map,
    buffer = .1,
    n.iter = 1,
    max.size2merge = 1e6)
  if(plotit)
    plot_blocksAndMapping(map = data.frame(merged.overlaps$map),
                          blk = data.frame(merged.overlaps$block),
                          ref.id = genomeIDs[1], altGenome2plot = genomeIDs[2])
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 4: Merging adjacent blocks...\n")

  merged.close = make_mergedBlocks(
    blk = data.table(merged.overlaps$block),
    map = data.table(merged.overlaps$map),
    buffer = (sqrt((min.block.size^2)+(min.block.size^2))-.1),
    n.iter = 5,
    max.size2merge = max.hits/8)

  if(plotit)
    plot_blocksAndMapping(map = data.frame(merged.overlaps$map),
                          blk = data.frame(merged.overlaps$block),
                          ref.id = genomeIDs[1], altGenome2plot = genomeIDs[2])

  #######################################################

  #######################################################
  # if (verbose)
  #   cat("##########\n# - Part 5: Re-running overlap-merged block-constrained orthofinder...\n")
  rerun.ovlp = rerun_orthofinderInBlk(
    genomeIDs = genomeIDs,
    blk = merged.overlaps$block,
    map = merged.overlaps$map,
    tmp.dir = tmp.dir,
    gff.dir = gff.dir,
    blast.dir = blast.dir,
    block.dir = block.dir,
    init.results = init.results,
    cull.blast.dir = cull.blast.dir)

  if(plotit)
    plot_blocksAndMapping(map = data.frame(rerun.ovlp$map),
                          blk = data.frame(rerun.ovlp$block),
                          ref.id = genomeIDs[1], altGenome2plot = genomeIDs[2])

  #######################################################

  rerun.close = rerun_orthofinderInBlk(
    genomeIDs = genomeIDs,
    blk = merged.close$block,
    map = merged.close$map,
    tmp.dir = tmp.dir,
    gff.dir = gff.dir,
    blast.dir = blast.dir,
    block.dir = block.dir,
    init.results = init.results,
    cull.blast.dir = cull.blast.dir)

  if(plotit)
    plot_blocksAndMapping(map = data.frame(rerun.close$map),
                          blk = data.frame(rerun.close$block),
                          ref.id = genomeIDs[1], altGenome2plot = genomeIDs[2])

  return(list(merged.close, init.results, rerun.ovlp, rerun.close))

  # #######################################################
  # if (verbose)
  #   cat("##########\n# - Part 7: Writing results and diagnostic plots...\n")
  #
  # write.csv(rerun.ovlp$block,
  #           file = file.path(results.dir,"finalBlocks.adjacentMerge.csv"),
  #           row.names = F)
  # write.csv(rerun.ovlp$map,
  #           file = file.path(results.dir,"finalMap.adjacentMerge.csv"),
  #           row.names = F)
  #
  # write.csv(rerun.close$block,
  #           file = file.path(results.dir,"finalBlocks.overlapMerge.csv"),
  #           row.names = F)
  # write.csv(rerun.close$map,
  #           file = file.path(results.dir,"finalMap.overlapMerge.csv"),
  #           row.names = F)
  #
  #
  # pdf(file.path(results.dir, "overallDiagnostics.OverlapMerge.pdf"),
  #     height = 11, width = 8.5, useDingbats = F)
  # comb = combn(genomeIDs,2,simplify = F)
  # for(i in 1:length(comb)){
  #   with(rerun.ovlp,
  #        plot_mapping(blk = block,
  #                     map = map,
  #                     genomes = comb[[i]]))
  # }
  # dev.off()
  #
  # pdf(file.path(results.dir, "overallDiagnostics.adjacentMerge.pdf"),
  #     height = 11, width = 8.5, useDingbats = F)
  #
  # for(i in 1:length(comb)){
  #   with(rerun.close,
  #        plot_mapping(blk = block,
  #                     map = map,
  #                     genomes = comb[[i]]))
  # }
  # dev.off()
  #
  # pdf(file.path(results.dir, "byBlockDiagnostics.OverlapMerge.pdf"),
  #     height = 11, width = 8.5, useDingbats = F)
  # comb = combn(genomeIDs,2,simplify = F)
  # for(i in 1:length(comb)){
  #   with(rerun.ovlp,
  #        plot_blocksAndMapping(blk = block,
  #                              map = map,
  #                              ref.id = comb[[i]][1],
  #                              altGenome2plot = comb[[i]][2],
  #                              main = paste(comb[[i]], collapse = " --> ")))
  # }
  # dev.off()
  #
  # pdf(file.path(results.dir, "byBlockDiagnostics.adjacentMerge.pdf"),
  #     height = 11, width = 8.5, useDingbats = F)
  #
  # for(i in 1:length(comb)){
  #   with(rerun.close,
  #        plot_blocksAndMapping(blk = block,
  #                              map = map,
  #                              ref.id = comb[[i]][1],
  #                              altGenome2plot = comb[[i]][2],
  #                              main = paste(comb[[i]], collapse = " --> ")))
  # }
  # dev.off()
  #
}
