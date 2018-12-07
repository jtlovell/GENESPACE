#' @title Run the orthofinder program
#'
#' @description
#' \code{run_orthofinder} A simple wrapper to run orthofinder from R.
#'
#' @param init.results Results from initial orthofinder run (via parse_orthofinder)
#' @param blk block file from make blocks or whatever
#' @param cull.blast.dir directory to write new output
#' @param block.dir directory to store the block output
#' @param blast.dir The path to the directory where the blast results should be stored
#' @param ... Additional arguments passed on to run_orthofinder
#' @details ...

#' @return Nothing, writes results to the cull.blast.dir directory
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
rerun_orthofinderInBlk = function(map, blk,
                                  gff.dir,
                                  blast.dir,
                                  init.results,
                                  cull.blast.dir,
                                  block.dir,
                                  tmp.dir,
                                  genomeIDs,
                                  ...){

  gff = init.results$gff
  sgff = split_gffByBlock(gff = gff, blk = blk)
  ogff = get_ofIDs(ogff = sgff,
                   of.geneIDs = init.results$ortho.info$gene.index,
                   of.speciesIDs = init.results$ortho.info$species.index)

  make_ofInputInBlk(blast.dir = blast.dir,
                    tmp.dir = tmp.dir,
                    ogff = ogff,
                    species.mappings = init.results$ortho.info$species.mappings,
                    of.speciesIDs = init.results$ortho.info$species.index)

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = tmp.dir,
    blast.dir = cull.blast.dir)

  culled.results = process_orthofinder(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    blast.dir = cull.blast.dir,
    mcscanx.input.dir = NULL,
    MCScanX.param = NULL,
    n.mappingWithinRadius = NULL)

  map = culled.results$blast
  gene.list = lapply(sgff, function(x) x$id)
  names(gene.list) = names(sgff)

  rr = rbindlist(lapply(names(gene.list), function(i){
    tmp = map[map$id1 %in% gene.list[[i]] &
                map$id2 %in% gene.list[[i]],]
    tmp$block.id = i
    return(tmp)
  }))

  rr$rank1 = frank(rr, chr1, start1, ties.method = "dense")
  rr$rank2 = frank(rr, chr2, start2, ties.method = "dense")
  bl = make_blocks(data.frame(rr),rename.blocks = T)

  return(bl)
}
