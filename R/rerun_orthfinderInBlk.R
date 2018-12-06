#' @title Run the orthofinder program
#'
#' @description
#' \code{run_orthofinder} A simple wrapper to run orthofinder from R.
#'
#' @param init.results Results from initial orthofinder run (via parse_orthofinder)
#' @param blk block file from make blocks or whatever
#' @param cull.blast.dir directory to write new output
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
rerun_orthofinderInBlk = function(blk,
                                  blast.dir,
                                  init.results,
                                  cull.blast.dir,
                                  ...){
  gff = init.results$gff
  sgff = with(blk, split_gffByBlock(gff = gff, blk = blk))
  ogff = get_ofIDs(ogff = sgff,
                   of.geneIDs = init.results$ortho.info$gene.index,
                   of.speciesIDs = init.results$ortho.info$species.index)

  make_ofInputInBlk(blast.dir = blast.dir,
                    out.dir = cull.blast.dir,
                    ogff = ogff,
                    species.mappings = init.results$ortho.info$species.mappings,
                    of.speciesIDs = init.results$ortho.info$species.index)

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = cull.blast.dir,
    blast.dir = block.dir,
    ...)

  culled.results = process_orthofinder(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    blast.dir = block.dir,
    mcscanx.input.dir = NULL,
    MCScanX.param = NULL,
    n.mappingWithinRadius = NULL)

  map = culled.results$blast
  gene.list = lapply(sgff, function(x) x$id)
  names(gene.list) = names(sgff)

  ml = rbindlist(lapply(names(gene.list), function(i){
    tmp = map[map$id1 %in% gene.list[[i]] &
                map$id2 %in% gene.list[[i]],]
    tmp$block.id = i
    return(tmp)
  }))

  bl = make_blocks(ml)
  return(bl)
}
