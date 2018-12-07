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
rerun_orthofinderInBlk = function(blk,
                                  blast.dir,
                                  init.results,
                                  cull.blast.dir,
                                  block.dir,
                                  ...){

  gff = init.results$gff
  sgff = with(merged.close, split_gffByBlock(gff = gff, blk = block))
  ogff = get_ofIDs(ogff = sgff,
                   of.geneIDs = init.results$ortho.info$gene.index,
                   of.speciesIDs = init.results$ortho.info$species.index)

  make_ofInputInBlk(blast.dir = dirs$blast,
                    tmp.dir = dirs$tmp,
                    ogff = ogff,
                    species.mappings = init.results$ortho.info$species.mappings,
                    of.speciesIDs = init.results$ortho.info$species.index)

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = dirs$tmp,
    blast.dir = dirs$block)

  culled.results = process_orthofinder(
    gff.dir = dirs$gff,
    genomeIDs = genomeIDs,
    blast.dir = dirs$block,
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

  mlo = ml[,colnames(synteny.results$map)[1:29],with = F]

  bl = make_blocks(ml,rename.blocks = T, rerank = T)

  return(bl)
}
