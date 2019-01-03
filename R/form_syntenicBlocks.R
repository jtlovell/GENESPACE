#' @title form_syntenicBlocks
#'
#' @description
#' \code{form_syntenicBlocks} pipe_syntenicBlocks
#'
#' @param genomeIDs Character vector of genome IDs to consider
#' @param dir.list Directory list, produced by `check_environment`.
#' @param min.block.size Numeric, minimum required block size
#' @param gap.multiplier Numeric, specification of the -m parameter in
#' MCScanX, as a function of the min.block.size.
#' @param mcscan.m.param Numeric, alternative to gap.multiplier to exactly
#' specify the MCScanX -m parameter.
#' @param merge.overlaps Logical, should overlapping blocks be merged?
#' @param clean.byDBscan Logical, should blocks be merged and cleaned by dbscan?
#' @param cull.byDBscan Logical, should blast hits be culled by dbscan?
#' @param cull.byMCscan Logical, should blast hits be culled by MCScanX?
#' @param max.size2merge the maximum block size to merge
#' @param return.ogblast Logical, should blast results of all orthogroups be
#' returned?
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
form_syntenicBlocks <- function(genomeIDs,
                                dir.list,
                                min.block.size = 3,
                                n.cores = 1,
                                gap.multiplier = 8,
                                mcscan.m.param = NULL,
                                cull.byDBscan = TRUE,
                                cull.byMCscan = TRUE,
                                return.ogblast = TRUE,
                                verbose = TRUE,
                                ...){




  #######################################################
  peptide.dir <- dir.list$peptide
  gff.dir <- dir.list$gff
  tmp.dir <- dir.list$tmp
  blast.dir <- dir.list$blast
  cull.blast.dir <- dir.list$cull.blast
  block.dir <- dir.list$block
  results.dir <- dir.list$results
  mcscan.dir <- dir.list$mcscan
  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 1: Parsing orthofinder ...\n#\tProgress:\n")

  if(!is.null(mcscan.m.param)){
    mcsp <- paste("-a -s", min.block.size,
                  "-m", mcscan.m.param,
                  "-w 2 -e 1")
  }else{
    mcsp <- paste("-a -s", min.block.size,
                  "-m", min.block.size*gap.multiplier,
                  "-w 2 -e 1")
  }

  init.results <- process_orthofinder(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    mcscan.dir = mcscan.dir,
    eps.radius = c(min.block.size*5*gap.multiplier,
                   min.block.size*2*gap.multiplier,
                   min.block.size*1*gap.multiplier),
    n.mappingWithinRadius = c(min.block.size,
                              min.block.size,
                              min.block.size),
    cull.byDBscan = cull.byDBscan,
    cull.byMCscan = cull.byMCscan,
    mcscan.param = mcsp,
    return.ogblast = return.ogblast)

  gff <- init.results$gff

  #######################################################

  #######################################################
  if (verbose)
    cat("##########\n# - Part 2: Building collinear blocks with MCScanX...\n")


  if(!is.null(mcscan.m.param)){
    mcsp <- paste("-a -s", min.block.size,
                  "-m", mcscan.m.param,
                  "-w 2")
  }else{
    mcsp <- paste("-a -s", min.block.size,
                  "-m", min.block.size*gap.multiplier,
                  "-w 2")
  }

  synteny.results <- pipe_mcs(blast = init.results$blast$map,
                              gff = gff,
                              mcscan.dir = mcscan.dir,
                              mcscan.param = mcsp)

  blk <- synteny.results$block
  map <- synteny.results$map

  synteny.results <- make_blocks(map)

  #######################################################
  #######################################################
  #######################################################
  if (verbose)
    cat("##########\n#\tDone!\n")
  out <- list(synteny.results = synteny.results,
              init.results = init.results)
  return(out)

}
