#' @title pipe_syntenicBlocks
#'
#' @description
#' \code{pipe_syntenicBlocks} pipe_syntenicBlocks
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
pipe_syntenicBlocks <- function(genomeIDs,
                                dir.list,
                                min.block.size = 3,
                                n.cores = 1,
                                max.size2merge = 1e6,
                                gap.multiplier = 8,
                                mcscan.m.param = NULL,
                                merge.overlaps = F,
                                clean.byDBscan = T,
                                cull.byDBscan = T,
                                cull.byMCscan = T,
                                return.ogblast = T,
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

  if(!clean.byDBscan){
    merged_blk <- synteny.results
    merged_blk.out <- NULL
  }else{
    if (verbose)
      cat("##########\n# - Part 3: Cleaning / merging via dbscan\n",
          "initial n blocks / mappings =",
          nrow(synteny.results$block),
          "/",
          nrow(synteny.results$map),"\n\t")
    tmp = synteny.results$map
    tmp$unique = with(tmp, paste(genome1, genome2, chr1, chr2))
    spl.map = split.data.table(tmp, "unique")

    min.rad = min.block.size
    merged_map = rbindlist(lapply(spl.map, function(tmp){
      x <- run_dbs(y = tmp[,c("rank1","rank2"),with = F],
                   eps.radius = min.rad+.1,
                   mappings = min.block.size)
      tmp$block.id = x$cluster
      cols = sample(rainbow(length(x$cluster)),size = length(x$cluster), replace = F)
      return(tmp)
    }))
    merged_map = merged_map[merged_map$block.id != 0,]

    merged_map$block.id = with(merged_map, as.numeric(as.factor(paste(unique, block.id))))
    merged_blk = make_blocks(merged_map,rename.blocks = F, rerank = T)
    merged_blk.out = merged_blk
    if (verbose)
      cat("cleaned n blocks / mappings =",
          nrow(merged_blk$block),
          "/",
          nrow(merged_blk$map),"\n")
  }

  #######################################################
  #######################################################
  #######################################################

  if(!merge.overlaps){
    merged.overlaps <- NULL
  }else{
    if (verbose)
      cat("##########\n# - Part 4: Merging overlapping blocks\n",
          "initial n blocks / mappings =",
          nrow(merged_blk$block),
          "/",
          nrow(merged_blk$map),"\n\t")

    merged.overlaps <- merge_blocks(
      blk = merged_blk$block,
      map = merged_blk$map,
      buffer = 0,
      n.cores =  n.cores,
      max.size2merge = 1e6)

    if (verbose)
      cat("merged to n blocks / mappings =",
          nrow(merged.overlaps$block),
          "/",
          nrow(merged.overlaps$map),"\n")

  }

  out <- list(synteny.results = synteny.results,
              init.results = init.results,
              merged.dbs = merged_blk.out,
              merged.results = merged.overlaps)
  return(out)

}
