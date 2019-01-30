#' @title build_syntenicBlocks
#'
#' @description
#' \code{build_syntenicBlocks} pipe_syntenicBlocks
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
build_syntenicBlocks <- function(genomeIDs,
                                 dir.list,
                                 min.block.size = 5,
                                 n.cores = 1,
                                 gap.multiplier = 8,
                                 mcscan.m.param = NULL,
                                 cull.byDBscan = TRUE,
                                 cull.byMCscan = TRUE,
                                 return.ogblast = TRUE,
                                 verbose = TRUE,
                                 str2drop = "Name=",
                                 str2parse = ";",
                                 whichAttr = 2,
                                 n.mappingWithinRadius = c(5,5),
                                 eps.radius = c(30,15),
                                 ...){

  #######################################################
  gff.dir <- dir.list$gff
  blast.dir <- dir.list$blast
  mcscan.dir <- dir.list$mcscan
  #######################################################

  #######################################################
  if (verbose)
    cat("Part #1 -- Importing and parsing gff3-formatted annotations ...\n")
  gff <- import_gff(
    gff.dir = gff.dir,
    genomeIDs = genomeIDs,
    verbose = verbose,
    str2drop = str2drop,
    str2parse = str2parse,
    whichAttr = whichAttr)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Part #2 -- Importing and parsing orthofinder output ...\n")
  of.results <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    verbose = verbose)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Part #3 -- Importing blast hits and merging with orthofinder results ...\n")
  blast <- with(of.results,
                import_ofBlast(
                  species.mappings = species.mappings,
                  genomeIDs = genomeIDs,
                  orthogroups = orthogroups,
                  gff = gff,
                  gene.index = gene.index,
                  verbose = verbose))
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if(cull.byDBscan){
    if (verbose)
      cat("Part #4 -- Culling blast results by 2d density ...\n")

    cull.tmp <- clean_blocks(
      map = blast,
      n.mappings = n.mappingWithinRadius,
      radius = eps.radius,
      n.cores = n.cores,
      clean.by.unique.genes = F,
      clean.by.og = F,
      clean.columns = F,
      verbose = T)

    cull.dbs<-cull.tmp$map

    if (verbose)
      cat("\tDone!\n")

    if (!cull.byMCScanX){
      if (verbose)
        cat("Skipping Part #5 (Culling blast results by multiple-collinearity)\n",
          "Returning clustering from dbscan\nDone!\n")
      map <- cull.tmp$map
      blk <- cull.tmp$block
    }
  }else{
    if (verbose)
      cat("Skipping Part #4 (Culling blast results by 2d density)\n")
    cull.dbs <- blast
  }
  #######################################################

  #######################################################
  if(cull.byMCScanX){
    if (verbose)
      cat("Part #5 -- Culling blast results by multiple-collinearity ...\n")

    if (!is.null(mcscan.m.param)) {
      mcsp <- paste("-a -s", min.block.size,
                    "-m", mcscan.m.param,
                    "-w 2 -e 1")
    } else {
      mcsp <- paste("-a -s", min.block.size,
                    "-m", min.block.size*gap.multiplier,
                    "-w 2 -e 1")
    }

    synteny.results <- pipe_mcscanx(blast = blast,
                                    gff = gff,
                                    mcscan.dir = mcscan.dir,
                                    mcscan.param = mcsp,
                                    verbose = T)
    map = synteny.results$map
    blk = synteny.results$block
    if (verbose)
      cat("\tDone!\n")

  }
  #######################################################
  if (verbose)
    cat("##########\n#\tDone!\n")
  out <- list(synteny.results = list(map = map,
                                     block = blk),
              init.results = init.results)
  return(out)
}
