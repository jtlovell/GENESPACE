#' @title build pwBlocks
#'
#' @description
#' \code{build_pwBlocks} build_pwBlocks
#'
#' @param dir.list dir.list
#' @param genomeIDs genomeIDs
#' @param mcscan.dir path to the directory where MCScanX will be run
#' @param gap.multiplier parameters to supply MCScanX
#' @param MCScanX.path the path to the the MCScanX program. If in the path,
#' just use "MCScanX".
#' @param min.blockSize min.blockSize
#' @param radius radius
#' @param clean.before.mcscanx clean.before.mcscanx
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
build_pwBlocks <- function(dir.list,
                           genomeIDs,
                           of.cores = 6,
                           min.blockSize = 5,
                           gap.multiplier = 8,
                           clean.before.mcscanx = T,
                           MCScanX.path,
                           radius = 100,
                           runFrom.pw.of = NULL,
                           verbose = T){

  if(verbose)
    cat("Importing gff annotations as data.tables ... \n")
  gff <- import_gff(gff.dir = dir.list$gff,
                    genomeIDs = genomeIDs)
  if(verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if(verbose)
    cat("Running pairwise orthofinder calls ... \n")
  comb <- combn(genomeIDs, 2, simplify = F)
  if(is.null(runFrom.pw.of)){
    pw.of <- lapply(comb, function(x){
      if(verbose)
        cat(paste0("\t",x[1]), "<-->", x[2],"... ")
      out <- rerun_pairwiseOF(dirs = dir.list,
                              genomeIDs = x,
                              of.cores = 6,
                              gff = gff,
                              verbose = F)
      if(verbose)
        cat("Done!\n")
      return(out)
    })
    if(verbose)
      cat("Saving pairwise gff data.tables as", file.path(getwd(),"pw.of.rda"),"... ")
    save(pw.of, file = "pw.of1.rda")

    if(verbose)
      cat("Done!\n")
  }else{
    load(runFrom.pw.of)
  }

  #######################################################

  #######################################################
  if(verbose)
    cat("Merging pairwise hits with gff annotations\n")
  pw.map <- merge_ofGff(comb = comb, pw.of = pw.of, gff = gff)
  if(verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if(clean.before.mcscanx){
    cull.dbs <- clean_blocks(
      map = pw.map,
      n.mappings = min.blockSize,
      radius = 100,
      n.cores = 1,
      verbose = T,
      clean.columns = F)
    cull.dbs <- cull.dbs$map
  }else{
    cull.dbs <- pw.map
  }
  #######################################################

  #######################################################
  if(verbose)
    cat("Forming syntenic blocks via MCScanX ... \n")
  mcsp <- paste("-a -s", min.blockSize,
                "-m", min.blockSize*gap.multiplier,
                "-w 2 -e 1")

  syn.blks <- pipe_mcscanx(blast = pw.map,
                           gff = gff,
                           MCScanX.path = MCScanX.path,
                           mcscan.dir = dir.list$mcscanx,
                           mcscan.param = mcsp,
                           verbose = T)
  if(verbose)
    cat("Done!")
  return(syn.blks)
}
