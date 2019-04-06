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
#' @importFrom utils combn
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
  #######################################################
  # -- Error checking
  stop_withMessage(c(dir.exists(unlist(dir.list))),
                   paste(dir.list[!dir.exists(unlist(dir.list))], "does not exist"))
  stop_withMessage(c(is.character(genomeIDs),
                     length(genomeIDs) > 1,
                     missing(genomeIDs)),
                   "genomeIDs must be a character vector of length > 1")
  stop_withMessage(c(is.numeric(of.cores),
                     is.numeric(min.blockSize),
                     is.numeric(gap.multiplier),
                     is.numeric(radius)),
                   "of.cores, min.blockSize, gap.multiplier, radius all must be numeric")
  stop_withMessage(c(dir.exists(MCScanX.path), file.exists(),
                     is.character(str2parse)),
                   "str2parse must be a single character")

  MCScanX.tool <- file.path(MCScanX.path,"MCScanX")
  stop_withMessage(c(is.logical(verbose),
                     is.logical(clean.before.mcscanx)),
                   "verbose and clean.before.mcscanx must be logical")
  stop_withMessage(is.null(runFrom.pw.of) | file.exists(runFrom.pw.of),
                   "verbose and clean.before.mcscanx must be logical")
  of.cores <- as.integer(of.cores)
  min.blockSize <- as.integer(min.blockSize)
  gap.multiplier <- as.integer(gap.multiplier)
  radius <- as.integer(radius)
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing gff annotations as data.tables ... \n")
  gff <- import_gff(gff.dir = dir.list$gff,
                    genomeIDs = genomeIDs)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Running pairwise orthofinder calls ... \n")
  comb <- combn(genomeIDs, 2, simplify = F)
  if (is.null(runFrom.pw.of)) {
    pw.of <- lapply(comb, function(x){
      if (verbose)
        cat(paste0("\t",x[1]), "<-->", x[2],"... ")
      out <- rerun_pairwiseOF(dirs = dir.list,
                              genomeIDs = x,
                              of.cores = 6,
                              gff = gff,
                              verbose = F)
      if (verbose)
        cat("Done!\n")
      return(out)
    })
    if (verbose)
      cat("Saving pairwise gff data.tables as",
          file.path(getwd(),"pw.of.rda"),"... ")
    save(pw.of, file = "pw.of1.rda")

    if (verbose)
      cat("Done!\n")
  }else{
    load(runFrom.pw.of)
  }
  #######################################################

  #######################################################
  if (verbose)
    cat("Merging pairwise hits with gff annotations\n")
  pw.map <- merge_ofGff(comb = comb,
                        pw.of = pw.of,
                        gff = gff,
                        dir.list = dir.list)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  genome.dt <- data.table(rbind(
    t(combn(genomeIDs,2)),
    cbind(genomeIDs,genomeIDs)))
  setnames(genome.dt, c("genome1","genome2"))
  pw.map <- merge(genome.dt, pw.map, by = c("genome1","genome2"))
  if (clean.before.mcscanx) {
    cull.dbs <- clean_blocks(
      map = pw.map,
      n.mappings = min.blockSize,
      radius = radius,
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
                           genomeIDs = genomeIDs,
                           MCScanX.path = MCScanX.path,
                           mcscan.dir = dir.list$mcscanx,
                           mcscan.param = mcsp,
                           verbose = T)
  if(verbose)
    cat("Done!")
  return(syn.blks)
}
