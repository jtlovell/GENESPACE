#' @title build_synBlocks
#'
#' @description
#' \code{build_synBlocks} build_synBlocks
#'
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param n.cores integer length 1, the number of parallel processes
#' to run.
#' @param method character length 1, specifying whether the 'pairwise' or
#' 'global' method should be employed.
#' @param use.topn logical, should the top n hit culled datasets be used.
#' @param only.orthogroups logical, should blast results be only those
#' found in the same orthogroups?
#' @param verbose logical, should updates be printed to the console?
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
rerun_orthofinder <- function(dir.list,
                              gff,
                              genomeIDs,
                              n.cores = 1,
                              use.topn = TRUE,
                              method = "pairwise",
                              only.orthogroups = T,
                              verbose = TRUE){
  #######################################################

  #######################################################
  rerun_pairwiseOF <- function(dir.list,
                               gff,
                               genomeIDs,
                               n.cores = 6,
                               verbose = T){

    gff <- subset(gff, genome %in% genomeIDs)

    tmp.dir <- dir.list$tmp
    pw.dir <- file.path(tmp.dir, "pw")
    cull.score.blast.dir <- dir.list$cull.score.blast
    peptide.dir <- dir.list$peptide

    if (file.exists(pw.dir))
      unlink(pw.dir, recursive = T)
    dir.create(pw.dir)
    on.exit(expr = unlink(pw.dir, recursive = T))

    if (verbose)
      cat("Preparing new orthofinder-formatted species ID database ... \n")
    map.db <- prep_ofDB(
      tmp.dir = tmp.dir,
      orig.dir = cull.score.blast.dir,
      peptide.dir = peptide.dir,
      output.dir = pw.dir,
      n.cores = n.cores,
      min.score = 0,
      verbose = FALSE,
      genomeIDs = genomeIDs)
    if (verbose)
      cat("\tDone!\n")

    com <- paste("orthofinder", "-b", pw.dir,
                 "-a", n.cores,
                 "-og 1>/dev/null 2>&1")
    system(com)

    gf <- read_ogs(
      of.dir = pw.dir,
      gff = gff)

    return(gf)
  }
  #######################################################

  #######################################################
  if(method == "pairwise"){
    if (verbose)
      cat("Running pairwise orthofinder calls ... \n")

    comb <- combn(genomeIDs, 2, simplify = F)
    combl <- rbindlist(lapply(1:length(comb), function(i)
      data.table(run = i,
                 expand.grid(run.genome1 = comb[[i]],
                             run.genome2 = comb[[i]]))))
    combl <- combl[!duplicated(combl[,c("run.genome1","run.genome2")]),]

    gfl <- rbindlist(lapply(1:length(comb), function(i){
      x <- comb[[i]]
      if (verbose)
        cat(paste0("\t",x[1]), "<-->", x[2],"... ")
      gf <- rerun_pairwiseOF(dir.list = dir.list,
                                genomeIDs = x,
                                n.cores = n.cores,
                                gff = gff,
                                verbose = F)
      gf[,run.genome1 := x[1]]
      gf[,run.genome2 := x[2]]
      if (verbose)
        cat("found", nrow(gf), "hits",
            with(gf,
                 sum(genome1 != genome2)),
            "in orthogroups\n")
      return(gf)
    }))

    if(verbose)
      cat("Processing results ... ")
    gfl <- merge(combl[,c("run.genome1","run.genome2")],
                 gfl,
                 by = c("run.genome1","run.genome2"))
    gfl[,run.genome1 := NULL]
    gfl[,run.genome2 := NULL]

    blast <- import_allBlasts(
      gff = gff,
      keep.geneNum = F,
      add.orthogroups = F,
      orthofinder.dir = dir.list$cull.score.blast,
      genomeIDs = genomeIDs,
      blast.files = list.files(
        path = dir.list$cull.score.blast,
        pattern = "^Blast",
        full.names = T))

    out <- merge(gfl,
                 blast,
                 by = c("genome1","genome2","id1","id2"))
  }else{
    tmp.dir <- dir.list$tmp
    if (use.topn) {
      blast.dir <- dir.list$cull.score.blast
    }else{
      blast.dir <- dir.locs$cull.blast
    }

    if (dir.exists(tmp.dir))
      un <- unlink(tmp.dir, recursive = T)
    dir.create(tmp.dir)

    fs <- list.files(blast.dir, full.names = T)
    fc <- file.copy(from = fs,
                    to = tmp.dir)

    com <- paste("orthofinder", "-b", tmp.dir,
                 "-a", n.cores,
                 "-og")
    system(com)

    blast <- import_allBlasts(
      gff = gff,
      keep.geneNum = F,
      add.orthogroups = T,
      orthofinder.dir = tmp.dir,
      genomeIDs = genomeIDs)

    out <- data.table(blast)
  }

  comb <- data.table(rbind(t(combn(genomeIDs, 2, simplify = T)),
                           t(sapply(genomeIDs, rep, 2))))
  setnames(comb, c("genome1", "genome2"))
  out <- merge(out, comb, by = c("genome1", "genome2"))

  if(verbose)
    cat("Done!\n")
  return(mirror_map(out))
}
