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
  mirror_blast <- function(blast,
                           w2ki,
                           genomes){
    g1 <- genomes[1]
    g2 <- genomes[2]
    bl <- data.table(subset(blast, genome1 != genome2))
    bl2 <- data.table(bl)
    setnames(bl2,
             grep("1$", colnames(bl2)),
             gsub("1$","xxxx",
                  colnames(bl2)[grep("1$", colnames(bl2))]))
    setnames(bl2,
             grep("2$", colnames(bl2)),
             gsub("2$","1",
                  colnames(bl2)[grep("2$", colnames(bl2))]))
    setnames(bl2,
             grep("xxxx", colnames(bl2)),
             gsub("xxxx","2", colnames(bl2)[grep("xxxx", colnames(bl2))]))

    blo <- rbind(bl, bl2)
    blo <- subset(blo, genome1 == g1 & genome2 == g2)
    blo <- rbind(blo, subset(blast, genome1 == genome2))
    blo[,ns := -score]
    setkey(blo, ns)
    blo <- blo[!duplicated(blo[,c("id1","id2")]),]
    blo$ns <- NULL

    if (all(genomes %in% names(w2ki)))
      out <- blo
    if (!any(genomes %in% names(w2ki)))
      out <- subset(blo, genome1 != genome2)
    if (g1 %in% names(w2ki) & !g2 %in% names(w2ki))
      out <- subset(blo, genome1 == g1 | genome2 == g1)
    if (g2 %in% names(w2ki) & !g1 %in% names(w2ki))
      out <- subset(blo, genome1 == g2 | genome2 == g2)
    return(out)
  }
  #######################################################

  #######################################################
  rerun_pairwiseOF <- function(dir.list,
                               gff,
                               genomeIDs,
                               n.cores = 6,
                               verbose = T){

    gff <- subset(gff, genome %in% genomeIDs)

    if (verbose)
      cat("Preparing new orthofinder-formatted species ID database ... \n")

    pw.dir <- file.path(dirname(dir.list$tmp), "pw")
    if (file.exists(pw.dir))
      unlink(pw.dir, recursive = T)
    dir.create(pw.dir)

    on.exit(expr = unlink(pw.dir, recursive = T))

    make_newOFdb(tmp.dir = dir.list$tmp,
                 cull.blast.dir = pw.dir,
                 peptide.dir = dir.list$peptide,
                 genomeIDs = genomeIDs,
                 verbose = verbose)
    if (verbose)
      cat("\tDone!\n")

    if (verbose)
      cat("Importing new and old orthofinder gene and species IDs ... ")
    old.ids <- read_speciesIDs(of.dir = dir.list$cull.score.blast,
                               genomeIDs = genomeIDs)
    new.ids <- read_speciesIDs(of.dir = pw.dir,
                               genomeIDs = genomeIDs)
    id.db <- merge(old.ids, new.ids, by = "genome")
    setnames(id.db,2:3,c("n.old","n.new"))

    map.db <- make_mapDB(id.db = id.db,
                         of.dir = dir.list$cull.score.blast,
                         cull.blast.dir = pw.dir)

    old.genes <- read_geneIDs(of.dir = dir.list$cull.score.blast,
                              gff = gff)
    new.genes <- read_geneIDs(of.dir = pw.dir,
                              gff = gff)
    genes <- merge(old.genes[,c("genome","id","gene.num")],
                   new.genes[,c("genome","id","gene.num")],
                   by = c("genome","id"))
    g1 <- with(genes,
               data.table(V1 = gene.num.x,
                          new1 = gene.num.y,
                          key = "V1"))
    g2 <- with(genes,
               data.table(V2 = gene.num.x,
                          new2 = gene.num.y,
                          key = "V2"))
    if (verbose)
      cat("\tDone!\n")

    if (verbose)
      cat("Reading blast file, replacing IDs and renaming ... \n")
    in.blast.files <- map.db$filename
    out.blast.files <- map.db$new.filename

    d <- lapply(1:nrow(map.db), function(i){
      if (verbose)
        cat(paste0("\t",map.db$genome1[i]),
            "-->",map.db$genome2[i],"... ")
      bl <- readRename_blastGenes(gene.dict1 = g1,
                                  gene.dict2 = g2,
                                  blast.file.in = in.blast.files[i],
                                  blast.file.out = out.blast.files[i],
                                  verbose = verbose)
      if (verbose)
        cat("Done!\n")
      return(bl)
    })

    com <- paste("orthofinder", "-b", pw.dir,
                 "-a", n.cores,
                 "-og 1>/dev/null 2>&1")
    system(com)


    gf <- read_ogs(of.dir = pw.dir, gff = gff)
    blast <- read_allBlasts(gff = gf,
                            keep.geneNum = F,
                            add.gff = T,
                            check.ogs = F,
                            of.dir = pw.dir,
                            genomeIDs = genomeIDs,
                            verbose = F)
    return(blast)
  }
  #######################################################

  #######################################################
  if(method == "pairwise"){
    if (verbose)
      cat("Running pairwise orthofinder calls ... \n")
    comb <- combn(genomeIDs, 2, simplify = F)
    wh.2keep <- sapply(genomeIDs, USE.NAMES = T, function(x)
      min(which(sapply(comb, function(y) x %in% y))))

    out <- rbindlist(lapply(1:length(comb), function(i){
      x <- comb[[i]]
      w2ki <- wh.2keep[wh.2keep == i]

      if (verbose)
        cat(paste0("\t",x[1]), "<-->", x[2],"... ")
      blast <- rerun_pairwiseOF(dir.list = dir.list,
                                genomeIDs = x,
                                n.cores = n.cores,
                                gff = gff,
                                verbose = F)
      bl <- mirror_blast(blast = blast,
                         w2ki = w2ki,
                         genomes = x)
      if (verbose)
        cat("found", nrow(bl), "hits",
            with(bl,
                 sum(og1 == og2 &
                       genome1 != genome2)),
            "in orthogroups\n")

      if (only.orthogroups) {
        bl <- subset(bl, og1 == og2)
      }
      bl[, og.id := og1]
      return(bl)
    }))
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
    fc <- file.copy(from = fs, to = tmp.dir)

    com <- paste("orthofinder", "-b", tmp.dir,
                 "-a", n.cores,
                 "-og")
    system(com)

    gf <- read_ogs(of.dir = tmp.dir,
                   gff = gff)
    out <- read_allBlasts(gff = gf,
                          keep.geneNum = F,
                          add.gff = T,
                          check.ogs = F,
                          of.dir = tmp.dir,
                          genomeIDs = genomeIDs,
                          verbose = F)

    if (only.orthogroups) {
      out <- subset(out, og1 == og2)
    }
    out[, og.id := og1]
  }

  comb <- data.table(rbind(t(combn(genomeIDs, 2, simplify = T)),
                           t(sapply(genomeIDs, rep, 2))))
  setnames(comb, c("genome1", "genome2"))
  out <- merge(out, comb, by = c("genome1", "genome2"))

  return(out)
}
