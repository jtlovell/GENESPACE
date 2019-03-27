#' @title Pseudogene utility functions
#' @description
#' \code{pwof_utils} Five utilities functions meant for internal calls in compareGeneSpace
#' @name pwof_utils
#'
#' @param map map results data.table
#' @param assembly.dir path to assembly fastas
#' @param tmp.dir path to temp directory
#' @param peptide.dir path to peptide directory
#' @param buffer numeric, the number of basepairs outside of the range to look at
#' @param genomeIDs character, indicating genomeIDs to consider.
#' @param clean logical, should the intermediate files be removed?
#' @param diamond.sensitive logical, should diamond blastx be run in the
#' sensitive mode?
#'
#' @param verbose logical, should updates be printed?
#' @param ... not currently in use
#'
#' @note \code{pwof_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{pwof_utils} returns its own arguments.
#'
#' @title read_ogs
#' @description
#' \code{read_ogs} fread_ogs
#' @rdname pwof_utils
#' @import data.table
#' @export
read_ogs <- function(of.dir, gff){
  og2 <- readLines(file.path(of.dir,
                             "Orthogroups.txt"))
  og2 <- lapply(og2, function(x) strsplit(x, " ")[[1]])
  og.name <- sapply(og2, function(x) x[1])
  og.length <- sapply(og2, length) - 1
  og.ids <- sapply(og2, function(x) x[-1])
  og.dt <- data.table(og = rep(og.name, og.length),
                      id = unlist(og.ids),
                      stringsAsFactors = F)
  og.gff <- merge(gff, og.dt, by = "id")
  return(og.gff)
}

#' @title readRename_blastGenes
#' @description
#' \code{readRename_blastGenes} readRename_blastGenes
#' @rdname pwof_utils
#' @import data.table
#' @export
readRename_blastGenes <- function(gene.dict1,
                                  gene.dict2,
                                  blast.file.in,
                                  blast.file.out,
                                  verbose = T){
  d <- fread(blast.file.in, key = "V2")
  if(verbose)
    cat("Processing", nrow(d), "raw hits ")
  d <- merge(gene.dict2, d)
  setkey(d, V1)
  d <- merge(gene.dict1, d)
  d$V1 <- NULL
  d$V2 <- NULL
  setnames(d,1:2,c("V1","V2"))

  d$V13 <- d$V12 * (-1)
  setkey(d, V13)
  d <- d[!duplicated(d[,c("V1","V2")])]
  d$V13 <- NULL
  if(verbose)
    cat(paste0("(",nrow(d)," unique) "))

  write.table(d, file = blast.file.out, sep = "\t",
              row.names = F, col.names =F, quote = F)
}


#' @title mirror_blast
#' @description
#' \code{mirror_blast} mirror_blast
#' @rdname pwof_utils
#' @import data.table
#' @export
mirror_blast <- function(blast.file1, blast.file2, verbose = T){
  d1 <- fread(blast.file1, key = "V12")
  d2 <- fread(blast.file2, key = "V12")
  d2 <- data.table(d2[,c(2,1,3:6,9:10,7:8,11:12)])
  setnames(d2, colnames(d1))
  d <- data.table(rbind(d1, d2), key = "V12")
  if(verbose)
    cat("Read in", nrow(d),"hits ... ")

  do <- d[,tail(.SD, 1), by = list(V1, V2)]
  if(verbose)
    cat("Writing", nrow(do), "unique best hits ... ")

  write.table(do, file = blast.file1, sep = "\t",
              row.names = F, col.names =F, quote = F)
  write.table(do[,c(2,1,3:6,9:10,7:8,11:12)], file = blast.file2, sep = "\t",
              row.names = F, col.names =F, quote = F)
}


#' @title make_mapDB
#' @description
#' \code{read_ogs} make_mapDB
#' @rdname pwof_utils
#' @import data.table
#' @export
make_mapDB <- function(id.db, blast.dir,cull.blast.dir){
  sm <- with(id.db, cbind(expand.grid(n.old,n.old, stringsAsFactors = F),
                          expand.grid(genome,genome, stringsAsFactors = F)))
  sm2 <- with(id.db, cbind(expand.grid(n.new,n.new, stringsAsFactors = F),
                           expand.grid(genome,genome, stringsAsFactors = F)))
  colnames(sm) <- c("n.old1","n.old2","genome1","genome2")
  colnames(sm2) <- c("n.new1","n.new2","genome1","genome2")
  sm <- merge(sm, sm2, by = c("genome1","genome2"))
  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n.old1, "_", sm$n.old2,
                                  ".txt"))
  sm$new.filename <- file.path(cull.blast.dir,
                               paste0("Blast",
                                      sm$n.new1, "_", sm$n.new2,
                                      ".txt"))
  return(sm)
}

#' @title read_geneIDs
#' @description
#' \code{read_geneIDs} read_geneIDs
#' @rdname pwof_utils
#' @import data.table
#' @export
read_geneIDs <- function(of.dir,
                         gff){
  gi <- fread(file.path(of.dir,
                        "SequenceIDs.txt"),
              sep = ":",
              stringsAsFactors = F,
              header = F,
              strip.white = T,
              col.names = c("gene.num", "id"))
  setkey(gff, id)
  setkey(gi, id)
  gi <- merge(gff, gi)
  return(gi)
}


#' @title read_speciesIDs
#' @description
#' \code{read_speciesIDs} read_speciesIDs
#' @rdname pwof_utils
#' @import data.table
#' @export
read_speciesIDs <- function(of.dir,
                            genomeIDs){
  si <- read.delim(file.path(of.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  si <- si[match(genomeIDs, si$genome),]
  rownames(si) <- si$genome
  return(si)
}

#' @title make_newOFdb
#' @description
#' \code{make_newOFdb} make_newOFdb
#' @rdname pwof_utils
#' @import data.table
#' @export
make_newOFdb <- function(tmp.dir,
                         cull.blast.dir,
                         peptide.dir,
                         genomeIDs,
                         verbose = T,
                         n.of.cores = 1){
  if (verbose)
    cat("\tCleaning tmp dir ... ")
  unlink(tmp.dir, recursive = T)
  dir.create(tmp.dir)

  unlink(cull.blast.dir, recursive = T)
  dir.create(cull.blast.dir)

  if (verbose)
    cat("Done!\n\tCopying peptide fastas ... ")
  cpd <- file.copy(from = file.path(peptide.dir,
                                    paste0(genomeIDs,".fa")),
                   to = tmp.dir)

  if (verbose)
    cat("Done!\n\tRe-making orthofinder input format ... ")
  system(paste("orthofinder", "-f", tmp.dir,
               "-a", n.of.cores, "-S diamond",
               "-op >/dev/null 2>&1"))

  if (verbose)
    cat("Done!\n\tMoving results to cull.blast directory ... ")
  blast.loc <- dirname(list.files(tmp.dir, pattern = "SequenceIDs",
                                  recursive = T, full.names = T)[1])
  fa.files <- list.files(blast.loc, pattern = "Species*",
                         full.names = T)
  fa.files <- fa.files[grep(".fa$", fa.files)]
  dmnd.files <- list.files(blast.loc, pattern = "diamondDBSpecies*",
                           full.names = T)

  sp.id.files <- file.path(blast.loc, "SpeciesIDs.txt")
  seq.id.files <- file.path(blast.loc, "SequenceIDs.txt")
  files <- c(fa.files, dmnd.files,
             sp.id.files, seq.id.files)
  nu <- file.copy(files, cull.blast.dir)
  if (verbose)
    cat("Done!\n")
}

#' @title cull_blastByScore
#' @description
#' \code{cull_blastByScore} cull_blastByScore
#' @rdname pwof_utils
#' @import data.table
#' @export
cull_blastByScore <- function(blast.file.in,
                              blast.file.out,
                              maxn,
                              verbose = T){
  d <- fread(blast.file.in, key = "V12")
  if(verbose)
    cat("Read in", nrow(d),"hits ... ")

  do <- d[,tail(.SD, maxn), by = list(V1)]
  if(verbose)
    cat("Writing", nrow(do), "hits in the top", maxn, "... ")

  write.table(do, file = blast.file.out, sep = "\t",
              row.names = F, col.names =F, quote = F)
}

#' @title merge_ofGff
#' @description
#' \code{read_ogs} merge_ofGff
#' @rdname pwof_utils
#' @import data.table
#' @export
merge_ofGff <- function(comb,
                        pw.of,
                        gff,
                        use.recip = T,
                        verbose = T){

  pw.of2 <- lapply(1:length(comb), function(i){
    x = pw.of[[i]]
    if(verbose)
      cat(paste0("\t",comb[[i]][1]),"<-->", comb[[i]][2],"... ")
    x[,n := length(unique(genome)), by = og]
    x <- subset(x, n > 1)

    x$og.id <- paste0(comb[[i]][1],"_",comb[[i]][2],"_", as.numeric(as.factor(x$og)))

    x$n <- NULL; x$order <- NULL; x$strand <- NULL; x$og <- NULL

    x1 <- data.table(x)
    x2 <- data.table(x)
    setnames(x1, paste0(colnames(x1), "1"))
    setnames(x2, paste0(colnames(x2), "2"))
    setkey(x1, id1)
    setkey(x2, id2)
    g1.ids <- x$id[x$genome == comb[[i]][1]]
    g2.ids <- x$id[x$genome == comb[[i]][2]]

    y <- read_allBlasts(gff = subset(gff, genome %in% comb[[i]]),
                        genomeIDs = comb[[i]],
                        blast.dir = dirs$blast,
                        verbose = F)
    y2 <- data.table(y[, c(2, 1, 3:6, 9:10, 7:8, 11:12)])
    setnames(y2, colnames(y))
    y0 <- rbind(y, y2)

    y <- subset(y0, (id1 %in% c(g2.ids, g1.ids) &
                      id2 %in% c(g2.ids, g1.ids)))

    y$neg.score <- y$score * (-1)
    setkey(y, neg.score)
    y <- y[!duplicated(y[,c("id1","id2")]),]
    y$neg.score <- NULL

    if(verbose)
      cat("initial hits =", nrow(y),"... ")
    setkey(y, id2)
    out <- merge(x2, y)
    setkey(out, id1)
    out <- merge(x1, out)
    out <- subset(out, og.id1 == og.id2)
    out$unique <- with(out, paste0(genome1, "_", genome2,".",comb[[i]][1], comb[[i]][2]))
    if(verbose)
      cat("found", sum(out$genome1 != out$genome2),"in orthogroups\n")
    return(out)
  })
  map <- rbindlist(pw.of2)
  map$og.id2 <- NULL
  setnames(map, "og.id1","og.id")
  map$unique.genome <- with(map, paste0(genome1, "_", genome2))
  pw.info <- map[,c("unique.genome","unique")]
  pw.info <- pw.info[!duplicated(pw.info$unique.genome),]
  pw.map2 <- map[map$unique %in% pw.info$unique,]
  return(pw.map2)
}

#' @title read_allBlasts
#' @description
#' \code{read_allBlasts} read_allBlasts
#' @rdname pwof_utils
#' @import data.table
#' @export
read_allBlasts <- function(gff,
                           genomeIDs,
                           blast.dir,
                           subset2gff = F,
                           verbose = T){
  if(verbose)
    cat("Importing original orthofinder database ... ")
  ofdat <- import_ofResults(
    gff = gff,
    genomeIDs = genomeIDs,
    blast.dir = blast.dir,
    verbose = F)
  if(verbose)
    cat("Done!\n")
  of.geneIndex <- ofdat$gene.index
  of.blastFiles <- ofdat$species.mappings
  of.blastFiles <- subset(of.blastFiles, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
  if(verbose)
    cat("Reading all blast files into memory ... ")
  all.blast <- rbindlist(lapply(of.blastFiles$filename, fread,
                                col.names = c("gn1", "gn2", "perc.iden", "align.length",
                                              "n.mismatch", "n.gapOpen", "q.start",
                                              "q.end", "s.start",
                                              "s.end", "eval", "score"),
                                key = "gn2"))
  if(verbose)
    cat("Done!\nParsing results and merging with geneIDs ... ")
  of.geneIndex1 <- data.table(of.geneIndex)
  setnames(of.geneIndex1, c("gn1","id1"))
  setkey(of.geneIndex1, "gn1")
  of.geneIndex2 <- data.table(of.geneIndex)
  setnames(of.geneIndex2, c("gn2","id2"))
  setkey(of.geneIndex2, "gn2")
  m <- merge(of.geneIndex2, all.blast)
  setkey(m, "gn1")
  m <- merge(of.geneIndex1, m)
  m$gn1 <- NULL
  m$gn2 <- NULL
  if(verbose)
    cat("Done!\n")
  return(m)
}

#' @title remake_ofInput
#' @description
#' \code{remake_ofInput} remake_ofInput
#' @rdname pwof_utils
#' @import data.table
#' @export
remake_ofInput <- function(dirs,
                           genomeIDs,
                           ploidy,
                           cull.blastByScore = T,
                           max.dup = 2,
                           verbose = T){
  # -- Step 1. Reformat peptides etc.
  if(dir.exists(dirs$cull.score.blast))
    unlink(dirs$cull.score.blast, recursive = T)
  dir.create(dirs$cull.score.blast)

  if(verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")
  make_newOFdb(tmp.dir = dirs$tmp,
               cull.blast.dir = dirs$cull.score.blast,
               peptide.dir = dirs$peptide,
               genomeIDs = genomeIDs)
  if(verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing gff annotations as a data.table ... \n")
  gff <- import_gff(gff.dir = dirs$gff,
                    genomeIDs = genomeIDs)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if(verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = dirs$blast, genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = dirs$cull.score.blast, genomeIDs = genomeIDs)
  id.db <- merge(old.ids, new.ids, by = "genome")
  setnames(id.db,2:3,c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       blast.dir = dirs$blast,
                       cull.blast.dir = dirs$cull.score.blast)
  #
  old.genes <- read_geneIDs(of.dir = dirs$blast, gff = gff)
  new.genes <- read_geneIDs(of.dir = dirs$cull.score.blast, gff = gff)
  genes <- merge(old.genes[,c("genome","id","gene.num")],
                 new.genes[,c("genome","id","gene.num")],
                 by = c("genome","id"))
  g1 <- with(genes,  data.table(V1 = gene.num.x, new1 = gene.num.y, key = "V1"))
  g2 <- with(genes,  data.table(V2 = gene.num.x, new2 = gene.num.y, key = "V2"))
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if(verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")

  d <- lapply(1:nrow(map.db), function(i){
    if(verbose)
      cat(paste0("\t",map.db$genome1[i]),"-->",map.db$genome2[i],"... ")
    bl <-readRename_blastGenes(gene.dict1 = g1,
                               gene.dict2 = g2,
                               blast.file.in = map.db$filename[i],
                               blast.file.out = map.db$new.filename[i])
    if(verbose)
      cat("Done!\n")
    return(bl)
  })

  #######################################################


  #######################################################

  #######################################################
  if(verbose)
    cat("Culling blast by score ... \n")

  map.db$score.filename <- file.path(dirs$cull.score.blast, basename(map.db$new.filename))
  bl <- lapply(1:nrow(map.db), function(i){
    s1 <- map.db$genome1[i]
    s2 <- map.db$genome2[i]
    maxn <- (ploidy[s1]*max.dup)/2
    if(verbose)
      cat(paste0("\t",s1),"-->",s2,"... ")
    blast.file.in = map.db$new.filename[i]
    blast.file.out = map.db$score.filename[i]
    cull_blastByScore(blast.file.in = blast.file.in,
                      blast.file.out = blast.file.out,
                      maxn = maxn,
                      verbose = T)
    if(verbose)
      cat("Done!\n")
  })
  if(verbose)
    cat("\tDone!\n")
  return(map.db)
}
