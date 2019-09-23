#' @title GENESPACE utility functions
#' @description
#' \code{genespace_utils} Functions for GENESPACE
#' @name genespace_utils
#'
#' @param add.gff logical, should
#' @param add.metadata logical, should
#' @param blast map-formatted data.table, without blockID info.
#' @param blast.dir file path to
#' @param blast.file.in file path to
#' @param blast.file.out file path to
#' @param blast.ids character, length
#' @param blk block-formatted data.table
#' @param blk.border numeric, length
#' @param blk.ids character, length
#' @param border.color character, length
#' @param check.ogs logical, should
#' @param chr character, length
#' @param clean.columns logical, should
#' @param comb genome combinations
#' @param cull.blast.dir file path to
#' @param dir.list list of file paths to
#' @param drop.NAs logical, should
#' @param end numeric, length
#' @param eps.radius numeric, length
#' @param fais fai-like data.table
#' @param fasta.dir file path to
#' @param gene.dict1 gene dictionary
#' @param gene.dict2 gene dictionary
#' @param geneID.abbrev.fun function to
#' @param geneid.cex numeric, length
#' @param geneid.offset numeric, length
#' @param genomeIDs character, length
#' @param genomes character, length
#' @param gff gff-like data.table
#' @param gff.dir file path to
#' @param gff.file file path to
#' @param id.db database of gene IDs / numbers
#' @param is.peptide logical, should
#' @param keep.best.id1.hit logical, should
#' @param keep.geneNum logical, should
#' @param m.param numeric, length
#' @param map map-formatted data.table
#' @param map.bychr map-formatted data.table split by chr
#' @param mappings numeric, length
#' @param maxn numeric, length
#' @param mcs.file file path to
#' @param mcscan.dir file path to
#' @param mcscan.param character, length
#' @param MCScanX.path file path to
#' @param min.blockSize numeric, length
#' @param min.dist2end numeric, length
#' @param min.perc.iden numeric, length
#' @param min.score numeric, length
#' @param n.cores numeric, length
#' @param n.mapping numeric, length
#' @param n.mappings numeric, length
#' @param n.of.cores numeric, length
#' @param n.out numeric, length
#' @param n.reps numeric, length
#' @param n.sample numeric, length
#' @param num numeric, length
#' @param of.blast blast file from orthofinder
#' @param of.dir file path to
#' @param of.ids character, length
#' @param only.orthogroups logical, should
#' @param ortho.col character, length
#' @param parse_fastaHeader.FUN function to
#' @param pattern character, length
#' @param peptide.dir file path to
#' @param prop.of.best numeric, length
#' @param pw.of logical, should
#' @param radius numeric, length
#' @param rank.buffer numeric, length
#' @param rename.blocks logical, should
#' @param rerank logical, should
#' @param return.start logical, should
#' @param silent.mcs logical, should
#' @param simplify.poly logical, should
#' @param start numeric, length
#' @param syn.blast syntenic blast, map format data.table
#' @param syn.ortho.map syntenic orthologous blast, map format data.table
#' @param ties.method character, length
#' @param tmp.dir file path to
#' @param use.rank logical, should
#' @param use.recip logical, should
#' @param use.topn logical, should
#' @param verbose logical, should
#' @param w2ki numeric vector specifying which genome comparisons to keep
#' @param which.in.blk numeric, length
#' @param whichAttr numeric, length
#' @param x numeric or character, length
#' @param y numeric, length
#'
#' @note \code{genespace_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespace_utils} returns its own arguments.
#'
#'


######################################################################
######################################################################
# Generic high-level functions
######################################################################

#' @title make_newOFdb
#' @description
#' \code{make_newOFdb} make_newOFdb
#' @rdname genespace_utils
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


#' @title readRename_blastGenes
#' @description
#' \code{readRename_blastGenes} readRename_blastGenes
#' @rdname genespace_utils
#' @import data.table
#' @export
readRename_blastGenes <- function(gene.dict1,
                                  gene.dict2,
                                  blast.file.in,
                                  blast.file.out,
                                  min.score = NULL,
                                  verbose = T){
  d <- fread(blast.file.in, key = "V2")
  if (verbose)
    cat("Processing", nrow(d), "raw hits ... ")
  if (!is.null(min.score)) {
    d <- subset(d, V12 >= min.score)
    if (verbose)
      cat(nrow(d),"with score >",min.score,"... ")
  }
  d <- merge(gene.dict2, d)
  setkey(d, V1)
  d <- merge(gene.dict1, d)
  d$V1 <- NULL
  d$V2 <- NULL
  setnames(d, 1:2, c("V1", "V2"))

  d$V13 <- d$V12 * (-1)
  setkey(d, V13)
  d <- d[!duplicated(d[, c("V1", "V2")])]
  d$V13 <- NULL
  if (verbose)
    cat(paste0("(", nrow(d), " unique) "))

  write.table(d,
              file = blast.file.out,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}


######################################################################
######################################################################
# read_allBlasts functions
######################################################################

#' @title read_speciesIDs
#' @description
#' \code{read_speciesIDs} read_speciesIDs
#' @rdname genespace_utils
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

#' @title make_mapDB
#' @description
#' \code{make_mapDB} make_mapDB
#' @rdname genespace_utils
#' @import data.table
#' @export
make_mapDB <- function(id.db,
                       blast.dir,
                       cull.blast.dir){
  sm <- with(id.db, cbind(expand.grid(n.old, n.old,
                                      stringsAsFactors = F),
                          expand.grid(genome, genome,
                                      stringsAsFactors = F)))
  sm2 <- with(id.db, cbind(expand.grid(n.new, n.new,
                                       stringsAsFactors = F),
                           expand.grid(genome, genome,
                                       stringsAsFactors = F)))
  colnames(sm) <- c("n.old1", "n.old2", "genome1", "genome2")
  colnames(sm2) <- c("n.new1", "n.new2", "genome1", "genome2")
  sm <- merge(sm, sm2, by = c("genome1", "genome2"))
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
#' @rdname genespace_utils
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

#' @title import_ofResults
#' @description
#' \code{import_ofResults} import_ofResults
#' @rdname genespace_utils
#' @import data.table
#' @export
import_ofResults <- function(gff,
                             blast.dir,
                             verbose = T,
                             genomeIDs){
  gff <- data.table(gff)
  #######################################################
  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("\tDecompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(blast.dir,
                           "*.gz")))
  }
  #######################################################
  if (verbose)
    cat("\tReading Species IDs\n")
  si <- read.delim(file.path(blast.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  rownames(si) <- si$genome
  #######################################################
  sm <- cbind(expand.grid(si$genome.num,
                          si$genome.num,
                          stringsAsFactors = F),
              expand.grid(si$genome,
                          si$genome,
                          stringsAsFactors = F))
  colnames(sm) <- c("n1","n2","genome1","genome2")
  for (i in rev(genomeIDs))
    sm$ref[sm$genome1 == i | sm$genome2 == i] <- i
  sm$alt <- ifelse(sm$genome1 == sm$ref, sm$genome2, sm$genome1)
  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = genomeIDs))
  #######################################################
  if (verbose)
    cat("\tReading gene IDs\n")
  sequence.index <- fread(file.path(blast.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))
  #######################################################
  if (verbose)
    cat("\tReading orthogroup networks\n")
  og <- readLines(file.path(blast.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])

  if (verbose)
    cat("\tBuilding orthogroup data.table\n")
  og2 <- readLines(file.path(blast.dir,
                             "Orthogroups.txt"))
  og2 <- lapply(og2, function(x) strsplit(x, " ")[[1]])
  og.name <- sapply(og2, function(x) x[1])
  og.length <- sapply(og2, length) - 1
  og.ids <- sapply(og2, function(x) x[-1])
  og2 <- data.table(block.id = NA,
                    og = rep(og.name, og.length),
                    id = unlist(og.ids),
                    stringsAsFactors = F)
  #######################################################
  if (verbose)
    cat("\tCompiling metadata\n")

  gffi <- gff[, c("id","genome")]
  setkey(gffi, "id")
  setkey(og2, "id")
  ogo <- merge(gffi, og2)
  setkey(ogo, "block.id", "genome", "og", "id")
  ogo[, og.n.genes := length(unique(id)), by = list(block.id, og)]
  ogo[, og.n.genomes := length(unique(genome)), by = list(block.id, og)]

  ogo.meta = ogo[!duplicated(ogo[, -c(1:2), with = F]),-c(1:2), with = F]
  #######################################################
  return(list(orthogroups = og,
              species.mappings = sm,
              species.index = si,
              orthogroup.metadata = ogo.meta,
              orthogroup.data.table = ogo,
              gene.index = sequence.index))
}

######################################################################
######################################################################
# remake_ofInput functions
######################################################################

#' @title cull_blastByScore
#' @description
#' \code{cull_blastByScore} cull_blastByScore
#' @rdname genespace_utils
#' @import data.table
#' @export
cull_blastByScore <- function(blast.file.in,
                              blast.file.out,
                              maxn,
                              verbose = T){
  d <- fread(blast.file.in, key = "V12")
  if (verbose)
    cat("Read in", nrow(d), "hits ... ")
  do <- d[,tail(.SD, maxn), by = list(V1)]
  if (verbose)
    cat(nrow(do),
        "top", maxn, "hits ... ")
  write.table(do,
              file = blast.file.out,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}

######################################################################
######################################################################
# build_synBlocks functions
######################################################################
#' @title read_ogs
#' @description
#' \code{read_ogs} fread_ogs
#' @rdname genespace_utils
#' @import data.table
#' @export
read_ogs <- function(of.dir,
                     gff){

  ortho.loc <- dirname(list.files(of.dir,
                                  pattern = "Orthogroups.txt",
                                  recursive = T,
                                  full.names = T)[1])
  og2 <- readLines(file.path(ortho.loc,
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

#' @title run_of
#' @description
#' \code{run_of} run_of
#' @rdname genespace_utils
#' @export
run_of <- function(tmp.dir,
                   blast.dir,
                   genomeIDs,
                   n.cores,
                   gff){

  if (dir.exists(tmp.dir))
    un <- unlink(tmp.dir, recursive = T)
  dir.create(tmp.dir)

  fs <- list.files(blast.dir, full.names = T)
  fc <- file.copy(from = fs, to = tmp.dir)

  com <- paste("orthofinder", "-b", tmp.dir,
               "-a", n.cores,
               "-og")
  system(com)

  og <- read_ogs(of.dir = tmp.dir, gff = gff)
  og1 <- data.table(id1 = og$id, og1 = og$og)
  og2 <- data.table(id2 = og$id, og2 = og$og)

  y <- read_allBlasts(gff = gff,
                      genomeIDs = genomeIDs,
                      blast.dir = tmp.dir,
                      verbose = T,
                      add.gff = T)

  out <- merge(merge(y, og1, by = "id1"), og2, by = "id2")
  return(out)
}

#' @title build_globalBlocks
#' @description
#' \code{build_globalBlocks} build_globalBlocks
#' @rdname genespace_utils
#' @export
build_globalBlocks <- function(dir.list,
                               gff,
                               genomeIDs,
                               MCScanX.path,
                               n.cores,
                               use.topn,
                               min.blockSize,
                               m.param,
                               only.orthogroups){
  tmp.dir <- dir.list$tmp
  mcscan.dir <- dir.list$mcscanx

  if (use.topn) {
    blast.dir <- dir.list$cull.score.blast
  }else{
    blast.dir <- dir.locs$cull.blast
  }

  bl <- run_of(tmp.dir = tmp.dir,
               blast.dir = blast.dir,
               n.cores = n.cores,
               gff = gff,
               genomeIDs = genomeIDs)

  if (only.orthogroups) {
    bl <- subset(bl, og1 == og2)
  }

  bl[, og.id := og1]
  bl$og1 <- NULL
  bl$og2 <- NULL

  MCScanX.tool <- file.path(MCScanX.path,"MCScanX")
  mcsp <- paste("-a -s", min.blockSize,
                "-m", m.param)

  geni <- t(combn(genomeIDs, 2))
  geni <- rbind(geni, t(sapply(genomeIDs, function(x) c(x,x))))
  geni <- data.table(geni)
  setnames(geni, c("genome1", "genome2"))
  blo <- merge(geni, bl, by = c("genome1", "genome2"))
  syn.blks <- pipe_mcscanx(blast = blo,
                           gff = gff,
                           genomeIDs = genomeIDs,
                           MCScanX.path = MCScanX.tool,
                           mcscan.dir = mcscan.dir,
                           mcscan.param = mcsp,
                           verbose = T)
  syn.blks$map$og.id <- syn.blks$map$og1


  return(list(blast = blo,
              map = syn.blks$map,
              blk = syn.blks$block))
}

######################################################################
######################################################################
# pipe_mcscanx functions
######################################################################

#' @title run_mcs
#' @description
#' \code{run_mcs} run_mcs
#' @rdname genespace_utils
#' @import data.table
#' @export
run_mcs <- function(blast,
                    gff,
                    genomeIDs,
                    mcscan.dir,
                    MCScanX.path,
                    mcscan.param,
                    silent.mcs){
  ######################################################################
  ######################################################################
  parse_mcs <- function(mcs.file){
    if (length(readLines(paste0(mcs.file, ".collinearity"))) < 12) {
      return(NULL)
    }else{
      mcscan.raw <- read.delim(paste0(mcs.file, ".collinearity"),
                               sep = "\t", header = F,
                               comment.char = "#", strip.white = T,
                               stringsAsFactors = F)

      fac <- as.numeric(as.factor(sapply(as.character(mcscan.raw$V1), function(x)
        strsplit(x, "-")[[1]][1])))
      genes <- with(mcscan.raw, c(V2, V3))
      out <- data.table(id1 = mcscan.raw$V2,
                        id2 = mcscan.raw$V3,
                        block.id = fac)
      out <- data.table(out)
      setkey(out, id1, id2)
      return(out)
    }
  }
  ######################################################################
  ######################################################################

  ######################################################################
  ######################################################################
  prep_mcs <- function(blast,
                       gff,
                       genomeIDs,
                       mcscan.dir,
                       mcscan.param,
                       MCScanX.path,
                       silent.mcs){

    if (length(unique(gff$genome)) == 1) {
      ga <- gff
      gb <- gff
      ga$genome <- "a"
      gb$genome <- "b"
      gff <- rbind(gb, ga)
      gff$genome <- factor(gff$genome,
                           levels = c("a", "b"))
    } else {
      gff$genome <- factor(gff$genome,
                           levels = genomeIDs)
    }

    setkey(gff, genome)

    gff[,chr.num := frank(chr, ties.method = "dense"),
        by = genome]
    gff$rank.start = frank(gff, genome, chr, start, ties.method = "dense")
    gff$rank.end = frank(gff, genome, chr, end, ties.method = "dense")
    gff[,genome.num := frank(genome, ties.method = "dense")]

    lets <- paste0(letters,letters)[1:length(unique(gff$genome))]
    gff$genome.abbrev <- paste0(lets[gff$genome.num],1)

    gff.in <- gff[,c("genome.abbrev", "id", "rank.start", "rank.end")]


    blast.in <- blast[,c("id1", "id2", "perc.iden", "align.length",
                         "n.mismatch", "n.gapOpen", "q.start", "q.end",
                         "s.start", "s.end", "eval", "score")]
    write.table(gff.in,
                file = file.path(mcscan.dir, "xyz.gff"),
                row.names = F,
                col.names = F,
                quote = F, sep = "\t")
    write.table(blast.in,
                file = file.path(mcscan.dir, "xyz.blast"),
                row.names = F,
                col.names = F,
                quote = F, sep = "\t")
    if (silent.mcs) {
      com <- paste(MCScanX.path, mcscan.param,
                   file.path(mcscan.dir, "xyz"), "&> /dev/null")
    } else {
      com <- paste(MCScanX.path, mcscan.param,
                   file.path(mcscan.dir, "xyz"))
    }

    system(com)
    return(file.path(mcscan.dir, "xyz"))
  }
  ######################################################################
  ######################################################################

  mcs.file <- prep_mcs(blast = blast,
                       gff = gff,
                       genomeIDs,
                       mcscan.dir = mcscan.dir,
                       mcscan.param = mcscan.param,
                       MCScanX.path = MCScanX.path,
                       silent.mcs = silent.mcs)
  mcs.parsed <- parse_mcs(mcs.file)
  if (!is.null(mcs.parsed)) {
    blast <- data.table(blast)
    setkey(blast, id1, id2)
    setkey(mcs.parsed, id1, id2)
    map.out <- data.table(merge(mcs.parsed, blast))
    setkey(map.out, block.id, chr1, start1)
    return(map.out)
  }else{
    return(blast[0,])
  }
}


######################################################################
######################################################################
# clean_blocks functions
######################################################################

#' @title simplify_map
#' @description
#' \code{simplify_map} simplify_map
#' @rdname genespace_utils
#' @import data.table
#' @export
simplify_map <- function(map,
                         gff,
                         genomeIDs,
                         mirror = F,
                         just.rank = F){
  tpg <- data.table(t(combn(genomeIDs, 2)))
  setnames(tpg, c("genome1","genome2"))
  tpg <- rbind(tpg,
               data.table(genome1 = genomeIDs,
                          genome2 = genomeIDs))
  map <- merge(map, tpg, by = c("genome1","genome2"))
  map[,block.id := ifelse(is.na(block.id), NA,
                          paste0("blk_",
                                 as.numeric(
                                   as.factor(
                                     paste(genome1, genome2,
                                           chr1, chr2,
                                           block.id)))))]

  gff <- subset(gff, genome %in% genomeIDs)
  spl.map <- split(map, by = c("genome1","genome2","chr1"))
  spl.gff <- split(gff, by = c("genome","chr"))
  mo1 <- rbindlist(lapply(spl.map, function(x){
    gc <- paste(x$genome1[1], x$chr1[1], sep = ".")
    y <- subset(spl.gff[[gc]], id %in% x$id1)
    x <- x[,c("genome1","genome2","chr1","chr2","id1","id2","block.id")]
    y[,rank := frankv(y, cols = c("start","end"), ties.method = "random")]
    o <- merge(x,
               with(y,
                    data.table(id1 = id,
                               start1 = start,
                               end1 = end,
                               rank1 = rank)),
               by = "id1")

  }))

  spl.map <- split(mo1, by = c("genome1","genome2","chr2"))
  map <- rbindlist(lapply(spl.map, function(x){
    gc <- paste(x$genome2[1], x$chr2[1], sep = ".")
    y <- subset(spl.gff[[gc]], id %in% x$id2)
    y[,rank := frankv(y, cols = c("start","end"), ties.method = "random")]
    o <- merge(x,
               with(y,
                    data.table(id2 = id,
                               start2 = start,
                               end2 = end,
                               rank2 = rank)),
               by = "id2")

  }))
  setcolorder(map, c(3:6,2,1,7:8,11,9,12,10,13))
  setkey(map, chr1, chr2, start1, start2)

  if (mirror) {
    m2 <- data.table(map)
    setcolorder(map, c("block.id","genome1","genome2","chr1","chr2",
                      "id1","id2","start1","start2",
                      "end1","end2","rank1","rank2"))
    setcolorder(m2, c("block.id","genome2","genome1","chr2","chr1",
                      "id2","id1","start2","start1",
                      "end2","end1","rank2","rank1"))
    setnames(m2, colnames(map))
    map <- rbind(map, m2)
    map <- map[!duplicated(map),]
  }

  map[,block.id := ifelse(is.na(block.id), NA,
                          paste0("blk_",
                                 as.numeric(
                                   as.factor(
                                     paste(genome1, genome2,
                                           chr1, chr2,
                                           block.id)))))]

  setkey(map, genome1, genome2, chr1, chr2, rank1, rank2)
  if(!just.rank){
    blk <- map[,list(start1 = min(start1),
                     end1 = max(end1),
                     start2 = min(start2),
                     end2 = max(end2),
                     rank.start1 = min(rank1),
                     rank.end1 = max(rank1),
                     rank.start2 = min(rank2),
                     rank.end2 = max(rank2),
                     nhits1 = length(unique(id1)),
                     nhits2 = length(unique(id2)),
                     nhits = length(unique(c(id1, id2)))),
               by = list(genome1, genome2, chr1, chr2, block.id)]

    return(list(map = map,
                blk = blk))
  }else{
    return(map[,c("genome1","genome2","chr1","chr2","id1","id2","rank1","rank2","block.id")])
  }
}

#' @title spl_ovlGap
#' @description
#' \code{spl_ovlGap} spl_ovlGap
#' @rdname genespace_utils
#' @import data.table
#' @export
spl_ovlGap <- function(map){
  map <- data.table(map)
  setkey(map, rank1)
  rl1 <- rle(map$block.id)
  map[,block.id := paste(block.id, rep(1:length(rl1$lengths), rl1$lengths), sep = ".")]
  setkey(map, rank2)
  rl2 <- rle(map$block.id)
  map[,block.id := paste(block.id, rep(1:length(rl2$lengths), rl2$lengths), sep = ".")]
  return(map)
}

#' @title any_dupInBlk
#' @description
#' \code{any_dupInBlk} any_dupInBlk
#' @rdname genespace_utils
#' @import data.table
#' @export
any_dupInBlk <- function(map, block.id1, block.id2){
  m1 <- subset(map, block.id == block.id1)
  m2 <- subset(map, block.id == block.id2)
  g1 <- unique(c(m1$id1, m1$id2))
  g2 <- unique(c(m2$id1, m2$id2))
  return(any(duplicated(c(g1, g2))))
}


#' @title chk_olv
#' @description
#' \code{chk_olv} chk_olv
#' @rdname genespace_utils
#' @import data.table
#' @export
chk_olv <- function(blk, type){
  b1 <- blk[,c("genome1","genome2","chr1","chr2","start1","end1","block.id")]
  b2 <- blk[,c("genome1","genome2","chr1","chr2","start2","end2","block.id")]
  setnames(b1, c("start1","end1"), c("start","end"))
  setnames(b2, c("start2","end2"), c("start","end"))

  spl1 <- split(b1, by = c("genome1","genome2","chr1","chr2"))
  spl2 <- split(b2, by = c("genome1","genome2","chr1","chr2"))

  ovl.dt <- rbindlist(lapply(names(spl1), function(i){
    x1 <- spl1[[i]]
    x2 <- spl2[[i]]
    setkey(x1, start, end)
    setkey(x2, start, end)
    ovl1 <- foverlaps(x1, x1, type = type, which=TRUE)
    ovl2 <- foverlaps(x2, x2, type = type, which=TRUE)
    ovl1[,inside := x1$block.id[xid]]
    ovl2[,inside := x2$block.id[xid]]
    ovl1[,outside := x1$block.id[yid]]
    ovl2[,outside := x2$block.id[yid]]
    ovl1 <- subset(ovl1, inside != outside)
    ovl2 <- subset(ovl2, inside != outside)
    if(nrow(ovl1) > 0 & nrow(ovl2) > 0){
      out <- rbind(ovl1[,c("inside","outside")],
                   ovl2[,c("inside","outside")])
      out <- out[duplicated(out),]
      if(nrow(out) > 0){
        return(out)
      }
    }
  }))
  if(nrow(ovl.dt) > 0){
    ovl.dt[,nhits.inside := blk$nhits[match(inside, blk$block.id)]]
    ovl.dt[,nhits.outside := blk$nhits[match(outside, blk$block.id)]]
  }
  return(ovl.dt)
}

#' @title merge_dupOvlBlks
#' @description
#' \code{merge_dupOvlBlks} merge_dupOvlBlks
#' @rdname genespace_utils
#' @import data.table
#' @export
merge_dupOvlBlks <- function(map, gff, genomeIDs,  max.iter = 10, verbose = T){


  nod <- 1
  iter <- 0
  while(nod > 0 & iter < max.iter){
    iter = iter + 1
    if(verbose)
      cat("\tIter",iter,"... ")
    simp <- simplify_map(map = map,
                         gff = gff,
                         genomeIDs = genomeIDs,
                         mirror = F)
    map <- data.table(simp$map)
    blk <- data.table(simp$blk)
    ovl.tot <- chk_olv(blk, type = "any")
    mtmp <- subset(map, block.id %in% unique(unlist(ovl.tot[,c("inside","outside")])))
    ovl.tot[,any.dup := sapply(1:nrow(ovl.tot), function(i)
      any_dupInBlk(map = mtmp,
                   block.id1 =fz inside[i],
                   block.id2 = outside[i]))]
    blk2merge <- subset(ovl.tot, any.dup)[,c("inside","outside")]
    nod <- nrow(blk2merge)
    if(nod == 0){
      if(verbose)
        cat("\tNo additional duplicate blocks ... Done!\n")
    }else{
      if(verbose)
        cat("\tFound", nod,"blocks with duplicate hits to be merged\n")
    }
    all.blkid <- data.table(block.id = unique(map$block.id))

    blk2merge2 <- rbind(data.table(block.id = blk2merge$inside,
                                   bid2 = blk2merge$outside),
                        data.table(block.id = all.blkid$block.id,
                                   bid2 = all.blkid$block.id))
    blk2merge2[,num1 := as.numeric(as.factor(block.id))]
    blk2merge2[,num2 := as.numeric(as.factor(bid2))]
    blk2merge2 <- subset(blk2merge2, num1 <= num2)
    blk2merge2 <- blk2merge2[!duplicated(blk2merge2),]
    blk2merge3 <- with(blk2merge2,
                       data.table(block.id = c(block.id, bid2),
                                  new.id = c(block.id, block.id)))
    setkey(blk2merge3, block.id, new.id)
    map <- merge(map,
                 blk2merge3[!duplicated(blk2merge3$block.id),],
                 by = "block.id")
    map[,block.id := new.id]
    map[,new.id := NULL]
  }
  return(simplify_map(map,
                      gff = gff,
                      genomeIDs = genomeIDs,
                      mirror = F))
}

#' @title run_dbs
#' @description
#' \code{run_dbs} run_dbs
#' @rdname genespace_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
run_dbs <- function(y,
                    eps.radius,
                    mappings){
  nn <- frNN(data.frame(y[, c("rank1", "rank2"), with = F]),
             eps = eps.radius)
  dbs <- dbscan(nn,
                minPts = mappings)
  y$cluster <- dbs$cluster
  return(y)
}

######################################################################
######################################################################
# cull_syntenicBlast functions
######################################################################

#' @title cull_blast2MapChr
#' @description
#' \code{cull_blast2MapChr} cull_blast2MapChr
#' @rdname genespace_utils
#' @import data.table
#' @export
cull_blast2MapChr <- function(map,
                              blast){
  map[, uniq := paste(genome1, genome2, chr1, chr2)]
  u <- unique(map$uniq)
  blast[, uniq := paste(genome1, genome2, chr1, chr2)]

  out <- subset(blast, uniq %in% u)
  return(out)
}

#' @title cull_blast2NewIds
#' @description
#' \code{cull_blast2NewIds} cull_blast2NewIds
#' @rdname genespace_utils
#' @import data.table
#' @export
cull_blast2NewIds <- function(blast,
                              map){
  setkey(map, id1, id2)
  setkey(blast, id1, id2)
  t.map <- map[,c("id1", "id2")]
  t.map$in.map <- TRUE

  blast <- merge(blast, t.map, all.x = T)
  blast <- blast[is.na(blast$in.map),colnames(blast), with = F]
  return(blast)
}

#' @title find_whichInBuffer
#' @description
#' \code{find_whichInBuffer} find_whichInBuffer
#' @rdname genespace_utils
#' @import data.table
#' @importFrom dbscan frNN
#' @export
find_whichInBuffer <- function(x,
                               y,
                               which.in.blk,
                               rank.buffer){

  nn <- frNN(x = data.frame(x, y),
             eps = rank.buffer)

  all.near.blk <- unique(unlist(nn$id[which.in.blk]))

  return(all.near.blk[order(all.near.blk)])
}

#' @title find_hitsInBuffer
#' @description
#' \code{find_hitsInBuffer} find_hitsInBuffer
#' @rdname genespace_utils
#' @import data.table
#' @export
find_hitsInBuffer <- function(map,
                              blast,
                              rank.buffer,
                              verbose){
  spl.map <- split(map, by = "unique.genome")
  spl.blast <- split(blast, by = "unique.genome")
  ns <- unique(names(spl.blast))
  ns <- ns[ns %in% unique(names(spl.map))]

  res.by.genome <- lapply(ns, function(i){
    i.map = spl.map[[i]]
    i.blast = spl.blast[[i]]
    if (verbose)
      cat(paste0("\t", i, " ... (new.hits = ",
                 nrow(i.blast), ", ", "map.size = ",
                 nrow(i.map), ")"))



    spl.i.map <- split(i.map, by = "unique.chr")
    spl.i.blast <- split(i.blast, by = "unique.chr")
    nis <- unique(names(spl.i.blast))
    nis <- nis[nis %in% unique(names(spl.i.map))]

    ids2keep <- rbindlist(lapply(nis, function(j){
      j.map <- spl.i.map[[j]]
      j.blast <- spl.i.blast[[j]]
      j.map.spl <- split(j.map, by = "block.id")

      ko <- lapply(names(j.map.spl), function(k){
        k.map <- j.map.spl[[k]]
        j.blast <- j.blast[with(j.blast,
                                rank1 >= (min(k.map$rank1) - (rank.buffer * 2)) &
                                  rank2 >= (min(k.map$rank2) - (rank.buffer * 2)) &
                                  rank1 <= (max(k.map$rank1) + (rank.buffer * 2)) &
                                  rank2 <= (max(k.map$rank2) + (rank.buffer * 2))), ]
        j.out <- rbind(k.map[ , c("id1", "id2")],
                       j.blast[ , c("id1", "id2")])
        x <- c(k.map$rank1, j.blast$rank1)
        y <- c(k.map$rank2, j.blast$rank2)
        wh <- 1:nrow(k.map)

        tokeep <- find_whichInBuffer(x = x,
                                     y = y,
                                     which.in.blk = wh,
                                     rank.buffer = rank.buffer)
        j.out$block.id <- k
        return(j.out[tokeep,])
      })

      return(rbindlist(ko))
    }))

    ids2keep <- rbind(ids2keep, i.map[,c("id1","id2","block.id")])
    ids2keep <- ids2keep[!duplicated(ids2keep),]

    if (verbose)
      cat(" returning", nrow(ids2keep), "\n")
    return(ids2keep)
  })
  out <- rbindlist(res.by.genome)
  setkey(out, id1, id2)
  return(out)
}

#' @title rerank_fromIDs
#' @description
#' \code{rerank_fromIDs} rerank_fromIDs
#' @rdname genespace_utils
#' @import data.table
#' @export
rerank_fromIDs <- function(map,
                           gff,
                           cull2genes.inmap = FALSE){

  m <- with(map, data.table(id = c(id1, id2), genome = c(genome1, genome2)))
  m <- m[!duplicated(m),]
  if(cull2genes.inmap){
    g <- merge(gff, m, by = c("id", "genome"))
  }else{
    g <- data.table(gff)
  }

  g[,rank := frank(start, ties.method = "random"),
     by = list(genome, chr)]

  if("rank1" %in% colnames(map))
    map[,rank1 := NULL]
  if("rank2" %in% colnames(map))
    map[,rank2 := NULL]

  g1 <- data.table(g[,c("genome","id","rank")])
  g2 <- data.table(g[,c("genome","id","rank")])
  setnames(g1, paste0(colnames(g1), "1"))
  setnames(g2, paste0(colnames(g2), "2"))
  out <- merge(g1,
               merge(g2, map, by = c("genome2","id2")),
               by = c("genome1","id1"))

  return(out)
}

######################################################################
######################################################################
# pull_orthonet functions
######################################################################

######################################################################
######################################################################
# extend_blocks functions
######################################################################

#' @title pull_synHomos
#' @description
#' \code{pull_synHomos} pull_synHomos
#' @rdname genespace_utils
#' @export
pull_synHomos <- function(syn.blast,
                          syn.ortho.map,
                          min.score,
                          prop.of.best,
                          min.perc.iden){

  if (is.na(min.score)) {
    sb <- subset(syn.blast, perc.iden >= min.perc.iden)

  }else{
    sb <- subset(syn.blast, score >= min.score)
  }

  sb[,propscore1 := score/max(score),
     by = list(genome2, chr2, id1)]
  sb[,propscore2 := score/max(score),
     by = list(genome1, chr1, id2)]
  sb <- subset(sb,
               propscore1 > prop.of.best |
                 propscore2 > prop.of.best)
  mu <- with(syn.ortho.map, unique(paste(id1, id2)))
  so <- subset(sb, !paste(id1, id2) %in% mu)
  return(so)
}

######################################################################
######################################################################
#  functions
######################################################################

#' @title cull_blast2blk
#' @description
#' \code{cull_blast2blk} cull_blast2blk
#' @rdname genespace_utils
#' @export
cull_blast2blk <- function(blast,
                           blk,
                           blk.rank.buffer){
  spl.blast <- split(blast, by = c("genome1","genome2","chr1","chr2"))
  spl.blk <- split(blk, by = c("genome1","genome2","chr1","chr2"))
  blast <- rbindlist(lapply(names(spl.blk), function(i){
    y <- spl.blk[[i]]
    x <- spl.blast[[i]]
    setkey(y, n.mapping)
    xc <- rbindlist(lapply(1:nrow(y), function(j){
      int1 <- findInterval(x$rank1,
                           c(y$rankstart1[j] - blk.rank.buffer,
                             y$rankend1[j] + blk.rank.buffer)) == 1
      int2 <- findInterval(x$rank2,
                           c(y$rankstart2[j] - blk.rank.buffer,
                             y$rankend2[j] + blk.rank.buffer)) == 1
      wh <- which(int1 & int2)
      if (length(wh) < 1) {
        return(NULL)
      }else{
        if ("block.id" %in% colnames(x))
          x[, block.id := NULL]

        return(data.table(x[wh,],
                          block.id = y$block.id[j]))
      }
    }))
    return(xc)
  }))
  return(blast)
}

#' @title cull_blast2blk
#' @description
#' \code{cull_blast2blk} cull_blast2blk
#' @rdname genespace_utils
#' @export
cull_synBlast <- function(map,
                          blast,
                          blast.in,
                          rank.buffer,
                          verbose){

  if("rank1" %in% colnames(map))
    map$rank1 <- NULL
  if("rank1" %in% colnames(map))
    map$rank2 <- NULL
  map$rank2 <- NULL
  map$what <- "map"
  blast$what <- "blast"

  spl.map <- split(map, by = "block.id")
  spl.blast <- split(blast, by = "block.id")

  blks <- table(map$block.id)
  blks <- names(blks)[order(blks)]
  idl <- rbindlist(lapply(1:length(blks), function(j){
    i = blks[j]
    if (verbose)
      if (j %% 100 == 0)
        cat("\tCompleted",j,"/",length(spl.map),"blocks\n")

    y <- spl.map[[i]]
    cn <- colnames(y)
    x <- spl.blast[[i]][,cn, with = F]
    z <- rbind(y,x)
    z <- z[!duplicated(z[,c("id1","id2")]),]

    wh <- find_whichInBuffer(x = frank(z$start1, ties.method = "dense"),
                             y = frank(z$start2, ties.method = "dense"),
                             which.in.blk = which(z$what == "map"),
                             rank.buffer = rank.buffer)
    return(data.table(id1 = z$id1[wh],
                      id2 = z$id2[wh],
                      block.id = i))
  }))
  idl <- idl[!duplicated(idl[,c("id1","id2")]),]

  blast.out <- merge(idl, blast.in, by = c("id1","id2"))
  return(blast.out)
}

#' @title cull_blast2blk
#' @description
#' \code{cull_blast2blk} cull_blast2blk
#' @rdname genespace_utils
#' @export
io_blast4of <- function(of.dir, gff, blast, genomeIDs){
  gi1 <- read_geneIDs(of.dir = of.dir, gff = gff)[,c("id","gene.num")]
  gi2 <- data.table(gi1)
  setnames(gi1, c("id1","gn1"))
  setnames(gi2, c("id2","gn2"))
  blast$gn2 <- blast$gn1 <- NULL
  blast <- merge(gi1, merge(gi2, blast, by = "id2"), by = "id1")

  #######################################################
  si <- read_speciesIDs(of.dir = of.dir, genomeIDs = genomeIDs)
  combs <- expand.grid(si$genome, si$genome)
  combs.n <- expand.grid(si$genome.num, si$genome.num)
  combs.file <- file.path(of.dir, paste0("Blast",combs.n[,1],"_", combs.n[,2],".txt"))
  cols2write <- c("gn1","gn2","perc.iden","align.length","n.mismatch","n.gapOpen",
                  "q.start","q.end","s.start","s.end","eval","score")
  for (i in 1:nrow(combs)) {
    tmp <- subset(blast,
                  genome1 == combs[i,1] &
                    genome2 == combs[i,2])[,cols2write, with = F]
    write.table(tmp, sep = "\t",
                file = combs.file[i],
                quote = F,
                col.names = F,
                row.names = F)
  }
}

#' @title proc_MCScanBlocks
#' @description
#' \code{proc_MCScanBlocks} proc_MCScanBlocks
#' @rdname genespace_utils
#' @export

#' @title extend_mcsBlks
#' @description
#' \code{extend_mcsBlks} extend_mcsBlks
#' @rdname genespace_utils
#' @export
extend_mcsBlks <- function(map,
                           dir.list,
                           gff,
                           genomeIDs,
                           radius = 100,
                           min.blockSize,
                           m.param,
                           verbose = T){
  MCScanX.tool <- file.path(MCScanX.path,"MCScanX")

  if (verbose)
    cat("Loading all blast files into memory ... ")
  blast <- read_allBlasts(gff = gff,
                          keep.geneNum = F,
                          add.gff = T,
                          check.ogs = F,
                          blast.dir = dir.list$cull.score.blast,
                          genomeIDs = genomeIDs,
                          verbose = F)
  if(verbose)
    cat("Done!\nMerging gff, blast and map data.tables ... ")
  # blast <- data.table(synblks$blast)
  map <- simplify_map(map,
                      gff = gff,
                      genomeIDs = genomeIDs,
                      mirror = F,
                      just.rank = F)$map

  bl <- merge(map[,c("block.id","genome1","genome2","id1","id2")],
              blast,
              all = T,
              by = c("genome1","genome2","id1","id2"))
  bl <- subset(bl, !is.na(chr1))
  blo <- simplify_map(map = bl,
                      gff = gff,
                      genomeIDs = genomeIDs,
                      mirror = F,
                      just.rank = T)
  newblk <- subset(blo, !is.na(block.id))
  newblk <- newblk[,list(rankstart1 = min(rank1, na.rm = T),
                         rankend1 = max(rank1, na.rm = T),
                         rankstart2 = min(rank2, na.rm = T),
                         rankend2 = max(rank2, na.rm = T),
                         n.mapping = .N),
                   by = list(genome1, genome2, chr1, chr2, block.id)]
  spl.blk <- split(newblk, by = c("genome1","genome2","chr1","chr2"))
  spl.blast <- split(blo, by = c("genome1","genome2","chr1","chr2"))

  if(verbose)
    cat("Done!\nCulling blast to syntenic hits ... ")

  blast.out <- rbindlist(lapply(names(spl.blk), function(i){
    y <- spl.blk[[i]]
    x <- spl.blast[[i]]
    setkey(y, n.mapping)
    xc <- rbindlist(lapply(1:nrow(y), function(j){

      int1 <- findInterval(x$rank1,
                           c(y$rankstart1[j] - radius,
                             y$rankend1[j] + radius)) == 1
      int2 <- findInterval(x$rank2,
                           c(y$rankstart2[j] - radius,
                             y$rankend2[j] + radius)) == 1
      wh <- which(int1 & int2)
      if (length(wh) < 1) {
        return(NULL)
      }else{

        tmp <- x[wh,]
        wh2 <- find_whichInBuffer(x = tmp$rank1,
                                  y = tmp$rank2,
                                  which.in.blk = which(!is.na(tmp$block.id)),
                                  rank.buffer = radius)
        out <- tmp[wh2,]
        if ("block.id" %in% colnames(out))
          out[, block.id := NULL]

        return(data.table(out,
                          block.id = y$block.id[j]))
      }
    }))
    return(xc)
  }))
  if(verbose)
    cat("Done!\nRe-forming syntenic blocks with MCScanX ... \n")
  blast.out <- blast.out[!duplicated(blast.out[,c("id1","id2")]),]
  bl.cull <- merge(blast, blast.out[,c("id1","id2")], by = c("id1","id2"))
  mcsp <- paste("-a -s", min.blockSize,
                "-m", m.param)
  tpg <- data.table(t(combn(genomeIDs, 2)))
  setnames(tpg, c("genome1","genome2"))
  tpg <- rbind(tpg,
               data.table(genome1 = genomeIDs,
                          genome2 = genomeIDs))
  bl.cull <- merge(bl.cull, tpg, by = c("genome1","genome2"))
  syn.out <- pipe_mcscanx(blast = bl.cull,
                          gff = gff,
                          dir.list = dir.list,
                          genomeIDs = genomeIDs,
                          MCScanX.path = MCScanX.tool,
                          mcscan.dir = dir.list$mcscanx,
                          mcscan.param = mcsp,
                          verbose = T)
  if(verbose)
    cat("\tDone!\nRe-formatting map and block data.tables ... ")
  simp.out <- simplify_map(map = syn.out$map, gff = gff, genomeIDs = genomeIDs, mirror = F)
  if(verbose)
    cat("Done!\n")
  return(simp.out)
}

######################################################################
######################################################################
#  functions
######################################################################

#' @title cull_blast2blk
#' @description
#' \code{cull_blast2blk} cull_blast2blk
#' @rdname genespace_utils
#' @export
find_dupBlks <- function(map, verbose){

  if(verbose)
    cat("Checking genome1 --> genome2 ... ")
  xb <- map[!duplicated(map[,c("id1","block.id")]),]
  dup.inblk <- subset(xb, id1 %in% xb$id1[duplicated(xb[,c("id1","genome2")])])
  db <- dup.inblk[,list(blks = list(unique(block.id))),
                  by = list(genome1, genome2, chr1, id1)]

  db <- db$blks[sapply(db$blks, length)>1]
  db <- sapply(db[!duplicated(db)], function(x) combn(x,2))
  db <- t(do.call(cbind, db))
  db1 <- db[!duplicated(db),]

  if(verbose)
    cat("Done!\nChecking genome2 --> genome1 ... ")
  xb <- map[!duplicated(map[,c("id2","block.id")]),]
  dup.inblk <- subset(xb, id2 %in% xb$id2[duplicated(xb[,c("id2","genome1")])])
  db <- dup.inblk[,list(blks = list(unique(block.id))),
                  by = list(genome1, genome2, chr2, id2)]

  db <- db$blks[sapply(db$blks, length)>1]
  db <- sapply(db[!duplicated(db)], function(x) combn(x,2))
  db <- t(do.call(cbind, db))
  db2 <- db[!duplicated(db),]

  out <- rbind(db1, db2, db1[,2:1], db2[,2:1])
  out <- out[!duplicated(out),]
  if(verbose)
    cat("Done!\nFound",nrow(out),"overlapping dupliate blocks\n")
  return(out)
}

######################################################################
######################################################################
#  functions
######################################################################

