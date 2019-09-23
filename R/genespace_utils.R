#' @title GENESPACE utility functions
#' @description
#' \code{genespace_utils} Functions for GENESPACE
#' @name genespace_utils

#' @param map data.table, containing the merged gff and blast results
#' @param n.mappings numeric, the number of mappings required within a
#' given radius. Length must match that of radius
#' @param radius numeric, the radius to search within for syntenic
#' mappings. Length must match that of n.mappings.
#' @param blast data.table, containing the merged gff and blast results.
#' Unlike the 'map' object, which really just needs merged gff coordinates,
#' this must have all the blast8 columns. See details.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param MCScanX.tool file.path, specifying the path to the
#' MCScanX program.
#' @param silent.mcs logical, should MCScanX progress be reported?
#' @param id.db genome ID dictionary
#' @param cull.blast.dir file.path, to the subdirectory
#' where the culled blast results should be written
#' @param blast.file.in file.path, to the blast file to read in
#' @param blast.file.out file.path, to the filename to which the
#' blast results should be written
#' @param maxn integer length 1, indicating the maximum number of
#' hits to be retained per query.
#' @param of.dir file.path, to the subdirectory containing
#' blast results
#' @param gene.dict1  data.table with column names V1, new1, keyed on
#' V1. V1 contains the original geneIDs, new1 contains the new geneIDs
#' @param gene.dict2  data.table with column names V2, new2, keyed on
#' V2. V2 contains the original geneIDs, new2 contains the new geneIDs
#' @param min.score numeric length 1, specifying the minimum blast
#' bit score for a hit to be retained.
#' @param n.cores integer length 1, the number of parallel processes
#' to run.
#' @param tmp.dir file.path, to the subdirectory
#' where the temporary results should be written
#' @param peptide.dir file.path, to the subdirectory containing
#' the parsed peptide files
#' @param mcscan.dir file.path, to the subdirectory containing
#' the mcscanx output, or where mcscanx results should be written.
#' @param mcscan.param character string, of all parameters to
#' pass to MCScanX
#'
#'
#' @note \code{genespace_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespace_utils} returns its own arguments.
#'
#'

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
                         n.cores = 1){
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
               "-a", n.cores, "-S diamond",
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
                       of.dir,
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
  sm$filename <- file.path(of.dir,
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
                             of.dir,
                             verbose = T,
                             genomeIDs){
  gff <- data.table(gff)
  #######################################################
  gz <- list.files(of.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("\tDecompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(of.dir,
                           "*.gz")))
  }
  #######################################################
  if (verbose)
    cat("\tReading Species IDs\n")
  si <- read.delim(file.path(of.dir,
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
  sm$filename <- file.path(of.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = genomeIDs))
  #######################################################
  if (verbose)
    cat("\tReading gene IDs\n")
  sequence.index <- fread(file.path(of.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))
  #######################################################
  if (verbose)
    cat("\tReading orthogroup networks\n")
  og <- readLines(file.path(of.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])

  if (verbose)
    cat("\tBuilding orthogroup data.table\n")
  og2 <- readLines(file.path(of.dir,
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

#' @title run_dbs
#' @description
#' \code{run_dbs} run_dbs
#' @rdname genespace_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
run_dbs <- function(map,
                    radius,
                    n.mappings){
  nn <- frNN(
    data.frame(map[, c("rank1", "rank2"), with = F]),
    eps = radius)
  dbs <- dbscan(nn,
                minPts = n.mappings)
  map[,cluster := dbs$cluster]
  return(map)
}

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
                    MCScanX.tool,
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
                       MCScanX.tool,
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
      com <- paste(MCScanX.tool, mcscan.param,
                   file.path(mcscan.dir, "xyz"), "&> /dev/null")
    } else {
      com <- paste(MCScanX.tool, mcscan.param,
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
                       MCScanX.tool = MCScanX.tool,
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

#' @title mirror_map
#' @description
#' \code{mirror_map} mirror_map
#' @rdname genespace_utils
#' @import data.table
#' @export
mirror_map <- function(map, keycol = NULL){
  cn <- colnames(map)
  cn1 <- cn[c(grep("1$",cn),grep("2$", cn),which(!grepl("1$|2$", cn)))]
  cn2 <- cn[c(grep("2$",cn),grep("1$", cn),which(!grepl("1$|2$", cn)))]
  out1 <- data.table(map[, cn1, with = F])
  out2 <- data.table(map[, cn2, with = F])
  setnames(out2, cn1)
  out3 <- rbind(out1, out2)
  if (!is.null(keycol))
    if (keycol %in% colnames(out3))
      setkeyv(out3, keycol)

  out3 <- out3[!duplicated(out3[, c("id1", "id2")]), ]
  return(out3)
}

#' @title clean_it
#' @description
#' \code{clean_it} clean_it
#' @rdname genespace_utils
#' @import data.table
#' @export
clean_it <- function(map,
                     genomeIDs,
                     rerank,
                     radius,
                     n.mappings,
                     verbose){

  map <- data.table(map)
  setkey(map, chr1, chr2, start1, start2)

  if (rerank) {
    map[,rank1 := frank(start1,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr2)]
  }

  if (!"unique.genome" %in% colnames(map)) {
    map$unique.genome <- with(map, paste(genome1, genome2))
  }
  if (!"unique" %in% colnames(map)) {
    map$unique <- with(map, paste(genome1, genome2, chr1, chr2))
  }

  if (is.null(genomeIDs))
    genomeIDs <- unique(gff$genome)

  map <- subset(map,
                genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  map[, genome1 := factor(genome1, levels = genomeIDs)]
  map[, genome2 := factor(genome2, levels = genomeIDs)]
  setkey(map, genome1, genome2)
  map[, genome1 := as.character(genome1)]
  map[, genome2 := as.character(genome2)]

  spl.gen <- split(map, by = c("genome1", "genome2"))
  merged_map <- rbindlist(lapply(spl.gen, function(x){
    g1 = x$genome1[1]
    g2 = x$genome2[1]

    if (verbose)
      cat(paste0("\t", g1), "-->", g2,
          paste0("(initial hits = ", nrow(x) , ") ... "))
    spl.map <- split(x, by = c("chr1", "chr2"))

    chr.map <- rbindlist(lapply(spl.map, function(tmp){
      x <- run_dbs(
        map = tmp[, c("rank1", "rank2"), with = F],
        radius = radius,
        n.mappings = n.mappings)

      tmp[, block.id := x$cluster]
      return(tmp)
    }))

    chr.map <- chr.map[chr.map$block.id != 0, ]
    if (verbose)
      cat(nrow(chr.map), "hits in",
          length(unique(paste(chr.map$unique,
                              chr.map$block.id))),
          "blocks\n")
    return(chr.map)
  }))

  merged_map[,block.id := as.numeric(
    as.factor(
      paste(genome1, genome2,
            chr1, chr2, block.id)))]

  merged_blk <- make_blocks(
    map = merged_map,
    rename.blocks = T,
    rerank = rerank,
    clean.columns = F,
    ties.method = "dense")
  return(merged_blk)
}
