#' @title selection utility functions
#' @description
#' \code{genespace_utils} Functions that allow selection stat calculation in GENESPACE
#' @name genespace_utils
#'
#' @param add.gff logical, should
#' @param add.metadata logical, should
#' @param alpha numeric, length
#' @param bg.col character, length
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
#' @param chr.abbrev.fun function to
#' @param chr.bg.cex numeric, length
#' @param chr.bg.col character, length
#' @param chr.bg.pch numeric, length
#' @param chr.buff numeric, length
#' @param chr.buffer numeric, length
#' @param chr.id.cex numeric, length
#' @param chr.id.col character, length
#' @param chr.lab.buff numeric, length
#' @param chr.list list of chromosomes to plot
#' @param chr.segm.col character, length
#' @param chr.segm.lwd numeric, length
#' @param clean.columns logical, should
#' @param col character, length
#' @param cols character, length
#' @param comb genome combinations
#' @param cull.blast.dir file path to
#' @param dir.list list of file paths to
#' @param do.cumulative logical, should
#' @param dodge.geneIDs numeric, length
#' @param dodge.x numeric, length
#' @param drop.NAs logical, should
#' @param e1 numeric, length
#' @param e2 numeric, length
#' @param end numeric, length
#' @param eps.radius numeric, length
#' @param fais fai-like data.table
#' @param fasta.dir file path to
#' @param fill.color character, length
#' @param forCircos logical, should
#' @param gap.prop numeric, length
#' @param gene.colors character, length
#' @param gene.dict1 gene dictionary
#' @param gene.dict2 gene dictionary
#' @param geneID.abbrev.fun function to
#' @param geneid.cex numeric, length
#' @param geneid.offset numeric, length
#' @param genes2plot character, length
#' @param genomeIDs character, length
#' @param genomes character, length
#' @param gff gff-like data.table
#' @param gff.dir file path to
#' @param gff.file file path to
#' @param id.db database of gene IDs / numbers
#' @param is.peptide logical, should
#' @param keep.best.id1.hit logical, should
#' @param keep.geneNum logical, should
#' @param lab.chr logical, should
#' @param lab.chr.1only logical, should
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
#' @param points.per.curve numeric, length
#' @param prop.of.best numeric, length
#' @param pw.of logical, should
#' @param radius numeric, length
#' @param rank.buffer numeric, length
#' @param rename.blocks logical, should
#' @param rerank logical, should
#' @param return.start logical, should
#' @param s1 numeric, length
#' @param s2 numeric, length
#' @param scale.it logical, should
#' @param scale2dodge numeric, length
#' @param silent.mcs logical, should
#' @param simplify.poly logical, should
#' @param start numeric, length
#' @param str2drop character, length
#' @param str2parse character, length
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
#' @param y.end numeric, length
#' @param y.start numeric, length
#'
#' @note \code{genespace_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespace_utils} returns its own arguments.
#'
#'
#' @title Parse fasta headers
#' @description
#' \code{parse_fastaHeader} Rename fasta header to match the gff.
#' @rdname genespace_utils
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
#' @export
parse_fastaHeader <- function(fasta.dir,
                              is.peptide = T,
                              pattern = "fa",
                              verbose = T,
                              parse_fastaHeader.FUN){

  files <- list.files(fasta.dir,
                      pattern = pattern,
                      full.names = T)

  if (verbose)
    cat("Renaming fasta headers ...\n")
  ss <- lapply(files, function(i){
    if (verbose)
      cat("...", i, "\n\t")
    if (is.peptide) {
      x <- readAAStringSet(i)
    }else{
      x <- readDNAStringSet(i)
    }
    if (verbose)
      cat("original names (e.g.):",
          names(x)[1])
    names(x) <- sapply(names(x), parse_fastaHeader.FUN)
    if (verbose)
      cat("\n\tparsed names (e.g.):",
          names(x)[1],"\n")
    writeXStringSet(x, filepath = i)
  })
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

#' @title parse_gff
#' @description
#' \code{parse_gff} parse_gff
#' @rdname genespace_utils
#' @import data.table
#' @export
parse_gff <- function(gff.file,
                      str2drop = "Name=",
                      str2parse = ";",
                      use = "gene",
                      whichAttr = 2){
  g <- suppressWarnings(
    fread(gff.file,
          showProgress = F,
          verbose = F))
  g <- g[g$V3 == use, c(9, 1, 4, 5, 7)]
  g$V9 <- sapply(g$V9, function(x)
    gsub(str2drop, "",
         strsplit(x, str2parse)[[1]][whichAttr]))
  setnames(g, c("id", "chr", "start", "end", "strand"))
  return(g)
}


#' @title import_gff
#' @description
#' \code{import_gff} import_gff
#' @rdname genespace_utils
#' @import data.table
#' @export
import_gff <- function(gff.dir,
                       genomeIDs,
                       use = "gene",
                       verbose = T,
                       str2drop = "Name=",
                       str2parse = ";",
                       whichAttr = 2){
  gff.files <- file.path(gff.dir,
                         paste0(genomeIDs, ".gff3"))
  names(gff.files) <- genomeIDs
  #######################################################
  gff <- rbindlist(lapply(names(gff.files), function(i){
    if (verbose)
      cat("\tReading",i,"... ")
    tmp <- parse_gff(gff.file = gff.files[[i]],
                     str2drop = str2drop,
                     str2parse = str2parse,
                     use = use,
                     whichAttr = whichAttr)

    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    if (verbose)
      cat("n. genes =", nrow(tmp),"\n")
    return(tmp)
  }))
  #######################################################
  setkey(gff, "genome", "id")
  return(gff)
}

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
                                  verbose = T){
  d <- fread(blast.file.in, key = "V2")
  if (verbose)
    cat("Processing", nrow(d), "raw hits ... ")
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
    cat("Writing", nrow(do),
        "hits in the top", maxn, "... ")

  write.table(do,
              file = blast.file.out,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}

#' @title read_ogs
#' @description
#' \code{read_ogs} fread_ogs
#' @rdname genespace_utils
#' @import data.table
#' @export
read_ogs <- function(of.dir,
                     gff){
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

#' @title make_blocks
#' @description
#' \code{make_blocks} make_blocks
#' @rdname genespace_utils
#' @import data.table
#' @export
make_blocks <- function(map,
                        rerank = T,
                        drop.NAs = F,
                        rename.blocks = F,
                        add.metadata = F,
                        clean.columns = T,
                        ties.method = "dense"){
  if (clean.columns) {
    cols2keep <- c("block.id", "orthogroup",
                   "genome1", "genome2",
                   "id1", "id2",
                   "og1", "og2",
                   "chr1", "start1", "end1", "strand1", "order1",
                   "chr2", "start2", "end2", "strand2", "order2",
                   "rank1", "rank2")
    cols2keep <- cols2keep[cols2keep %in% colnames(map)]
    map <- map[,cols2keep, with = F]
  }
  map <- data.table(map)
  if (rename.blocks) {
    map$block.id <- as.numeric(as.factor(with(map, paste(genome1, genome2, block.id))))
  }

  setkey(map, chr1, chr2, start1, start2)
  if (rerank) {
    map[,rank1 := frank(start1,
                        ties.method = ties.method),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = ties.method),
        by = list(genome1, genome2, chr2)]
  }
  if (drop.NAs) {
    map <- map[complete.cases(map),]
  }

  if (!add.metadata) {
    out.blk <- map[,list(chr1 = chr1[1],
                         chr2 = chr2[1],
                         start1 = min(start1),
                         start2 = min(start2),
                         end1 = max(end1),
                         end2 = max(end2),
                         rankstart1 = min(rank1),
                         rankstart2 = min(rank2),
                         rankend1 = max(rank1),
                         rankend2 = max(rank2),
                         n.mapping = length(start1),
                         orient = ifelse(length(start1) <= 1, "+",
                                         ifelse(cor(jitter(start1),
                                                    jitter(start2)) > 0,"+", "-"))),
                   by = list(block.id, genome1, genome2)]
  }else{
    out.blk <- map[,list(chr1 = chr1[1],
                         chr2 = chr2[1],
                         start1 = min(start1),
                         start2 = min(start2),
                         end1 = max(end1),
                         end2 = max(end2),
                         rankstart1 = min(rank1),
                         rankstart2 = min(rank2),
                         rankend1 = max(rank1),
                         rankend2 = max(rank2),
                         n.mapping = length(start1),
                         n.unique.map1 = length(unique(id1)),
                         n.unique.map2 = length(unique(id2)),
                         orient = ifelse(length(start1) <= 1, "+",
                                         ifelse(cor(jitter(start1),
                                                    jitter(start2)) > 0,"+", "-"))),
                   by = list(block.id, genome1, genome2)]
    out.blk$width1 <- with(out.blk, rankend1 - rankstart1) + 1
    out.blk$width2 <- with(out.blk, rankend2 - rankstart2) + 1
    out.blk$prop.map1 <- with(out.blk, n.unique.map1 / width1)
    out.blk$prop.map2 <- with(out.blk, n.unique.map2 / width2)
    out.blk$prop.ratio1 <- with(out.blk, prop.map1 / prop.map2)
    out.blk$prop.ratio2 <- with(out.blk, prop.map2 / prop.map1)
  }


  map <- data.table(map,
                    stringsAsFactors = F)
  blk <- data.table(out.blk,
                    stringsAsFactors = F)

  return(list(block = blk,
              map = map))
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
                     by.blk = F,
                     verbose){

  map <- data.table(map)
  setkey(map, chr1, chr2, start1, start2)
  if(by.blk & "block.id" %in% colnames(map)){
    if (rerank) {
      map[,rank1 := frank(start1,
                          ties.method = "dense"),
          by = list(genome1, genome2, chr1, block.id)]
      map[,rank2 := frank(start2,
                          ties.method = "dense"),
          by = list(genome1, genome2, chr2, block.id)]
    }

    map$unique.genome <- with(map, paste(genome1, genome2, block.id))
    map$unique <- with(map, paste(genome1, genome2, chr1, chr2, block.id))

  }else{
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
  }


  map[, genome1 := factor(genome1, levels = genomeIDs)]
  map[, genome2 := factor(genome2, levels = genomeIDs)]
  setkey(map, genome1, genome2)
  map[, genome1 := as.character(genome1)]
  map[, genome2 := as.character(genome2)]

  spl.gen <- split(map, by = "unique.genome")
  merged_map <- rbindlist(lapply(spl.gen, function(x){
    g1 = x$genome1[1]
    g2 = x$genome2[1]
    if (verbose)
      cat(paste0("\t",g1),"-->",g2,paste0("(initial hits = ",nrow(x),") ... "))
    spl.map <- split(x, by = "unique")
    chr.map <- rbindlist(lapply(spl.map, function(tmp){
      x <- run_dbs(y = tmp[, c("rank1", "rank2"), with = F],
                   eps.radius = radius,
                   mappings = n.mappings)
      tmp$block.id <- x$cluster
      return(tmp)
    }))
    chr.map <- chr.map[chr.map$block.id != 0, ]
    if (verbose)
      cat(nrow(chr.map), "hits in",
          length(unique(paste(chr.map$unique, chr.map$block.id))), "blocks\n")
    return(chr.map)
  }))

  merged_map$block.id <- with(merged_map,
                              as.numeric(as.factor(paste(unique, block.id))))

  merged_blk <- make_blocks(map = merged_map,
                            rename.blocks = T,
                            rerank = T,
                            clean.columns = F,
                            ties.method = "dense")
  return(merged_blk)
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

#' @title read_allBlasts
#' @description
#' \code{read_allBlasts} read_allBlasts
#' @rdname genespace_utils
#' @import data.table
#' @export
read_allBlasts <- function(gff,
                           genomeIDs,
                           blast.dir,
                           check.ogs = T,
                           add.gff = F,
                           keep.geneNum = F,
                           verbose = T){
  if (check.ogs) {
    if (verbose)
      cat("Importing orthofinder database ... ")
    ofdat <- import_ofResults(
      gff = gff,
      genomeIDs = genomeIDs,
      blast.dir = blast.dir,
      verbose = F)
    if (verbose)
      cat("Done!\n")
    of.geneIndex <- ofdat$gene.index
    of.blastFiles <- ofdat$species.mappings
    of.blastFiles <- subset(of.blastFiles, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
  }else{
    of.blastFiles <- list(filename = list.files(blast.dir, pattern = "^Blast", full.names = T))
  }

  if (verbose)
    cat("Reading all blast files into memory ... ")
  all.blast <- rbindlist(lapply(of.blastFiles$filename, fread,
                                col.names = c("gn1", "gn2", "perc.iden", "align.length",
                                              "n.mismatch", "n.gapOpen", "q.start",
                                              "q.end", "s.start",
                                              "s.end", "eval", "score"),
                                key = "gn2"))
  of.geneIndex <- read_geneIDs(of.dir = blast.dir,
                               gff = gff)
  setnames(of.geneIndex, "gene.num","gn")
  if (verbose)
    cat("Done!\nParsing results and merging with geneIDs ... ")
  if (!add.gff){
    of.geneIndex1 <- data.table(of.geneIndex)
    setnames(of.geneIndex1, c("gn1","id1"))
    setkey(of.geneIndex1, "gn1")
    of.geneIndex2 <- data.table(of.geneIndex)
    setnames(of.geneIndex2, c("gn2","id2"))
    setkey(of.geneIndex2, "gn2")
    m <- merge(of.geneIndex2, all.blast)
    setkey(m, "gn1")
    m <- merge(of.geneIndex1, m)

    if (!keep.geneNum) {
      m$gn1 <- NULL
      m$gn2 <- NULL
    }
  }else{
    if (add.gff) {
      gff1 <- data.table(of.geneIndex)
      gff2 <- data.table(of.geneIndex)
      setnames(gff1, paste0(colnames(gff1), "1"))
      setnames(gff2, paste0(colnames(gff2), "2"))
      setkey(gff2, gn2)
      setkey(all.blast, gn2)
      m <- merge(gff2, all.blast)
      setkey(gff1, gn1)
      setkey(m, gn1)
      m <- merge(gff1, m)
      if (!keep.geneNum) {
        m$gn1 <- NULL
        m$gn2 <- NULL
      }
    }
  }

  if (verbose)
    cat("Done!\n")
  return(m)
}


#' @title merge_ofGff
#' @description
#' \code{read_ogs} merge_ofGff
#' @rdname genespace_utils
#' @import data.table
#' @export
merge_ofGff <- function(comb,
                        pw.of,
                        gff,
                        dir.list,
                        use.recip = T,
                        verbose = T){

  pw.of2 <- lapply(1:length(comb), function(i){
    x = pw.of[[i]]
    if (verbose)
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
                        blast.dir = dir.list$blast,
                        verbose = F, check.ogs = F, add.gff = T,keep.geneNum = F)
    y2 <- data.table(y[, c(2, 1, 3:6, 9:10, 7:8, 11:12)])
    setnames(y2, colnames(y))
    y0 <- rbind(y, y2)

    y <- subset(y0, (id1 %in% c(g2.ids, g1.ids) &
                       id2 %in% c(g2.ids, g1.ids)))

    y$neg.score <- y$score * (-1)
    setkey(y, neg.score)
    y <- y[!duplicated(y[,c("id1","id2")]),]
    y$neg.score <- NULL

    if (verbose)
      cat("initial hits =", nrow(y),"... ")
    setkey(y, id2)
    out <- merge(x2, y)
    setkey(out, id1)
    out <- merge(x1, out)
    out <- subset(out, og.id1 == og.id2)
    out$unique <- with(out, paste0(genome1, "_", genome2,".",comb[[i]][1], comb[[i]][2]))
    if (verbose)
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

#' @title rerun_pairwiseOF
#' @description
#' \code{rerun_pairwiseOF} rerun_pairwiseOF
#' @rdname genespace_utils
#' @import data.table
#' @export
rerun_pairwiseOF <- function(dir.list,
                             gff,
                             genomeIDs,
                             n.cores = 6,
                             verbose = T){

  gff <- subset(gff, genome %in% genomeIDs)
  # -- Step 1. Reformat peptides etc.
  if(verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")

  pw.dir <- file.path(dirname(dir.list$tmp), "pw")
  if(file.exists(pw.dir))
    unlink(pw.dir, recursive = T)
  dir.create(pw.dir)

  on.exit(expr =  unlink(pw.dir, recursive = T))

  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = pw.dir,
               peptide.dir = dir.list$peptide,
               genomeIDs = genomeIDs,
               verbose = verbose)
  if (verbose)
    cat("\tDone!\n")
  #######################################################
  if (verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = dir.list$cull.score.blast, genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = pw.dir, genomeIDs = genomeIDs)
  id.db <- merge(old.ids, new.ids, by = "genome")
  setnames(id.db,2:3,c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       blast.dir = dir.list$cull.score.blast,
                       cull.blast.dir = pw.dir)
  #
  old.genes <- read_geneIDs(of.dir = dir.list$cull.score.blast, gff = gff)
  new.genes <- read_geneIDs(of.dir = pw.dir, gff = gff)
  genes <- merge(old.genes[,c("genome","id","gene.num")],
                 new.genes[,c("genome","id","gene.num")],
                 by = c("genome","id"))
  g1 <- with(genes,  data.table(V1 = gene.num.x, new1 = gene.num.y, key = "V1"))
  g2 <- with(genes,  data.table(V2 = gene.num.x, new2 = gene.num.y, key = "V2"))
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")
  in.blast.files <- map.db$filename
  out.blast.files <- map.db$new.filename

  d <- lapply(1:nrow(map.db), function(i){
    if (verbose)
      cat(paste0("\t",map.db$genome1[i]),"-->",map.db$genome2[i],"... ")
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
                          blast.dir = pw.dir,
                          genomeIDs = genomeIDs,
                          verbose = F)
  return(blast)
}

#' @title prep_mcs
#' @description
#' \code{prep_mcs} prep_mcs
#' @rdname genespace_utils
#' @import data.table
#' @export
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

#' @title parse_mcs
#' @description
#' \code{parse_mcs} parse_mcs
#' @rdname genespace_utils
#' @import data.table
#' @export
parse_mcs <- function(mcs.file){
  if(length(readLines(paste0(mcs.file, ".collinearity"))) < 12){
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

  mcs.file <- prep_mcs(blast = blast,
                       gff = gff,
                       genomeIDs,
                       mcscan.dir = mcscan.dir,
                       mcscan.param = mcscan.param,
                       MCScanX.path = MCScanX.path,
                       silent.mcs = silent.mcs)
  mcs.parsed <- parse_mcs(mcs.file)
  if(!is.null(mcs.parsed)){
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

#' @title pipe_mcscanx
#' @description
#' \code{pipe_mcscanx} pipe_mcscanx
#' @rdname genespace_utils
#' @import data.table
#' @export
pipe_mcscanx <- function(blast,
                         gff,
                         genomeIDs,
                         mcscan.dir,
                         mcscan.param,
                         MCScanX.path,
                         split.by.blk = F,
                         silent.mcs = T,
                         verbose = T){
  #######################################################
  gff.tmp <- gff
  gff.tmp$genome <- paste0(gff.tmp$genome,"xxxx")
  gff.tmp$id <- paste0(gff.tmp$id,"xxxx")
  gff <- data.table(rbind(gff, gff.tmp))

  genomeIDs <- c(genomeIDs, paste0(genomeIDs, "xxxx"))

  bl.dif <- blast[blast$genome1 != blast$genome2,]
  bl.same <- blast[blast$genome1 == blast$genome2,]
  bl.same$id2 <- paste0(bl.same$id2, "xxxx")
  bl.same$genome2 <- paste0(bl.same$genome2, "xxxx")
  blast <- rbind(bl.dif, bl.same)
  #######################################################
  blast <- blast[!duplicated(blast),]
  if(split.by.blk){
    blast$unique <- with(blast, paste(genome1, genome2, block.id))
  }else{
    blast$unique <- with(blast, paste(genome1, genome2))
  }

  if ("block.id" %in% colnames(blast))
    blast$block.id <- NULL

  spl <- split(blast, by = "unique")
  #######################################################
  out <- rbindlist(lapply(spl, function(x){
    genomes <- c(x$genome1[1], x$genome2[1])

    if (verbose)
      cat(paste0("\t", gsub("xxxx", "", genomes[1])),
          "-->", gsub("xxxx", "", genomes[2]),
          paste0("(initial hits = ", nrow(x), ")"))

    gff.x <- gff[gff$genome %in% genomes,]

    tmp <- run_mcs(blast = x,
                   gff = gff.x,
                   genomeIDs = genomeIDs,
                   MCScanX.path = MCScanX.path,
                   mcscan.dir = mcscan.dir,
                   mcscan.param = mcscan.param,
                   silent.mcs = silent.mcs)
    if(nrow(tmp) == 0){
      if (verbose)
        cat(" culled hits = 0!\n")
      return(tmp)
    }else{

      tmp$block.id <- with(tmp, paste0(unique, block.id))

      if (verbose)
        cat(" culled hits =", nrow(tmp), "\n")

      return(tmp)
    }
  }), fill = T)
  #######################################################
  out$block.id <- as.numeric(as.factor(out$block.id))
  out$genome2 <- gsub("xxxx", "",  out$genome2)
  out$id2 <- gsub("xxxx", "",  out$id2)
  #######################################################
  return(make_blocks(out, clean.columns = F))
}


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

#' @title cull_syntenicBlast
#' @description
#' \code{cull_syntenicBlast} cull_syntenicBlast
#' @rdname genespace_utils
#' @import data.table
#' @export
cull_syntenicBlast <- function(map,
                               blast,
                               gff,
                               rank.buffer,
                               verbose = T){
  if (verbose)
    cat("Dropping chromosome combinations in blast not found in map... ")
  blast.in <- cull_blast2MapChr(blast = blast,
                                map = map)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Dropping blast hits in the map already ... ")
  t.blast <- cull_blast2NewIds(blast = blast.in,
                               map = map)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Making new blast and map with gff-based ranks ... ")
  id.names <- c("genome1", "genome2",
                "id1", "id2",
                "chr1", "chr2",
                "start1", "start2",
                "end1", "end2")
  all.ids <- rbind(map[, id.names, with = F],
                   t.blast[, id.names, with = F])
  all.ids <- all.ids[!duplicated(all.ids), ]
  all.ids[ ,rank1 := frank(start1, ties.method = "dense"),
           by = list(genome1, chr1)]
  all.ids[ ,rank2 := frank(start2, ties.method = "dense"),
           by = list(genome2, chr2)]
  setkey(all.ids, id1, id2)
  setkey(map, id1, id2)
  setkey(blast, id1, id2)
  r.map <- merge(all.ids, map[ ,c("id1", "id2")])
  r.blast <- merge(all.ids, t.blast[ ,c("id1", "id2")])

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Finding neighbors to mappings in blast hits ... completed:\n")
  r.map[, unique.genome := paste(genome1, genome2)]
  r.blast[, unique.genome := paste(genome1, genome2)]
  r.map[, unique.chr := paste(unique.genome, chr1, chr2)]
  r.blast[, unique.chr := paste(unique.genome, chr1, chr2)]

  r.map <- merge(r.map, map[,c("id1","id2","block.id")], by = c("id1","id2"))

  ids2keep <- find_hitsInBuffer(map = r.map,
                                blast = r.blast,
                                rank.buffer = rank.buffer,
                                verbose = verbose)

  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reformatting blast output ... ")

  setkey(blast, id1, id2)

  out.blast <- merge(ids2keep, blast)
  re.out <- rerank_fromIDs(map = out.blast,
                           gff = gff)

  if (verbose)
    #######################################################
  cat("Done!\n")
  return(re.out)
}


#' @title rerank_fromIDs
#' @description
#' \code{rerank_fromIDs} rerank_fromIDs
#' @rdname genespace_utils
#' @import data.table
#' @export
rerank_fromIDs <- function(map,
                           gff){

  m.gff1 <- map[, c("genome1", "id1", "chr1", "start1","end1")]
  m.gff2 <- map[, c("genome2", "id2", "chr2", "start2","end2")]
  gff <- gff[, c("genome", "id", "chr", "start", "end")]
  setnames(m.gff1, colnames(gff))
  setnames(m.gff2, colnames(gff))
  g <- rbind(gff, m.gff1, m.gff2)
  g <- g[!duplicated(g), ]

  g[, unique.chr := paste(genome, chr)]
  spl <- split(g, by = "unique.chr")
  g <- rbindlist(lapply(spl, function(x){
    x$rank <- frank(x, start, end, id, ties.method = "dense")
    return(x)
  }))
  g$unique.chr <- NULL

  g1 <- data.table(g)
  g2 <- data.table(g)

  setnames(g1, paste0(colnames(g1), "1"))
  setnames(g2, paste0(colnames(g2), "2"))
  setkeyv(g1, cols = colnames(g1)[-6])
  setkeyv(g2, cols = colnames(g2)[-6])

  setkeyv(map, cols = colnames(g2)[-6])
  if ("rank1" %in% colnames(map))
    map$rank1 <- NULL
  if ("rank2" %in% colnames(map))
    map$rank2 <- NULL
  if ("rank" %in% colnames(map))
    map$rank <- NULL
  if ("n" %in% colnames(map))
    map$n <- NULL

  m1 <- merge(g2, map)
  setkeyv(m1, cols = colnames(g1)[-6])
  out <- merge(g1, m1)

  return(out)
}


#' @title extract_seqIDs
#' @description
#' \code{extract_seqIDs} extract_seqIDs
#' @rdname genespace_utils
#' @import data.table
#' @export
extract_seqIDs <- function(blast,
                           genomeIDs){
  ug <- with(blast,
             data.table(genome = c(genome1, genome2),
                        chr = c(chr1,chr2),
                        pos = c(start1, start2),
                        id = c(id1, id2),
                        stringsAsFactors = F))
  ug$genome <- factor(ug$genome, levels = genomeIDs)
  setkey(ug, genome, chr, pos)
  ug <- ug[!duplicated(ug$id),]
  ug <- ug[ug$genome %in% genomeIDs,]
  ug <- ug[with(ug, order(genome, chr, pos)),]
  ugenes <- data.table(rbindlist(lapply(split(ug, "genome"), function(x){
    x$num <- 0:(nrow(x) - 1)
    return(x)
  })))

  ugenes$genome.num <- as.numeric(factor(ugenes$genome, levels = genomeIDs)) - 1
  ugenome <- ugenes[!duplicated(ugenes[,c("genome","genome.num")]),c("genome","genome.num")]

  seqs <- with(ugenes, paste0(genome.num, "_", num, ": ", id))
  spes <- with(ugenome, paste0(genome.num, ": ", genome, ".fa"))

  dmnd.ids <- with(ugenome, paste0("diamondDBSpecies", genome.num))
  names(dmnd.ids) <- ugenome$genome

  fa.ids <- with(ugenome, paste0("Species", genome.num, ".fa"))
  names(fa.ids) <- ugenome$genome

  ugenome <- data.frame(ugenome)
  rownames(ugenome) <- ugenome$genome
  eg <- as.matrix(expand.grid(genomeIDs, genomeIDs))
  eg <- data.table(genome1 = eg[,1],
                   genome2 = eg[,2],
                   genome.num1 = ugenome[eg[,1], "genome.num"],
                   genome.num2 = ugenome[eg[,2], "genome.num"])
  return(list(species.ids = spes,
              sequence.ids = seqs,
              dmnd.ids = dmnd.ids,
              fa.ids = fa.ids,
              blast.ids = data.frame(eg),
              gene.dict = with(ugenes,
                               data.table(id = id,
                                          num = paste0(genome.num, "_", num))),
              id.list = with(ugenes, split(id, genome)),
              num.list = with(ugenes, split(paste0(genome.num, "_", num), genome))))
}

#' @title write_ofBlast
#' @description
#' \code{write_ofBlast} write_ofBlast
#' @rdname genespace_utils
#' @import data.table
#' @export
write_ofBlast <- function(blast,
                          of.ids,
                          blast.ids,
                          of.dir){
  bl <- blast[,c("id1", "id2",
                 "perc.iden", "align.length",
                 "n.mismatch","n.gapOpen",
                 "q.start", "q.end",
                 "s.start", "s.end",
                 "eval", "score",
                 "genome1","genome2")]
  bl1 <- data.table(bl)
  bl2 <- data.table(bl[,c(2,1,3:6,8,7,10,9,11:12,14,13)])
  setnames(bl2, colnames(bl1))
  bl <- rbind(bl1, bl2)
  bl$negscore <- bl$score * (-1)
  setkey(bl, negscore)
  bl <- bl[!duplicated(bl[,c("id1", "id2")]),]
  bl$negscore <- NULL

  d1 <- with(of.ids$gene.dict, data.table(id1 = id, num1 = num))
  d2 <- with(of.ids$gene.dict, data.table(id2 = id, num2 = num))
  setkey(bl, id2)
  setkey(d2, id2)
  setkey(d1, id1)
  blt <-  merge(d2, bl)
  setkey(blt, id1)
  bl <- merge(d1, blt)[,c(2, 4, 5:16)]

  eg <- data.frame(blast.ids)
  for (i in 1:nrow(eg)) {
    x1 = eg[i, 1]
    x2 = eg[i, 2]
    n1 = eg[i, 3]
    n2 = eg[i, 4]
    tmp.bl <- bl[with(bl, genome1 == x1 & genome2 == x2), -c(13:14)]
    write.table(tmp.bl,
                file = file.path(of.dir,
                                 paste0("Blast", n1, "_", n2, ".txt")),
                quote = F, col.names = F,
                row.names = F, sep = "\t")
  }
}

#' @title write_ofData
#' @description
#' \code{write_ofData} write_ofData
#' @rdname genespace_utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet writeXStringSet
#' @export
write_ofData <- function(blast,
                         genomeIDs,
                         of.dir,
                         peptide.dir,
                         verbose = T){
  if (verbose)
    cat("Extracting sequence and species IDs ... ")
  of.ids <- extract_seqIDs(blast = blast,
                           genomeIDs = genomeIDs)
  cat(of.ids$species.ids,
      sep = "\n",
      file = file.path(of.dir, "SpeciesIDs.txt"))
  cat(of.ids$sequence.ids,
      sep = "\n",
      file = file.path(of.dir, "SequenceIDs.txt"))
  if (verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if (verbose)
    cat("Parsing and databasing peptide sequences ... ")
  pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readAAStringSet(file.path(peptide.dir,
                              paste0(i, ".fa")))))
  pep.spl <- sapply(names(of.ids$id.list), USE.NAMES = T, simplify = F, function(i){
    x <- pep.fastas[of.ids$id.list[[i]]]
    names(x) <- of.ids$num.list[[i]]
    return(x)
  })

  pep.files <- sapply(names(pep.spl), function(i){
    outf <- file.path(of.dir, of.ids$fa.ids[i])
    outdb <- file.path(of.dir, of.ids$dmnd.ids[i])
    writeXStringSet(pep.spl[[i]], filepath = outf)
    system(paste("diamond makedb --quiet",
                 "--in", outf,
                 "-d", outdb))
    return(list(fa = outf,
                db = outdb))
  })
  if (verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if (verbose)
    cat("Converting blast results to orthofinder-formatted text files ... ")
  write_ofBlast(blast = blast,
                of.ids = of.ids,
                blast.ids = of.ids$blast.ids,
                of.dir = of.dir)
  if (verbose)
    cat("Done!\n")
  ########################################################
}

#' @title format_gffChrlist
#' @description
#' \code{format_gffChrlist} rename and simplify gff
#' @rdname genespace_utils
#' @export
format_gffChrlist <- function(gff,
                              genomes,
                              chr.list,
                              use.rank,
                              gap.prop,
                              do.cumulative = T){
  gff <- rbindlist(lapply(1:length(genomes), function(i){
    tmp <- subset(gff, genome == genomes[i] &
                    chr %in% chr.list[[i]])
    tmp$genome <- names(chr.list)[i]
    return(tmp)
  }))

  if (use.rank) {
    gff[,start := frank(start, ties.method = "random"),
        by = list(genome, chr)]
    gff[,end := frank(end, ties.method = "random"),
        by = list(genome, chr)]
  }
  gff$chr <- as.character(gff$chr)
  if (do.cumulative) {
    gff <- convert_gff2coords(gff = gff,
                              chr.list,
                              gap.prop = gap.prop)
  }else{
    gff$start.p <- gff$start
    gff$end.p <- gff$end
  }

  return(gff)
}

#' @title format_mapChrlist
#' @description
#' \code{format_mapChrlist} rename and simplify map
#' @rdname genespace_utils
#' @export
format_mapChrlist <- function(map,
                              gff,
                              chr.list,
                              genomes,
                              forCircos = F,
                              do.cumulative = T){
  map <- rbindlist(lapply(1:(length(genomes) - 1), function(i){
    tmp <- subset(map, genome1 == genomes[i] &
                    chr1 %in% chr.list[[i]] &
                    genome2 == genomes[i + 1] &
                    chr2 %in% chr.list[[i + 1]])
    tmp$genome1 <- names(chr.list)[i]
    tmp$genome2 <- names(chr.list)[i + 1]
    return(tmp)
  }))

  gi <- genomes
  genomes <- as.character(names(chr.list))

  gf1 <- with(gff,
              data.table(id1 = id,
                         chr1 = chr,
                         pos1 = start.p))
  gf2 <- with(gff,
              data.table(id2 = id,
                         chr2 = chr,
                         pos2 = start.p))

  map <- map[!duplicated(map[,c("block.id","og.id",
                                "genome1","genome2",
                                "id1","id2")]),]
  gf1 <- gf1[!duplicated(gf1), ]
  gf2 <- gf2[!duplicated(gf2), ]

  map <- merge(gf2,
               map[,c("block.id","og.id",
                      "genome1","genome2",
                      "id1","id2")],
               by = c("id2"))
  map <- merge(gf1,
               map, by = c("id1"))

  if (forCircos) {
    combns <- data.table(t(combn(genomes,2)))
    setnames(combns, c("genome1", "genome2"))
    map <- merge(combns,
                 map,
                 by = c("genome1", "genome2"))

  }else{
    combs <- data.table(y1 = 0:(length(genomes) - 2),
                        y2 = 1:(length(genomes) - 1),
                        genome1 = genomes[1:(length(genomes) - 1)],
                        genome2 = genomes[2:length(genomes)],
                        stringsAsFactors = F)
    map <- rbind(map,
                 with(map,
                      data.table(block.id = block.id,
                                 og.id = og.id,
                                 genome1 = genome2,
                                 genome2 = genome1,
                                 id1 = id2,
                                 id2 = id1,
                                 chr1 = chr2,
                                 chr2 = chr1,
                                 pos1 = pos2,
                                 pos2 = pos1,
                                 stringsAsFactors = F)))
    map <- map[!duplicated(map),]
    map <- merge(combs,
                 map,
                 by = c("genome1", "genome2"))
  }

  return(map)
}

#' @title convert_gff2coords
#' @description
#' \code{convert_gff2coords} makes blocks and merges by block.id
#' @rdname genespace_utils
#' @export
convert_gff2coords <- function(gff,
                               chr.list,
                               gap.prop){

  spl <- split(gff, by = "genome")
  out <- rbindlist(lapply(names(spl), function(i){
    x <- spl[[i]]
    y <- chr.list[[i]]
    x$chr <- factor(x$chr, levels = y)
    setkey(x, chr)
    x[, start.l := calc_linearCoord(chr = chr,
                                    start = start,
                                    end = end,
                                    gap.prop = gap.prop,
                                    scale.it = F,
                                    return.start = T)]
    x[, end.l := calc_linearCoord(chr = chr,
                                  start = start,
                                  end = end,
                                  gap.prop = gap.prop,
                                  scale.it = F,
                                  return.start = F)]
    x[, start.p := calc_linearCoord(chr = chr,
                                    start = start,
                                    end = end,
                                    gap.prop = gap.prop,
                                    scale.it = T,
                                    return.start = T)]
    x[, end.p := calc_linearCoord(chr = chr,
                                  start = start,
                                  end = end,
                                  gap.prop = gap.prop,
                                  scale.it = T,
                                  return.start = F)]
    x$chr <- as.character(x$chr)
    return(x)
  }))
  return(out)
}


#' @title calc_linearCoord
#' @description
#' \code{calc_linearCoord} convert start / end positions into
#' linear coordinates..
#' @rdname genespace_utils
#' @export
calc_linearCoord <- function(chr,
                             start,
                             end,
                             gap.prop,
                             scale.it,
                             return.start){
  ml <- data.table(chr = chr,
                   start = start,
                   end = end,
                   stringsAsFactors = F)
  fais <- ml[, list(start = 0,
                    end = max(end)),
             by = list(chr)]
  setkey(fais, chr)
  fais <- add_gap(fais = fais,
                  gap.prop = gap.prop)
  maxf <- max(fais$endl)
  if (!scale.it)
    maxf <- 1

  bt <- merge(fais[, c("chr","startl")], ml, by = "chr")

  if (return.start) {
    return((bt$start + bt$startl) / maxf)
  }else{
    return((bt$end + bt$startl) / maxf)
  }
}

#' @title add_gap
#' @description
#' \code{add_gap} add gap to linear chromosome x axis.
#' @rdname genespace_utils
#' @export
add_gap <- function(fais,
                    gap.prop){
  gap <- sum(fais$end) * gap.prop
  fais$endl <- (cumsum(fais$end + gap)) - gap
  fais$startl <- fais$endl - fais$end
  return(fais)
}

#' @title make_fais
#' @description
#' \code{make_fais} make data.table akin to fasta index.
#' @rdname genespace_utils
#' @export
make_fais <- function(genomes,
                      gff){
  ypos <- data.table(y = 0:(length(genomes) - 1),
                     genome = genomes)

  fais <- gff[,list(start = min(start.p),
                    end = max(start.p)),
              by = list(genome, chr)]
  fais <- merge(fais,
                ypos,
                by = "genome")
  fais$genome <- factor(fais$genome, levels = genomes)
  setkey(fais, genome)
  fais$genome <- as.character(fais$genome)
  return(fais)
}

#' @title draw_blkPolygon
#' @description
#' \code{draw_blkPolygon} draw blk Polygon
#' @rdname plot_utils
#' @export
draw_blkPolygon <- function(blk,
                            chr.buffer,
                            blk.border,
                            simplify.poly,
                            points.per.curve){
  blk[,unique := paste(genome1, genome2)]
  spl <- split(blk, by = "unique")
  for (j in names(spl)) {
    tmp <- spl[[j]]
    ys <- chr.buffer + tmp$y1[1]
    ye <- (1 - chr.buffer) + tmp$y1[1]
    y <- cos_y(y.start = ys,
               y.end = ye,
               n.out = points.per.curve)

    for (i in 1:nrow(tmp)) {
      with(tmp[i,], make_polygon(
        s1 = start1, e1 = end1,
        s2 = start2, e2 = end2,  y = y,
        fill.color = blk.col,
        simplify.poly = simplify.poly,
        border.color = blk.border))
    }
  }
}

#' @title draw_genePolygon
#' @description
#' \code{draw_genePolygon} draw gene Polygon
#' @rdname genespace_utils
#' @export
draw_genePolygon <- function(map,
                             genes2plot,
                             chr.buff,
                             ortho.col,
                             simplify.poly,
                             points.per.curve){
  ogs2plot <- map$og.id[map$id1 %in% genes2plot |
                          map$id2 %in% genes2plot]
  mapt <- subset(map,og.id %in% ogs2plot)
  mapt[,unique := paste(genome1, genome2)]
  spl <- split(mapt, by = "unique")
  for (j in names(spl)) {
    tmp <- spl[[j]]
    ys <- chr.buff + tmp$y1[1]
    ye <- (1 - chr.buff) + tmp$y1[1]
    y <- cos_y(y.start = ys,
               y.end = ye,
               n.out = points.per.curve)

    for (i in 1:nrow(tmp)) {
      with(tmp[i,], make_polygon(
        s1 = pos1,
        e1 = pos1,
        s2 = pos2,
        e2 = pos2,
        y = y,
        fill.color = ortho.col,
        simplify.poly = simplify.poly,
        border.color = ortho.col))
    }
  }
}

#' @title annotate_riparian
#' @description
#' \code{annotate_riparian} annotate_riparian
#' @rdname genespace_utils
#' @export
annotate_riparian <- function(map,
                              genomes,
                              chr.buff,
                              dodge.geneIDs,
                              scale2dodge,
                              geneid.cex,
                              geneid.offset,
                              gene.colors,
                              geneID.abbrev.fun,
                              genes2plot){
  mapg1 <-  with(subset(map,
                        genome1 == genomes[1] &
                          !duplicated(id1) &
                          id1 %in% genes2plot),
                 data.table(genome = genome1,
                            id = id1,
                            y = y1,
                            x = pos1))
  mapg2 <-  with(subset(map,
                        genome2 == genomes[length(genomes)] &
                          !duplicated(id2) &
                          id2 %in% genes2plot),
                 data.table(genome = genome2,
                            id = id2,
                            y = y2,
                            x = pos2))
  mapg1$x1 <- get_xClus(x = mapg1$x,
                        n.reps = 5,
                        dodge.x = dodge.geneIDs,
                        scale2dodge = scale2dodge)
  mapg2$x1 <- get_xClus(x = mapg2$x,
                        n.reps = 5,
                        dodge.x = dodge.geneIDs,
                        scale2dodge = scale2dodge)

  mapg1$y1 <- mapg1$y - geneid.offset
  mapg1$y <- mapg1$y - chr.buff

  mapg2$y1 <- mapg2$y + geneid.offset
  mapg2$y <- mapg2$y + chr.buff

  with(rbind(mapg1, mapg2),
       segments(y0 = y,
                y1 = y1,
                x0 = x,
                x1 = x1))


  mapg1$lab1 <- geneID.abbrev.fun(mapg1$id)
  for (i in 1:nrow(mapg1)) {
    if (mapg1$id[i] %in% genes2plot) {
      if (is.null(names(genes2plot))) {
        mapg1$lab1[i] <- mapg1$id[i]
      }else{
        mapg1$lab1[i] <- names(genes2plot)[genes2plot == mapg1$id[i]]
      }
    }
  }
  mapg2$lab2 <- geneID.abbrev.fun(mapg2$id)
  for (i in 1:nrow(mapg2)) {
    if (mapg2$id[i] %in% genes2plot) {
      if (is.null(names(genes2plot))) {
        mapg2$lab2[i] <- mapg2$id[i]
      }else{
        mapg2$lab2[i] <- names(genes2plot)[genes2plot == mapg2$id[i]]
      }
    }
  }

  with(mapg1,
       text(lab1,
            y = y1,  x = x1, srt = 90,
            cex = geneid.cex, adj = c(1.05,.5),
            col = gene.colors[id]))

  with(mapg2,
       text(lab2,
            y = y1,  x = x1, srt = 90,
            cex = geneid.cex, adj = c(-0.05,.5),
            col = gene.colors[id]))

}

#' @title label_riparian
#' @description
#' \code{label_riparian} label_riparian
#' @rdname genespace_utils
#' @export
label_riparian <- function(fais,
                           chr.segm.lwd,
                           chr.segm.col,
                           genomes,
                           chr.lab.buff,
                           chr.bg.col,
                           chr.bg.cex,
                           chr.bg.pch,
                           chr.id.col,
                           chr.id.cex,
                           lab.chr,
                           lab.chr.1only,
                           chr.abbrev.fun){

  if (length(chr.segm.col) == length(genomes)) {
    names(chr.segm.col) <- genomes
    chr.segm.col <- chr.segm.col[match(fais$genome, names(chr.segm.col))]
  }else{
    chr.segm.col <- chr.segm.col[1]
  }

  if (length(chr.segm.lwd) == length(genomes)) {
    names(chr.segm.lwd) <- genomes
    chr.segm.lwd <- chr.segm.lwd[match(fais$genome, names(chr.segm.lwd))]
  }else{
    chr.segm.lwd <- chr.segm.lwd[1]
  }

  # add in chromosomes and IDs
  with(fais, segments(x0 = start,
                      x1 = end,
                      y0 = y,
                      y1 = y,
                      lwd = chr.segm.lwd,
                      col = chr.segm.col))

  fais$class <- with(fais,
                     ifelse(genome == genomes[1],
                            "bottom",
                            ifelse(genome == genomes[length(genomes)],
                                   "top","middle")))
  fais$lab.y <- with(fais,
                     ifelse(class == "middle", y,
                            ifelse(class == "top", y - chr.lab.buff,
                                   y + chr.lab.buff)))

  if (!is.null(chr.bg.col) & lab.chr) {
    points(x = rowMeans(fais[,c("start","end")]),
           y = fais$lab.y,
           col = chr.bg.col,
           pch = chr.bg.pch,
           cex = chr.bg.cex)
  }
  if (lab.chr) {
    if (lab.chr.1only) {
      fais2 <- fais[fais$genome == genomes[1],]
    }else{
      fais2 <- fais
    }
    text(chr.abbrev.fun(fais2$chr),
         x = rowMeans(fais2[,c("start","end")]),
         y = fais2$lab.y,
         col = chr.id.col,
         cex = chr.id.cex)
  }
}


#' @title add_alpha
#' @description
#' \code{add_alpha} add transparency to a color.
#' @rdname genespace_utils
#' @importFrom grDevices col2rgb rgb
#' @export
add_alpha <- function(col,
                      alpha = 1){
  if (missing(col))
    stop("Please provide a vector of colours.")

  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1],
              x[2],
              x[3],
              alpha = alpha))
}


#' @title make_polygon
#' @description
#' \code{make_polygon} given a y string and start/end coords, draw a polygon
#' @rdname genespace_utils
#' @importFrom graphics polygon
#' @export
make_polygon <- function(s1,
                         s2,
                         e1,
                         e2,
                         y,
                         fill.color,
                         border.color,
                         simplify.poly){
  y <- y[!is.na(y)]
  ye <- max(y)
  ys <- min(y)

  x1 <- seq(from = s1,
            to = s2,
            length.out = length(y))
  x3 <- seq(from = e2,
            to = e1,
            length.out = length(y))

  x2 <- c(s2,e2)
  x4 <- c(e1,s1)

  y1 <- y
  y2 <- c(ye,ye)
  y3 <- rev(y)
  y4 <- c(ys,ys)
  if (!is.null(simplify.poly)) {
    qs <- qnorm(seq(.0001,
                    .9999,
                    length.out = length(y)))
    wh <- which(!duplicated(round_toAny(qs, simplify.poly)))
    y1 <- y1[wh]
    y3 <- y3[wh]
    x1 <- x1[wh]
    x3 <- x3[wh]
  }
  polygon(c(x1,x2,x3,x4),
          c(y1,y2,y3,y4),
          col = fill.color,
          border = border.color)
}

#' @title round_toAny
#' @description
#' \code{round_toAny} round_toAny
#' @rdname genespace_utils
#' @importFrom dbscan frNN dbscan
#' @export
round_toAny <- function(x,
                        num){
  num * round(x / num)
}

#' @title get_xClus
#' @description
#' \code{get_xClus} for a set of x,
#' @rdname genespace_utils
#' @importFrom dbscan frNN dbscan
#' @export
get_xClus <- function(x,
                      n.reps,
                      dodge.x,
                      scale2dodge = 2){
  for (i in n.reps) {
    clus <- run_dbs(y = data.table(rank1 = as.numeric(x),
                                   rank2 = as.numeric(x)),
                    eps.radius = dodge.x*scale2dodge,
                    mappings = 1)$cluster
    out <- data.table(x = x,
                      clus = clus)
    out[,round := mean(x), by = list(clus)]
    out[, rank := as.numeric(frank(x, ties.method = "random")),
        by = list(clus)]
    out[,n.inclus := .N,
        by = list(clus)]
    out[,n.inclus := .N,
        by = list(clus)]
    out[, rank := ifelse(n.inclus == 1, 1,(rank - mean(rank)) + 1),
        by = list(clus)]
    out[, x1 := (rank * dodge.x + round) - dodge.x,
        by = list(clus)]
    x <- out$x1
  }
  return(x)
}

#' @title cos_y
#' @description
#' \code{cos_y} given a range of x, return a cosine curve.
#' @rdname genespace_utils
#' @export
cos_y <- function(y.start,
                  y.end,
                  n.out,
                  n.sample = n.out){
  rs <- seq(from = 0,
            to = pi,
            length.out = n.sample)
  xt <- 1 - cos(rs)
  xo <- seq(from = min(xt),
            to = max(xt),
            length.out = n.out)
  yvals <- sapply(xo, function(x)
    which.min(abs(x - xt)))
  y <- y.start + (yvals - min(yvals)) *
    ((y.end - y.start) / (max(yvals) - min(yvals)))
  return(y)
}

#' @title track_blocks
#' @description
#' \code{track_blocks} track blocks
#' @rdname genespace_utils
#' @export
track_blocks <- function(map,
                         blk.ids,
                         cols,
                         bg.col){
  ulist <- unlist(lapply(blk.ids, function(x)
    unique(subset(map, block.id == x)$og.id)))
  if (any(duplicated(ulist)))
    ulist <- ulist[!duplicated(ulist)]

  if (length(cols) != length(blk.ids))
    stop("cols and blk.ids must be of the same length")

  blk.cols <- bg.col
  for (i in 1:length(blk.ids)) {
    bi <- blk.ids[i]
    ci <- cols[i]

    mb <- subset(map, block.id == bi)
    ob <- unique(mb$og.id)

    map$block.id <- with(map, ifelse(og.id %in% ob, paste0(block.id, "_inblk",i,"col"), block.id))
    bo <- rep(ci, length(grep(paste0("_inblk",i,"col"), map$block.id)))
    names(bo) <- map$block.id[grep(paste0("_inblk",i,"col"), map$block.id)]
    bo <- bo[!duplicated(names(bo))]
    blk.cols <- c(blk.cols, bo)
  }
  return(list(map = map,
              blk.cols = blk.cols))
}

#' @title mirror_blast
#' @description
#' \code{mirror_blast} mirror_blast
#' @rdname genespace_utils
#' @export
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

#' @title run_of
#' @description
#' \code{run_of} run_of
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


#' @title pull_nonSynOrthos
#' @description
#' \code{pull_nonSynOrthos} pull_nonSynOrthos
#' @rdname genespace_utils
#' @export
pull_nonSynOrthos <- function(of.blast,
                              map,
                              gff,
                              rank.buffer,
                              verbose = T){

  ob <- data.table(subset(of.blast, genome1 != genome2))
  ob2 <- data.table(ob)
  setnames(ob2,
           grep("1$", colnames(ob2)),
           gsub("1$","xxxx",
                colnames(ob2)[grep("1$", colnames(ob2))]))
  setnames(ob2,
           grep("2$", colnames(ob2)),
           gsub("2$","1",
                colnames(ob2)[grep("2$", colnames(ob2))]))
  setnames(ob2,
           grep("xxxx", colnames(ob2)),
           gsub("xxxx","2", colnames(ob2)[grep("xxxx", colnames(ob2))]))

  blo <- rbind(ob, ob2)
  blo <- rbind(blo, subset(of.blast, genome1 == genome2))


  ob <- data.table(subset(map, genome1 != genome2))
  ob2 <- data.table(ob)
  setnames(ob2,
           grep("1$", colnames(ob2)),
           gsub("1$","xxxx",
                colnames(ob2)[grep("1$", colnames(ob2))]))
  setnames(ob2,
           grep("2$", colnames(ob2)),
           gsub("2$","1",
                colnames(ob2)[grep("2$", colnames(ob2))]))
  setnames(ob2,
           grep("xxxx", colnames(ob2)),
           gsub("xxxx","2", colnames(ob2)[grep("xxxx", colnames(ob2))]))

  mapo <- rbind(ob, ob2)
  mapo <- rbind(mapo, subset(map, genome1 == genome2))
  mapo <- mapo[!duplicated(mapo[,c("id1","id2")])]
  blo <- blo[!duplicated(blo[,c("id1","id2")])]
  blast.in.map <- extend_blocks(
    gff = gff,
    map = mapo,
    blast = blo,
    n.iter = 1,
    rank.buffer = rank.buffer,
    clean.it = F,
    verbose = verbose)$map

  ui.inmap <- with(blast.in.map, unique(paste(id1, id2)))
  blo2 <- subset(blo, !paste(id1, id2) %in% ui.inmap)
  if ("uniq" %in% colnames(blo2))
    blo2$uniq <- NULL
  return(blo2)
}


#' @title simplify_map
#' @description
#' \code{simplify_map} simplify_map
#' @rdname genespace_utils
#' @export
simplify_map <- function(map,
                         n.mapping,
                         radius,
                         genomeIDs,
                         keep.best.id1.hit){

  mapc <- with(map,
               data.table(genome1 = genome2, genome2 = genome1,
                          id1 = id2,  id2 = id1,
                          chr1 = chr2, chr2 = chr1,
                          start1 = start2, start2 = start1,
                          end1 = end2, end2 = end1,
                          score = score))

  mapc <- data.table(rbind(map[,colnames(mapc), with = F],
                           mapc))
  mapc <- mapc[!duplicated(mapc),]

  geni <- t(combn(genomeIDs, 2))
  geni <- rbind(geni, t(sapply(genomeIDs, function(x) c(x,x))))
  geni <- data.table(geni)
  setnames(geni, c("genome1", "genome2"))
  mapm <- merge(geni, mapc, by = c("genome1", "genome2"))

  mapm <- clean_blocks(
    map = mapm,
    radius = radius,
    n.mappings = n.mapping,
    genomeIDs = genomeIDs)$map


  mapo <- with(mapm,
               data.table(block.id = block.id,
                          genome1 = genome2, genome2 = genome1,
                          id1 = id2,  id2 = id1,
                          chr1 = chr2, chr2 = chr1,
                          start1 = start2, start2 = start1,
                          end1 = end2, end2 = end1,
                          score = score))
  out <- data.table(rbind(mapm[,colnames(mapo), with = F],
                           mapo))

  out <- out[!duplicated(out),]
  if (keep.best.id1.hit)
    out <- drop_mapDup(out)

  return(out)
}

#' @title drop_mapDup
#' @description
#' \code{drop_mapDup} drop_mapDup
#' @rdname genespace_utils
#' @export
drop_mapDup <- function(map){
  map[,ns := -score]
  setkey(map, ns)
  out <- map[, head(.SD, 1), by = list(block.id, id1)]
  out$ns <- NULL
  return(out)
}

#' @title extend_blk2chrend
#' @description
#' \code{extend_blk2chrend} extend_blk2chrend
#' @rdname genespace_utils
#' @export
extend_blk2chrend <- function(blk,
                              gff,
                              map,
                              genomeIDs,
                              min.dist2end){

  gl <- gff[,list(min = min(start),
                  max = max(end)),
            by = list(genome, chr)]

  spl.gff <- split(gff, by = "genome")
  spl.map <- split(map, by = c("genome1","genome2"))

  for (x in spl.map) {
    g1 <- subset(spl.gff[[x$genome1[1]]], id %in% x$id1)
    g2 <- subset(spl.gff[[x$genome2[1]]], id %in% x$id2)
    gm1 <- g1[,list(g.min = id[frank(start, ties.method = "random") <= min.dist2end],
                    g.max = id[frank(-end, ties.method = "random") <= min.dist2end]),
              by = list(genome, chr)]
    gm2 <- g2[,list(g.min = id[frank(start, ties.method = "random") <= min.dist2end],
                    g.max = id[frank(-end, ties.method = "random") <= min.dist2end]),
              by = list(genome, chr)]

    gm1 <- merge(gl, gm1, by = c("genome", "chr"))
    gm2 <- merge(gl, gm2, by = c("genome", "chr"))

    for (j in 1:nrow(gm1)) {
      wh.blks <- unique(x$block.id[x$id1 %in% gm1$g.min[j]])
      blk$start1[with(blk, block.id %in% wh.blks)] <- gm1$min[j]

      wh.blks <- unique(x$block.id[x$id1 %in% gm1$g.max[j]])
      blk$end1[with(blk, block.id %in% wh.blks)] <- gm1$max[j]
    }

    for (j in 1:nrow(gm2)) {
      wh.blks <- unique(x$block.id[x$id2 %in% gm2$g.min[j]])
      blk$start2[with(blk, block.id %in% wh.blks)] <- gm2$min[j]

      wh.blks <- unique(x$block.id[x$id2 %in% gm2$g.max[j]])
      blk$end2[with(blk, block.id %in% wh.blks)] <- gm2$max[j]
    }
  }

  return(blk)
}

#' @title pull_duplicatesInChr
#' @description
#' \code{pull_duplicatesInChr} pull_duplicatesInChr
#' @rdname genespace_utils
#' @export
pull_duplicatesInChr <- function(map.bychr,
                                 gff){
  dup.blks <- character()
  b <- map.bychr[,list(start1 = min(start1),
                       end1 = max(end1),
                       start2 = min(start2),
                       end2 = max(end2)),
                 by = list(genome1,genome2,chr1,chr2,block.id)]
  if (nrow(b) > 1) {
    b$block.id <- as.character(b$block.id)

    ovl <- 1
    nb <- nrow(b)
    while (ovl > .1 & nb > 1) {
      genes.byBlk <- lapply(1:nrow(b), function(i){
        subset(gff, genome == b$genome1[i] & chr == b$chr1[i] &
                 start >= b$start1[i] & end <= b$end1[i])$id
      })
      names(genes.byBlk) <- as.character(b$block.id)

      eg <- data.table(expand.grid(as.character(b$block.id),
                                   as.character(b$block.id),
                                   stringsAsFactors = F)[,c(2,1)])
      eg <- subset(eg, Var1 != Var2)

      eg$prop.overlap <- apply(eg, 1, function(x)
        length(intersect(genes.byBlk[[x[1]]],
                         genes.byBlk[[x[2]]]))/length(genes.byBlk[[x[1]]]))
      eg <- eg[order(-eg$prop.overlap),]
      ovl <- eg$prop.overlap[1]
      if (eg$prop.overlap[1] > 0.1) {
        dup.blks <- c(dup.blks, eg$Var2[1])
        b <- b[!b$block.id %in% dup.blks]
        nb <- nrow(b)
      }
    }
  }
  return(dup.blks)
}

#' @title find_blocksToExtend
#' @description
#' \code{find_blocksToExtend} find_blocksToExtend
#' @rdname genespace_utils
#' @export
find_blocksToExtend <- function(blk,
                                chr.end.buffer = 5,
                                ovl.buffer = 3,
                                same.chr.buffer = 500,
                                verbose = T){

  if(verbose)
    cat("Finding adjacent blocks ...")
  adj.blk1 <- blk[,find_adjBlks1(.SD,
                                 ovl.buffer = ovl.buffer,
                                 same.chr.buffer = same.chr.buffer),
                  by = list(genome1, genome2, chr1)]
  adj.blk2 <- blk[,find_adjBlks2(.SD,
                                 ovl.buffer = ovl.buffer,
                                 same.chr.buffer = same.chr.buffer),
                  by = list(genome1, genome2, chr2)]
  adjb <- merge(adj.blk1, adj.blk2, by = c("genome1","genome2","block.id"))
  adjb <- subset(adjb, block.id != blk.left | block.id != blk.right |
                   block.id != blk.down | block.id != blk.up)
  rs <- rowSums(is.na(adjb[,c("blk.left","blk.right","blk.up","blk.down")]))
  blks.adj <- adjb[rs <= 3,]
  if(verbose)
    cat("Done\n")

  m <- melt(blks.adj, id.vars = c("genome1","genome2","chr1","chr2"),
            measure.vars = c("block.id","blk.left","blk.right","blk.up","blk.down"))
  m[,variable:= NULL]
  m <- m[!duplicated(m),]
  spl.adj <- split(m, by = c("genome1","genome2","chr1","chr2"))
  blk[,chr.rankend1 := max(rankend1),
      by = list(genome1, genome2, chr1)]
  blk[,chr.rankend2 := max(rankend2),
      by = list(genome1, genome2, chr2)]
  blk[,n.inchr := .N,
      by = list(genome1, genome2, chr1, chr2)]
  blks.solo <- subset(blk, n.inchr == 1 &
                        rankstart1 > chr.end.buffer &
                        (rankend1 - chr.end.buffer) < chr.rankend1 &
                        rankstart2 > chr.end.buffer &
                        (rankend2 - chr.end.buffer) < chr.rankend2)

  blks2test <- unique(c(blks.solo$block.id, blks.adj$block.id))
  return(blks2test)
}

#' @title find_leftBlk
#' @description
#' \code{find_leftBlk} find_leftBlk
#' @rdname genespace_utils
#' @export
find_leftBlk <- function(st.rank, end.ranks, chrs, chr, ovl.buffer, same.chr.buffer){

  is.left.onchr <- st.rank > (end.ranks - ovl.buffer) & chrs == chr
  is.left.offchr <- st.rank > (end.ranks - ovl.buffer) & chrs != chr
  if(sum(is.left.onchr | is.left.offchr) == 0){
    out <- NA
  }else{
    if(sum(is.left.offchr) == 0 | sum(is.left.onchr) == 0){
      ens <- end.ranks[is.left.onchr | is.left.offchr]
      out <- which(end.ranks == ens[which.min(st.rank - ens)])[1]
    }else{
      min.dist.onchr <- min(st.rank - end.ranks[is.left.onchr])
      min.dist.offchr <- min(st.rank - end.ranks[is.left.offchr])
      if(min.dist.onchr < min.dist.offchr + same.chr.buffer){
        ens <- end.ranks[is.left.onchr]
        out <- which(end.ranks == ens[which.min(st.rank - ens)] &
                       is.left.onchr)[1]
      }else{
        ens <- end.ranks[is.left.offchr]
        out <- which(end.ranks == ens[which.min(st.rank - ens)] &
                       is.left.offchr)[1]
      }
    }
  }
  return(out)
}

#' @title find_rightBlk
#' @description
#' \code{find_rightBlk} find_rightBlk
#' @rdname genespace_utils
#' @export
find_rightBlk <- function(st.ranks, end.rank, chrs, chr, ovl.buffer, same.chr.buffer){

  is.right.onchr <- st.ranks > (end.rank - ovl.buffer) & chrs == chr
  is.right.offchr <- st.ranks > (end.rank - ovl.buffer) & chrs != chr
  if(sum(is.right.onchr | is.right.offchr) == 0){
    out <- NA
  }else{
    if(sum(is.right.offchr) == 0 | sum(is.right.onchr) == 0){
      sts <- st.ranks[is.right.onchr | is.right.offchr]
      out <- which(st.ranks == sts[which.min(sts - end.rank)])[1]
    }else{
      min.dist.onchr <- min(st.ranks[is.right.onchr] - end.rank)
      min.dist.offchr <- min(st.ranks[is.right.offchr] - end.rank)
      if(min.dist.onchr < min.dist.offchr + same.chr.buffer){
        sts <- st.ranks[is.right.onchr]
        out <-  which(st.ranks == sts[which.min(sts - end.rank)] &
                        is.right.onchr)[1]
      }else{
        sts <- st.ranks[is.right.offchr]
        out <-  which(st.ranks == sts[which.min(sts - end.rank)] &
                        is.right.offchr)[1]
      }
    }
  }
  return(out)
}

find_adjBlks1 <- function(b, ovl.buffer, same.chr.buffer){
  which.left <- sapply(1:nrow(b), function(i)
    find_leftBlk(st.rank = b$rankstart1[i], end.ranks = b$rankend1,
                 chr = b$chr2[i], chrs = b$chr2,
                 ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  which.right <- sapply(1:nrow(b), function(i)
    find_rightBlk(st.ranks = b$rankstart1, end.rank = b$rankend1[i],
                  chr = b$chr2[i], chrs = b$chr2,
                  ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))

  return(data.table(block.id = b$block.id,
                    start1 = b$start1,
                    end1 = b$end1,
                    blk.left = b$block.id[which.left],
                    blk.right = b$block.id[which.right],
                    blk.left.end = b$end1[which.left],
                    blk.right.start = b$start1[which.right]))
}

find_adjBlks2 <- function(b, ovl.buffer, same.chr.buffer){
  which.left <- sapply(1:nrow(b), function(i)
    find_leftBlk(st.rank = b$rankstart2[i], end.ranks = b$rankend2,
                 chr = b$chr1[i], chrs = b$chr1,
                 ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  which.right <- sapply(1:nrow(b), function(i)
    find_rightBlk(st.ranks = b$rankstart2, end.rank = b$rankend2[i],
                  chr = b$chr1[i], chrs = b$chr1,
                  ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  return(data.table(block.id = b$block.id,
                    start2 = b$start2,
                    end2 = b$end2,
                    blk.down = b$block.id[which.left],
                    blk.up = b$block.id[which.right],
                    blk.down.end = b$end1[which.left],
                    blk.up.start = b$start1[which.right]))
}


