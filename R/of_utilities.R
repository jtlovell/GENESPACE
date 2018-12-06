#' @title Orthofinder utility functions
#' @name of_utilities
#' @aliases split.data.table
#' @aliases get_ofIDs
#' @aliases make_ofInputInBlk
#' @aliases split_gffByBlock
#' @aliases make_mergedBlocks
#' @aliases make_MCSBlocks
#' @aliases cull_blastByMCS
#' @aliases cull_blastByDBS
#' @aliases cull_blast2og
#' @aliases merge_gffBlast
#' @aliases import_blast
#' @aliases import_ofResults
#' @aliases import_gff
#' @aliases parse_gff
#'
#' @description
#' \code{of_utilities} Several utilities functions meant for internal calls in compareGeneSpace
#' @name of_utilities
#' @param gff.dir character, directory containing the gff3 formatted annotation files
#' @param blast.dir character, directory containing the orthofinder output
#' @param gff.file character, gff file name and path
#' @param genomeIDs character, genome identifiers
#' @param abbrevs character, genome abbreviations
#' @param blk data.table containing the block information
#' @param map data.table containing the map information
#' @param blast data.table containing blast hits
#' @param gff data.table containing the parsed gff annotation data
#' @param ogff list of data.tables, split gff by block.
#' @param orthogroups orthogroup object from import_ofResults
#' @param gene.index gene index from import_ofResults
#' @param MCScanX.params character, parameters to be passed to MCScanX
#' @param mcscanx.input.dir directory for mcscan temporary files to be stored
#' @param n.mappingWithinRadius numeric, number of hits required to be in the radius
#' @param eps.radius numeric, size of the radius
#' @param pairs.only logical, should only pairs of hits in orthofinder output be retained
#' @param min.propMax numeric, minimum proportion of max score for a gene
#' @param min.score numeric, minimum score for a hit to be retained
#' @param max.hitsPerGene numeric, maximum number of hits to be retained per gene
#' @param str2drop character, string in attribute column of gff file to be dropped
#' @param str2parse character, string in attribute column of gff file to use as the separator
#' @param whichAttr numeric, which attribute should be returned in the
#' gff attribute column
#' @param out.dir directory where should output be written.
#' @param species.mappings the species mapping from parse_orthofinder
#' @param of.speciesIDs orthofinder species ids, from parse_orthofinder
#' @param of.geneIDs orthofinder gene ids, from parse_orthofinder
#' @param x data.table
#' @param f factor
#' @param by factor
#' @param drop drop
#' @param flatten non-recursive unlisting
#' @param verbose logical, should updates be reported?
#' @param ... additional arguments passed to data.table
#' @note \code{of_utilities} is a generic name for the functions documented.
#' \cr
#' If called, \code{of_utilities} returns its own arguments.
#'

#' @title Fast split of data.table
#' @description
#' \code{split.data.table} Much faster than base split.
#' @rdname utilities
#' @import data.table
#' @export
split.data.table <- function(x,
                             f,
                             drop = FALSE,
                             by,
                             flatten = FALSE,
                             ...){
  if(missing(by) && !missing(f)) by = f
  stopifnot(!missing(by), is.character(by), is.logical(drop), is.logical(flatten), !".ll" %in% names(x), by %in% names(x), !"nm" %in% by)
  if(!flatten){
    .by = by[1L]
    tmp = x[, list(.ll=list(.SD)), by = .by, .SDcols = if(drop) setdiff(names(x), .by) else names(x)]
    setattr(ll <- tmp$.ll, "names", tmp[[.by]])
    if(length(by) > 1L) return(lapply(ll, split.data.table, drop = drop, by = by[-1L])) else return(ll)
  } else {
    tmp = x[, list(.ll=list(.SD)), by=by, .SDcols = if(drop) setdiff(names(x), by) else names(x)]
    setattr(ll <- tmp$.ll, 'names', tmp[, .(nm = paste(.SD, collapse = ".")), by = by, .SDcols = by]$nm)
    return(ll)
  }
}

#' @title Fast reading of gff3 files
#' @description
#' \code{parse_gff} Much faster parsing and reading of gff3 files than
#' bioconductor txdb methods.
#' @rdname utilities
#' @import data.table
#' @export
parse_gff <- function(gff.file,
                      str2drop = "Name=",
                      str2parse = ";",
                      whichAttr = 2){
  g <- suppressWarnings(
    data.table::fread(gff.file,
                      showProgress = F,
                      verbose = F))
  g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
  g$V9 <- sapply(g$V9, function(x)
    gsub(str2drop, "",
         strsplit(x, str2parse)[[1]][whichAttr]))
  data.table::setnames(g, c("id", "chr", "start", "end", "strand"))
  return(g)
}

#' @title import gff annotation
#' @description
#' \code{import_gff} wrapper for parse_gff
#' @rdname utilities
#' @import data.table
#' @export
import_gff <- function(gff.dir,
                       genomeIDs,
                       verbose = T,
                       str2drop = "Name=",
                       str2parse = ";",
                       whichAttr = 2){

  if (verbose)
    cat("Importing gff3 annotation files\n")
  gff.files <- file.path(gff.dir,
                         paste0(genomeIDs, ".gff3"))
  names(gff.files) <- genomeIDs

  gff <- rbindlist(lapply(names(gff.files), function(i){
    if (verbose)
      cat("\tReading",i,"... ")
    tmp <- parse_gff(gff.file = gff.files[[i]],
                     str2drop = str2drop,
                     str2parse = str2parse,
                     whichAttr = whichAttr)
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    if (verbose)
      cat("n. genes =", nrow(tmp),"\n")
    return(tmp)
  }))

  setkey(gff, "genome", "id")
  if (verbose)
    cat("\tDone!\n")
  return(gff)
}

#' @title Import orthofinder results
#' @description
#' \code{import_ofResults} Primary processing and import of orthofinder metadata
#' @rdname of_utilities
#' @import data.table
#' @export
import_ofResults <- function(gff,
                             blast.dir,
                             verbose = T,
                             genomeIDs){
  if (verbose)
    cat("Importing orthofinder results:\n\t")
  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("Decompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(blast.dir,
                           "*.gz")))
  }

  #######################################################
  if (verbose)
    cat("Reading Species IDs\n\t")
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
  sm <- sm[sm$genome1 != sm$genome2,]
  sm$ref <- NA
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
    cat("Reading gene IDs\n\t")
  sequence.index <- fread(file.path(blast.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))

  #######################################################
  if (verbose)
    cat("Reading orthogroup networks\n\t")
  og <- readLines(file.path(blast.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])

  if (verbose)
    cat("Building orthogroup data.table\n\t")
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


  if (verbose)
    cat("Compiling metadata\n")

  gffi = gff[,c("id","genome")]
  setkey(gffi, "id")
  setkey(og2, "id")
  ogo = merge(gffi, og2)
  setkey(ogo, "block.id", "genome", "og", "id")
  ogo[, og.n.genes := length(unique(id)), by = list(block.id, og)]
  ogo[, og.n.genomes := length(unique(genome)), by = list(block.id, og)]

  ogo.meta = ogo[!duplicated(ogo[, -c(1:2), with = F]),-c(1:2), with = F]

  #######################################################
  if (verbose)
    cat("Done!\n")
  return(list(orthogroups = og,
              species.mappings = sm,
              species.index = si,
              orthogroup.metadata = ogo.meta,
              orthogroup.data.table = ogo,
              gene.index = sequence.index))
}

#' @title Import orthofinder blast results
#' @description
#' \code{import_blast} Simple import and processing of orthofinder blast results
#' @rdname of_utilities
#' @import data.table
#' @export
import_blast <- function(species.mappings,
                         genomeIDs,
                         verbose = T,
                         orthogroups,
                         gene.index,
                         gff,
                         mcscanx.input.dir,
                         MCScanX.param = "-a -s 10 -m 25 -w 2 -e 1",
                         n.mappingWithinRadius = c(5,5,5),
                         eps.radius = c(100,50,25),
                         pairs.only = T){

  if (verbose)
    cat("Parsing gff annotations\n")

  gff1 <- data.table(gff)
  gff2 <- data.table(gff)
  setnames(gff1,c("id1", "chr1", "start1", "end1",
                  "strand1", "genome1", "order1"))
  setnames(gff2,c("id2", "chr2", "start2", "end2",
                  "strand2", "genome2", "order2"))
  setkey(gff1, "id1")
  setkey(gff2, "id2")

  spl.gff1 <- split.data.table(gff1, "genome1")
  spl.gff2 <- split.data.table(gff2, "genome2")

  if (verbose)
    cat("Parsing orthogroups\n")
  wh2 <- which(sapply(orthogroups, length) > 1)
  g2 <- orthogroups[wh2]
  l2 <- sapply(g2, length)
  n2 <- rep(names(l2), l2)

  g1 <- data.table(og = n2,
                   id = unlist(g2),
                   stringsAsFactors = F)
  setkey(g1, "id")
  setkey(gene.index, "id")
  g1 <- merge(g1, gene.index)
  g2 <- data.table(g1)

  setnames(g1, c("id1","og1", "gn1"))
  setnames(g2, c("id2", "og2", "gn2"))
  setkey(g1, "gn1")
  setkey(g2, "gn2")


  if (verbose)
    cat("Importing BLAST results\n")
  sm <- species.mappings
  comb <- combn(genomeIDs, 2, simplify = F)
  blast <- rbindlist(lapply(comb, function(x){
    if (verbose)
      cat(paste0("\t" ,x[1]), "-->",
          x[2], "... retaining hits:\n\t\t")
    smo <- sm[sm$ref == x[1] & sm$alt == x[2],]
    smo <- smo[order(smo$map.rank),]
    suppressWarnings(b1 <- fread(smo$filename[1],
                                 showProgress = F))
    suppressWarnings(b2 <- fread(smo$filename[2],
                                 showProgress = F))
    b2 <- data.table(b2[, c(2, 1, 3:6, 9:10, 7:8, 11:12)])
    setnames(b2, colnames(b1))

    blast.in <- rbind(b1, b2)
    setnames(blast.in, c("gn1", "gn2", "perc.iden",
                         "align.length", "n.mismatch",
                         "n.gapOpen", "q.start", "q.end",
                         "s.start", "s.end",
                         "eval", "score"))
    blast.in$neg.score <- blast.in$score * (-1)
    setkey(blast.in, "gn1", "gn2", "neg.score")
    blast.out <- blast.in[!duplicated(blast.in[, c("gn1","gn2"), with = F])]

    if (verbose)
      cat(nrow(blast.out), "(unique), ")

    blast.cull <- blast.out

    setkey(g1, "gn1")
    setkey(g2, "gn2")

    setkey(blast.cull, "gn2")
    bl1 <- merge(g2, blast.cull)
    setkey(bl1, "gn1")
    bl2 <- merge(g1, bl1)

    if (pairs.only) {
      bl2 <- data.table(bl2[with(bl2, og1 == og2),])
      if (verbose)
        cat(nrow(bl2), "(in orthogroups), ")
    }else{
      if (verbose)
        cat(nrow(bl2), "(in orthogroups), ")
    }

    gf1 <- spl.gff1[[x[1]]]
    gf2 <- spl.gff2[[x[2]]]

    setkey(bl2, "id2")
    blo <- merge(gf2, bl2)
    setkey(blo, "id1")
    blo2 <- merge(gf1, blo)

    if (!is.null(n.mappingWithinRadius) &
        !is.null(eps.radius)) {
      blo3 <- cull_blastByDBS(
        blast = blo2,
        n.mappingWithinRadius = n.mappingWithinRadius,
        eps.radius = eps.radius,
        verbose = F)

      if (verbose)
        cat(nrow(blo3), "(DBS), ")
    }else{
      if (verbose)
        cat("\n")
      blo3 <- blo2
    }

    if (!is.null(MCScanX.param) &
        !is.null(mcscanx.input.dir)) {
      blo4 <- cull_blastByMCS(
        blast = blo3,
        MCScanX.param = MCScanX.param,
        mcscanx.input.dir = mcscanx.input.dir,
        verbose = F)

      if (verbose)
        cat(nrow(blo4), "(MCScanX)\n")
    }else{
      if (verbose)
        cat("\n")
      blo4 <- blo3
    }

    return(blo4)
  }))
  if (verbose)
    cat("\tDone!\n")
  return(blast)
}

#' @title Combine gff and blast files
#' @description
#' \code{merge_gffBlast} Simple merging of imported gff and blast data.tables from orthofinder.
#' @rdname of_utilities
#' @import data.table
#' @export
merge_gffBlast <- function(gff,
                           blast,
                           verbose = T){
  if (verbose)
    cat("Merging", nrow(blast),"BLAST hits with",
        nrow(gff), "gene annotation entries")
  bl <- blast

  gff1 <- data.table(gff)
  gff2 <- data.table(gff)
  setnames(gff1, c("id1", "chr1", "start1", "end1",
                   "strand1", "genome1", "order1"))
  setnames(gff2, c("id2", "chr2", "start2", "end2",
                   "strand2", "genome2", "order2"))
  setkey(gff1, "id1")
  setkey(gff2, "id2")

  bl$uniq.n = paste(sapply(bl$gn1, function(x) strsplit(x, "_")[[1]][1]),
                    sapply(bl$gn1, function(x) strsplit(x, "_")[[1]][2]))


  setkey(bl, "id2")
  mg1 <- merge(gff2, bl)
  setkey(mg1, "id1")
  mg2 <- merge(gff1, mg1)
  mg2$gn1 <- NULL
  mg2$gn2 <- NULL
  mg2$neg.score <- NULL
  if (verbose)
    cat("\n\tDone!\n")
  return(mg2)
}


#' @title Reduce blast results to orthogroups
#' @description
#' \code{cull_blast2og} Combine orthogroup data with blast results.
#' Drop to only ghits within orthogroups if pairs.only = T.
#' @rdname of_utilities
#' @import data.table
#' @export
cull_blast2og <- function(blast,
                          orthogroups,
                          pairs.only = T,
                          verbose = T){
  blast.cull <- blast
  if (verbose)
    cat("Culling", nrow(blast),
        "BLAST hits to genes within orthogroups\n\t")
  wh2 <- which(sapply(orthogroups, length) > 1)
  g2 <- orthogroups[wh2]
  l2 <- sapply(g2, length)
  n2 <- rep(names(l2), l2)

  g1 <- data.table(og = n2,
                   id = unlist(g2),
                   stringsAsFactors = F)
  g2 <- data.table(g1)
  setnames(g1, c("og1", "id1"))
  setnames(g2, c("og2", "id2"))
  setkey(g1, "id1")
  setkey(g2, "id2")
  setkey(blast.cull, "id2")
  bl1 <- merge(g2, blast.cull)
  setkey(bl1, "id1")
  bl2 <- merge(g1, bl1)

  if (verbose)
    cat("retained", nrow(bl2),
        "hits where either hit is in an orthogroup\n")

  if (pairs.only){
    bl2 <- data.table(bl2[with(bl2, og1 == og2),])
    if(verbose)
      cat("\tretained", nrow(bl2),
          "hits where both hits are in the same orthogroup\n")
  }
  if (verbose)
    cat("\tDone!\n")
  return(bl2)
}

#' @title Reduce blast results by DBS
#' @description
#' \code{cull_blastByDBS} Search pairwise gene-orders for regions of dense hits.
#' Speeds up and improves block construction accuracy.
#' @rdname of_utilities
#' @import data.table
#' @export
cull_blastByDBS <- function(blast,
                            n.mappingWithinRadius = c(5,5,5),
                            eps.radius = c(100,50,25),
                            verbose = T){
  blast.cull <- blast
  if(verbose)
    cat("Culling", nrow(blast.cull),
        "BLAST hits by 2d Density\n")
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

  if(length(n.mappingWithinRadius) != length(eps.radius)){
    warning("eps.radius and n.mappingWithinRadius not of same length")
    if(length(n.mappingWithinRadius) > length(eps.radius)){
      n.mappingWithinRadius <- n.mappingWithinRadius[1:length(eps.radius)]
    }else{
      eps.radius <- eps.radius[1:length(n.mappingWithinRadius)]
    }
  }
  map <- blast.cull
  map$unique <- with(map, paste(genome1, genome2))
  spl <- split.data.table(map, "unique")
  out <- rbindlist(lapply(spl, function(x){
    if (verbose)
      cat(paste0("\t", x$genome1[1]), "-->", x$genome2[1],
          paste0("(initial hits = ", nrow(x),") "))
    for (i in 1:length(eps.radius)) {
      x$rank1 <- frank(x, "chr1", "start1", ties.method = "dense")
      x$rank2 <- frank(x, "chr2", "start2", ties.method = "dense")
      x <- run_dbs(y = x,
                   eps.radius = eps.radius[i],
                   mappings = n.mappingWithinRadius[i])
      x <- x[x$cluster != 0,]
      if (nrow(x) < min(n.mappingWithinRadius)) {
        break
      }
    }
    if (verbose)
      cat("culled hits =", nrow(x), "\n")
    return(x)
  }))
  if (verbose)
    cat("\tDone!\n")
  return(out)
}

#' @title Reduce blast results by MCScanX
#' @description
#' \code{cull_blastByMCS} Search pairwise gene-orders for runs of hits,
#' via MCScanX.
#' @rdname of_utilities
#' @import data.table
#' @export
cull_blastByMCS <- function(blast,
                            MCScanX.param = "-a -s 5 -m 50 -w 5 -e 1",
                            mcscanx.input.dir,
                            verbose = T){

  blast.results <- blast
  if (verbose)
    cat("Culling BLAST results by MCScanX clusters\n")
  blast.results$unique = with(blast.results,
                              paste(genome1, genome2))
  mcscan.cull <- rbindlist(
    lapply(split.data.table(blast.results, "unique"), function(x){
      if (verbose)
        cat(paste0("\t", x$genome1[1]), "-->", x$genome2[1],
            paste0("(initial hits = ", nrow(x), ") "))

      abbrevs = paste0(LETTERS, letters)[1:2]
      names(abbrevs) = c(x$genome1[1], x$genome1[2])

      out = run_MCScanX(blast.results = x,
                        abbrevs = abbrevs,
                        mcscanx.input.dir = mcscanx.input.dir,
                        MCScanX.params = MCScanX.param,
                        verbose = FALSE)
      if (verbose)
        cat("culled to", nrow(out), "\n")
      return(out)
    }))

  blast.results.mcs <- merge(blast.results,
                             mcscan.cull[, 1:2,with = F])
  if (verbose)
    cat("\tDone!\n")
  return(blast.results.mcs)
}

#' @title Build collinear blocks
#' @description
#' \code{cull_blastByMCS} Call MCScanX to build collinear blocks.
#' @rdname of_utilities
#' @import data.table
#' @export
make_MCSBlocks <- function(blast,
                           genomeIDs,
                           mcscanx.input.dir,
                           verbose = T,
                           MCScanX.params = "-a -s 10 -m 10 -w 2"){

  abbrevs = paste0(LETTERS, letters)[length(genomeIDs)]

  blast.results <- blast
  spl = split.data.table(blast.results, "unique")
  comb <- combn(genomeIDs, 2, simplify = F)

  mcscan.list = lapply(1:length(comb), function(i){
    if (verbose)
      cat(comb[[i]][1], "-->", comb[[i]][2])
    x = spl[[paste(comb[[i]], collapse = " ")]]
    mctmp = run_MCScanX(
      blast.results = x,
      abbrevs = abbrevs,
      mcscanx.input.dir = mcscanx.input.dir,
      MCScanX.params = MCScanX.params,
      verbose = FALSE)
    tmp = make_blocks(mctmp)
    if (verbose)
      cat(" ...", nrow(tmp$map), "(of",
          paste0(nrow(x), ")"),
          "hits in", nrow(tmp$block), "blocks\n")
    return(tmp)
  })

  mcscan.blk <- rbindlist(lapply(mcscan.list, function(x) x$block))
  mcscan.map <- rbindlist(lapply(mcscan.list, function(x) x$map))


  return(list(block = mcscan.blk,
              map = mcscan.map))
}

#' @title Merge collinear blocks
#' @description
#' \code{make_mergedBlocks} Merge collinear blocks using the
#' R geometry package.
#' @rdname of_utilities
#' @import data.table
#' @importFrom geometry mesh.drectangle
#' @export
make_mergedBlocks <- function(blk,
                              map,
                              buffer = -1,
                              max.iter = 10,
                              max.size2merge = 200,
                              genomeIDs,
                              verbose = T){
  map$unique <- paste(map$genome1, map$genome2)
  blk$unique <- paste(blk$genome1, blk$genome2)
  spl.blk <- split.data.table(blk, "unique")
  spl.map <- split.data.table(map, "unique")
  comb <- combn(genomeIDs, 2, simplify = F)
  merge.list <- lapply(1:length(comb), function(i){
    if (verbose)
      cat(comb[[i]][1], "-->", comb[[i]][2])

    xblk <- spl.blk[[paste(comb[[i]], collapse = " ")]]
    xmap <- spl.map[[paste(comb[[i]], collapse = " ")]]
    if (verbose)
      cat(" ...", nrow(xmap), "hits in", nrow(xblk), "blocks\n")

    small.blocks <- merge_blocks(blk = xblk,
                                 map = xmap,
                                 buffer = buffer,
                                 max.iter = max.iter,
                                 max.size2merge = max.size2merge)
    return(small.blocks)
  })

  merge.blk <- rbindlist(lapply(merge.list, function(x) x$block))
  merge.map <- rbindlist(lapply(merge.list, function(x) x$map))

  if (verbose)
    cat("Done!")
  return(list(block = merge.blk,
              map = merge.map))
}

#' @title Split gff by block
#' @description
#' \code{split_gffByBlock} Split gff data.tables into blocks.
#' @rdname of_utilities
#' @import data.table
#' @export
split_gffByBlock <- function(blk,
                             gff,
                             verbose = T){
  gff$unique <- with(gff, paste(genome, chr, sep = "_"))
  blk$unique1 <- with(blk, paste(genome1, chr1, sep = "_"))
  blk$unique2 <- with(blk, paste(genome2, chr2, sep = "_"))

  if (verbose)
    cat("Splitting gff into", nrow(blk), "blocks ... ")
  sgff <- split.data.table(gff, "unique")
  gffo <- lapply(1:nrow(blk), function(i){
    btmp <- blk[i,]
    g1 <- sgff[[btmp$unique1]]
    g2 <- sgff[[btmp$unique2]]
    g1 <- g1[g1$start <= btmp$end1 & g1$end >= btmp$start1,]
    g2 <- g2[g2$start <= btmp$end2 & g2$end >= btmp$start2,]
    go <- rbind(g1, g2)
    return(go)
  })

  names(gffo) <- with(blk, paste(genome1, genome2, block.id))
  if (verbose)
    cat("Done!\n")
  return(gffo)
}


#' @title Make orthofinder input
#' @description
#' \code{make_ofInputInBlk} Prepare orthofinder input from initial run and
#' block object.
#' @rdname of_utilities
#' @import data.table
#' @export
make_ofInputInBlk <- function(blast.dir,
                              out.dir,
                              ogff,
                              species.mappings,
                              of.speciesIDs,
                              verbose = T){
  if (verbose)
    cat("Copying orthofinder output to", out.dir,"\n")
  tmp.blast.dir = file.path(out.dir,"tmp")
  tmp.blast.dir.files = file.path(out.dir,"tmp","blasts")

  if (file.exists(tmp.blast.dir)) {
    system(paste("rm -rf", tmp.blast.dir))
  }
  system(paste("mkdir", tmp.blast.dir))
  system(paste("cp -r", blast.dir, out.dir))
  ogn1 <- sapply(names(ogff), function(x)
    strsplit(x, " ")[[1]][1])
  ogn2 <- sapply(names(ogff), function(x)
    strsplit(x, " ")[[1]][2])

  if (verbose)
    cat("Culling blast hits ...\n\t")

  for (i in 1:nrow(species.mappings)) {
    x <- species.mappings[i,]
    g <- ogff[ogn1 == x$ref & ogn2 == x$alt]
    xfile <- basename(x$filename)
    if (verbose)
      cat(xfile)
    f <- fread(file.path(blast.dir,xfile),
               header = F,
               stringsAsFactors = F,
               check.names = F)
    if (verbose)
      cat(paste0(" (", nrow(f), " hits) ... "))
    setkey(f, V1, V2)
    fl <- rbindlist(lapply(1:length(g), function(j){
      y <- g[[j]]
      genes <- y$gene.num
      return(f[f$V1  %in% genes & f$V2 %in% genes,])
    }))

    if (verbose)
      cat("retained", nrow(fl), "hits\n\t")
    write.table(fl,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = file.path(out.dir,xfile))
  }

  for (i in genomeIDs){
    ogl <- ogff[ogn1 == i | ogn2 == i]
    ogs <- lapply(ogl, function(x) x$gene.num[x$genome == i])
    ognum <- of.speciesIDs$genome.num[of.speciesIDs$genome == i]
    xfile <- paste0("Blast", ognum,
                    "_", ognum, ".txt")

    if (verbose)
      cat(xfile)
    f <- fread(file.path(blast.dir,xfile), header = F, stringsAsFactors = F, check.names = F)

    if (verbose)
      cat(paste0(" (", nrow(f), " hits) ... "))
    setkey(f, V1, V2)

    fl <- rbindlist(lapply(1:length(ogs), function(j){
      genes <- ogs[[j]]
      return(f[f$V1  %in% genes & f$V2 %in% genes,])
    }))

    if (verbose)
      cat("retained", nrow(fl), "hits\n\t")
    write.table(fl,
                sep = "\t",
                row.names = F,
                col.names = F,
                quote = F,
                file = file.path(out.dir,xfile))
  }

  system(paste("cp",
               file.path(blast.dir, "SequenceIDs.txt"),
               out.dir))
  system(paste("cp",
               file.path(blast.dir, "SpeciesIDs.txt"),
               out.dir))

  system(paste("cp",
               file.path(blast.dir, "diamondDB*.dmnd"),
               out.dir))
  system(paste("cp",
               file.path(blast.dir, "Species*.fa"),
               out.dir))

  if (verbose)
    cat("Done!\n\t")
}

#' @title Add orthofinder info to gff
#' @description
#' \code{get_ofIDs} Add orthofinder gene ID numbers to gff.
#' @rdname of_utilities
#' @import data.table
#' @export
get_ofIDs <- function(ogff,
                      of.geneIDs,
                      of.speciesIDs,
                      verbose = T){
  if (verbose)
    cat("Adding orthofinder IDs to gffs ... ")
  of.geneIDs <- data.table(of.geneIDs)
  setkey(of.geneIDs, "id")
  of.speciesIDs <- data.table(of.speciesIDs)
  setkey(of.speciesIDs, "genome")

  ofg.out <- lapply(ogff, function(x){
    setkey(x, "genome")
    x1 <- merge(of.speciesIDs, x)
    setkey(x1, "id")
    x2 <- merge(of.geneIDs, x1)
    return(x2)
  })
  names(ofg.out) <- names(ogff)
  if (verbose)
    cat("Done!\n")
  return(ofg.out)
}
