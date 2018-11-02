#' @title Process orthofinder results
#'
#' @description
#' \code{parse_orthofinder} Cull orthofinder BLAST hits to those that are highest
#' confidence
#'
#' @param blast.dir The path to the directory where the blast results should be stored
#' @param gff.dir The path to the directory containing the annotations in gff3 format.
#' @param min.propMax Numeric, A dynamic blast score threshold to retain mappings, where
#' only mappings with scores > min.propMax of the highest score for each gene
#' @param min.score Numeric, the hard threshold, below which all mappings are removed.
#' @param nmapsPerHaplotype Numeric, the number of mappings retained per haplotype. For
#' example, if set to 1 (default), two mappings are retained for each gene in a diploid
#' genome, provided the scores satisfy the two score thresholds.
#' @param eps.radius Numeric, the gene position rank radius within which to count the number
#' of mappings. If length > 1, must be of same length as n.mappingWithinRadius
#' @param n.mappingWithinRadius Numeric, see eps.radius.
#' @param ploidy Named numeric vector, where names match genome ids. Informs the number
#' of mappings retained for each gene.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Runs the following pipeline:
#' \enumerate{
#'   \item Read in reciprocal BLASTs via Diamond within orthofinder
#'   \item Find hierarchical order of maps, from order of genomes in `ploidy`
#'   \item Parse blast hits to genes within multi-genome orthogroups
#'   \item Parse BLAST hits by score
#'   \item Parse BLAST hits by x-y positional density, where x and y gene order
#' }
#' @return A culled data.table with culled blast hits.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @export
parse_orthofinder <- function(blast.dir,
                              gff.dir,
                              verbose = TRUE,
                              ploidy,
                              min.propMax = .5,
                              min.score = 50,
                              nmapsPerHaplotype = 1,
                              eps.radius = c(100, 50, 20),
                              n.mappingWithinRadius = c(10, 10, 10),
                              ...){


  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if(length(gz) > 0){
    if(verbose)
      cat("Decompressing blast results\n")
    system(paste("gunzip",
                 file.path(blast.dir,
                           "*.gz")))
  }

  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Parsing species IDs\n")

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
  si <- si[genomeIDs,]

  sm <- cbind(expand.grid(si$genome.num,
                          si$genome.num,
                          stringsAsFactors = F),
              expand.grid(si$genome,
                          si$genome,
                          stringsAsFactors = F))
  colnames(sm) <- c("n1","n2","genome1","genome2")
  sm <- sm[sm$genome1 != sm$genome2,]
  sm$ref <- NA
  for(i in rev(names(ploidy)))
    sm$ref[sm$genome1 == i | sm$genome2 == i] <- i
  sm$alt <- ifelse(sm$genome1 == sm$ref, sm$genome2, sm$genome1)


  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = genomeIDs))
  sm <- sm[order(sm$map.rank),]
  sl <- split(sm,sm$unique)

  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Parsing gene IDs\n")
  sequence.index <- fread(file.path(blast.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))

  og <- readLines(file.path(blast.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])
  ng <- sapply(og, length)


  genes.in.orthogroups <- unlist(og[ng > 1])
  sequence.index <- sequence.index[sequence.index$id %in% genes.in.orthogroups,]
  ogs <- data.table(id = unlist(og),
                    og = rep(ons, ng),
                    stringsAsFactors = F)
  setkey(sequence.index, id)
  setkey(ogs, id)
  si1 <- merge(sequence.index, ogs)
  si2 <- data.table(si1)
  setnames(si1, c("id", "gn", "og"))
  si1$genome.num <- as.numeric(sapply(si1$gn, function(x)
    strsplit(x, "_")[[1]][1]))
  setkey(si1, genome.num)
  si.in <- data.table(si)
  setkey(si.in, genome.num)
  si2 <- merge(si1, si.in)
  setkey(si2, genome, id)

  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Parsing gff files\n")
  gff.files <- list.files(gff.dir,
                          full.names = T)
  names(gff.files) <- gsub(".gff3$", "",
                           basename(gff.files))

  parse_gff <- function(gff){
    g <- suppressWarnings(
      data.table::fread(gff,
                        showProgress = F,
                        verbose = F))
    g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
    g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                          strsplit(x, ";")[[1]][2]))
    data.table::setnames(g, c("id", "chr", "start", "end", "strand"))
    return(g)
  }

  gff <- rbindlist(lapply(names(gff.files), function(i){
    tmp <- parse_gff(gff.files[[i]])
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "random")
    return(tmp)
  }))
  setkey(gff, genome, id)

  gff.in = merge(gff, si2)


  loop_dbs <- function(map,
                       radii,
                       mappings,
                       verbose = T){
    run_dbs <- function(map,
                        eps_radius,
                        mappings){
      x <- data.frame(map)
      x$x_a <- frank(x[,c("chr1", "start1")],
                     ties.method = "dense")
      x$x_b <- frank(x[,c("chr2", "start2")],
                     ties.method = "dense")

      nn <- frNN(x[,c("x_a", "x_b")],
                 eps = eps_radius)
      dbs <- dbscan(nn,
                    minPts = mappings)
      return(data.frame(rank1 = x$x_a,
                        rank2 = x$x_b,
                        cluster = dbs$cluster,
                        stringsAsFactors = F))
    }
    if(length(radii) != length(mappings))
      stop("radii and mappings must be same length\n")
    for(i in 1:length(radii)){
      dclus <- run_dbs(map, eps = radii[i], mappings = mappings[i])
      mo <- cbind(map, data.table(dclus))
      map <- map[mo$cluster != 0,]
      if(verbose)
        cat("\t... retained",
            nrow(map),
            "where",
            mappings[i],
            "hits were found within",
            radii[i],
            "genes\n")
    }
    return(map)
  }

  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Parsing blast files ... \n")

  blast.out <- lapply(1:length(sl), function(i){
    x <- sl[[i]]
    if(verbose)
      cat("\t",
          x$genome1[which.min(x$n1)],
          "mapped against",
          x$genome1[which.max(x$n1)],
          ":\n")
    x <- x[order(x$map.rank),]
    suppressWarnings(b1 <- fread(x$filename[1],
                                 showProgress = F))
    suppressWarnings(b2 <- fread(x$filename[2],
                                 showProgress = F))
    b2 <- data.table(b2[,c(2, 1, 3:6, 9:10, 7:8, 11:12)])
    setnames(b2, colnames(b1))
    blast.in <- rbind(b1, b2)
    setnames(blast.in, c("gn1", "gn2", "perc.iden",
                         "align.length", "n.mismatch",
                         "n.gapOpen", "q.start", "q.end",
                         "s.start", "s.end",
                         "eval", "score"))

    if(verbose)
      cat("\t... read in",
          nrow(blast.in),
          "total mappings\n")

    gff1 <- data.table(gff.in)
    gff2 <- data.table(gff.in)
    setnames(gff1, paste0(colnames(gff1), "1"))
    setnames(gff2, paste0(colnames(gff2), "2"))
    setkey(gff1, gn1)
    setkey(blast.in, gn1)
    m1 <- merge(gff1, blast.in)
    setkey(gff2, gn2)
    setkey(m1, gn2)
    merged <- merge(gff2, m1)

    d <- data.table(merged[with(merged, og1 == og2),])
    if(verbose)
      cat("\t... culled to",
          nrow(d),
          "mappings in orthogroups\n")

    d$negscore <- d$score * (-1)
    setkey(d, id1, id2, negscore)
    d <- data.table(d[!duplicated(d[,c("id1", "id2"), with = F]),])


    if(verbose)
      cat("\t... culled to",
          nrow(d),
          "unique mappings\n")

    ploidy2 <- ploidy[d$genome2[1]]
    ploidy1 <- ploidy[d$genome1[1]]

    d[, rank1 := frank(score, ties.method = "dense"),
      by = list(id1)]
    d[, rank2 := frank(score, ties.method = "dense"),
      by = list(id2)]
    cullrank <- d[d$rank1 <= ploidy1*nmapsPerHaplotype |
                    d$rank2 <= ploidy2*nmapsPerHaplotype ,]

    cullrank[, prop2 := score/max(score),
             by = list(id2)]
    cullrank[, prop1 := score/max(score),
             by = list(id1)]

    cullscore <- cullrank[cullrank$prop1>=min.propMax |
                            cullrank$prop2>=min.propMax,]

    cullscore$rank1 <- NULL
    cullscore$prop1 <- NULL
    cullscore$rank2 <- NULL
    cullscore$prop2 <- NULL

    culldbs <- loop_dbs(map = cullscore,
                        radii = eps.radius,
                        mappings = n.mappingWithinRadius)


    return(culldbs)
  })

  if(verbose)
    cat("Parsing orthofinder results ... DONE!\n")
  ret <- rbindlist(blast.out)
  ret$rank1 <- frank(ret, chr1, start1)
  ret$rank2 <- frank(ret, chr2, start2)
  ret <- ret[,c("genome1", "id1", "chr1", "start1", "end1", "rank1",
                "genome2", "id2", "chr2", "start2", "end2", "rank2",
                "perc.iden", "align.length", "n.mismatch", "n.gapOpen",
                "q.start", "q.end", "s.start", "s.end",
                "eval", "score"), with = F]
  return(ret)
}
