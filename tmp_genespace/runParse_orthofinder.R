#' @title Pipe for initial orthofinder run
#'
#' @description
#' \code{runParse_orthofinder} Builds a table of best hits from reciprocal blasts
#'
#' @param peptide.dir The path to the directory containing the peptide fasta sequence files
#' @param blast.dir The path to the directory where the blast results should be stored
#' @param tmp.dir The path to the directory where temporary files will be stored then deleted
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
#' @param runOF Logical, should the full orthofinder pipe be run,
#' or should blast results just be parsed?
#' @param fasta.pattern Character, the character string to grep for in the fasta directory.
#' Useful if the fasta files are indexed.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Runs the following pipeline:
#' \enumerate{
#'   \item Run reciprocal BLASTs via Diamond within orthofinder (if `runOrtho`)
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
runParse_orthofinder<-function(peptide.dir,
                               tmp.dir,
                               blast.dir,
                               gff.dir,
                               ploidy,
                               min.propMax = .5,
                               min.score = 50,
                               nmapsPerHaplotype = 1,
                               eps.radius = c(100,50,20),
                               n.mappingWithinRadius = c(10,10,10),
                               runOF = TRUE,
                               fasta.pattern = "*.fa",
                               verbose = T,
                               returnAllBlast = F,
                               ...){
  if(runOF){
    ################   ################   ################
    ################   ################   ################
    if(verbose)
      cat("Copying peptide fasta files to", tmp.dir,"\n")
    if(file.exists(tmp.dir)){
      system(paste("rm -r", tmp.dir))
    }
    system(paste("mkdir", tmp.dir))
    system(paste("cp",
                 file.path(peptide.dir,
                           fasta.pattern),
                 tmp.dir))

    ################   ################   ################
    ################   ################   ################
    if(verbose)
      cat("Running blasts within Orthofinder\n")
    system(paste("orthofinder -f",
                 tmp.dir,
                 "-S diamond -og"))

    ################   ################   ################
    ################   ################   ################
    if(verbose)
      cat("Moving blasts results to", blast.dir,"\n")
    if(file.exists(blast.dir)){
      system(paste("rm -r", blast.dir))
    }
    system(paste("mkdir", blast.dir))

    blast.loc <- dirname(list.files(tmp.dir,
                                    pattern = "SequenceIDs",
                                    recursive = T,
                                    full.names = T)[1])
    ortho.loc <- dirname(list.files(tmp.dir,
                                    pattern = "Orthogroups.txt",
                                    recursive = T,
                                    full.names = T)[1])

    system(paste("cp",
                 file.path(blast.loc, "SequenceIDs.txt"),
                 blast.dir))
    system(paste("cp",
                 file.path(blast.loc, "SpeciesIDs.txt"),
                 blast.dir))
    system(paste("cp",
                 file.path(blast.loc, "Blast*"),
                 blast.dir))
    system(paste("cp",
                 file.path(ortho.loc, "Orthogroups.txt"),
                 blast.dir))
  }


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
    g <- g[g$V3 == "gene",c(9,1,4,5,7)]
    g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                          strsplit(x,";")[[1]][2]))
    data.table::setnames(g, c("id","chr","start","end","strand"))
    return(g)
  }

  gff <- rbindlist(lapply(names(gff.files), function(i){
    tmp <- parse_gff(gff.files[[i]])
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr","start")],
                       ties.method = "random")
    return(tmp)
  }))

  ################   ################   ################
  ################   ################   ################
  if(verbose)
    cat("Parsing orthofinder blast results\n")

  si <- read.delim(file.path(blast.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num", "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  rownames(si) <- si$genome
  si <- si[names(ploidy),]
  si$ploidy <- ploidy

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

  gz <- list.files(blast.dir,
                   pattern = ".gz$")
  if(length(gz)>0){
    system(paste("gunzip",file.path(blast.dir,"*.gz")))
  }
  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = names(ploidy)))
  sl <- split(sm,sm$unique)

  sequence.index <- fread(file.path(blast.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))

  og <- readLines(file.path(blast.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x," ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])
  ng <- sapply(og, length)

  if(any(file.exists(file.path(blast.dir,"*gz"))))
    system(paste("gunzip",
                 file.path(blast.dir,"*gz")))

  genes.in.orthogroups <- unlist(og[ng>1])
  sequence.index <- sequence.index[sequence.index$id %in% genes.in.orthogroups,]
  si1 <- data.table(sequence.index)
  si2 <- data.table(sequence.index)
  setnames(si1, c("gn1", "id1"))
  setnames(si2, c("gn2", "id2"))

  a <- og
  al <- unlist(a)
  g <- rep(seq_along(a),
           sapply(a, length))

  gff.of <- gff
  b <- gff.of$id
  gff.of$og <- g[match(b, al)]
  gff.o <- gff.of[,c(1:4,6,8)]
  setnames(gff1 <- data.table(gff.o), paste0(names(gff.o),"1"))
  setnames(gff2 <- data.table(gff.o), paste0(names(gff.o),"2"))

  loop_dbs <- function(map,
                       radii,
                       mappings){
    run_dbs <- function(map,
                        eps_radius,
                        mappings){
      x <- data.frame(map)
      x$x_a <- frank(x[,c("chr1","start1")],
                     ties.method = "dense")
      x$x_b <- frank(x[,c("chr2","start2")],
                     ties.method = "dense")

      nn <- frNN(x[,c("x_a","x_b")],
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
    }
    return(map)
  }

  bl <- lapply(sl, function(i){
    i <- i[order(i$map.rank),]
    suppressWarnings(b1 <- fread(i$filename[1], showProgress = F))
    suppressWarnings(b2 <- fread(i$filename[2], showProgress = F))
    b2 <- data.table(b2[,c(2,1,3:6,9:10,7:8,11:12)])
    setnames(b2, colnames(b1))
    d <- rbind(b1, b2)
    n <- nrow(d)
    setnames(d, c("gn1","gn2","perc.iden","align.length","n.mismatch",
                  "n.gapOpen","q.start","q.end","s.start","s.end","eval","score"))
    d <- merge(si1,
               merge(si2, d, by = "gn2"),
               by = "gn1")
    d$gn1 <- NULL
    d$gn2 <- NULL
    d <- merge(gff1,
               merge(gff2,  d, by = "id2"),
               by = "id1")
    dallout = d
    d <- d[d$og1 == d$og2,]

    d <- d[d$score >= min.score,]
    setkey(d, score)
    d <- d[!duplicated(d[,c("id1","id2"), with = F]),]

    d <- d[order(d$id1, d$id2, -d$score), ]

    subject.ploidy <- ploidy[d$genome2[1]]
    query.ploidy <- ploidy[d$genome1[1]]

    d[, rank := frank(score, ties.method = "dense"),
      by = list(id1)]
    d <- d[d$rank <= subject.ploidy*nmapsPerHaplotype,]

    d[, rank := frank(score, ties.method = "dense"),
      by = list(id2)]
    d <- d[d$rank <= query.ploidy*nmapsPerHaplotype,]

    d[, prop := score/max(score),
      by = list(id2)]
    d <- d[d$prop>=min.propMax,]
    d[, prop := score/max(score),
      by = list(id1)]
    d <- d[d$prop>=min.propMax,]

    d$rank <- NULL
    d$prop <- NULL

    d <- loop_dbs(map = d,
                  radii = eps.radius,
                  mappings = n.mappingWithinRadius)


    if(returnAllBlast){
      if(verbose)
        cat("Done",
            d$genome1[1],"vs.",d$genome2[1],
            ": initial hits =",n, "culled to",
            nrow(dallout),"\n")
      return(dallout)
    }else{
      if(verbose)
        cat("Done",
            d$genome1[1],"vs.",d$genome2[1],
            ": initial hits =",n, "culled to",
            nrow(d),"\n")
      return(d)
    }

  })
  out <- rbindlist(bl)

  if(verbose)
    cat("Writing culled blast output to",
        results.dir, "\n")
  if(!file.exists(results.dir))
    system(paste("mkdir", results.dir))
  write.csv(out,
            file = file.path(results.dir, "culledblasts.csv"),
            row.names = F)

  return(out)
}
