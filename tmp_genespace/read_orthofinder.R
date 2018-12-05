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
read_orthofinder <- function(blast.dir,
                             gff,
                             parse_genomeID.FUN = function(x) strsplit(x,"_")[[1]][3],
                             verbose = TRUE,
                             ...){

  if(verbose)
    cat("Reading blasts from", blast.dir,"... ")
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


  si <- read.delim(file.path(blast.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- sapply(gsub(".fa$", "", si$genome),parse_genomeID.FUN)

  sm <- cbind(expand.grid(si$genome.num,
                          si$genome.num,
                          stringsAsFactors = F),
              expand.grid(si$genome,
                          si$genome,
                          stringsAsFactors = F))
  colnames(sm) <- c("n1","n2","genome1","genome2")



  sm$filename <- file.path(blast.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$genome1, sm$genome2)
  sl <- split(sm,sm$unique)

  ################   ################   ################
  ################   ################   ################

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
  setkey(gff, genome,id)
  gff.in = merge(gff, si2)

  og.info = gff.in[,list(n.unique.genes = length(unique(id)),
               n.hits = length(id),
               n.genomes = length(unique(genome))),
         by = list(og)]

  gff.in = merge(gff.in, og.info, by = "og")


  ################   ################   ################
  ################   ################   ################


  ################   ################   ################
  ################   ################   ################

  blast.out <- lapply(1:length(sl), function(i){
    x <- sl[[i]]
    suppressWarnings(b1 <- fread(x$filename[1],
                                 showProgress = F))
    if(file.exists(x$filename[2])){
      suppressWarnings(b2 <- fread(x$filename[2],
                                   showProgress = F))
      b2 <- data.table(b2[,c(2, 1, 3:6, 9:10, 7:8, 11:12)])
      setnames(b2, colnames(b1))
      blast.in <- rbind(b1, b2)
    }else{
      blast.in <- b1
    }
    setnames(blast.in, c("gn1", "gn2", "perc.iden",
                         "align.length", "n.mismatch",
                         "n.gapOpen", "q.start", "q.end",
                         "s.start", "s.end",
                         "eval", "score"))

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
    d<-merged
    d$negscore <- d$score * (-1)
    setkey(d, id1, id2, negscore)
    d <- data.table(d[d$id1 != d$id2,])
    d <- data.table(d[!duplicated(d[,c("id1", "id2"), with = F]),])

    d$negscore <- d$score * (-1)
    setkey(d, id1, id2, negscore)
    d <- data.table(d[!duplicated(d[,c("id1", "id2"), with = F]),])

    return(d)
  })
  ret <- rbindlist(blast.out)
  ret$rank1 <- frank(ret, chr1, start1, ties.method = "dense")
  ret$rank2 <- frank(ret, chr2, start2, ties.method = "dense")
  ret <- ret[,c("genome1", "id1", "chr1", "start1", "end1", "rank1",
                "genome2", "id2", "chr2", "start2", "end2", "rank2",
                "og1","n.unique.genes1", "n.hits1", "n.genomes1",
                "perc.iden", "align.length", "n.mismatch", "n.gapOpen",
                "q.start", "q.end", "s.start", "s.end",
                "eval", "score"), with = F]
  if(verbose)
    cat("found the following number of hits:\n",
        paste(with(ret, table(genome1)),collapse = ", "),
        paste0("(",paste(names(with(ret, table(genome1))),collapse = ", "),")"),"\n")

  ret$unique = paste(ret$genome1, ret$genome2)
  same = ret[ret$genome1 == ret$genome2,]
  dif = ret[ret$genome1 != ret$genome2,]
  if(nrow(dif)>0){
    spl  =split(dif, dif$unique)
    difo<-rbindlist(lapply(spl, function(y){
      wh1 = which(genomeIDs == y$genome1[1])
      wh2 = which(genomeIDs == y$genome2[1])
      if(wh1>wh2){
        y<-y[,c("genome2", "id2", "chr2", "start2", "end2", "rank2",
                "genome1", "id1", "chr1", "start1", "end1", "rank1",
                "og1","n.unique.genes1", "n.hits1", "n.genomes1",
                "perc.iden", "align.length", "n.mismatch", "n.gapOpen",
                "q.start", "q.end", "s.start", "s.end",
                "eval", "score","unique")]
        setnames(y,colnames(dif))
      }
      return(y)
    }))
  }
  out = rbind(same, difo)
  setnames(out,c("genome1", "id1", "chr1", "start1", "end1", "rank1",
                 "genome2", "id2", "chr2", "start2", "end2", "rank2","og",
                 "n.unique.genes", "n.hits", "n.genomes",
                 "perc.iden", "align.length", "n.mismatch", "n.gapOpen",
                 "q.start", "q.end", "s.start", "s.end",
                 "eval", "score","unique"))
  out$negscore = out$score * (-1)
  setkey(out, id1, id2, negscore)
  out = out[!duplicated(out[,c("id1","id2"), with = F]),]
  return(list(blast = out, ortho = gff.in))
}
