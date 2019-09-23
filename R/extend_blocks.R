#' @title Extend and complete syntenic mappings
#'
#' @description
#' \code{extend_blocks} Extend and complete syntenic mappings.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL default, taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param use.score.cull.blast logical, should the blasts be from the
#' score.cull directory, or just the genome.cull directory?
#' @param ... Not currently in use
#' @param rank.buffer numeric, the radius to search within for syntenic
#' mappings.
#'
#' @details ...
#'
#' @return A blast data.table, without block.id information
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
extend_blocks <- function(gff,
                          genomeIDs,
                          dir.list,
                          map,
                          verbose = T,
                          use.score.cull.blast = T,
                          rank.buffer){
  rerank_mapFromGff <- function(gff, map, genomeIDs){
    gff <- subset(gff, genome %in% genomeIDs)
    spl.map <- split(map, by = c("genome1","genome2","chr1"))
    spl.gff <- split(gff, by = c("genome","chr"))
    mo1 <- rbindlist(lapply(spl.map, function(x){
      gc <- paste(x$genome1[1], x$chr1[1], sep = ".")
      y <- subset(spl.gff[[gc]], id %in% x$id1)
      x <- x[,c("genome1","genome2","chr1","chr2","id1","id2")]
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
    setcolorder(map, c(3:6,2,1,7:8,11,9,12,10))
    setkey(map, chr1, chr2, start1, start2)
    return(map)
  }

  if (use.score.cull.blast) {
    of.dir <- dir.list$cull.score.blast
  }else{
    of.dir <- dir.list$cull.blast
  }

  all.blast <- read_allBlasts(
    gff = gff,
    genomeIDs = genomeIDs,
    of.dir = of.dir,
    check.ogs = F,
    add.gff = T,
    keep.geneNum = F,
    verbose = verbose)

  comb <- data.table(rbind(t(combn(genomeIDs, 2, simplify = T)),
                           t(sapply(genomeIDs, rep, 2))))
  setnames(comb, c("genome1", "genome2"))
  map <- merge(map, comb, by = c("genome1", "genome2"))

  if(verbose)
    cat("Subsetting blast hits to chromosomes found in syntenic blocks ... ")
  all.blast[,what := "blast"]
  map[,what := "map"]
  cn <- colnames(all.blast)[colnames(all.blast) %in% colnames(map)]
  spl.map <- split(map, by = c("genome1","genome2","chr1","chr2"))
  spl.blast <- split(all.blast, by = c("genome1","genome2","chr1","chr2"))
  blast <- rbindlist(spl.blast[names(spl.map)])

  z <- rbind(map[, cn, with = F],
             blast[, cn, with = F])
  z <- z[!duplicated(z[,c("id1","id2")]),]

  if(verbose)
    cat("Done!\nReranking blast hits by position in gff ... ")
  zo <- rerank_mapFromGff(
    map = z,
    gff = gff,
    genomeIDs = genomeIDs)

  zo <- merge(zo, z[,c("id1","id2","what")], by = c("id1","id2"))
  if(verbose)
    cat("Done!\nSearching for syntenic hits within", length(spl.map),"pairwise chromosome combinations ...\n")
  spl.map <- split(zo, by = c("genome1","genome2","chr1","chr2"))
  out <- rbindlist(lapply(names(spl.map), function(i){
    if(which(names(spl.map) == i) %% 100 == 0)
      cat("\tCompleted",which(names(spl.map) == i),"/",length(spl.map),"\n")
    z <- spl.map[[i]]
    wh <- find_whichInBuffer(x = z$rank1,
                             y = z$rank2,
                             which.in.blk = which(z$what == "map"),
                             rank.buffer = rank.buffer)
    return(data.table(z[wh,]))
  }))
  if(verbose)
    cat("\tDone!\n")
  out <- merge(all.blast, out[,c("genome1","genome2","id1","id2")], by = c("genome1","genome2","id1","id2"))
  out[,what := NULL]
  return(out)
}
