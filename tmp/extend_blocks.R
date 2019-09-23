#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{extend_blocks} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param map The map data.frame or data.table
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id':gene identifier, 'chr',
#' 'start', 'end', 'strand', 'genome': matching an element in genomeIDs,
#' 'order': gene order within that genome.
#' @param blast the blast dataset to screen for syntenic hits
#' @param n.iter Number of iterations to run
#' @param rank.buffer The buffer, in gene rank order.
#' @param radius Numeric, length of 1, specifiying the
#' radius to search for hits. Passed to clean_blocks
#' @param n.mappings Numeric, length of 1, specifiying the
#' minimum number of hits in a radius. Passed to clean_blocks
#' @param clean.it Should the cleaning step be run?
#' @param verbose Logical, should updates be printed
#' @param ... Additional arguments passed to cull_syntenicBlast
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
extend_blocks <- function(map,
                          dir.list,
                          gff,
                          genomeIDs,
                          radius = 100,
                          min.blockSize = MCScanX.s.param,
                          MCScanX.s.param,
                          MCScanX.m.param,
                          MCScanX.path,
                          verbose = T){

  if (verbose)
    cat("Loading all blast files into memory ... ")
  blast <- read_allBlasts(
    gff = gff,
    keep.geneNum = F,
    add.gff = T,
    check.ogs = F,
    blast.dir = dir.list$cull.score.blast,
    genomeIDs = genomeIDs,
    verbose = F)
  if(verbose)
    cat("Done!\nMerging gff, blast and map data.tables ... ")
  # blast <- data.table(synblks$blast)
  map <- simplify_map(
    map,
    gff = gff,
    genomeIDs = genomeIDs,
    mirror = F,
    just.rank = F)$map

  bl <- merge(map[,c("block.id","genome1","genome2","id1","id2")],
              blast,
              all = T,
              by = c("genome1","genome2","id1","id2"))
  bl <- subset(bl, !is.na(chr1))
  blo <- simplify_map(
    map = bl,
    gff = gff,
    genomeIDs = genomeIDs,
    mirror = F,
    just.rank = T)
  newblk <- subset(blo, !is.na(block.id))
  newblk <- newblk[,list(
    rankstart1 = min(rank1, na.rm = T),
    rankend1 = max(rank1, na.rm = T),
    rankstart2 = min(rank2, na.rm = T),
    rankend2 = max(rank2, na.rm = T),
    n.mapping = .N),
    by = list(genome1, genome2, chr1, chr2, block.id)]
  spl.blk <- split(newblk, by = c("genome1","genome2","chr1","chr2"))
  spl.blast <- split(blo, by = c("genome1","genome2","chr1","chr2"))

  if (verbose)
    cat("Done!\nCulling blast to syntenic hits ... ")

  blast.out <- rbindlist(lapply(names(spl.blk), function(i){
    y <- spl.blk[[i]]
    x <- spl.blast[[i]]
    setkey(y, n.mapping)
    xc <- rbindlist(lapply(1:nrow(y), function(j){

      int1 <- findInterval(
        x$rank1,
        c(y$rankstart1[j] - radius,
          y$rankend1[j] + radius)) == 1
      int2 <- findInterval(
        x$rank2,
        c(y$rankstart2[j] - radius,
          y$rankend2[j] + radius)) == 1
      wh <- which(int1 & int2)
      if (length(wh) < 1) {
        return(NULL)
      }else{

        tmp <- x[wh,]
        wh2 <- find_whichInBuffer(
          x = tmp$rank1,
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
  if (verbose)
    cat("Done!\nRe-forming syntenic blocks with MCScanX ... \n")
  blast.out <- blast.out[!duplicated(blast.out[,c("id1","id2")]),]
  bl.cull <- merge(blast, blast.out[,c("id1","id2")], by = c("id1","id2"))

  tpg <- data.table(t(combn(genomeIDs, 2)))
  setnames(tpg, c("genome1","genome2"))
  tpg <- rbind(tpg,
               data.table(genome1 = genomeIDs,
                          genome2 = genomeIDs))
  bl.cull <- merge(bl.cull, tpg, by = c("genome1","genome2"))
  syn.out <- pipe_mcscanx(
    blast = bl.cull,
    gff = gff,
    dir.list = dir.list,
    genomeIDs = genomeIDs,
    MCScanX.path = MCScanX.path,
    MCScanX.s.param = MCScanX.s.param,
    MCScanX.m.param = MCScanX.m.param,
    silent.mcs = T,
    verbose = T)
  if (verbose)
    cat("Re-formatting map and block data.tables ... ")
  simp.out <- simplify_map(
    map = syn.out$map,
    gff = gff,
    genomeIDs = genomeIDs,
    mirror = F)
  if (verbose)
    cat("Done!\n")
  return(simp.out)
}

