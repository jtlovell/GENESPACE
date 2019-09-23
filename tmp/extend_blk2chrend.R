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


check_brkpts <- function(map,
                         gff,
                         verbose = T,
                         genomeIDs,
                         ploidy,
                         n.cores = 1){

  map <- simplify_map(map = map,
                      gff = gff,
                      genomeIDs = genomeIDs,
                      mirror = F)

  gffu <- rbindlist(lapply(genomeIDs, function(i) {
    tmp <- data.table(gff)
    tmp[,genome1 := genome]
    tmp[,genome2 := i]
    return(tmp)
  }))

  blk.assign <- assign_gffBlks(map = map$map,
                               gff = gff,
                               verbose = verbose,
                               n.cores = n.cores)
  cat("\tDone!\nParsing block assignments and merging with gff ...")

  blk.assign2 <- blk.assign[,c("block.id","genome2","genome1",
                               "chr2","chr1","start2","end2","start1","end1",
                               "id","chr","start","end","strand","genome","order","inmap")]
  setnames(blk.assign2, colnames(blk.assign))
  blk.assign <- rbind(blk.assign, blk.assign2)
  blk.assign <- subset(blk.assign, genome == genome1 & chr == chr1)
  blk.assign <- blk.assign[!duplicated(blk.assign[,c("id","genome1","genome2","chr1","chr2","start1","start2")])]
  blk.assign[,n.blks := length(unique(block.id)),
             by = list(genome1, genome2, genome, chr1, id)]
  ba <- blk.assign[,c("genome","id","genome1","genome2","n.blks")]
  ba <- ba[!duplicated(ba),]
  go <- merge(gffu, ba, all.x = T, by = c("genome","genome1","genome2","id"))
  go$n.blks[is.na(go$n.blks)] <- 0
  if(verbose)
    cat("\tDone!")

  blk.assign[,n.exp.hits := ploidy[genome]/2]
  with(blk.assign, table(n.blks, n.exp.hits))
}





