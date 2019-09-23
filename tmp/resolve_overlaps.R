split_blksByRef <- function(map, genomeIDs){
  mi1 <- subset(map, genome1 == genomeIDs[1])
  mi2 <- subset(map, genome2 == genomeIDs[1])
  mio1 <- mi1[,list(u = unique(id2)),
              by = list(genome1, block.id, id1)]
  mio2 <- mi2[,list(u = unique(id1)),
              by = list(genome2, block.id, id2)]

  tpg <- data.table(t(combn(genomeIDs, 2)))
  setnames(tpg, c("genome1","genome2"))
  tpg <- tpg[!duplicated(tpg$genome1),]
  tpm <- merge(map, tpg, by = c("genome1","genome2"))

  spl <- split(tpm, by = c("genome1","genome2"))
  for(i in 1:nrow(tpg)){
    y <- spl.gff[[paste0(x$genome1[1], ".", x$chr1[1])]]
  }

}

resolve_overlaps <- function(map, genomeIDs, gff, max.gap){

  # 1. Split map by genome and chr1
  spl.map <- split(map, by = c("genome1","genome2","chr1"))
  spl.gff <- split(gff, by = c("genome","chr"))

  x <- spl.map[["PhalliiHAL.Sbicolor.Chr03"]]
  y <- spl.gff[[paste0(x$genome1[1], ".", x$chr1[1])]]
  y <- subset(y, y$id %in% x$id1)
  y$rank <- frankv(y, cols = c("start", "end"), ties.method = "dense")

  z <- merge(y[,c("id","rank")],
             data.table(id = x$id1,
                        block.id = x$block.id),
             by = "id")

  zb <- z[,list(start = min(rank),
                end = max(rank)),
          by = "block.id"]
  setkey(zb, start, end)

  ovl <- foverlaps(zb, zb, type="any", which=TRUE)
  ovl <- subset(ovl, xid < yid)
  if(nrow(ovl) > 0){
    novl <- apply(ovl, 1, function(j){
      b1 <- zb$block.id[j[1]]
      b2 <- zb$block.id[j[2]]
      inter <- intersect(z$id[z$block.id == b1],
                         z$id[z$block.id == b2])
      return(length(inter))
    })

  }


  # 2. Check if there are any duplicate hits between blocks
  x.nodup <- x[!duplicated(x[,c("id1","block.id")]),]
  x.dup <- x.nodup$id1[duplicated(x.nodup$id1)]
  dup.map <- subset(x.nodup, id1 %in% x.dup)
  dup.blk.list <- split(dup.map$block.id, dup.map$id1)

  # 3. Check if any overlap occurs, rank from smallest to largest
  int1 <- findInterval(x$rank1,
                       c(y$rankstart1[j] - blk.rank.buffer,
                         y$rankend1[j] + blk.rank.buffer)) == 1
  int2 <- findInterval(x$rank2,
                       c(y$rankstart2[j] - blk.rank.buffer,
                         y$rankend2[j] + blk.rank.buffer)) == 1

  # ... For non-duplicated overlapping blocks
  # 4. Calculate RLE
  # 5. Any RLE < max.gap -> drop hits
  # 6. Re-calculate RLE
  # 7. Any RLE >= max.gap, split into new blocks
}




# 1. Resolve overlapping blocks with interleaved hits
# -- for each genome and
# 2. For a given block, what are all the genes adjacent doing?
# -- if going to a different block, stop
# -- if no strong blast hit, drop from dataset, re-check
# -- if has a strong blast hit near given block, add to block, re-check
# -- repeat until bounding hits go to a different block

# 2.
