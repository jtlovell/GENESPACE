simplify_map <- function(map, gff, genomeIDs){
  map <- subset(map, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
  gff <- subset(gff, genome %in% genomeIDs)
  spl.map <- split(map, by = c("genome1","genome2","chr1"))
  spl.gff <- split(gff, by = c("genome","chr"))
  mo1 <- rbindlist(lapply(spl.map, function(x){
    gc <- paste(x$genome1[1], x$chr1[1], sep = ".")
    y <- subset(spl.gff[[gc]], id %in% x$id1)
    x <- x[,c("genome1","genome2","chr1","chr2","id1","id2","block.id")]
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
  setcolorder(map, c(3:6,2,1,7:8,11,9,12,10,13))
  setkey(map, chr1, chr2, start1, start2)
  return(map)
}
