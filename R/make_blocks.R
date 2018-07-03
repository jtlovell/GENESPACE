make_blocks<-function(map){
  out.blk = with(map,
                 data.table(
                   block.id = as.character(tapply(block.id, block.id, function(x) x[1])),
                   genome1 = tapply(genome1, block.id, function(x) x[1]),
                   genome2 = tapply(genome2, block.id, function(x) x[1]),
                   chr1 = tapply(chr1, block.id, function(x) x[1]),
                   chr2 = tapply(chr2, block.id, function(x) x[1]),
                   start1 = tapply(start1, block.id, min),
                   start2 = tapply(start2, block.id, min),
                   end1 = tapply(end1, block.id, max),
                   end2 = tapply(end2, block.id, max),
                   rankstart1 = tapply(rank1, block.id, min),
                   rankstart2 = tapply(rank2, block.id, min),
                   rankend1 = tapply(rank1, block.id, max),
                   rankend2 = tapply(rank2, block.id, max),
                   medianscore = tapply(score, block.id, median),
                   n.mapping = tapply(score, block.id, length),
                   stringsAsFactors = F))
  out.blk$density = with(out.blk,
                         n.mapping/((abs(rankend1-rankstart1)+abs(rankend2-rankstart2))/2))

  orient = sapply(split(map, map$block.id), function(x) cor(x$start1, x$start2))
  orient = data.table(block.id = names(orient), orient = ifelse(orient>0,"+","-"))
  out.blk = merge(out.blk, orient, by = "block.id")

  map = data.frame(map[order(map$genome1, map$chr1, map$start1),], stringsAsFactors = F)
  blk = data.frame(out.blk, stringsAsFactors = F)
  return(list(block = blk, map = map))
}
