merge_blocks = function(blk, map, buffer = 3, verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")


  blk1 = blk[,c(1:5,10:13)]
  blk2 = blk[,c(1:5,10:13)]
  colnames(blk1)[c(1,6:9)]<-paste0(colnames(blk1)[c(1,6:9)],"_a")
  colnames(blk2)[c(1,6:9)]<-paste0(colnames(blk2)[c(1,6:9)],"_b")
  m = merge(blk1, blk2,by = colnames(blk1)[2:5])
  m = m[!duplicated(m),]
  m = m[m$block.id_a != m$block.id_b,]
  m$unique = apply(m[,c("block.id_a","block.id_b")],1,function(x) paste(x[order(x)], collapse = "_"))
  m = m[!duplicated(m$unique),]
  m$leftin1 = with(m,
                   (rankstart1_a >= rankstart1_b - buffer &
                      rankstart1_a <= rankend1_b + buffer))
  m$bottomin1 = with(m,
                     (rankstart2_a >= rankstart2_b - buffer &
                        rankstart2_a <= rankend2_b + buffer))
  m$leftin2 = with(m,
                   (rankstart1_b >= rankstart1_a - buffer &
                      rankstart1_b <= rankend1_a + buffer))
  m$bottomin2 = with(m,
                     (rankstart2_b >= rankstart2_a - buffer &
                        rankstart2_b <= rankend2_a + buffer))
  m$rightin1 = with(m,
                    (rankend1_a <= rankend1_b + buffer &
                       rankend1_a >= rankstart1_b - buffer))
  m$topin1 = with(m,
                  (rankend2_a <= rankend2_b + buffer &
                     rankend2_a >= rankstart2_b - buffer))
  m$rightin2 = with(m,
                    (rankend1_b <= rankend1_a + buffer &
                       rankend1_b >= rankstart1_a - buffer))
  m$topin2 = with(m,
                  (rankend2_b <= rankend2_a + buffer &
                     rankend2_b >= rankstart2_a - buffer))
  tomerge = m[rowSums(m[16:23]) >= 4,]
  tomerge$lowblock = apply(tomerge[,c("block.id_a","block.id_b")],1,min)
  tomerge$highblock = apply(tomerge[,c("block.id_a","block.id_b")],1,max)
  for(i in 1:nrow(tomerge))
    map$block.id[map$block.id == tomerge$highblock[i]]<-tomerge$lowblock[i]

  map$mapping = paste(map$genome1, map$genome2)
  map$rank1 = frank(map[,c("mapping","chr1","start1")], ties.method = "dense")
  map$rank2 = frank(map[,c("mapping","chr2","start2")], ties.method = "dense")

  out = make_blocks(map)

  map = data.frame(out[["map"]], stringsAsFactors = F)
  blk = data.frame(out[["block"]], stringsAsFactors = F)
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}
