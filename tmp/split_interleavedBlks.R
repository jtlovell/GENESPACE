split_interleavedBlks <- function(map, n.iter = 3){
  chk_it<- function(joind, what = "1"){
    if(what == "1"){
      chk <- data.table(joind$ovl.1.index)
    }else{
      chk <- data.table(joind$ovl.2.index)
    }
    if(nrow(chk) >= 1){
      ck <- apply(chk, 1, paste, collapse = "_")
      dup <- unique(c(apply(joind$dup.index, 1, paste, collapse = "_"),
                      apply(joind$dup.index[,c(2,1)], 1, paste, collapse = "_")))
      chk <- data.table(t(apply(chk[!ck %in% dup],1,function(x) x[order(x)])))
      smap <- split(joind$map, by = "block.id")
      spl <- apply(chk, 1, function(x)
        rbindlist(smap[x]))
      tmap <- joind$map
      for(i in 1:nrow(chk)){
        blks <- unlist(chk[i,])
        if(all(blks %in% tmap$block.id)){
          x <- subset(tmap, block.id %in% blks)
          if(nrow(x) > 1){
            if(what == "1"){
              setkey(x, start1)
            }else{
              setkey(x, start2)
            }
            x[,rl:=add_rlelen(block.id)]
            x[,rlid:=add_rlelenid(block.id)]
            x <- subset(x, rl >= 3)
            x[,block.id := paste0(block.id, ".", rlid)]
            x[,rl := NULL]
            x[,rlid := NULL]
            tmap <- rbind(x, subset(tmap, !block.id %in% blks))
          }
        }
      }
    }
    out <- make_blocks(tmap)
    return(out)
  }

  for(i in 1:n.iter){
    cat("Iteration",i,": ")

    map <- join_ovlpBlks(map = map, buffer = 0, max.iter = 1, verbose = F)


    map <- chk_it(map)
    map <- join_ovlpBlks(map = map$map, buffer = 0, max.iter = 1, verbose = F)
    map <- chk_it(map, what = "2")
    cat("split into", nrow(map$block),"\n")
    map <- map$map
  }
  return(make_blocks(map))
}
