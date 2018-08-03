find_newBlocks = function(blk, map, blast.results,
                          min.genes = 10, verbose = T,
                          min.unique.genes = 5){

  run_dbs = function(map, eps_radius, mappings){
    x=data.frame(map)

    nn = dbscan::frNN(x[,c("rank1","rank2")], eps = eps_radius)
    dbs = dbscan::dbscan(nn, minPts = mappings)
    return(data.frame(rank1 = x$rank1,
                      rank2 = x$rank2,
                      cluster = dbs$cluster,
                      stringsAsFactors = F))
  }

  blast.results$mapping =  paste(blast.results$genome1, blast.results$genome2)
  all_ass = merge(blast.results,
                  map[,c("id1","id2","block.id")],
                  all = T,
                  by = c("id1","id2"))

  spl = split(all_ass, paste(all_ass$genome1, all_ass$genome2))
  all_ass = rbindlist(lapply(spl, function(x){
    x = data.frame(x)
    x$rank1 = frank(x[,c("chr1","start1")], ties.method = "dense")
    x$rank2 = frank(x[,c("chr2","start2")], ties.method = "dense")
    return(x)
  }))
  no_ass = all_ass[is.na(all_ass$block.id),]

  no_assi = no_ass[!duplicated(no_ass[,c("genome1","genome2","chr1","chr2")]),]
  tab = data.frame(table(no_ass[,c("genome1","genome2","chr1","chr2")]),
                   stringsAsFactors = F)
  x = tab[tab$Freq > min.genes,]
  for(i in 1:4) x[,i]<-as.character(x[,i])
  for(i in 1:nrow(x)){

    to_ass = no_ass[with(no_ass,
                         genome1 == x$genome1[i] &
                           genome2 == x$genome2[i] &
                           chr1 == x$chr1[i] &
                           chr2 == x$chr2[i]),]


    to_ass$clus = run_dbs(to_ass, eps_radius = 20, mappings = 10)$cluster
    if(sum(to_ass$clus!=0)>=min.genes &
       length(unique(to_ass$id1))>=min.unique.genes &
       length(unique(to_ass$id2))>=min.unique.genes ){
      if(verbose) cat("added", sum(to_ass$clus!=0),"genes in a block for",
                      paste(x[i,1:4], collapse = ", "),"\n")
      to_ass$block.id = paste0("new",i)
      to_ass = to_ass[,colnames(map), with = F]
      map = rbind(map, to_ass)
    }
  }
  map$mapping = paste(map$genome1, map$genome2)
  map$rank1 = frank(map[,c("mapping","chr1","start1")], ties.method = "dense")
  map$rank2 = frank(map[,c("mapping","chr2","start2")], ties.method = "dense")
  map$block.id = as.character(as.numeric(as.factor(map$block.id)))
  out = make_blocks(map)
  map = out$map
  blk = out$block
  return(list(map = map, block = blk))
}
