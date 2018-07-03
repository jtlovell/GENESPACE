plot_mapping = function(cols = NULL, blk, map, genomes = NULL){
  if(is.null(cols))
    cols = rep_len(c("red3","salmon","darkorange3","gold",
                     "grey50","lightgreen","forestgreen","darkgreen",
                     "cyan","dodgerblue3","violet","purple"), length.out = nrow(blk))
  map = map[order(map$block.id),]
  combn = map[!duplicated(map[,c("genome1","genome2")]),c("genome1","genome2")]
  if(!is.null(genomes)){
    is.in = colSums(apply(combn,1,function(x) x %in% genomes))
    combn = combn[is.in == 2,]
  }
  trsh = apply(combn, 1, function(x){
    tmp = map[map$genome1 == x[1] &
                map$genome2 == x[2],]
    tmp$rank1 = frank(tmp[,c("chr1","start1")], ties.method = "dense")
    tmp$rank2 = frank(tmp[,c("chr2","start2")], ties.method = "dense")
    tmp$col = cols[as.numeric(as.factor(tmp$block.id))]
    with(tmp,
         plot(rank1, rank2,
              col = col,
              pch = ".",
              ylab = paste(x[2],"gene order"),
              xlab = paste(x[1],"gene order"),
              axes = F, bty = "n", asp = 1,
              main = paste(x[1],"vs", x[2])))

    lab1 = tapply(tmp$rank1,tmp$chr1, mean)
    end1 = tapply(tmp$rank1,tmp$chr1, max)[-length(lab1)]

    lab2 = tapply(tmp$rank2,tmp$chr2, mean)
    end2 = tapply(tmp$rank2,tmp$chr2, max)[-length(lab2)]

    axis(1, at = lab1, labels = names(lab1))
    axis(2, at = lab2, labels = names(lab2))
    segments(x0 = end1, x1 = end1, y0 = min(tmp$rank2), y1 = max(tmp$rank2),
             col = "lightgrey", lty = 2)
    segments(y0 = end2, y1 = end2, x0 = min(tmp$rank1), x1 = max(tmp$rank1),
             col = "lightgrey", lty = 2)
  })
}
