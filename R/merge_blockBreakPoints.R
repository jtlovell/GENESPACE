merge_blockBreakPoints = function(map,
                                  genomeIDs,
                                  checkOvl.bp = 100000,
                                  checkOvl.rank = 3){

  blk.bygenome = do.call(rbind,lapply(genomeIDs, function(x){
    wh1 = which(map$genome1 == x)
    wh2 = which(map$genome2 == x)
    m = with(map,
             data.frame(genome = x,
                        alt.genome = c(genome2[wh1], genome1[wh2]),
                        chr = c(chr1[wh1], chr2[wh2]),
                        block.id = c(block.id[wh1],
                                     block.id[wh2]),
                        start = c(start1[wh1], start2[wh2]),
                        end = c(end1[wh1], end2[wh2]),
                        stringsAsFactors = F))
    m = data.table(m)

    setkey(m, chr, start)
    m[, rankstart := frank(start, ties.method = "dense"),
      by = list(chr)]
    setkey(m, chr, start)
    m[, rankend := frank(end, ties.method = "dense"),
      by = list(chr)]
    setkey(m, chr, block.id)

    b = m[,list(start=min(start),
                end=max(end),
                rankstart=min(rankstart),
                rankend=max(rankend)),
          by=list(genome, alt.genome, chr, block.id)]
    b$min.bp = ave(b$start, b$chr, FUN = min)
    b$max.bp = ave(b$end, b$chr, FUN = max)
    b$min.rank = ave(b$rankstart, b$chr, FUN = min)
    b$max.rank = ave(b$rankend, b$chr, FUN = max)

    setkey(b, chr, start)

    df = with(b, data.table(genome = c(genome, genome),
                            alt.genome = c(alt.genome, alt.genome),
                            chr = c(chr, chr),
                            block.id = c(block.id, block.id),
                            type = c(rep("start", length(start)),
                                     rep("end",length(end))),
                            bp = c(start, end),
                            rank = c(rankstart, rankend),
                            min.bp = c(min.bp, min.bp),
                            max.bp = c(max.bp,max.bp),
                            min.rank = c(min.rank,min.rank),
                            max.rank = c(max.rank,max.rank),
                            stringsAsFactors = F))
    df$chr.ends = with(df,
                       ifelse(bp == min.bp | rank == min.rank, "start",
                              ifelse(bp == max.bp | rank == max.rank, "end","middle")))
    return(data.frame(df))
  }))

  spl = split(blk.bygenome, paste(blk.bygenome$chr, blk.bygenome$genome))

  new.coords = rbindlist(lapply(spl, function(y){
    nn = with(y, dbscan::frNN(cbind(rank,rep(0,length(rank))), eps = checkOvl.rank))
    y$rank.group = dbscan::dbscan(nn, minPts = 2)$cluster
    nn = with(y, dbscan::frNN(cbind(bp,rep(0,length(bp))), eps = checkOvl.bp))
    y$bp.group = dbscan::dbscan(nn, minPts = 2)$cluster
    out = y[y$bp.group == 0 & y$rank.group ==0,]
    y = y[y$bp.group != 0 | y$rank.group != 0,]
    y$rank.group[y$rank.group==0]<-NA
    y$bp.group[y$bp.group==0]<-NA
    if(nrow(out)>0){
      out$final.group = 1:nrow(out)
      i=max(out$final.group)
    }else{
      i=0
    }
    y$final.group = NA

    while(nrow(y)>0){
      i = i+1
      wh.bp = which(y$bp.group == y$bp.group[1])
      wh.rank = which(y$rank.group == y$rank.group[1])
      n.bp = length(wh.bp)
      n.rank = length(wh.rank)

      if(n.bp > n.rank){
        wh.choose = wh.bp
      }else{
        wh.choose = wh.rank
      }
      y$final.group[wh.choose]<-i
      out<-rbind(out, y[wh.choose,])
      y = y[-wh.choose,]
    }
    out$bp.low = ave(out$bp, out$final.group, FUN = min)
    out$bp.high = ave(out$bp, out$final.group, FUN = max)
    out$bp.med = round(ave(out$bp, out$final.group, FUN = median),0)
    out = do.call(rbind, lapply(split(out, out$final.group), function(z){
      z$new.bp = ifelse(any(z$chr.ends == "start"),z$bp.low,
                        ifelse(any(z$chr.ends == "end"), z$bp.high,z$bp.med))
      return(z)
    }))
    out = out[,c("genome","alt.genome","chr","block.id","type","bp","rank","final.group","new.bp")]

    return(out)
  }))

  nblk = data.frame(new.coords, stringsAsFactors = F)
  blk = make_blocks(map)$block
  blko = blk
  for(i in unique(blk$block.id)){
    y = nblk[nblk$block.id == i,]
    tmp = blk[blk$block.id ==i,]
    se = y$new.bp[with(y, c(which(genome == tmp$genome1 & type == "start"),
                            which(genome == tmp$genome2 & type == "start"),
                            which(genome == tmp$genome1 & type == "end"),
                            which(genome == tmp$genome2 & type == "end")))]
    blk[blk$block.id ==i,c("start1","start2","end1","end2")]<-se
  }
  return(blk)
}
