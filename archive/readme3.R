parse_blast2blocks = function(blast.results,
                              abbrevs,
                              mcscanx.input.dir,
                              MCScanX.path,
                              MCScanX.params = "-a -s 5 -m 25",
                              merge.buffer = 1.5,
                              plotDiagnostics = T,
                              plotPairwise = T,
                              verbose = T){

  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######
  plot_blocksAndMapping = function(map,
                                   blk,
                                   ref.id,
                                   altGenome2plot,
                                   chr1toplot,
                                   chr2toplot,
                                   main = NULL){
    tpb = blk[blk$genome1 == ref.id & blk$genome2 == altGenome2plot,]
    tp = map[map$genome1 == ref.id & map$genome2 == altGenome2plot,]
    tpb$s1 = with(tpb, ifelse(orient == "+", rankstart1, rankend1))
    tpb$e1 = with(tpb, ifelse(orient == "+", rankend1, rankstart1))
    tpb$rankstart1 = tpb$s1
    tpb$rankend1 = tpb$e1
    tpb$s1 = NULL
    tpb$e1 = NULL
    if(!is.null(chr1toplot)){
      tpb = tpb[tpb$chr1 %in% chr1toplot,]
      tp = tp[tp$chr1 %in% chr1toplot,]
    }
    if(!is.null(chr2toplot)){
      tpb = tpb[tpb$chr2 %in% chr2toplot,]
      tp = tp[tp$chr2 %in% chr2toplot,]
    }
    sb1 = split(tpb, tpb$chr1)
    st1 = split(tp, tp$chr1)
    for(i in names(sb1)){
      sb2 = split(sb1[[i]], sb1[[i]]$chr2)
      st2 = split(st1[[i]], st1[[i]]$chr2)
      for(j in names(sb2)){
        t2 = st2[[j]]
        b2 = sb2[[j]]
        cols = rep_len(c("red3","salmon","darkorange3","gold",
                         "grey50","lightgreen","forestgreen","darkgreen",
                         "cyan","dodgerblue3","violet","purple"), length.out = nrow(b2))
        with(t2, plot(rank1, rank2,
                      col = cols[as.numeric(as.factor(block.id))],
                      pch = 16, cex = .5,
                      ylab = paste(altGenome2plot,j,"gene order"),
                      xlab = paste(ref.id,i,"gene order"),
                      main = main))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 =rankend2,
                          col = "black", lwd = 1.5))
        with(b2, text(x = rowMeans(b2[,c("rankstart1","rankend1")]),
                      y = rowMeans(b2[,c("rankstart2","rankend2")]),
                      labels = block.id,
                      col = "black", cex = .5, adj = c(1,-1)))
      }
    }
  }

  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######
  merge_adjacentBlocks = function(map, blk, buffer = 1.5,
                                    verbose = T){
    if(verbose)
      cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
    if(verbose)
      cat("Looking for overlapping blocks ...\n")

    whichSegmentsIntersect <- function(segment1, segment2) {
      getBoundingBox <- function(P0, P1) {
        llx <- apply(cbind(P0[,1], P1[,1]),1, min)
        lly <- apply(cbind(P0[,2], P1[,2]),1, min)
        urx <- apply(cbind(P0[,1], P1[,1]),1, max)
        ury <- apply(cbind(P0[,2], P1[,2]),1, max)
        return(cbind(llx, lly, urx, ury))
      }
      box1 <- getBoundingBox(segment1[,1:2], segment1[,3:4])
      box2 <- getBoundingBox(segment2[,1:2], segment2[,3:4])
      out = apply(box2, 1, function(k){
        chk1 <- box1[1] <= k[3]
        chk2 <- box1[3] >= k[1]
        chk3 <- box1[2] <= k[4]
        chk4 <- box1[4] >= k[2]
        return(chk1 && chk2 && chk3 && chk4)
      })
      return(out)
    }

    blk$block.id = as.character(blk$block.id)
    blk$uniq = paste0(blk$genome1,"_",blk$genome2,"_",blk$chr1, "_", blk$chr2)
    spl = split(blk, blk$uniq)
    nr = sapply(spl, nrow)
    spl = spl[nr>1]
    nm = 1
    for(i in names(spl)){
      x = spl[[i]]
      x = x[order(x$start1, x$start2),]

      has.ovl = sapply(1:nrow(x), function(y){
        k = x[,c("rankstart1","rankstart2","rankend1","rankend2")]
        z = x[y,c("rankstart1","rankstart2","rankend1","rankend2")]
        z[,1]<-max(z[,1]-buffer,0)
        z[,2]<-max(z[,2]-buffer,0)
        z[,3]<-max(z[,3]+buffer,0)
        z[,4]<-max(z[,4]+buffer,0)
        return(whichSegmentsIntersect(z, k))
      })

      colnames(has.ovl)<-x$block.id
      rownames(has.ovl)<-x$block.id
      al = any(apply(has.ovl,1, all))
      diag(has.ovl)<-FALSE
      an = any(has.ovl)
      if(!al & an){
        has.ovl = has.ovl[rowSums(has.ovl)>0, colSums(has.ovl)>0]
        t1 = apply(has.ovl, 1, function(x) rownames(has.ovl)[which(x)])
        t2 = apply(has.ovl,2, function(x) rownames(has.ovl)[which(x)])
        mlist = list()
        for(j in colnames(has.ovl)){
          tmp = unique(unlist(c(j, t1[sapply(t1, function(x) any(grepl(j,x)))],
                                t2[sapply(t2, function(x) any(grepl(j,x)))],
                                names(t1)[sapply(t1, function(x) any(grepl(j,x)))],
                                names(t2)[sapply(t2, function(x) any(grepl(j,x)))])))
          mlist[[j]] = tmp[order(tmp)]
        }
        mlist = mlist[!duplicated(mlist)]

        for(j in 1:length(mlist)){
          map$block.id[map$block.id %in% mlist[[j]]]<-mlist[[j]][1]
        }
      }else{
        if(al){
          map$block.id[map$block.id %in% x$block.id]<-x$block.id[1]
        }
      }
    }
    out = make_blocks(map)
    map = data.frame(out[["map"]], stringsAsFactors = F)
    blk = data.frame(out[["block"]], stringsAsFactors = F)
    if(verbose)
      cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
    return(list(block = blk, map = map))
  }

  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######
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

  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######


  run_MCScanX = function(blast.results,
                         abbrevs,
                         mcscanx.input.dir,
                         MCScanX.params = "-a -s 5 -m 25"){

    if(file.exists(mcscanx.input.dir)){
      system(paste("rm -r",mcscanx.input.dir))
    }
    system(paste("mkdir",mcscanx.input.dir))

    br = data.frame(blast.results, stringsAsFactors = F)
    bs = br[!duplicated(br[,c("genome1","genome2")]),c("genome1","genome2")]
    g1 = names(table(bs$genome1)[order(table(bs$genome1), decreasing = T)])
    glast = unique(bs$genome2[!bs$genome2 %in% g1])
    g = 1:length(c(g1,glast))
    names(g) = c(g1,glast)

    bs$index1 = NA
    for(i in names(g)) bs$index1[bs$genome1 == i]<-g[i]
    bs$index2 = NA
    for(i in names(g)) bs$index2[bs$genome2 == i]<-g[i]
    bs$abbrev1 = NA
    for(i in names(g)) bs$abbrev1[bs$genome1 == i]<-abbrevs[i]
    bs$abbrev2 = NA
    for(i in names(g)) bs$abbrev2[bs$genome2 == i]<-abbrevs[i]

    bs = bs[order(bs$index1, bs$index2),]

    spl = split(blast.results, paste(blast.results$genome1, blast.results$genome2))

    out = lapply(1:nrow(bs), function(i){
      id = paste(bs[i,1:2], collapse = " ")
      x = spl[[id]]
      gff1 = x[,c("chr1","id1","start1","end1")]
      chr1 =  as.numeric(as.factor(gff1$chr1))
      id1 = x$genome1[1]
      abbrev1 = abbrevs[id1]
      gff1$chr1 = paste0(abbrev1,chr1)

      gff2 = x[,c("chr2","id2","start2","end2")]
      chr2 =  as.numeric(as.factor(gff2$chr2))
      id2 = x$genome2[1]
      abbrev2 = abbrevs[id2]
      gff2$chr2 = paste0(abbrev2,chr2)
      names(gff2) = names(gff1)
      gff = data.table(rbind(gff1, gff2))
      gff = gff[order(gff$chr1, gff$start1),]

      blast = x[,c("id1","id2","perc.iden","align.length",
                   "n.mismatch", "n.gapOpen", "q.start", "q.end",
                   "s.start", "s.end", "eval", "score")]

      return(list(gff = gff, blast = blast))
    })

    gff.o = rbindlist(lapply(out, function(x) x$gff))
    blast.o = rbindlist(lapply(out, function(x) x$blast))

    write.table(gff.o, file = file.path(mcscanx.input.dir,"all.gff"),
                quote = F, sep = "\t", row.names = F, col.names = F)
    write.table(blast.o, file = file.path(mcscanx.input.dir, "all.blast"),
                quote = F, sep = "\t", row.names = F, col.names = F)

    mcscan.input = file.path(mcscanx.input.dir,"all")

    if(is.null(MCScanX.params)){
      com = paste("MCScanX",mcscan.input)
    }else{
      com = paste("MCScanX",MCScanX.params,mcscan.input)
    }
    system(com)
    mcscan.raw = read.delim(paste0(mcscan.input,".collinearity"),
                     sep = "\t", header  =F,
                     comment.char = "#", strip.white = T,
                     stringsAsFactors = F)

    fac = sapply(as.character(mcscan.raw$V1),function(x) strsplit(x,"-")[[1]][1])
    m = data.table(mcscan.raw[,2:3])
    setnames(m, c("id1","id2"))
    m2 = data.table(m)
    setnames(m, c("id2","id1"))
    m$block.id = as.numeric(as.factor(fac))
    m2$block.id = as.numeric(as.factor(fac))
    out1 = merge(m, blast.results, by = c("id1","id2"))
    out2 = merge(m2, blast.results, by = c("id1","id2"))
    out = rbind(out1, out2)
    out$rank1 = frank(out[,c("chr1","start1")], ties.method = "dense")
    out$rank2 = frank(out[,c("chr2","start2")], ties.method = "dense")
    return(out)
  }

  merge_overlappingBlocks = function(map, blk, buffer = 1.5, verbose = T){
    if(verbose)
      cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
    if(verbose)
      cat("Looking for overlapping blocks ...\n")
    map$block.id = paste(map$block.id,map$genome1, map$genome2, map$chr1, map$chr2)
    o = make_blocks(map)
    map = o$map
    blk = o$block
    spl = split(blk, paste(blk$genome1, blk$genome2, blk$chr1, blk$chr2))
    spl2 = split(map, paste(map$genome1, map$genome2, map$chr1, map$chr2))
    test = lapply(names(spl), function(i){
      x = spl[[i]]
      x = x[order(-x$n.mapping),]

      y = spl2[[i]]

      for(j in x$block.id){
        if(j %in% y$block.id){
          wh = with(x, which(
            ((rankstart1 + buffer) >= rankstart1[block.id == j]  &
               (rankstart2 + buffer) >= rankstart2[block.id == j]  &
               (rankstart1 - buffer) <= rankend1[block.id == j]  &
               (rankstart2 - buffer)  <= rankend2[block.id == j]) |
              ((rankend1 + buffer) >= rankstart1[block.id == j] &
                 (rankend2 + buffer) >= rankstart2[block.id == j]  &
                 (rankend1 - buffer)  <= rankend1[block.id == j] &
                 (rankend2 - buffer)  <= rankend2[block.id == j])))
          if(length(wh)>0){
            tomerge = x$block.id[wh]
            y$block.id[y$block.id %in% tomerge]<-j
            o = make_blocks(y)
            y = o$map
            x = o$block
          }
        }
      }
      return(y)
    })
    tmp = rbindlist(test)
    tmp$block.id = as.numeric(as.factor(tmp$block.id))
    out = make_blocks(tmp)

    map = data.frame(out[["map"]], stringsAsFactors = F)
    blk = data.frame(out[["block"]], stringsAsFactors = F)
    if(verbose)
      cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
    return(list(block = blk, map = map))
  }


  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######


  split_byDensity<-function(map, max.dist = 5,
                            quantile = 0.999, verbose = T){

    split_it = function(map, radius,quantile){
      spl = split(map, paste(map$genome1, map$genome2))
      map = rbindlist(lapply(spl, function(x){
        x$x_a = frank(x$start1,  ties.method = "dense")
        x$x_b = frank(x$start2, ties.method = "dense")
        return(x)
      }))

      spl = split(map, map$block.id)
      out = do.call(rbind, lapply(names(spl), function(i){
        x = spl[[i]]

        x$was_split = F
        if(nrow(x)>radius*2){
          nn = dbscan::frNN(x[,c("x_a","x_b")], eps = radius)

          xbar = mean(unlist(nn$dist))
          xsd = sd(unlist(nn$dist))
          qthr = qnorm(quantile, mean = xbar, sd = xsd)

          dbs = dbscan::dbscan(nn, minPts = qthr)
          x$newclust = dbs$cluster
          if(all(x$newclust==0)){
            x$newclust<-1
          }
          x = x[x$newclust!=0,]
          x$block.id = paste0(x$block.id, "_",x$newclust)
          x$was_split = length(unique(dbs$cluster))>1
          x$newclust = NULL
        }
        return(x)
      }))
    }

    ntospl = 1
    map$was_split = TRUE
    while(any(map$was_split) & ntospl <= 3){
      ntospl = ntospl+1
      if(verbose)
        cat("checking", sum(map$was_split[!duplicated(map$block.id)]),"blocks\n")
      map = split_it(map,
                     radius = sqrt((max.dist^2)*2),
                     quantile = quantile)
    }
    map$block.id = as.character(as.numeric(as.factor(map$block.id)))
    out = make_blocks(map)
    if(verbose)
      cat("split into",nrow(out$block),"... Done")
    return(list(map = out$map, block = out$block))
  }


  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######

  mcscan.raw = run_MCScanX(blast.results = blast.results,
                           abbrevs = abbrevs,
                           mcscanx.input.dir = mcscanx.input.dir,
                           MCScanX.path = MCScanX.path,
                           MCScanX.params = MCScanX.params)
  mcscan.blocks = make_blocks(mcscan.raw)
  old = mcscan.blocks
  new = old
  ndiff = 1000
  if(verbose)
    cat("Merging", nrow(old$block),"Blocks, if they overlap ...\n")
  while(ndiff>2){
    old = new
    new = merge_overlappingBlocks(map = old$map,
                                     blk = old$block,
                                     buffer = merge.buffer,
                                     verbose = F)
    ndiff = nrow(old$block) - nrow(new$block)
    if(verbose)
      cat("\tTo", nrow(new$block),"\n")
  }
  merged = new

  blk = merged$block
  if(plotDiagnostics){
    ui = unique(with(blk[blk$n.mapping > 20,], paste(genome1, genome2, chr1, chr2)))
    pdf("~/Downloads/test.pdf", height = 5, width = 5)
    for(i in ui){
      plot_blocksAndMapping(map = merged$map,
                            blk= merged$block,
                            ref.id = strsplit(i," ")[[1]][1],
                            altGenome2plot = strsplit(i," ")[[1]][2],
                            chr1toplot = strsplit(i," ")[[1]][3],
                            chr2toplot =  strsplit(i," ")[[1]][4],
                            main = "Raw MCScanX blocks")
    }
    dev.off()


  ###### ###### ###### ###### ###### ######
  ###### ###### ###### ###### ###### ######


}



  ################
  ################




  ###
  # 1. Prep the blast scheme from the results


  prep_MCScanX



}
