finalize_mappings = function(map,
                             blk,
                             gff.dir,
                             genomeIDs,
                             prop.overlap= .25,
                             verbose = T){
  if(verbose)
    cat("Parsing blocks and mappings\n")
  bg = rbindlist(lapply(genomeIDs, function(x){
    wh1 = which(blk$genome1 == x)
    wh2 = which(blk$genome2 == x)
    b = with(blk,
             data.table(genome = x,
                        alt.genome = c(genome2[wh1], genome1[wh2]),
                        chr = c(chr1[wh1], chr2[wh2]),
                        block.id = c(block.id[wh1],
                                     block.id[wh2]),
                        start = c(start1[wh1], start2[wh2]),
                        end = c(end1[wh1], end2[wh2]),
                        stringsAsFactors = F))

    df = with(b, data.table(genome = c(genome, genome),
                            chr = c(chr, chr),
                            bp = c(start, end),
                            stringsAsFactors = F))
    df = df[!duplicated(df),]
    setkey(df, genome, chr, bp)
    return(df)
  }))

  mg = rbindlist(lapply(genomeIDs, function(x){
    wh1 = which(map$genome1 == x)
    wh2 = which(map$genome2 == x)
    b = with(map,
             data.table(genome = x,
                        alt.genome = c(genome2[wh1], genome1[wh2]),
                        id = c(id1[wh1], id2[wh2]),
                        alt.id = c(id2[wh1], id1[wh2]),
                        stringsAsFactors = F))
    setkey(b, alt.genome)
    return(b)
  }))

  if(verbose)
    cat("Parsing gff files\n")

  gff.files = list.files(gff.dir, full.names = T)
  names(gff.files) = gsub(".gff3$","",basename(gff.files))

  parse_gff = function(gff){
    g = suppressWarnings(
      data.table::fread(gff,showProgress = F, verbose = F))
    g = g[g$V3 == "gene",c(9,1,4,5,7)]
    g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
    data.table::setnames(g, c("id","chr","start","end","strand"))
    return(g)
  }

  gff = rbindlist(lapply(names(gff.files), function(i){
    tmp = parse_gff(gff.files[[i]])
    tmp$genome = i
    tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
    return(tmp)
  }))

  setkey(gff,genome,chr, end)

  if(verbose)
    cat("Finding gene overlaps between blocks\n")

  tmp.gff = gff
  spl.gff = list()
  for(i in 1:nrow(br)){
    x = bg[i,]
    wh = tmp.gff$genome == x$genome &
      tmp.gff$chr == x$chr &
      tmp.gff$end <= x$bp
    out = tmp.gff[wh,]
    out$new.block = paste0(x$genome,"_",i)
    spl.gff[[i]]<-out
    tmp.gff<-tmp.gff[!wh,]
  }

  spl.gff = spl.gff[sapply(spl.gff, nrow)>0]
  spl.ids = lapply(spl.gff, function(x) {
    tmp = x[,c(6,8,1,2)]
    tmp$start = min(x$start)
    tmp$end = max(x$end)

    setkey(tmp, genome, id)
    return(tmp)
  })

  if(verbose)
    cat("Merging blocks by percent gene overlap\n")

  index = rbindlist(spl.ids)
  setnames(index,c("alt.genome","alt.block","alt.id","alt.chr","alt.start","alt.end"))
  setkey(mg, genome, id)
  blk.tab = rbindlist(lapply(spl.ids, function(x){
    len = nrow(x)
    me = merge(merge(x, mg), index, by = c("alt.genome","alt.id"))
    if(nrow(me)>1){
      tab = table(me$alt.block)

      out = data.table(ref.block = x$new.block[1],
                       block.id = names(tab),
                       n.hits = as.numeric(tab),
                       n.genes = len)
    }else{
      out = data.table(ref.block = x$new.block[1],
                       block.id = NA,
                       n.hits = 0,
                       n.genes = len)
    }

    return(out)
  }))
  blk.tab$prop = blk.tab$n.hits / blk.tab$n.genes
  blk.table = blk.tab[blk.tab$prop>=prop.overlap,]
  if(verbose)
    cat("Initial n. unique blocks = ",nrow(blk.table),"\n")
  bm = lapply(unique(c(blk.table$ref.block, blk.table$block.id)), function(i){
    unique(unlist(blk.table[blk.table$ref.block == i | blk.table$block.id == i,1:2]))
  })
  bm = bm[!duplicated(bm)]
  for(i in 1:8){
    tmp = lapply(unlist(bm), function(x){
      unique(unlist(bm[sapply(bm,function(i) any(x %in% i))]))
    })
    bm = tmp[!duplicated(tmp)]
  }
  if(verbose)
    cat("Combined into",length(bm),"blocks\n")

  blk.map = bm
  outindex = index
  setnames(outindex, gsub("alt.","", colnames(outindex)))

  names(blk.map)<-paste0("block_",1:length(blk.map))
  blk.genes = sapply(names(blk.map), USE.NAMES = T, simplify = F, function(x){
    outindex$id[outindex$block %in% blk.map[[x]]]
  })

  if(verbose)
    cat("Checking for breaks within blocks\n")

  blk.gff = sapply(names(blk.genes), USE.NAMES = T, simplify = F, function(x){
    gff[gff$id %in% blk.genes[[x]]]
  })
  blk.out = rbindlist(lapply(names(blk.gff), function(x){
    tmp = blk.gff[[x]]
    tmp$block.id = x
    u = split(tmp, paste(tmp$chr,tmp$genome, sep = "_"))
    out = rbindlist(lapply(u, function(y){
      nn = dbscan::frNN(cbind(y$order,0), eps = 3)
      dbs = dbscan::dbscan(nn, minPts = 0)
      y$cluster = dbs$cluster
      y$sub.block = paste0(y$block.id, "_",dbs$cluster)
      return(y)
    }))
    return(out)
  }))

  bedlist = split(blk.out, blk.out$block.id)
  bedout = rbindlist(lapply(bedlist, function(x){
    tmp = split(x, x$sub.block)
    out = rbindlist(lapply(tmp, function(y){
      b = y[,list(start=min(start),
                  end=max(end),
                  ngenes = length(start)),
            by=list(genome, chr, block.id, sub.block)]
    }))
    return(out)
  }))

  if(verbose)
    cat("Done\n")

  genelist = sapply(names(blk.gff), USE.NAMES = T, simplify = F, function(x){
    blk.gff[[x]]$id
  })
  return(list(genelist = genelist, bed = bedout))
}
