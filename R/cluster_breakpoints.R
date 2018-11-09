#' @title Find breakpoint clusters
#'
#' @description
#' \code{cluster_breakpoints} Merge cluster breakpoints based on positions.
#'
#' @param genomeIDs genome identifiers. The first one will be used as the reference
#' coordinate system
#' @param blk The block object
#' @param checkOvl.rank How many positions (gene rank) should be allowed for the
#' breakpoints to move to be merged
#' @param verbose should updates be printed?
#' @param ... Not currently in use
#' @details Nothing yet
#' @return nothing
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @importFrom dbscan frNN dbscan
#' @import data.table
#' @export
cluster_breakpoints = function(blk,
                               map,
                               genomeIDs,
                               gff.dir,
                               verbose = T,
                               min.genesInBlock = 10,
                               min.blockBp = 5e4,
                               checkOvl.rank = 0, ...){

  if(verbose)
    cat("Using", genomeIDs[1],"as the reference coordinate system\n\t")

  if(verbose)
    cat("Parsing gff files\n")
  gff.files <- list.files(gff.dir,
                          full.names = T)
  names(gff.files) <- gsub(".gff3$", "",
                           basename(gff.files))

  parse_gff <- function(gff){
    g <- suppressWarnings(
      data.table::fread(gff,
                        showProgress = F,
                        verbose = F))
    g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
    g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                          strsplit(x, ";")[[1]][2]))
    data.table::setnames(g, c("id", "chr", "start", "end", "strand"))
    return(g)
  }

  gff <- rbindlist(lapply(names(gff.files), function(i){
    tmp <- parse_gff(gff.files[[i]])
    tmp$genome <- i
    tmp$order <- frank(tmp[,c("chr", "start")],
                       ties.method = "dense")
    return(tmp)
  }))
  setkey(gff, genome, id)

  cbp = function(b, genomeIDs, verbose, checkOvl.rank){
    r = with(b, cbind(c(rankstart1,rankend1),
                      c(rankstart1,rankend1)))
    nn = frNN(r, eps = checkOvl.rank)
    dbs = dbscan(nn, minPts = 2)$cluster
    b$start.clus = dbs[1:nrow(b)]
    b$end.clus = dbs[(nrow(b)+1):length(dbs)]
    merged<-b

    for(i in unique(dbs)){
      if(i!=0){
        wh.x = b$start.clus == i | b$end.clus == i
        x = b[wh.x,]
        if(all(x$start.clus == i)){
          merged$start1[wh.x] <- min(x$start1)
        }else{
          if(all(x$end.clus == i)){
            merged$end1[wh.x] <- max(x$end1)
          }else{
            wh.start = which(b$start.clus == i)
            wh.end = which(b$end.clus == i)
            whx.start = which(x$start.clus == i)
            whx.end = which(x$end.clus == i)
            mean.pos = mean(c(x$start1[whx.start], x$end1[whx.end]))
            merged$start1[wh.start]<-mean.pos
            merged$end1[wh.end]<-mean.pos
          }
        }
      }
    }

    b <- merged
    r = with(b, cbind(c(rankstart1,rankend1),
                      c(rankstart1,rankend1)))
    nn = frNN(r, eps = checkOvl.rank)
    dbs = dbscan(nn, minPts = 2)$cluster

    b$start.clus = dbs[1:nrow(b)]
    b$end.clus = dbs[(nrow(b)+1):length(dbs)]
    merged<-b
    for(i in unique(dbs)){
      if(i!=0){
        wh.x = b$start.clus == i | b$end.clus == i
        x = b[wh.x,]
        if(all(x$start.clus == i)){
          merged$start1[wh.x] <- min(x$start1)
        }else{
          if(all(x$end.clus == i)){
            merged$end1[wh.x] <- max(x$end1)
          }else{
            wh.start = which(b$start.clus == i)
            wh.end = which(b$end.clus == i)
            whx.start = which(x$start.clus == i)
            whx.end = which(x$end.clus == i)
            mean.pos = mean(c(x$start1[whx.start], x$end1[whx.end]))
            merged$start1[wh.start]<-mean.pos
            merged$end1[wh.end]<-mean.pos
          }
        }
      }
    }
    return(merged)
  }

  b = blk[blk$genome1 == genomeIDs[1],]
  if(verbose)
    cat("Searching", nrow(b)*2, "blocks ...\n\t")
  spl = split(b, b$chr1)
  out = lapply(spl, function(x){
    o = cbp(b = x, genomeIDs = genomeIDs,verbose = F,checkOvl.rank = checkOvl.rank)
    return(o)
  })

  bo = rbindlist(out)


  rmap = bo[bo$genome1 == genomeIDs[1],]
  rblk = blk[blk$genome1 == genomeIDs[1],]

  brpts = rbind(data.frame(chr = rblk$chr1, pos = rblk$start1,
                           stringsAsFactors = F),
                data.frame(chr = rblk$chr1, pos = rblk$end1,
                           stringsAsFactors = F))

  if(verbose)
    cat("Initial clustering into",nrow(brpts),"breakpoints\n\t")
  brpts = brpts[!duplicated(brpts),]
  brpts = brpts[order(brpts$chr, brpts$pos),]
  brpts = rbindlist(lapply(split(brpts, brpts$chr), function(x){
    x$start = c(x$pos[-nrow(x)],NA)
    x$end = c(x$pos[-1],NA)
    return(x)
  }))
  brpts = data.frame(brpts[complete.cases(brpts),])
  brpts$width = with(brpts, end - start)

  gff.ref = gff[gff$genome == genomeIDs[1],]
  gff.ref.spl = split(gff.ref, gff.ref$chr)
  brpts$n.genes <- sapply(1:nrow(brpts),function(i)
    with(gff.ref.spl[[brpts$chr[i]]],
         sum(chr == brpts$chr[i] & start <= brpts$end[i] & end >= brpts$start[i]))
  )
  wh0 = which(brpts$n.genes == 0)
  brpts$ext.start <- sapply(1:nrow(brpts),function(i){
    if(i %in% wh0){
      brpts$start[i]
    }else{
      with(gff.ref.spl[[brpts$chr[i]]],
           min(start[(chr == brpts$chr[i] & start <= brpts$end[i] & end >= brpts$start[i])]))
    }
  })
  brpts$ext.end <- sapply(1:nrow(brpts),function(i){
    if(i %in% wh0){
      brpts$end[i]
    }else{
      with(gff.ref.spl[[brpts$chr[i]]],
           max(end[(chr == brpts$chr[i] & start <= brpts$end[i] & end >= brpts$start[i])]))
    }
  })
  brpts$brk.id = 1:nrow(brpts)

  wh = brpts$width >=min.blockBp & brpts$n.genes >= min.genesInBlock
  good = brpts[wh, ]
  bad = brpts[!wh,]

  if(verbose)
    cat("Re-assigning",nrow(bad),"blocks that were smaller than",min.blockBp,"bp or", min.genesInBlock,"genes\n\t")

  spl.good = split(good, good$chr)
  fix.bad = rbindlist(lapply(1:nrow(bad), function(i){
    x = bad[i,]
    y = spl.good[[x$chr[1]]]
    close.start = x$start - y$end
    close.end = y$start - x$end
    wh.min = apply(cbind(close.start,close.end),1,function(x) min(abs(x)))
    x$brk.id<-y$brk.id[which.min(wh.min)]
    return(x)
  }))

  r.brpts = rbind(fix.bad, good)

  brpts.o = r.brpts[,list(chr = chr[1],
                          start = min(start),
                          end = max(end),
                          n.genes = sum(n.genes)),
                    by = list(brk.id)]
  brpts.o$width = with(brpts.o, end-start)

  if(verbose)
    cat("Generated a block database with",nrow(brpts.o),"unique breakpoints\n\t")


  rmap = data.table(map[map$genome1 == genomeIDs[1],])
  n.perblock = data.frame(table(block.id= rmap$block.id))
  colnames(n.perblock)[2]<-"n.total"
  rmap.spl = lapply(1:nrow(brpts.o),function(i){
    x = brpts.o[i,]
    y = rmap[with(rmap, chr1 == x$chr & end1 >= x$start & start1 <= x$end),]
    o = y[,list(genome = genome2[1],
                chr = chr2[1],
                start = min(start2),
                end = max(end2)),
          by = list(block.id)]
    o = rbind(data.table(block.id = NA,genome = genomeIDs[1],
                         chr = x$chr,
                         start = x$start, end = x$end),
              o)
    o$breakpoint.id = x$brk.id
    return(o)
  })
  brk.out = rbindlist(rmap.spl)

  if(verbose)
    cat("Done - ",nrow(blk), "blocks have been assigned to",
        length(unique(brk.out$breakpoint.id)),"breakpoints\n")
  if(verbose)
    cat("Building a new map object against the", genomeIDs[1],"reference coordinate system...\n\t",
        "Original map contained", nrow(rmap),"blast hits\n\t")
  test = lapply(rmap.spl, function(x){
    xr = x[x$genome == genomeIDs[1],]
    y = rmap[rmap$chr1 == xr$chr & rmap$end1 >= xr$start & rmap$start1 <= xr$end,]
    yspl = rbindlist(lapply(1:nrow(x), function(i){
      if(x$genome[i] == genomeIDs[1]){
        z = y[y$genome1 == x$genome[i],]
        z = z[z$chr1 == x$chr[i],]
        z = z[z$end1 >= x$start[i] & z$start1 <= x$end[i],]
        z$breakpoint.id <- x$breakpoint.id[1]
      }else{
        z = y[y$genome2 == x$genome[i],]
        z = z[z$chr2 == x$chr[i],]
        z = z[z$end2 >= x$start[i] & z$start2 <= x$end[i],]
        z$breakpoint.id <- x$breakpoint.id[1]
      }
      return(z)
    }))
  })

  if(verbose)
    cat("\tSplitting gff files by breakpoints\n")
  bp = brk.out
  gff.spl = sapply(unique(bp$genome), simplify = F, USE.NAMES = T, function(i){
    gg = gff[gff$genome == i,]
    out = sapply(unique(gg$chr), simplify = F, USE.NAMES = T, function(j){
      gg[gg$chr == j,]
    })
    return(out)
  })

  bp.spl = split(bp, bp$breakpoint.id)

  if(verbose)
    cat("\tExpanding breakpoints by blast hits in mapping\n")

  gff.bp = rbindlist(lapply(names(bp.spl), function(i){
    x = bp.spl[[i]]
    z = rbindlist(lapply(1:nrow(x), function(j){
      y = x[j,]
      go = gff.spl[[y$genome]][[y$chr]]
      go = go[go$end >= y$start & go$start <= y$end,]
      go$block.id = x$block.id[j]
      return(go)
    }))
    mo1 = map[map$id1 %in% z$id,]
    mo2 = map[map$id2 %in% z$id,]

    mo = rbind(mo1, mo2)
    mi = make_blocks(mo)
    mm = mi$map
    mm$block.id = paste0(x$breakpoint.id,"_", mm$block.id)
    mm$breakpoint.id = x$breakpoint.id[1]
    return(mm)
  }))

  blk = make_blocks(gff.bp, rename.blocks = F)
  map  = blk$map
  blk = blk$block
  if(verbose)
    cat("Returning",nrow(map),"mappings within", length(unique(blk$block.id)),"blocks across",length(unique(brk.out$breakpoint.id)),"unique breakpoints\n")
  return(list(map = map, block= blk, breakpoints = brk.out))
}
