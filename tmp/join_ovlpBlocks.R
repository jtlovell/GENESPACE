join_ovlpBlks <- function(map, buffer = 50, verbose = T, max.iter = 10){

  index_ovlpBlocks <- function(blk, use.id1 = TRUE, buffer){

    if(use.id1){
      x1 <- blk[,list(start = rankstart1 - buffer,
                      end = rankend1 + buffer),
                by = list(block.id, genome1, genome2, chr1)]
      setnames(x1, "chr1", "chr")
    }else{
      x1 <- blk[,list(start = rankstart2 - buffer,
                      end = rankend2 + buffer),
                by = list(block.id, genome1, genome2, chr2)]
      setnames(x1, "chr2", "chr")
    }
    spl <- split(x1, by = c("genome1","genome2","chr"))
    spl <- spl[sapply(spl, nrow) > 1]
    fo <- rbindlist(lapply(spl, function(x){
      setkey(x, start, end)
      ovl <- foverlaps(x, x, type = "any", which=TRUE)
      ovl[,V1 := x$block.id[xid]]
      ovl[,V2 := x$block.id[yid]]
      return(subset(ovl, xid < yid)[,c("V1","V2")])
    }))
    return(fo)
  }

  index_dupBlocks <- function(map){
    db <- map[,list(id = unique(c(id1, id2))),
              by = list(block.id, genome1, genome2)]
    spl <- split(db, by = c("genome1", "genome2"))
    dbd <- rbindlist(lapply(spl, function(x)
      subset(x, id %in% x$id[duplicated(x$id)])))
    dbdd <- dbd[,list(uniq = paste(unique(block.id), collapse = ",")),
                by = list(id, genome1, genome2)]
    dbdd <- dbdd[!duplicated(dbdd$uniq),]
    dup.comb <- data.table(
      do.call(rbind, lapply(dbdd$uniq, function(x)
        expand.grid(strsplit(x,",")[[1]],strsplit(x,",")[[1]]))))
    dup.comb <- data.table(t(apply(dup.comb, 1, function(x) x[order(x)])))
    dup.comb <- dup.comb[!duplicated(dup.comb),]
    dup.comb <- subset(dup.comb, V1 != V2)
    return(dup.comb)
  }

  blk <- make_blocks(map)$block
  i <- 0
  nb <- nrow(blk)-1
  while(i <= max.iter & nb < nrow(blk)){
    i <- i + 1
    if(verbose){
      cat(paste0("Iter ", i, ": Initial n.blks = "), ifelse(i == 1, nrow(blk), nb),"\n")
      cat("\tChecking for blocks with duplicate hits ... ")
    }
    dup.index <- index_dupBlocks(map)
    if(verbose)
      cat("Done!\n\tChecking for blocks with adjacent coordinates ... ")
    blk <- make_blocks(map)$block
    ovl1 <- index_ovlpBlocks(blk,  buffer = buffer)
    ovl2 <- index_ovlpBlocks(blk, buffer = buffer, use.id1 = F)
    ovl.index <- merge(ovl1, ovl2)

    ovl.nodup <- rbind(data.table(is.dup = T, dup.index),
                       data.table(is.dup = F, ovl.blks))
    ovl.nodup <- subset(ovl.nodup[!duplicated(ovl.nodup[,2:3]), ], !is.dup)[,c(2:3)]
    if(verbose)
      cat("Done!\n\tJoining adjacent blocks ... ")
    ovl.tm <- subset(data.table(block.id = unique(map$block.id)))
    rename.index <- merge(ovl.tm,
                          with(ovl.nodup,
                               data.table(block.id = V2,
                                          new.id = V1)),
                          by = "block.id", all.x = T)
    rename.index[,new.id := ifelse(is.na(new.id), block.id, new.id)]
    out <- merge(map, rename.index, by = "block.id")
    out[,block.id := NULL]
    setnames(out, "new.id", "block.id")
    out <- make_blocks(out)
    nb <- nrow(out$block)
    map <- data.table(out$map)
    if(verbose)
      cat("Done!\n\tReturning merged dataset with", nb,"blocks\n")
  }
  if(verbose)
    cat("Done!\n")
  return(list(ovl.nodup.index = ovl.nodup,
              dup.index = dup.index,
              ovl.index = ovl.index,
              ovl.1.index = ovl1,
              ovl.2.index = ovl2,
              ovl.nodup.index = chk,
              map = out$map, blk = out$block))
}
