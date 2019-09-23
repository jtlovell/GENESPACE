proc_MCScanBlocks <- function(map,
                              genomeIDs,
                              gff,
                              max.blockSize2merge,
                              verbose = T){
  ############################################################
  if(verbose)
    cat("Merging overlapping and duplicated blocks ... \n")
  movlp <- merge_dupOvlBlks(map = map,
                            genomeIDs = genomeIDs,
                            gff = gff)

  ############################################################
  if(verbose)
    cat("Splitting overlapping, non-duplicated blocks ... \n")
  map <- data.table(movlp$map)
  blk <- data.table(movlp$blk)
  ovl.tot <- chk_olv(blk, type = "any")
  mtmp <- subset(map, block.id %in% unique(unlist(ovl.tot[,c("inside","outside")])))
  ovl.tot[,any.dup := sapply(1:nrow(ovl.tot), function(i)
    any_dupInBlk(map = mtmp,
                 block.id1 = inside[i],
                 block.id2 = outside[i]))]
  blk2chk <- subset(ovl.tot, !any.dup)[, c("inside", "outside")]
  blk2chk <- data.table(t(apply(blk2chk, 1, function(x)
    x[order(as.factor(x))])))
  blk2chk <- blk2chk[!duplicated(blk2chk), ]
  ############################################################
  if (verbose)
    cat("\tFound",
        nrow(blk2chk),
        "overlapping blocks that need to be split\n")
  mlist <- apply(blk2chk, 1, function(x)
    subset(mtmp, block.id %in% x))
  mout <- rbindlist(lapply(mlist, spl_ovlGap))
  map <- rbind(mout, map)
  map <- map[!duplicated(map[,c("id1","id2")]),]

  if (verbose)
    cat("\tDone!\nJoining adjacent blocks ... \n\tStarting with",
        length(unique(map$block.id)), "blocks ...\n")
  spl <- split(map, by = c("genome1","genome2","chr1","chr2"))

  comb.out <- rbindlist(lapply(spl, function(x){
    db <- run_dbs(x,
                  eps.radius = max.blockSize2merge,
                  mappings = 1)
    dbl <- tapply(db$block.id, db$cluster, unique)
    if (length(dbl) > 1) {
      for (i in 1:length(dbl))
        x$block.id[x$block.id %in% dbl[[i]]] <- i
    }
    return(x)
  }))
  if (verbose)
    cat("\tFound",
        length(unique(comb.out$block.id)),
        "blocks after joining\n")
  out <- simplify_map(map = comb.out,
                      gff = gff,
                      genomeIDs = genomeIDs)
  if (verbose)
    cat("\tDone!\n")
  return(out)
}
