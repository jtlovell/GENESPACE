check_blkOverlap = function(blk, buffer = 1e5, verbose = T, n.cores = 1){
  if(verbose)
    cat("Building block database ... ")
  blk.long = with(blk, rbind(
    data.table(block.id = block.id,
               genome = genome1,
               chr = chr1,
               start = start1,
               end = end1),
    data.table(block.id = block.id,
               genome = genome2,
               chr = chr2,
               start = start2,
               end = end2)))
  blk.long$round.start = round(blk.long$start,-log10(buffer))
  blk.long$round.end = round(blk.long$end,-log10(buffer))
  blk.str = rbindlist(apply(blk.long,1,function(x)
    data.table(block.id = x[1], genome = x[2], chr = x[3],
               seqs = seq(from = x[6], to = x[7], by = buffer))))
  if(verbose)
    cat("Done!\nBuilding database of overlapping blocks:\n")
  spl.genome = split(blk.str,"genome")
  out.ovl <- rbindlist(lapply(names(spl.genome), function(i){
    if(verbose)
      cat("\tRunning",i,"... ")
    x = spl.genome[[i]]
    spl.chr = split(x, "chr")
    out.chr <- rbindlist(lapply(names(spl.chr), function(j){
      y = spl.chr[[j]]
      spl.blk = split(y, "block.id")
      out.blk = rbindlist(mclapply(names(spl.blk), mc.cores = 6, function(k){
        z = spl.blk[[1]]
        u = unique(y$block.id[y$seqs %in% z$seqs])
        return(data.table(genome = i, chr = j, block = k, ovl.blks =  u))
      }))
      return(out.blk)
    }))
    if(verbose)
      cat("Done!\n")
    return(out.chr)
  }))
  return(out.ovl)
}
