complete_ogNet = function(orthogroup.list,
                          gene.track.net,
                          verbose = T){
  og.list = orthogroup.list
  out.net = gene.track.net
  if(verbose)
    cat("Merging orthogroup list with gene-tracking network ... ")
  setkey(out.net, id)
  og.dt = data.table(og.id = rep(names(og.list),sapply(og.list,length)),
                     id = unlist(og.list))
  setkey(og.dt, id)
  og.blk = merge(og.dt, out.net)

  if(verbose)
    cat("Done!\nCulling database to only multi-genome orthogroups ... ")
  og.blk[,unique.chr := paste0(genome, chr, og.id)]
  og.blk[,n.genomes := length(unique(genome)),
         by = og.id]
  og.blk = og.blk[og.blk$n.genomes > 1,]
  if(verbose)
    cat("Done!\nMaking granges object from data.table ... ")
  gr <- makeGRangesFromDataFrame(
    og.blk,
    seqnames.field = "unique.chr",
    ignore.strand = T,
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE)

  if(verbose)
    cat("Done!\nFinding tracking overlaps ... ")
  og.blk$cluster <- findOverlaps(
    gr,
    ignore.strand = TRUE,
    maxgap = 1e7,
    select = "first",
    type = "any",
    drop.self = F,
    drop.redundant = F)

  if(verbose)
    cat("Done!\nSummarizing by overlap ... ")
  out <- og.blk[,list(start = median(start),
                      end = median(end)),
                by = list(og.id, genome, chr, cluster)]
  out$start <- round(out$start)
  out$end <- round(out$end)
  out$length = with(out, end - start)
  setkey(out, genome, chr)
  if(verbose)
    cat("Done!\n")
  return(out[,c("og.id","genome","chr","start","end")])
}

