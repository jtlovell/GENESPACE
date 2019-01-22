track_genes <- function(gff,
                        orthogroup.list,
                        buffer = 10,
                        merge.hit.dist = 1e7,
                        min.prp.hits = .75,
                        n.cores = 6,
                        verbose = T){
  if(verbose)
    cat("Merging orthogroup networks with gff ... ")
  og.list <- orthogroup.list
  og.lens = sapply(og.list, length)

  og.dt <- data.table(og.id = rep(names(og.list), og.lens),
                      id = unlist(og.list))
  gff.og = merge(gff, og.dt, by = "id")

  if(verbose)
    cat("Done!\nSummarizing gff/orthogroup ... ")
  gog<-gff.og[,list(order = frank(start, ties.method = "random"),
                    start = start,
                    end = end,
                    id = id,
                    og.id = og.id,
                    n.genes = length(start)),
              by = list(genome, chr)]
  gog<-gog[gog$n.genes>5,]
  setkey(gog, genome, chr, order)

  spl = split(gog, "genome")
  if(verbose)
    cat("Done!\nBuilding gene-tracking network for each genome ... \n")
  map.net <- lapply(names(spl),  function(i){
    cat(i,"\n\tBuilding in gene order buffer ...")
    # - pull the focal chromosome
    tmp.chr = spl[[i]]
    # - pull all other chromosomes
    tmp.no = rbindlist(spl[-which(names(spl)==i)])
    # - drop any with the same genome
    tmp.no = tmp.no[tmp.no$genome != tmp.chr$genome[1],]
    # for each gene in chr, find a the nearby ones
    index.list = lapply(1:nrow(tmp.chr),function(i){
      ibuf = (i-buffer):(i+buffer)
      ibuf = ibuf[ibuf>0 & ibuf < nrow(tmp.chr)]
      return(ibuf)
    })

    buf.out = rbindlist(mclapply(index.list,mc.cores = n.cores, function(x)
      tmp.chr[x,]))
    chr.start = rep(tmp.chr$chr, sapply(index.list, length))
    buf.out$focal.id = rep(tmp.chr$id, sapply(index.list, length))
    buf.out <- buf.out[buf.out$chr == chr.start,]

    # buf.out[,min.hits:=length(start)*min.prp.hits,
    #      by = list(focal.id)]
    cat("Done!\n\tSplitting database into ")
    setkey(buf.out, id)
    xlist = split(buf.out, "focal.id")
    xlist.chunk = split(xlist, ceiling(1:length(xlist)/1000))
    if(verbose)
      cat(length(xlist.chunk),"blocks ... ")
    x.gog = rbindlist(lapply(1:length(xlist.chunk), function(i){
      xl = xlist.chunk[[i]]
      all.chunk  = rbindlist(xl)
      og.ids.inchr = gog$og.id[gog$id %in% all.chunk$id]
      all.og.map = gog[gog$og.id %in% og.ids.inchr,]

      x.gog.tmp<-rbindlist(mclapply(names(xl), mc.cores = n.cores, function(j){
        z = xl[[j]]
        og.inbuf = all.og.map$og.id[all.og.map$id %in% z$id]
        map.inbuf = all.og.map[all.og.map$og.id %in% og.inbuf,]
        map.inbuf$focal.id = j
        return(map.inbuf)
      }))
      return(x.gog.tmp)
    }))
    cat("Done!\n\tConverting network to granges ... ")
    x.gog$unique.chr = with(x.gog, paste0(genome, chr))
    x.gog$unique = with(x.gog, paste0(unique.chr, focal.id))
    gr <- makeGRangesFromDataFrame(x.gog,
                                   seqnames.field = "unique",
                                   ignore.strand = T,
                                   start.field = "start",
                                   end.field = "end",
                                   keep.extra.columns = FALSE)
    cat("Done!\n\tFinding overlapping hits ... ")
    x.gog$cluster <- frank(findOverlaps(gr,
                                        ignore.strand=TRUE,
                                        maxgap = merge.hit.dist,
                                        select = "first",
                                        type = "any",
                                        drop.self=F,
                                        drop.redundant=F),
                           ties.method = "dense")

    out = x.gog[order(x.gog$cluster),]
    cat("Done!\n")
    out[,n.hit.by.genome := length(unique(id)),
        by = list(focal.id, genome)]
    out2 <- out[,list(start = min(start),
                     end = max(end),
                     n.hits = length(start)),
               by = list(focal.id, genome, chr, cluster, n.hit.by.genome)]
    out2$length <- with(out2, end - start)
    out2$prop.map = with(out2, n.hits / n.hit.by.genome)
    out2 = out2[out2$prop.map > min.prp.hits,]
    setnames(out2,1,"id")
    return(out2)
  })
  out.net <- rbindlist(map.net)
  return(out.net)
}
