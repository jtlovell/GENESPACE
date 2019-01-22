make_completeBlocks = function(og.ids = NULL,
                               orthogroup.list,
                               blk,
                               gene.blks,
                               chunk.size = 1000,
                               n.cores = 6,
                               verbose = T){

  ########################################################
  ########################################################
  get_completeBlk = function(og.genes,
                             gene.blks,
                             gb){
    og.block.ids = gene.blks[gene.blks$id %in% og.genes,]
    og.blk = data.frame(gb[gb$block.id %in% og.block.ids$block.id,])
    if(nrow(og.blk) == 0){
      return(data.table(genome = NA, chr = NA, start = NA, end = NA))
    }else{
      og.blk$unique.chr = paste0(og.blk$genome, og.blk$chr)
      og.blk.spl = split(og.blk, og.blk$unique.chr)
      og.intersect = rbindlist(lapply(og.blk.spl, function(x){
        if(nrow(x) == 1){
          x$cluster = 1
        }else{
          gr <- makeGRangesFromDataFrame(x,
                                         seqnames.field = "unique.chr",
                                         keep.extra.columns = TRUE)
          x$cluster <- findOverlaps(gr,
                                    ignore.strand=TRUE,
                                    maxgap = 1e5,
                                    select = "first",
                                    type = "any",
                                    drop.self=F,
                                    drop.redundant=F)
        }
        return(x)
      }))
      og.blk <- og.intersect[,list(start = median(start),
                                   end = median(end)),
                             by = list(genome, chr,cluster)]
      return(og.blk[,c("genome","chr","start","end")])
    }
  }
  ########################################################
  ########################################################
  # --  Make a long-formatted block file
  gb = with(blk, rbind(
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
  setkey(gb, genome, block.id)

  # -- load the orthogroup data
  if(is.null(og.ids)){
    og.ids <- names(orthogroup.list)
  }
  og.list = orthogroup.list[og.ids]

  # -- split into batches
  if(length(og.ids)>chunk.size){
    og.list.split = split(og.list, ceiling(1:length(og.list)/chunk.size))
  }else{
    og.list.split = list(og.list)
  }

  if(verbose)
    cat("Total orthogroups =",length(og.ids),"\n\t")
  out.all<-lapply(1:length(og.list.split), function(i){
    if(verbose)
      cat("Running:", 1+((i-1)*chunk.size),"\n\t")
    x = og.list.split[[i]]
    out = mclapply(names(x), mc.cores = n.cores, function(j)
      get_completeBlk(og.genes = x[[j]],
                      gene.blks = gene.blks,
                      gb = gb))
    names(out)<-names(x)
    return(out)
  })
  if(verbose)
    cat("Done!\n")
  return(unlist(out.all, recursive = F))
}

