make_geneBlockMapping <- function(blk, gff, verbose = T){
  if(verbose)
    cat("Preparing block and gff databases ... ")
  gf = gff
  bl = with(blk, rbind(data.table(block.id = block.id,
                                  genome1 = genome1,
                                  genome2 = genome2,
                                  chr1 = chr1,
                                  chr2 = chr2,
                                  start1 = start1,
                                  start2 = start2,
                                  end1 = end1,
                                  end2 = end2),
                       data.table(block.id = block.id,
                                  genome1 = genome2,
                                  genome2 = genome1,
                                  chr1 = chr2,
                                  chr2 = chr1,
                                  start1 = start2,
                                  start2 = start1,
                                  end1 = end2,
                                  end2 = end1)))
  bl$unique.chr = paste(bl$genome1, bl$chr1)
  spl = split(bl, "unique.chr")
  gf$unique.chr = paste(gf$genome, gf$chr)
  gf = gf[gf$unique.chr %in% bl$unique.chr,]
  spl.gf = split(gf, "unique.chr")
  genomes = sapply(names(spl.gf), function(x) strsplit(x, " ")[[1]][1])
  if(verbose)
    cat("Done!\n")
  if(verbose)
    cat("Pulling all blocks for each gene by genome:\n")
  out.comb <- list()
  for(i in unique(genomes)){
    if(verbose)
      cat("\t",i,"... ")
    brks = spl.gf[genomes == i]
    out.all <- unlist(lapply(names(brks), function(i){
      x.bl = data.table(spl[[i]])
      x.gf = data.table(spl.gf[[i]])
      out = apply(x.gf[,3:4],1,function(x)
        x.bl$block.id[x.bl$start1 <= x[1] & x.bl$end1 >= x[2]])
      names(out)<-x.gf$id
      return(out)
    }),recursive = F)
    out.comb[[i]]<-out.all
    if(verbose)
      cat("Done!\n")
  }
  if(verbose)
    cat("Combining blocks into a data.table:\n")
  gb.dt = rbindlist(lapply(names(out.comb), function(i){
    if(verbose)
      cat(i,"... ")
    x = out.comb[[i]]
    x = x[sapply(x, length)>0]
    lengths = sapply(x, length)
    ns = names(x)
    out <- data.table(genome = i, id = rep(ns, lengths), block.id = unlist(x))
    if(verbose)
      cat("Done!\n\t")
    return(out)
  }))
  return(gb.dt)
}
