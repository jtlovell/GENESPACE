find_regions2blast <- function(complete.blks.gff,
                               chunk.size = 1000,
                               n.cores = 1,
                               genomeID.rank.multiplier = 2,
                               verbose = T){
  if(verbose)
    cat("Finding best genes to blast for",
        length(unique(complete.blks.gff$og.id)),
        "orthogroups\n\t")

  spl = split(complete.blks.gff, "og.id")
  chunks = split(spl, ceiling(1:length(spl)/chunk.size))

  cl3.out<-lapply(1:length(chunks), function(ch.id){
    if(verbose)
      cat("Running chunk", ch.id, "/", length(chunks),"\n\t")
    chunk = chunks[[ch.id]]
    closest3 <- rbindlist(mclapply(chunk, mc.cores = n.cores, function(x){
      xn = x[is.na(x$id),]
      ng = x[!is.na(x$id),]
      if(nrow(xn) == 0 | nrow(ng) == 0){
        tmp <- x[1,]
        for(i in c(colnames(tmp), "best1","best2","best3")) tmp[,i] <- NA
        return(tmp)
      }else{
        ng$length.rank <- frank(ng, -length, ties.method = "random")

        tmp = data.table(t(sapply(1:length(xn$genome), function(i){
          y = xn$genome[i]
          wh = which(genomeIDs == y)
          dist2genome <-sapply(ng$genome, function(z)
            if(z == y){
              100
            }else{
              abs(which(genomeIDs == z) - wh)
            })
          g.rank = frank(dist2genome, ties.method = "dense")
          l.rank = ng$length.rank
          best3 = which(frank(g.rank*genomeID.rank.multiplier + l.rank, ties.method = "random")<=3)
          if(length(best3)==0){
            best3 <- rep(NA,3)
          }
          if(length(best3)==1){
            best3<-c(best3,NA,NA)
          }
          if(length(best3)==2){
            best3<-c(best3,NA)
          }
          return(ng$id[best3])
        })))
        setnames(tmp,c("best1","best2","best3"))
        xn <- cbind(xn, tmp)
        return(xn)
      }
    }))
    return(closest3)
  })
  return(cl3.out)
}
