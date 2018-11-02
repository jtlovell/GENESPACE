#' @title Conduct reciprocal protein alignments between a pair of genomes
#'
#' @description
#' \code{align_peptideByDiamond} Uses the Diamond Blast program to generate
#' databases, then reciprocally align.
#'
#' @param map pairwise mapping file
#' @param tsp.method Method passed to solve.TSP. Should be Concorde, since it
#' is by far the best.
#' @param Concorde.path If tsp.method = Concorde, the path the the program.
#' Should be something like `../concorde/TSP`
#' @param max.jump The maximum size gap between markers in the path
#' for which a block is not split. If a mapping is this far from
#' any other mapping, it is dropped. This is calculated within block.
#' @param plotit Should plots be returned for all multiple-block
#' blocks.
#' @param verbose Logical, should status updates be printed?
#' @param ... Additional arguements passed on to solve.TSP.
#' @details More here soon.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import TSP
#' @import data.table
#' @export
split_blocksByTSP<-function(map,
                            tsp.method = "Concorde",
                            Concorde.path = NULL,
                            max.jump = 5,
                            plotit = T,
                            verbose = T,
                            min.block.size = 5,
                            ...){

  if(tsp.method == "Concorde"){
    if(is.null(Concorde.path))
      stop("If using Concorde TSP, must supply path to program.\n")
    concorde_path(Concorde.path)
  }

  spl = split(map, map$block.id)
  spl = spl[order(-sapply(spl, nrow))]
  if(verbose)
    cat("Read",nrow(map),"pairwise mappings across",length(spl),
        "blocks\nSplitting within blocks via TSP:\n")
  n = ifelse(length(spl)<50,10,
             ifelse(length(spl)<200,40,
                    ifelse(length(spl)<1000,200,
                           ifelse(length(spl)<2000,400,1000))))
  brk = lapply(names(spl), function(i){
    if(which(names(spl) == i) %% n ==0)
      cat("completed",which(names(spl) == i),"/",length(spl),"\n")
    x = spl[[i]]
    r1<-frank(x$rank1, ties.method = "dense")
    r2<-frank(x$rank2, ties.method = "dense")
    d1 = diff(r1)
    d2 = diff(r2)
    if(nrow(x)>((max.jump*2)+1) & max(c(d1,d2)>max.jump)){
      nna = 1
      while(nna>0){
        tsp.out = run_TSP(x = x$rank1, y = x$rank2,
                          tsp.method = tsp.method,
                          unique = paste0("r_",1:nrow(x)),
                          rerank = TRUE,
                          max.jump = max.jump)
        xr= frank(x$rank1, ties.method = "dense")
        yr= frank(x$rank2, ties.method = "dense")
        if(plotit & length(unique(tsp.out))>1){
          plot(xr, yr, type = "n",
               ylab = x$genome2[1],
               xlab = x$genome1[1],
               main = paste(x$chr1[1], x$chr2[1], x$block.id[1]))
          points(xr, yr, col = ifelse(is.na(tsp.out),"black",tsp.out),
                 pch = ifelse(is.na(tsp.out),4,1),
                 cex = .3)
          text(tapply(xr,tsp.out,mean),
               tapply(yr,tsp.out,mean),
               labels = names(tapply(yr,tsp.out,mean)),
               col = as.numeric(names(tapply(yr,tsp.out,mean))),
               cex = .5, adj = c(2,-2))
        }

        x$tmp = x$block.id
        x$block.id= paste0(x$block.id,"_", tsp.out)
        nna = sum(is.na(tsp.out))
        x = x[!is.na(tsp.out),]
        if(nna>0){
          x$block.id= x$tmp
        }
        x$tmp = NULL
      }
    }
    return(x)
  })

  tmp = rbindlist(brk)
  tokeep = table(tmp$block.id)
  tokeep = names(tokeep)[tokeep>=min.block.size]
  tmp = tmp[tmp$block.id %in% tokeep,]

  tmp = make_blocks(tmp)

  if(verbose)
    cat("Done! Split",nrow(tmp$map),"pairwise mappings into",
        nrow(tmp$block),
        "blocks\n")

  return(list(map = tmp$map, block = tmp$block))
}
