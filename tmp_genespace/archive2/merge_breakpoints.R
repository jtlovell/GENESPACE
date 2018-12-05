#' @title Merge block breakpoints
#'
#' @description
#' \code{merge_breakpoints} Reduces the total number of block breakpoints, so that, when
#' blocks are concatenated, there are fewer small orphan blocks.
#'
#' @param blk The block object (data.frame/data.table)
#' @param map The map object (data.frame/data.table)
#' @param genomeIDs Character vector indicating the genome IDs to consider
#' @param max.bp Numeric, the physical distance (bp) that a breakpoint is permitted to move
#' @param checkOvl.rank Number, the radius (gene order rank) of 2D scan to look for nearby
#' block breakpoints
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details ...
#' @return A new block dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
merge_breakpoints = function(blk,
                             map,
                             genomeIDs,
                             checkOvl.rank = 4,
                             verbose = T,
                             max.bp = 1e5,
                             ...){
  blk = data.table(blk)
  map = data.table(map)
  if(verbose)
    cat("Starting with", length(unique(unlist(
      blk[,c("start1","start2","end1","end2"),with = F]))),
      "unique block breakpoints\n")
  blk_tomerge = rbindlist(lapply(genomeIDs,  function(x){
    y = map[map$genome1 == x | map$genome2 == x,]
    y$wh = ifelse(y$genome1 == x ,1,2)
    if(all(y$wh == 2)){
      tmp = y[,c(2,1,3,9:13,4:8,14:ncol(y)), with = F]
      setnames(tmp,colnames(y))
    }else{
      if(any(y$wh == 2) & any(y$wh == 1)){
        y1 = y[y$wh == 1,]
        y2 = y[y$wh == 2,]
        tmp = y2[,c(2,1,3,9:13,4:8,14:ncol(y2)), with = F]
        setnames(tmp,colnames(y))
        tmp = rbind(tmp,y1)
      }else{
        tmp = y
      }
    }
    z = make_blocks(tmp)$block
    zspl = split(z, z$chr1)
    merge_db = rbindlist(lapply(zspl, function(y){
      long = with(y,data.frame(genome = y$genome1[1],
                               block.id = c(block.id, block.id),
                               initial.chr = c(chr1, chr1),
                               initial.coord = c(start1, end1),
                               type = rep(c("start","end"), each = nrow(y)),
                               coord = c(start1, end1),
                               rank = c(rankstart1, rankend1),
                               y = 1,
                               stringsAsFactors = F))
      nn = frNN(long[,c("rank","y")], eps = checkOvl.rank)
      long$clus = dbscan(nn, minPts = 2)$cluster
      long = long[long$clus!=0,]
      out = rbindlist(lapply(split(long, long$clus), function(k){
        if(all(k$type == "start")){
          k$coord = min(k$coord)
        }else{
          if(all(k$type == "end")){
            k$coord = max(k$coord)
          }else{
            k$coord = mean(k$coord)
          }
        }
        return(k)
      }))
    }))
    return(merge_db)
  }))

  blk_tomerge$diff = with(blk_tomerge,initial.coord - coord)
  blk_tomerge = blk_tomerge[abs(blk_tomerge$diff)<=max.bp,]
  blk.out = data.frame(blk)
  for(i in 1:nrow(blk_tomerge)){
    x = blk_tomerge[i,]
    wh1 = with(blk.out,which(genome1 == x$genome &
                               chr1 == x$initial.chr &
                               start1 == x$initial.coord &
                               x$type == "start"))
    wh2 = with(blk.out,which(genome1 == x$genome &
                               chr1 == x$initial.chr &
                               end1 == x$initial.coord &
                               x$type == "end"))
    wh3 = with(blk.out,which(genome2 == x$genome &
                               chr2 == x$initial.chr &
                               start2 == x$initial.coord &
                               x$type == "start"))
    wh4 = with(blk.out,which(genome2 == x$genome &
                               chr2 == x$initial.chr &
                               end2 == x$initial.coord &
                               x$type == "end"))

    if(length(wh1)>0){
      blk.out$start1[wh1]<-x$coord
    }
    if(length(wh2)>0){
      blk.out$end1[wh2]<-x$coord
    }
    if(length(wh3)>0){
      blk.out$start2[wh3]<-x$coord
    }
    if(length(wh4)>0){
      blk.out$end2[wh4]<-x$coord
    }
  }
  if(verbose)
    cat("Returning with", length(unique(unlist(
      blk.out[,c("start1","start2","end1","end2")]))),
      "unique block breakpoints\n")
  return(blk.out)
}

