#' @title Find breakpoint clusters
#'
#' @description
#' \code{cluster_breakpoints} Merge cluster breakpoints based on positions.
#'
#' @param genomeIDs genome identifiers. The first one will be used as the reference
#' coordinate system
#' @param blk The block object
#' @param checkOvl.rank How many positions (gene rank) should be allowed for the
#' breakpoints to move to be merged
#' @param verbose should updates be printed?
#' @param ... Not currently in use
#' @details Nothing yet
#' @return nothing
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @importFrom dbscan frNN dbscan
#' @import data.table
#' @export
cluster_breakpoints = function(blk,
                               genomeIDs,
                               verbose = T,
                               checkOvl.rank = 5, ...){
  if(verbose)
    cat("Using", genomeIDs[1],"as the reference coordinate system\n\t")
  b = blk[blk$genome1 == genomeIDs[1],]

  r = with(b, cbind(c(rankstart1,rankend1),
                    c(rankstart1,rankend1)))
  nn = frNN(r, eps = checkOvl.rank)
  dbs = dbscan(nn, minPts = 2)$cluster
  if(verbose)
    cat("Clustered",length(unique(c(b$start1, b$end1))),
        "unique breakpoints within",checkOvl.rank,"hits\n\t")
  b$start.clus = dbs[1:nrow(b)]
  b$end.clus = dbs[(nrow(b)+1):length(dbs)]
  merged<-b

  for(i in unique(dbs)){
    if(i!=0){
      wh.x = b$start.clus == i | b$end.clus == i
      x = b[wh.x,]
      if(all(x$start.clus == i)){
        merged$start1[wh.x] <- min(x$start1)
      }else{
        if(all(x$end.clus == i)){
          merged$end1[wh.x] <- max(x$end1)
        }else{
          wh.start = which(b$start.clus == i)
          wh.end = which(b$end.clus == i)
          whx.start = which(x$start.clus == i)
          whx.end = which(x$end.clus == i)
          mean.pos = mean(c(x$start1[whx.start], x$end1[whx.end]))
          merged$start1[wh.start]<-mean.pos
          merged$end1[wh.end]<-mean.pos
        }
      }
    }
  }
  if(verbose)
    cat("Returning",length(unique(c(merged$start1, merged$end1))),
        "merged breakpoints\n")

  b <- merged
  r = with(b, cbind(c(rankstart1,rankend1),
                    c(rankstart1,rankend1)))
  nn = frNN(r, eps = checkOvl.rank)
  dbs = dbscan(nn, minPts = 2)$cluster
  if(verbose)
    cat("Clustered",length(unique(c(b$start1, b$end1))),
        "unique breakpoints within",checkOvl.rank,"hits\n\t")
  b$start.clus = dbs[1:nrow(b)]
  b$end.clus = dbs[(nrow(b)+1):length(dbs)]
  merged<-b
  for(i in unique(dbs)){
    if(i!=0){
      wh.x = b$start.clus == i | b$end.clus == i
      x = b[wh.x,]
      if(all(x$start.clus == i)){
        merged$start1[wh.x] <- min(x$start1)
      }else{
        if(all(x$end.clus == i)){
          merged$end1[wh.x] <- max(x$end1)
        }else{
          wh.start = which(b$start.clus == i)
          wh.end = which(b$end.clus == i)
          whx.start = which(x$start.clus == i)
          whx.end = which(x$end.clus == i)
          mean.pos = mean(c(x$start1[whx.start], x$end1[whx.end]))
          merged$start1[wh.start]<-mean.pos
          merged$end1[wh.end]<-mean.pos
        }
      }
    }
  }
  if(verbose)
    cat("Returning",length(unique(c(merged$start1, merged$end1))),
        "merged breakpoints\n")
  return(merged)
}
