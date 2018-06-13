#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{make_inputFileMatrix} Utility function to build metadata
#'
#' @param map The map object (data.frame or data.table)
#' @param blk The block object (data.frame)
#' @param buffer Numeric, the overlapping distance between two blocks.
#' 0 indicates that blocks that overlap by >=0 should be merged.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Primarily used in the run_MCScanX pipeline.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
merge_overlappingBlocks = function(map, blk, buffer = 1.5,
                                   verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  if(verbose)
    cat("Looking for overlapping blocks ...\n")

  inrect = function(left,top,right, bottom, points, buffer = 1){
    apply(points,1, function(k){
      k[1]>(left+buffer) &
        k[1]<(right-buffer) &
        k[2]>(bottom+buffer) &
        k[2]<(top-buffer)
    })
  }

  blk$block.id = as.character(blk$block.id)
  blk$uniq = paste0(blk$genome1,"_",blk$genome2,"_",blk$chr1, "_", blk$chr2)
  spl = split(blk, blk$uniq)
  nr = sapply(spl, nrow)
  spl = spl[nr>1]
  for(i in names(spl)){
    x = spl[[i]]
    x = x[order(x$start1, x$start2),]

    has.ovl = sapply(1:nrow(x), function(y)
      rowSums(cbind(inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("start1","start2")]),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("end1","start2")]),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("start1","end2")]),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("end1","end2")]))))
    colnames(has.ovl)<-x$block.id
    rownames(has.ovl)<-x$block.id
    has.ovl = has.ovl>1

    if(sum(has.ovl)>0){
      if(sum(colSums(has.ovl)>0)==1){
        cn = colnames(has.ovl)[colSums(has.ovl)>0]
        wh = rownames(has.ovl)[which(has.ovl[,colSums(has.ovl)>0])]
        mlist = list(wh)
        names(mlist) = cn
      }else{
        mlist = apply(has.ovl[,colSums(has.ovl)>0],2,function(y) names(which(y)))
      }
      for(j in 1:length(mlist)){
        map$block.id[map$block.id %in% mlist[[j]]]<-names(mlist)[j]
      }
    }
  }
  out = make_blocks(map)
  map = data.frame(out[["map"]], stringsAsFactors = F)
  blk = data.frame(out[["block"]], stringsAsFactors = F)
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}
