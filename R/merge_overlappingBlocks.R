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
merge_overlappingBlocks = function(map, blk,
                                   verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  if(verbose)
    cat("Looking for overlapping blocks ... ")

  run.it = function(map,blk){
    blk$block.id = as.character(blk$block.id)
    blk$uniq = paste0(blk$genome1,"_",blk$genome2,"_",blk$chr1, "_", blk$chr2)
    spl = split(blk, blk$uniq)

    for(i in names(spl)){
      x = spl[[i]]
      x = x[order(x$start1, x$start2),]
      x$s1 = frank(x$start1, ties.method = "dense")
      x$s2 = frank(x$start2, ties.method = "dense")
      x$e1 = frank(x$end1, ties.method = "dense")
      x$e2 = frank(x$end2, ties.method = "dense")

      x$tomerge = apply(x[,c("s1","e1")],1,function(y) length(unique(y))>1) &
        apply(x[,c("s2","e2")],1,function(y) length(unique(y))>1)
      mlist = lapply(which(x$tomerge), function(i)
        unique(as.numeric(x[i,c("s1","e1")]))[order(unique(as.numeric(x[i,c("s1","e1")])))])
      mlist = mlist[!duplicated(mlist)]
      if(length(mlist)>0){
        mo = lapply(mlist, function(y) x$block.id[y])
        for(j in 1:length(mo)){
          map$block.id[map$block.id %in% mo[[j]]]<-mo[[j]][1]
        }
      }
    }
    out = make_blocks(map)
    map = data.frame(out[["map"]], stringsAsFactors = F)
    blk = data.frame(out[["block"]], stringsAsFactors = F)
    return(list(block = blk, map = map, novl = novl))
  }

  novl = 1
  iter = 1
  while(novl>0){
    ninit = nrow(blk)
    if(verbose)
      cat("Pass", iter,"(",ninit,")\t")
    tmp = run.it(map = map, blk = blk)
    nafter = nrow(tmp$block)
    novl= ninit-nafter
    if(verbose)
      cat("n overlaps = ", novl,"\n")
    if(novl>0){
      map = tmp$map
      blk = tmp$block
    }
  }
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}
