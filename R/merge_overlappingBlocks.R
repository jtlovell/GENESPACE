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
                                   verbose = T,
                                   buffer = 1){
  blk$uniq = paste0(blk$genome1,"_",blk$genome2)
  map$uniq = paste0(map$genome1,"_",map$genome2)

  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  tomergelist = lapply(1:nrow(blk), function(i){
    x = blk[i,]
    y = blk[-i,]
    yt = y[y$uniq == x$uniq &
             y$chr1 == x$chr1 &
             y$chr2 == x$chr2,]
    in1e = yt$rankend1>=x$rankstart1-buffer & yt$rankend1<=x$rankend1+buffer
    in1s = yt$rankstart1>=x$rankstart1-buffer & yt$rankstart1<=x$rankend1+buffer
    in2s = yt$rankend2>=x$rankstart2-buffer & yt$rankend2<=x$rankend2+buffer
    in2e = yt$rankstart2>=x$rankstart2-buffer & yt$rankstart2<=x$rankend2+buffer

    tobind = which((in1e | in1s) & (in2e | in2s))

    return(c(x$block.id,yt$block.id[tobind]))
  })
  names(tomergelist)<-blk$block.id
  if(verbose)
    cat("Merging adjacent blocks ... ")

  check.overlap = sapply(blk$block.id, function(i){
    x = tomergelist[[i]]
    g = sapply(tomergelist, function(j) any(j %in% x))
    return(unique(unlist(tomergelist[g])))
  })

  uoverl = check.overlap[!duplicated(check.overlap)]
  mapo = map
  for(i in names(uoverl)){
    mapo$block.id[map$block.id %in% uoverl[[i]]]<-i
  }

  out.blk = make_blocks(mapo)
  map = data.frame(mapo[["map"]], stringsAsFactors = F)
  blk = data.frame(out.blk[["block"]], stringsAsFactors = F)

  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = out.blk, map = map))
}
