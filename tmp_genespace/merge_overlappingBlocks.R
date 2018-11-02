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
merge_overlappingBlocks = function(map, blk, buffer = 1.5, verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  if(verbose)
    cat("Looking for overlapping blocks ...\n")
  map$block.id = paste(map$block.id,map$genome1, map$genome2, map$chr1, map$chr2)
  o = make_blocks(map)
  map = o$map
  blk = o$block
  spl = split(blk, paste(blk$genome1, blk$genome2, blk$chr1, blk$chr2))
  spl2 = split(map, paste(map$genome1, map$genome2, map$chr1, map$chr2))
  test = lapply(names(spl), function(i){
    x = spl[[i]]
    x = x[order(-x$n.mapping),]

    y = spl2[[i]]

    for(j in x$block.id){
      if(j %in% y$block.id){
        wh = with(x, which(
          ((rankstart1 + buffer) >= rankstart1[block.id == j]  &
             (rankstart2 + buffer) >= rankstart2[block.id == j]  &
             (rankstart1 - buffer) <= rankend1[block.id == j]  &
             (rankstart2 - buffer)  <= rankend2[block.id == j]) |
            ((rankend1 + buffer) >= rankstart1[block.id == j] &
               (rankend2 + buffer) >= rankstart2[block.id == j]  &
               (rankend1 - buffer)  <= rankend1[block.id == j] &
               (rankend2 - buffer)  <= rankend2[block.id == j])))
        if(length(wh)>0){
          tomerge = x$block.id[wh]
          y$block.id[y$block.id %in% tomerge]<-j
          o = make_blocks(y)
          y = o$map
          x = o$block
        }
      }
    }
    return(y)
  })
  tmp = rbindlist(test)
  tmp$block.id = as.numeric(as.factor(tmp$block.id))
  out = make_blocks(tmp)

  map = data.frame(out[["map"]], stringsAsFactors = F)
  blk = data.frame(out[["block"]], stringsAsFactors = F)
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}
