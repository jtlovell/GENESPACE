#' @title Merge close blocks
#'
#' @description
#' \code{merge_blocks} Checks for overlaps between blocks and merges accordingly.
#'
#' @param map the map object
#' @param blk the block object
#' @param verbose logical, should updates be printed?
#' @param reciprocal logical, If TRUE, do not merge blocks entirely within another
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
merge_blocks = function(blk,
                        map,
                        buffer = 0,
                        verbose = T,
                        reciprocal = FALSE){

  find_2merge = function(x, buffer, n.match){
    if(nrow(x)==1){
      return(NULL)
    }else{
      in.end = sapply(1:nrow(x), function(i){
        (x$rankend1[i] + buffer) > (x$rankstart1) &
          (x$rankend1[i] - buffer) < (x$rankend1) &
          (x$rankend2[i] + buffer) > (x$rankstart2) &
          (x$rankend2[i] - buffer) < (x$rankend2)
      })
      in.start = sapply(1:nrow(x), function(i){
        (x$rankstart1[i] + buffer) > (x$rankstart1) &
          (x$rankstart1[i] - buffer) < (x$rankend1) &
          (x$rankstart2[i] + buffer) > (x$rankstart2) &
          (x$rankstart2[i] - buffer) < (x$rankend2)
      })
      on.end = sapply(1:nrow(x), function(i){
        (x$rankend1[i] + buffer) > (x$rankstart1) &
          (x$rankend1[i] - buffer) < (x$rankend1) &
          (x$rankstart2[i] + buffer) > (x$rankstart2) &
          (x$rankstart2[i] - buffer) < (x$rankend2)
      })
      on.start = sapply(1:nrow(x), function(i){
        (x$rankstart1[i] + buffer) > (x$rankstart1) &
          (x$rankstart1[i] - buffer) < (x$rankend1) &
          (x$rankend2[i] + buffer) > (x$rankstart2) &
          (x$rankend2[i] - buffer) < (x$rankend2)
      })
      in.both = in.end + t(in.start)
      on.both = on.end + t(on.start)
      diag(in.both)<-0
      diag(on.both)<-0
      in.both <- in.both >= n.match
      on.both <- on.both >= n.match
      if(!any(in.both) & !any(on.both)){
        return(NULL)
      }else{
        if(!any(in.both)){
          wh = which(on.both, arr.ind = TRUE)[1,]
          to.merge = x$block.id[unique(as.numeric(wh))]
          return(to.merge)
        }else{
          wh = which(in.both, arr.ind = TRUE)[1,]
          to.merge = x$block.id[unique(as.numeric(wh))]
          return(to.merge)
        }
      }
    }
  }

  remap_merge = function(blk, map, verbose = TRUE,
                         buffer, n.match){
    if(verbose)
      cat("Checking", nrow(blk), "block coordinates\n\t")
    spl = split(blk, with(blk, paste(genome1, genome2, chr1, chr2)))

    merge.list = lapply(spl, function(x) find_2merge(x = x,
                                                     buffer = buffer,
                                                     n.match = n.match))
    merge.list <- merge.list[!sapply(merge.list, is.null)]
    if(length(merge.list)>1){
      merge.list = merge.list[sapply(merge.list, length)>1]

      for(i in 1:length(merge.list)){
        map$block.id[map$block.id %in% merge.list[[i]]]<-merge.list[[i]][1]
      }

      blk = make_blocks(map)
      map = blk$map
      blk = blk$block

      if(verbose)
        cat("Returning", nrow(blk), "merged blocks\n")
    }else{
      if(verbose)
        cat("Found no blocks to merge\n")
    }

    return(list(block = blk, map = map))
  }


  if(reciprocal){
    n.match = 2
  }else{
    n.match = 1
  }

  nblk = nrow(blk)+1
  i = 0
  while(i <= max.iter &
        nblk > nrow(blk)){
    i = i + 1
    if(verbose)
      cat("Iteration",i,"... ")
    nblk = nrow(blk)
    tmp = remap_merge(blk = blk,
                      map = map,
                      n.match = n.match,
                      buffer = buffer)

    blk = tmp$block
    map = tmp$map
  }


  if(verbose)
    cat("Done! Returning a dataset with", nrow(blk), "blocks\n")

  return(list(block = blk, map = map))
}
