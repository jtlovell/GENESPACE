#' @title Merge close blocks
#'
#' @description
#' \code{merge_blocks} Checks for overlaps between blocks and merges accordingly.
#'
#' @param map the map object
#' @param blk the block object
#' @param verbose logical, should updates be printed?
#' @param max.iter numeric, the maximum number of iterations to perform to look for
#' blocks to merge.
#' @param max.size2merge numeric the maximum sized block to be merged
#' @param ... Not currently in use
#' @details Needs to be run prior to the pipeline. Makes some objects that are required.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom geometry mesh.drectangle
#' @export
merge_blocks <- function(blk,
                        map,
                        buffer = 0,
                        verbose = T,
                        max.size2merge = 1e4,
                        max.iter = 2){

  find_2merge <- function(x, buffer){
    if(nrow(x) == 1){
      return(NULL)
    }else{
      mat <- matrix(NA,
                    nrow = nrow(x),
                    ncol = nrow(x))
      rownames(mat) <- colnames(mat) <- x$block.id
      for(i in rownames(mat)){
        for(j in rownames(mat)){
          p <- x[x$block.id == i,]
          o1 <- min(mesh.drectangle(p = rbind(
            as.numeric(x[x$block.id == j,c("rankstart1","rankstart2")]),
            as.numeric(x[x$block.id == j,c("rankend1","rankend2")]),
            as.numeric(x[x$block.id == j,c("rankend1","rankstart2")]),
            as.numeric(x[x$block.id == j,c("rankstart1","rankend2")]),
            as.numeric(x[x$block.id == j,c("rankstart1","rankend2")]),
            as.numeric(x[x$block.id == j,c("rankend1","rankstart2")])),
            x1 = p$rankstart1,
            x2 = p$rankend1,
            y1 = p$rankstart2,
            y2 = p$rankend2))
          mat[i,j] <- o1
        }
      }
      diag(mat) <- buffer + 1
      mat[lower.tri(mat)] <- buffer + 1
      if(all(as.numeric(mat) > buffer)){
        return(NULL)
      }else{
        wh <- which(mat == min(mat),
                    arr.ind = T)
        if(nrow(wh) > 1){
          totsize <- apply(wh, 1, function(z)
            min(x$n.mapping[z]))
          wh <- wh[which.min(totsize), ]
        }
        return(colnames(mat)[as.vector(wh)])
      }
    }
  }

  remap_merge <- function(blk,
                          map,
                          verbose = TRUE,
                          buffer,
                          n.match){
    if(verbose)
      cat("n. blocks:", nrow(blk), "--> ")
    blk.merge <- blk[blk$n.mapping <= max.size2merge,]
    blk.merge = data.table(blk.merge)
    blk.merge$uniq = with(blk.merge, paste(genome1, genome2, chr1, chr2))
    spl <- split.data.table(blk.merge, "uniq")
    merge.list <- lapply(spl, function(x)
      find_2merge(x = x,
                  buffer = buffer))
    merge.list <- merge.list[!sapply(merge.list, is.null)]
    if(length(merge.list) > 1){
      merge.list <- merge.list[sapply(merge.list, length) > 1]

      for(i in 1:length(merge.list)){
        ml <- merge.list[[i]]
        map$block.id[map$block.id %in% ml] <- ml[1]
      }

      blk <- make_blocks(map)
      map <- blk$map
      blk <- blk$block

      if(verbose)
        cat(nrow(blk), "\n")
    }else{
      if(verbose)
        cat("Found no blocks to merge\n")
    }

    return(list(block = blk, map = map))
  }

  smap = map
  sblk = blk

  sblk$unique <- with(sblk, paste(genome1, genome2))
  smap$unique <- with(smap, paste(genome1, genome2))
  splsmap <- split.data.table(smap, "unique")
  splsblk <- split.data.table(sblk, "unique")

  out = lapply(names(splsmap), function(j){
    blk = splsblk[[j]]
    map = splsmap[[j]]

    nblk <- nrow(blk) + 1
    i <- 0
    unmerge_unique <- vector()

    while(i < max.iter &
          nblk > nrow(blk)){
      i <- i + 1
      if(verbose)
        cat("\tIteration", i," ")
      nblk <- nrow(blk)
      blk$unique <- NULL
      map$unique <- NULL
      if(nrow(blk) == 1){
        tmp <- list(block = blk,
                    map = map)
      }else{
        tmp <- remap_merge(blk = blk,
                           map = map,
                           n.match = n.match,
                           buffer = buffer)
      }
      blk <- tmp$block
      map <- tmp$map
    }

    return(list(block = blk, map = map))
  })

  blk = rbindlist(lapply(out, function(x) x$block))
  map = rbindlist(lapply(out, function(x) x$map))

  return(list(block = blk, map = map))
}
