#' @title Merge close blocks
#'
#' @description
#' \code{merge_blocks} Checks for overlaps between blocks and merges accordingly.
#'
#' @param map the map object
#' @param blk the block object
#' @param verbose logical, should updates be printed?
#' @param max.size2merge numeric the maximum sized block to be merged
#' @param buffer the number of gene hits away from the bound of the block to
#' consider an overlap. For example, -1 means that the blocks must overlap by
#' at least 1 gene, and 9 means that blocks that are up to 9 genes apart will be
#' merged. See details
#' @param ... Not currently in use
#' @details This step is crucial to ensure that the MCScanX blocks ar in order.
#' By default, MCScanX will often join parts of inverted blocks, leading to
#' apparently duplicated regions, which are actually just overlapping
#' erroneously called block breakpoints. Because of this, we recommend running
#' merge_blocks with very large max.size2merge and max.iter. This will
#' ensure that all overlapping neighboring blocks are joined.
#'
#' In structurally diverged species comparisons, there may be many-many small
#' blocks that are adjacent, but broken becuase of gaps in BLAST hits or
#' real small inversions. To simplify the block structure, we also recommend
#' merging 'close' blocks that are separated by no more than the minimum
#' block size -1. Therefore, no merging over existing blocks will occur.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom geometry mesh.drectangle
#' @export
#'
#'
merge_blocks <- function(blk,
                         map,
                        buffer = 0,
                        verbose = T,
                        max.size2merge = 1e4){

  make_ovlMatrix = function(x, buffer){
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
          as.numeric(x[x$block.id == j,c("rankstart1","rankstart2")]),
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
    return(mat)
  }


  find_wh2merge = function(x, mat, buffer){

    wh <- which(mat <= buffer,
                arr.ind = T)

    if(nrow(wh) > 1){
      totsize <- apply(wh, 1, function(z)
        min(x$n.mapping[z]))
      totovl <- mat[mat<=buffer]
      wh = wh[order(totovl, totsize),]
      u = as.numeric(t(wh))
      u = unique(u[duplicated(u)])
      for(k in u){
        u2 = as.numeric(t(wh))
        u2 = unique(u2[duplicated(u2)])
        if(k %in% u2){
          rows = which(wh == k, arr.ind = T)[,1]
          rows = rows[rows != min(rows)]
          wh = wh[-rows,]
        }
      }
      if(length(wh) == 2){
        out = cbind(colnames(mat)[wh[1]],
                    colnames(mat)[wh[2]])

      }else{
        out = do.call(rbind,lapply(1:nrow(wh), function(k){
          colnames(mat)[wh[k,]]
        }))
      }
    }else{
      out = cbind(colnames(mat)[wh[1]],
                  colnames(mat)[wh[2]])
    }
    return(out)
  }


  blk$unique <- with(blk, paste(genome1, genome2, chr1, chr2))
  map$unique <- with(map, paste(genome1, genome2, chr1, chr2))
  splsmap <- split.data.table(map, "unique")
  splsblk <- split.data.table(blk, "unique")

  nmaps = sapply(splsblk, nrow)
  blk1 = unlist(sapply(splsblk, function(x) unique(x$block.id))[nmaps == 1])
  out.blk = blk[blk$block.id  %in% blk1,]
  out.map = map[map$block.id  %in% blk1,]

  in.blk = splsblk[names(nmaps)[nmaps > 1]]
  in.map = splsmap[names(nmaps)[nmaps > 1]]

  for(i in names(in.blk)){
    tmp.map = in.map[[i]]
    tmp.blk = in.blk[[i]]
    merge.it = TRUE
    if(verbose)
      cat(i,"\t merging",nrow(tmp.blk),"blocks ")
    while(merge.it){
      tblk = tmp.blk[tmp.blk$n.mapping <= max.size2merge, ]
      mat = make_ovlMatrix(tblk, buffer = buffer)
      merge.mat = find_wh2merge(x = tblk,  mat, buffer = buffer)
      if(is.null(merge.mat) | sum(complete.cases(merge.mat)) == 0){
        merge.it <- FALSE
      }else{
        for(j in 1:nrow(merge.mat)){
          mlj = merge.mat[j,]
          tmp.map$block.id[tmp.map$block.id %in% mlj] <- mlj[1]
          tmp.blk<-tmp.blk[tmp.blk$block.id != mlj[2],]
        }
      }
    }
    out.blk<-rbind(out.blk, tmp.blk)
    out.map<-rbind(out.map,tmp.map)
    if(verbose)
      cat("to",nrow(tmp.blk), "\n")
  }

  blk <- make_blocks(out.map)
  map <- data.table(blk$map)
  blk <- data.table(blk$block)

  return(list(block = blk, map = map))
}
