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
                        n.cores = 1,
                        max.size2merge = 1e6){

  make_ovlMatrix = function(x, buffer, n.cores){

    mat = do.call(cbind,mclapply(x$block.id, mc.cores = n.cores, function(i){
      sapply(x$block.id, function(j){
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
        return(o1)
      })
    }))
    colnames(mat)<-rownames(mat) <-x$block.id
    diag(mat) <- buffer + 1
    mat1 = mat
    mat2 = t(mat)
    mato = mat1
    mato[mato>mat2]<-mat2[mato>mat2]
    mato[lower.tri(mato)] <- buffer + 1
    return(mato)
  }

  find_all2keep = function(mat, buffer){
    wh = unique(as.numeric(which(mat<buffer,arr.ind = T)))
    return(unique(rownames(mat)[wh]))
  }



  find_wh2merge = function(mat, buffer){

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


  blk$unique = with(blk, paste(genome1, genome2, chr1, chr2))
  map$unique = with(map, paste(genome1, genome2, chr1, chr2))

  blk2check = blk$block.id[blk$n.mapping <= max.size2merge]
  n.inu = table(blk$unique)
  single = names(n.inu)[n.inu == 1]
  drop.these.blocks = blk$block.id[blk$unique %in% single]
  blk2check = blk2check[!blk2check %in% drop.these.blocks]
  k = 1
  while(length(blk2check > 0)){
    cat("Iteration",k,"... ",nrow(blk),"--> ")
    k = k+1
    sblk = split.data.table(blk[blk$block.id %in% blk2check,], "unique")
    smap = split.data.table(map[map$block.id %in% blk2check,], "unique")

    ovl.list = sapply(names(sblk),USE.NAMES = T, simplify = F, function(i){
      return(make_ovlMatrix(x = sblk[[i]], buffer = buffer, n.cores = 8))
    })

    do.merge = names(ovl.list[sapply(ovl.list,min) <= buffer])

    wh.list = sapply(do.merge, USE.NAMES = T, simplify = F, function(i){
        find_wh2merge(mat = ovl.list[[i]],
                      buffer = buffer)
    })

    blk2check = unique(unlist(lapply(do.merge, function(i)
      find_all2keep(mat = ovl.list[[i]],
                    buffer = buffer))))

    wh.out = do.call(rbind, wh.list)
    if(!is.null(wh.out)){
      for(i in 1:nrow(wh.out)){
        map$block.id[map$block.id %in% wh.out[i,]]<-wh.out[i,1]
      }
      out = make_blocks(map, rename.blocks = F, rerank = T)
      map = data.table(out$map)
      blk = data.table(out$block)
      blk$unique = with(blk, paste(genome1, genome2, chr1, chr2))
      map$unique = with(map, paste(genome1, genome2, chr1, chr2))

      n.inu = table(blk$unique)
      single = names(n.inu)[n.inu == 1]
      drop.these.blocks = blk$block.id[blk$unique %in% single]
      cat(nrow(blk),"blocks\n")
    }else{
      cat("no merging possible\n")
    }
  }
  return(list(map = map, block = blk))
}
