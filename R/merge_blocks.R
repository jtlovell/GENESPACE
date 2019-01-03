#' @title Merge close blocks
#'
#' @description
#' \code{merge_blocks} Checks for overlaps between blocks and merges accordingly.
#'
#' @param map the map data.table or data.frame
#' @param blk the block data.table or data.frame
#' @param max.size2merge numeric the maximum sized block to be merged
#' @param buffer the number of gene hits away from the bound of the block to
#' consider an overlap. For example, -1 means that the blocks must overlap by
#' at least 1 gene, and 9 means that blocks that are up to 9 genes apart will be
#' merged. See details
#' @param n.iter the number of iterations to run
#' @param n.cores The number of parallel processes to run.
#' @param ignore.orient Logical, should orientation of blocks be ignored?
#' @param verbose logical, should updates be printed?
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
#' @import parallel
#' @export
merge_blocks <- function(blk,
                         map,
                         max.size2merge = 1e6,
                         buffer = 0,
                         n.iter = 5,
                         verbose = T,
                         ignore.orient = T,
                         n.cores = 1){

  blk <- data.table(blk)
  map <- data.table(map)

  buffer.b <- buffer.a <- buffer


  check_ifBlksOvlp <- function(a, b, buffer.a, buffer.b){

    check_ovlp <- function(a0, a1, b0, b1) {
      (a0 <= b1 & b1 <= a1) ||
        (a0 <= b0 & b0 <= a1) ||
        (a0 >= b0 & a1 <= b1) | (a0 <= b0 & a1 >= b1)
    }

    ovlp.x <- check_ovlp(a[1] - buffer.a,
                         a[3] + buffer.a,
                         b[1] - buffer.b,
                         b[3] + buffer.b)

    ovlp.y <- check_ovlp(a[2] - buffer.a,
                         a[4] + buffer.a,
                         b[2] - buffer.b,
                         b[4] + buffer.b)

    return(ovlp.x & ovlp.y)
  }


  make_ovlMatrix <- function(block,
                             buffer.a,
                             buffer.b,
                             n.cores = 1,
                             verbose = T){

    mbuff <- max(c(buffer.a,
                   buffer.b))

    if (verbose)
      cat("Generating list of block coordinates\n")

    block$unique.chr <- with(block,
                             paste(genome1, genome2,
                                   chr1, chr2))

    spl.block <- split(block, "unique.chr")

    spl.coord <- mclapply(spl.block, mc.cores = n.cores, function(x){
      tmp <- data.matrix(x[,c("rankstart1","rankstart2",
                              "rankend1","rankend2")])
      rownames(tmp) <- x$block.id
      return(tmp)
    })

    if (verbose)
      cat("Generating distance matrix between blocks\n")

    ovl.out <- mclapply(spl.coord, mc.cores = n.cores, function(x){
      if (nrow(x) == 1) {
        out <- matrix(FALSE)
        rownames(out) <- colnames(out) <- rownames(x)
        return(out)
      }else{
        spl.mat <- sapply(rownames(x), function(i){
          sapply(rownames(x), function(j){
            check_ifBlksOvlp(a = x[i,],
                             b = x[j,],
                             buffer.a = buffer.a,
                             buffer.b = buffer.b)
          })
        })

        rownames(spl.mat) <- colnames(spl.mat) <- rownames(x)
        mat1 <- spl.mat
        mat2 <- t(spl.mat)
        mato <- mat1
        mato[mat2] <- TRUE
        diag(mato) <- FALSE
        return(mato)
      }
    })
    return(ovl.out)
  }

  find_wh2merge <- function(mat.list,
                            n.cores = 1){
    mclapply(mat.list,  mc.cores = n.cores, function(mat){

      wh <- which(mat,
                  arr.ind = T)

      if (nrow(wh) > 1) {
        u <- as.numeric(t(wh))
        u <- unique(u[duplicated(u)])
        for (k in u) {
          u2 <- as.numeric(t(wh))
          u2 <- unique(u2[duplicated(u2)])
          if (k %in% u2) {
            rows <- which(wh == k,
                          arr.ind = T)[,1]
            rows <- rows[rows != min(rows)]
            wh <- wh[-rows,]
          }
        }
        if (length(wh) == 2) {
          out <- cbind(colnames(mat)[wh[1]],
                      colnames(mat)[wh[2]])

        }else{
          out <- do.call(rbind,lapply(1:nrow(wh), function(k){
            colnames(mat)[wh[k,]]
          }))
        }
      }else{
        out <- cbind(colnames(mat)[wh[1]],
                    colnames(mat)[wh[2]])
      }
      return(out)
    })
  }

  process_whlistMap <- function(wh.list,
                                map,
                                blk){
    wh.out <- data.table(do.call(rbind, wh.list))
    wh.out <- wh.out[complete.cases(wh.out), ]
    wh.tmp <- data.table(wh.out)
    wh.tmp$V1 <- wh.tmp$V2
    wh.out <- rbind(wh.out, wh.tmp)
    setnames(wh.out, c("block.id", "new.block"))

    ob <- blk[!blk$block.id %in% c(wh.out$block.id, wh.out$new.block),]
    wh.out <- rbind(data.frame(block.id = ob$block.id,
                               new.block = ob$block.id),
                    wh.out)

    wh.out$block.id <- as.character(wh.out$block.id)
    map$block.id <- as.character(map$block.id)
    setkey(map, block.id)
    setkey(wh.out, block.id)

    m <- merge(wh.out, map)
    m$block.id <- m$new.block
    m$new.block <- NULL
    out <- make_blocks(m,
                       rename.blocks = T,
                       rerank = T)
    return(out)
  }

  if (verbose)
    cat("Merging blocks with an gene rank overlap of",
        buffer.a,
        "...\n\t")
  for (i in 1:n.iter) {
    if(verbose)
      cat("Iteration", i, "...",
          nrow(blk), "--> ")
    if (ignore.orient) {
      ovl.out <- make_ovlMatrix(block = blk[blk$n.mapping <= max.size2merge,],
                                n.cores = n.cores,
                                buffer.a = buffer.a,
                                buffer.b = buffer.b,
                                verbose = F)
      wh.list <- find_wh2merge(ovl.out,
                               n.cores = n.cores)
    }else{
      ovl.out <- make_ovlMatrix(block = blk[blk$n.mapping <= max.size2merge &
                                              blk$orient == "+",],
                                n.cores = n.cores,
                                buffer.a = buffer.a,
                                buffer.b = buffer.b,
                                verbose = F)
      wh.list1 <- find_wh2merge(ovl.out,
                                n.cores = n.cores)
      ovl.out <- make_ovlMatrix(block = blk[blk$n.mapping <= max.size2merge &
                                              blk$orient == "-",],
                                n.cores = n.cores,
                                buffer.a = buffer.a,
                                buffer.b = buffer.b,
                                verbose = F)
      wh.list2 <- find_wh2merge(ovl.out,
                                n.cores = n.cores)
      for (i in unique(names(wh.list1), names(wh.list2))) {
        wh.list[[i]] <- rbind(wh.list1[[i]], wh.list2[[i]])
      }
    }

    mb <- process_whlistMap(wh.list = wh.list,
                            map = map, blk = blk)
    map <- mb$map
    blk <- mb$block

    if (verbose)
      cat(nrow(blk),"\n\t")
  }

  if (verbose)
    cat("Done!\n")
  return(mb)
}
