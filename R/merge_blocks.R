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
#' @param clean.columns Should column names be cleaned out?
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
#' @importFrom parallel mclapply
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom compiler cmpfun
#' @export
merge_blocks = function(map,
                        blk,
                        verbose = T,
                        buffer = -1,
                        clean.columns = T,
                        max.iter = 10){
  #######################################################
  #######################################################
  merge_clusters = function(map, blk, buffer){
    b = blk
    m = map

    b$unique = with(b, paste0(genome1,genome2, chr1, chr2))
    b$cluster1 = find_blkCluster(blk = b,
                                 seqnames.field = "unique",
                                 start.field = "rankstart1",
                                 end.field = "rankend1",
                                 buffer = buffer)
    b$cluster2 = find_blkCluster(blk = b,
                                 seqnames.field = "unique",
                                 start.field = "rankstart2",
                                 end.field = "rankend2",
                                 buffer = buffer)

    b[,new.block.id := head(block.id,1),
      by = list(cluster1, cluster2)]
    b.clus = b[,c("new.block.id","block.id")]
    b.clus = b.clus[!duplicated(b.clus),]
    setkey(m, block.id)
    setkey(b.clus, block.id)
    m.clus = merge(b.clus, m)
    m.clus$block.id <- m.clus$new.block.id
    m.clus$new.block.id <- NULL
    out = make_blocks(map = m.clus,
                      rename.blocks = F,
                      rerank = T,
                      clean.columns = clean.columns,
                      ties.method = "dense")
    return(out)
  }
  #######################################################
  #######################################################
  find_blkCluster = function(blk,
                             seqnames.field = "chr",
                             start.field = "start",
                             end.field = "end",
                             buffer = 0){

    gr <- makeGRangesFromDataFrame(blk,
                                   seqnames.field = seqnames.field,
                                   ignore.strand = T,
                                   start.field = start.field,
                                   end.field = end.field,
                                   keep.extra.columns = FALSE)
    ovl = findOverlaps(gr, gr,
                       ignore.strand = T,
                       maxgap = buffer,
                       select = "arbitrary",
                       type = "any")

    cluster <- frank(ovl,
                     ties.method = "dense")
    return(cluster)
  }
  #######################################################
  #######################################################
  find_blkCluster <- cmpfun(find_blkCluster)
  merge_clusters <- cmpfun(merge_clusters)
  #######################################################
  #######################################################

  old.blk = blk
  n.old.blk <- nrow(blk)
  n.new.blk <- 0
  n.iter = 0


  while((n.new.blk<n.old.blk | n.iter == 0) & n.iter < max.iter){
    n.iter = n.iter + 1
    cat("\tIteration", paste0(n.iter,":"),
        nrow(blk),"blocks ")
    n.old.blk = nrow(blk)
    rund = merge_clusters(map = map,
                          blk = blk, buffer = buffer)
    blk <- rund$block
    map <- rund$map
    n.new.blk = nrow(blk)
    cat("-->", nrow(blk),"\n")
  }
  cat("\tDone!\n")

  return(list(map = map, block = blk))

}
