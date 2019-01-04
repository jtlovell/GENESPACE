#' @title Make detailed dotplot of mapping and blocks
#'
#' @description
#' \code{plot_blocksAndMapping} For each pairwise comparison between chromosomes
#' between two genomes, make a dotplot, with line segments connecting mappings
#' in the same block
#'
#' @param map The map data.table or data.frame
#' @param blk The block data.table or data.frame
#' @param ref.id Character of length one, indicating the genome to be plotted
#' on the x axis
#' @param altGenome2plot Character of length one, indicating the genome to be plotted
#' on the y axis
#' @param chr1toplot Chromosomes in ref.id to be plotted
#' @param chr2toplot Chromosomes in altGenome2plot to be plotted
#' @param main The title of the plot
#' @param colorSegment Logical, should the segments be colored (TRUE) or the
#' points (FALSE).
#' @param ... Other arguments passed to plot
#' @details More here
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
plot_blocksAndMapping <- function(map,
                                  blk,
                                  ref.id,
                                  altGenome2plot,
                                  chr1toplot = NULL,
                                  chr2toplot = NULL,
                                  main = NULL,
                                  colorSegment = T,
                                  ...){

  blk <- data.frame(blk)
  map <- data.frame(map)

  tpb <- blk[blk$genome1 == ref.id &
               blk$genome2 == altGenome2plot,]
  tp <- map[map$genome1 == ref.id &
              map$genome2 == altGenome2plot,]
  tpb$s1 <- with(tpb,
                 ifelse(orient == "+", rankstart1, rankend1))
  tpb$e1 <- with(tpb,
                 ifelse(orient == "+", rankend1, rankstart1))
  tpb$rankstart1 <- tpb$s1
  tpb$rankend1 <- tpb$e1
  tpb$s1 <- NULL
  tpb$e1 <- NULL
  if (!is.null(chr1toplot)) {
    tpb <- tpb[tpb$chr1 %in% chr1toplot,]
    tp <- tp[tp$chr1 %in% chr1toplot,]
  }
  if (!is.null(chr2toplot)) {
    tpb <- tpb[tpb$chr2 %in% chr2toplot,]
    tp <- tp[tp$chr2 %in% chr2toplot,]
  }
  sb1 <- split(tpb, tpb$chr1)
  st1 <- split(tp, tp$chr1)
  for (i in names(sb1)) {
    sb2 <- split(sb1[[i]], sb1[[i]]$chr2)
    st2 <- split(st1[[i]], st1[[i]]$chr2)
    for (j in names(sb2)) {
      t2 <- st2[[j]]
      b2 <- sb2[[j]]
      cols <- rep_len(c("red3", "salmon", "darkorange3", "gold",
                        "grey50", "lightgreen", "forestgreen", "darkgreen",
                        "cyan", "dodgerblue3", "violet", "purple"),
                      length.out = nrow(b2))
      if (colorSegment) {
        with(t2, plot(rank1, rank2,
                      col = rgb(0,0,0,.5),
                      pch = 16,
                      cex = .5,
                      ylab = paste(altGenome2plot, j, "gene order"),
                      xlab = paste(ref.id, i, "gene order"),
                      main = main, ...))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 = rankend2,
                          col = cols[as.numeric(as.factor(block.id))],
                          lwd = 2))
      }else{
        with(t2, plot(rank1, rank2,
                      col = cols[as.numeric(as.factor(block.id))],
                      pch = 16,
                      cex = .5,
                      ylab = paste(altGenome2plot, j, "gene order"),
                      xlab = paste(ref.id, i, "gene order"),
                      main = main, ...))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 = rankend2,
                          col = "black",
                          lwd = 1.5))
      }

      with(b2, text(x = rowMeans(b2[, c("rankstart1", "rankend1")]),
                    y = rowMeans(b2[, c("rankstart2", "rankend2")]),
                    labels = block.id,
                    col = "black",
                    cex = .5,
                    adj = c(1, -1)))
    }
  }
}

