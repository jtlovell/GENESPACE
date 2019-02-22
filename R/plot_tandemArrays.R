#' @title plot tandem arrays
#'
#' @description
#' \code{plot_tandemArrays} plot tandem arrays and nearby
#' syntenic orthologous blast hits
#'
#' @param map the map data.table, with an additional column 'array.id'
#' @param array.id a vector of tandem array ids to plot
#' @param buffer Regions surrounding array to plot
#' @param what Character, either 'bp' for phyical position plot, or
#' 'rank' for gene rank space
#' @param col Color for the points
#' @param pch Shape for the points
#' @param array.pch passed to pch in points for hits in the array
#' @param array.col passed to col in points for hits in the array
#' @param array.cex passed to cex in points for hits in the array
#' @param ... Additional argument assed on to plot
#'
#' @details A way to look and and define tandem array regions.
#'
#' @return The coordinates of the tandem array,
#' and, if plotit, a zoomed-in plot.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
plot_tandemArrays <- function(array.id,
                              map,
                              buffer = 1e5,
                              pch = 16,
                              col = "black",
                              array.pch = 16,
                              array.col = "red",
                              array.cex = .6,
                              what = "bp",
                              ...){

  if (!what %in% c("bp","rank"))
    stop("what must be either bp or rank\n")

  if (!"array.id" %in% colnames(map))
    stop("could not find a column id named array.id, are you sure you rank find_tandemArrays?\n")

  if (!array.id %in% map$array.id)
    stop("could not find", array.id, "in the array.id column of the map data.table\n")

  wh <- which(map$array.id == array.id)
  ta.hits <- map[wh, ]

  if (length(unique(ta.hits$chr1)) > 1) {
    wh.most <- table(ta.hits$chr1)
    wh.most <- wh.most[order(-wh.most)]
    wh.most <- names(wh.most)[1]
    ta.hits <- ta.hits[ta.hits$chr1 == wh.most,]
  }
  if (length(unique(ta.hits$chr2)) > 1) {
    wh.most <- table(ta.hits$chr2)
    wh.most <- wh.most[order(-wh.most)]
    wh.most <- names(wh.most)[1]
    ta.hits <- ta.hits[ta.hits$chr2 == wh.most,]
  }


  c1 <- with(ta.hits, chr1[1])
  c2 <- with(ta.hits, chr2[1])
  g1 <- with(ta.hits, genome1[1])
  g2 <- with(ta.hits, genome2[1])
  og <- with(ta.hits, og1[1])


  if (what == "bp") {
    pos1.start <- with(ta.hits, min(start1))
    pos2.start <- with(ta.hits, min(start2))
    pos1.end <- with(ta.hits, max(end1))
    pos2.end <- with(ta.hits, max(end2))

    buf.hits <- map[with(map,
                        genome1 == g1 & genome2 == g2 &
                          chr1 == c1 & chr2 == c2 &
                          end1 >= (pos1.start - buffer) &
                          end2 >= (pos2.start - buffer) &
                          start1 <= (pos1.end + buffer) &
                          start2 <= (pos2.end + buffer)),]
  } else {
    order1.start <- with(ta.hits, min(rank1))
    order2.start <- with(ta.hits, min(rank2))
    order1.end <- with(ta.hits, max(rank1))
    order2.end <- with(ta.hits, max(rank2))

    buf.hits <- map[with(map,
                        genome1 == g1 & genome2 == g2 &
                          chr1 == c1 & chr2 == c2 &
                          rank1 >= (order1.start - buffer) &
                          rank2 >= (order2.start - buffer) &
                          rank1 <= (order1.end + buffer) &
                          rank2 <= (order2.end + buffer)),]
  }

  n.ta.1 <- length(unique(ta.hits$id1))
  n.ta.2 <- length(unique(ta.hits$id2))

  if (what == "bp") {
    with(buf.hits, plot(start1, start2,
                        pch = pch,
                        col = col,
                        xlab = paste0(g1, " (", c1, ") gene position (n = ", n.ta.1,")"),
                        ylab = paste0(g2, " (", c1, ") gene position (n = ", n.ta.2,")"),
                        main = array.id[1],
                        ...))
    with(ta.hits, points(start1, start2,
                         pch = array.pch,
                         col = array.col,
                         cex = array.cex))
  } else {
    with(buf.hits, plot(order1, order2,
                        pch = pch,
                        col = col,
                        xlab = paste0(g1, " (", c1, ") gene order (n = ", n.ta.1, ")"),
                        ylab = paste0(g2, " (", c1, ") gene order (n = ", n.ta.2, ")"),
                        main = array.id[1],
                        ...))
    with(ta.hits, points(order1, order2,
                         pch = array.pch,
                         col = array.col,
                         cex = array.cex))
  }
}
