#' @title Run MCScanX's detect_collinear_tandem_arrays
#'
#' @description
#' \code{find_pairwiseTandemArrays} Run MCScanX's detect_collinear_tandem_arrays
#'
#'
#' @param blast the blast dataset to screen for syntenic hits
#' @param ta.id tandem array id to plot
#' @param buffer Regions surrounding array to plot
#' @param plotit logical, should plots be made?
#' @param ... Passed on to plot
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
zoom_tandemArrays <- function(ta.id,
                              blast,
                              buffer = 20,
                              plotit = T,
                              ...){

  ta.hits <- blast[blast$tandemarray.id == ta.id,]
  if(length(unique(ta.hits$chr1)) > 1){
    wh.most = table(ta.hits$chr1)
    wh.most = wh.most[order(-wh.most)]
    wh.most = names(wh.most)[1]
    ta.hits = ta.hits[ta.hits$chr1 == wh.most,]
  }
  if(length(unique(ta.hits$chr2)) > 1){
    wh.most = table(ta.hits$chr2)
    wh.most = wh.most[order(-wh.most)]
    wh.most = names(wh.most)[1]
    ta.hits = ta.hits[ta.hist$chr2 == wh.most,]
  }
  order1.start = with(ta.hits, min(order1))
  order2.start = with(ta.hits, min(order2))
  order1.end = with(ta.hits, max(order1))
  order2.end = with(ta.hits, max(order2))
  pos1.start = with(ta.hits, min(start1))
  pos2.start = with(ta.hits, min(start2))
  pos1.end = with(ta.hits, max(end1))
  pos2.end = with(ta.hits, max(end2))
  c1 = with(ta.hits, chr1[1])
  c2 = with(ta.hits, chr2[1])
  g1 = with(ta.hits, genome1[1])
  g2 = with(ta.hits, genome2[1])
  og = with(ta.hits, og1[1])
  buf.hits = blast[with(blast,
                        genome1 == g1 & genome2 == g2 &
                          chr1 == c1 & chr2 == c2 &
                          order1 >= (order1.start - buffer) &
                          order2 >= (order2.start - buffer) &
                          order1 <= (order1.end + buffer) &
                          order2 <= (order2.end + buffer)),]
  n.ta.1 <- length(unique(ta.hits$id1))
  n.ta.2 <- length(unique(ta.hits$id2))
  out = data.table(ta.id = ta.id,
                   genome1 = g1,
                   genome2 = g2,
                   chr1 = c1,
                   chr2 = c2,
                   start1 = pos1.start,
                   start2 = pos2.start,
                   end1 = pos1.end,
                   end2 = pos2.end,
                   rankstart1 = order1.start,
                   rankstart2 = order2.start,
                   rankend1 = order1.end,
                   rankend2 = order2.end,
                   n.genes1 = n.ta.1,
                   n.genes2 = n.ta.2,
                   og = og,
                   stringsAsFactors = F)
  if(plotit){
    with(buf.hits, plot(order1, order2,
                        pch = 16, col = "black",
                        xlab = paste0(g1, " gene order (n = ",n.ta.1,")"),
                        ylab = paste0(g2, " gene order (n = ",n.ta.2,")"),
                        main = ta.id, asp = 1, ...))
    with(ta.hits, points(order1, order2, pch = 16,
                         col = "red", cex = .6))
  }
  return(out)
}
