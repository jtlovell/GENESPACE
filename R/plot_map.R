#' @title build pwBlocks
#'
#' @description
#' \code{build_pwBlocks} build_pwBlocks
#'
#' @param map data.table, format similar to blast
#' @param gff data.table, format similar to gff
#' @param genomeIDs character, genomes to plot
#' @param chr.list list, named by genome labels (label will
#' replace the genomes element). Must be same length as genomes.
#' Character vector must contain chromosome IDs within the genomes.
#' @param gap.prop numeric, proportion of teh genome to use for gaps
#' between chromosomes
#' @param palette function, palette function for colors
#' @param use.rank logical, should ranks or bp coordinates be plotted?
#' @param col.byBlk logical, should the blocks be colored independently?
#' @param label.blocks logical, should text labeling the blocks
#' be added to the plots?
#' @param pch numeric, passed on to plot.
#' @param text.cex numeric, size of the block labels, passed on to plot.
#' @param cex.axis numeric, size of axis text, passed on to plot.
#' @param plot.self logical, should intra-genomic plots be made?
#' @param alpha numeric (0-1) specifying transparency of points.
#' @param chr.abbrev.fun function, to strip characters off of chromosome IDs
#' @param ... additional arguments passed to plot
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom utils combn
#' @importFrom grDevices rainbow
#' @export
plot_map <- function(palette = pal_genespace,
                     chr.list,
                     map,
                     gff,
                     genomeIDs,
                     gap.prop = 0.01,
                     use.rank = T,
                     col.byBlk = T,
                     label.blocks = F,
                     chr.abbrev.fun = function(x) x,
                     pch = ".",
                     alpha = .5,
                     text.cex = .5,
                     plot.self = F,
                     cex.axis = .8,
                     ...){

  if (!"block.id" %in% colnames(map))
    map[,block.id := 1]

  if(length(genomeIDs) != length(chr.list))
    stop("GenomeIDs must be of same length as names of chr.list\n")

  gff <- format_gffChrlist(gff = gff,
                           chr.list = chr.list,
                           genomes = genomeIDs,
                           use.rank = use.rank,
                           gap.prop = gap.prop)

  chr.coords <- gff[,list(start = min(start.p),
                        end = max(end.p)),
                  by = list(genome, chr)]
  chr.coords[,mean := (start+end) / 2]

  if(col.byBlk){
    m <- map[,c("id1","id2","genome1","genome2","block.id")]
    m <- rbind(m, data.table(id1 = m$id2,
                             id2 = m$id1,
                             genome1 = m$genome2,
                             genome2 = m$genome1,
                             block.id = m$block.id))
  }else{
    m <- map[,c("id1","id2","genome1","genome2")]
    m <- rbind(m, data.table(id1 = m$id2,
                             id2 = m$id1,
                             genome1 = m$genome2,
                             genome2 = m$genome1))
  }
  m <- m[!duplicated(m),]

  g1 <- data.table(gff)
  g2 <- data.table(gff)
  setnames(g1, paste0(colnames(gff), "1"))
  setnames(g2, paste0(colnames(gff), "2"))
  m <- merge(merge(m, g1, by = c("genome1","id1")), g2, by= c("genome2","id2"))


  if(length(genomeIDs) == 1){
    cmb <- rbind(genomeIDs, genomeIDs)
  }else{
    cmb <- combn(genomeIDs,2)
    if(plot.self){
      cmb <- cbind(cmb, sapply(genomeIDs, function(x) c(x,x)))
    }
  }

  for(i in 1:ncol(cmb)){
    g1 = cmb[1,i]
    g2 = cmb[2,i]
    tchr1 <- subset(chr.coords, genome == g1)
    tchr2 <- subset(chr.coords, genome == g2)
    tmap <- subset(m, genome1 == g1 & genome2 == g2)
    lab <- ifelse(use.rank, "gene-rank", "physical position")
    if(!is.function(palette)){
      tmap$col <- palette[1]
    }else{
      if(!col.byBlk){
        pal <- palette(length(unique(tmap$chr1)))
        tmap$col <- add_alpha(pal[match(tmap$chr1,chrl[[g1]])],alpha)
      }else{
        pal <- sample(palette(length(unique(tmap$block.id))))
        tmap$col <- add_alpha(pal[match(tmap$block.id,unique(tmap$block.id))],alpha)
      }
    }


    plot(c(0,1),c(0,1), type = "n", bty = "n", axes = F,
         xlab = paste(g1,"chromosome"),
         ylab = paste(g2,"chromosome"),
         main = paste("mapping by", lab), ...)
    for(j in 1:nrow(tchr1)){
      for(k in 1:nrow(tchr2)){
        rect(xleft = tchr1$start[j],
             xright = tchr1$end[j],
             ybottom = tchr2$start[k],
             ytop = tchr2$end[k],
             col = "grey95", border = NA)
        pts <- subset(tmap, chr1 == tchr1$chr[j] &
                        chr2 == tchr2$chr[k])
        if(nrow(pts) > 0){
          with(pts, points(start.p1, start.p2, pch =pch, col = col))
          if(label.blocks){
            bt <- pts[, list(min1 = min(start.p1),
                             min2 = min(start.p2),
                             max1 = max(end.p1),
                             max2 = max(end.p2)),
                      by = list(block.id)]
            bt[,pos1 := (min1 + max1) / 2]
            bt[,mean2 := (min2 + max2) / 2]
            bt[,pos2 := (max2+mean2)/2]
            with(bt, text(x = pos1, y = pos2, labels = block.id,
                          cex = text.cex))
          }
        }
      }
    }
    axis(1, at = tchr1$mean, labels = chr.abbrev.fun(tchr1$chr), cex.axis = cex.axis)
    axis(2, at = tchr2$mean, labels = chr.abbrev.fun(tchr2$chr), cex.axis = cex.axis)
  }
  return(list(map = m))
}
