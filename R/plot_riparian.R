#' @title plot_riparian
#'
#' @description
#' \code{plot_riparian} plot_riparian
#'
#' @param map data.table, format similar to blast
#' @param gff data.table, format similar to gff
#' @param genomes character, genomes to plot
#' @param chr.list list, named by genome labels (label will
#' replace the genomes element). Must be same length as genomes.
#' Character vector must contain chromosome IDs within the genomes.
#' @param use.rank logical, should ranks or bp coordinates be plotted?
#' @param chr.abbrev.fun function, to strip characters off of chromosome IDs
#' @param vertical.plot.buffer numeric, buffer above / below riparian plot
#' @param horiz.plot.buffer numeric, buffer left / right of riparian plot
#' @param chr.bg.col character, to be coerced to a color of the point behind
#' chromosome IDs in riparian plot
#' @param chr.bg.pch numeric, shape of point behind chromosomes in riparian plot
#' @param chr.bg.cex numeric, character expansion factor of points behind
#' chromosomes in riparian plot
#' @param chr.id.col character, to be coerced to a color of the
#' chromosome IDs in riparian plot
#' @param chr.id.cex  numeric, character expansion factor of
#' chromosome IDs in riparian plot
#' @param gap.prop numeric, proportion of total genome size represented
#' by gaps in riparian plot
#' @param chr.segm.lwd numeric, chromosome line segment thickness
#' in riparian plot
#' @param chr.buff numeric, how far should all line segments stop away
#' from chromosome line segments in riparian plot
#' @param blk.col character, to be coerced to a color of the polygons that
#' represent syntenic blocks in riparian plot
#' @param blk.border logical (or character coerced to color), of the block
#' borders in riparian plot
#' @param chr.lab.buff chr.lab.buff
#' @param ... not currently in use
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
plot_riparian <- function(map,
                          chr.list,
                          genomeIDs,
                          gff,
                          points.per.curve = 1000,
                          vertical.plot.buffer = 1,
                          horiz.plot.buffer = .01,
                          lab.chr = T,
                          lab.chr.1only = F,
                          use.rank = T,
                          chr.bg.col = rgb(1,1,1,.8),
                          chr.bg.pch = 16,
                          chr.bg.cex = 2.5,
                          chr.id.col = "black",
                          chr.id.cex = .8,
                          chr.lab.buff = chr.buff/2,
                          chr.abbrev.fun = function(x)
                            gsub("^0", "",
                                 gsub("_", "",
                                      gsub("scaffold", "",
                                           gsub("chr", "", tolower(x))))),
                          gap.prop = .02,
                          chr.segm.lwd = 4,
                          chr.segm.col = "black",
                          chr.buff = 0.05,
                          blk.col = rgb(0,0,1,.1),
                          blk.border = FALSE,
                          simplify.poly = .25){
  ####################################################################
  # 1. Set up plotting window
  omar <- par()$mar
  par(mar = c(1,1,1,1))
  on.exit(par(mar = omar))
  if (length(genomeIDs) != length(chr.list)) {
    stop("genomes must be the same length as chr.list ... \n")
  }
  ####################################################################
  # 2. Get genomic coordinates from the gff
  gff <- format_gffChrlist(
    gff = gff,
    chr.list = chr.list,
    genomes = genomeIDs,
    use.rank = use.rank,
    gap.prop = gap.prop)
  ####################################################################
  # 3. Ensure proper ordering of genomes
  combs <- data.table(t(combn(genomeIDs,2)))
  setnames(combs, c("genome1","genome2"))
  map <- merge(map, combs, by = c("genome1","genome2"))
  ####################################################################
  # 4. Project genomic coords onto the map
  map <- format_mapChrlist(
    map = map,
    gff = gff,
    chr.list = chr.list,
    genomes = genomeIDs)
  setkey(map, pos1)
  ####################################################################
  # 5. Make the blk data.table
  fais <- make_fais(genomes = genomeIDs, gff = gff)
  blk <- map[,list(start1 = pos1[1],
                   start2 = pos2[1],
                   end1 = pos1[length(pos1)],
                   end2 = pos2[length(pos2)]),
             by = list(block.id,
                       y1, y2,
                       genome1, genome2,
                       chr1, chr2)]
  ####################################################################
  # 6. make the plot window and y axis
  plot(range(fais$y),
       range(fais$y),
       type = "n",
       xlim = c(0 - horiz.plot.buffer*5,
                1 + horiz.plot.buffer),
       ylim = c(0 - vertical.plot.buffer,
                max(fais$y) + vertical.plot.buffer),
       axes = F,
       bty = "n",
       ylab = "genomes",
       xlab = "relative physical position")
  ####################################################################
  # 7. add genome IDs to the y axis
  genomes <- names(chr.list)
  text(labels = genomes,
       x = rep(-.01, length(genomes)),
       y = 0:(length(genomes) - 1),
       adj = c(1, .5))
  ####################################################################
  # 8. for each genome comparison, plot blocks
  draw_blkPolygon(
    blk = blk,
    chr.buffer = chr.buff,
    blk.border = blk.border,
    blk.col = blk.col,
    simplify.poly = simplify.poly,
    points.per.curve = points.per.curve)

  ####################################################################
  # 9. plot and label the chromosomes
  label_riparian(
    fais = fais,
    chr.segm.col = chr.segm.col,
    chr.segm.lwd = chr.segm.lwd,
    genomes = genomes,
    chr.lab.buff = chr.lab.buff,
    chr.bg.col = chr.bg.col,
    chr.bg.cex = chr.bg.cex,
    chr.bg.pch = chr.bg.pch,
    chr.id.col = chr.id.col,
    chr.id.cex = chr.id.cex,
    lab.chr = lab.chr,
    lab.chr.1only = lab.chr.1only,
    chr.abbrev.fun = chr.abbrev.fun)
  ####################################################################
  # 10. Return data to add genes, highlight blocks, etc.
  return(list(
    gff = gff,
    map = map,
    blk = blk,
    fais = fais,
    rip.param = list(
      points.per.curve = points.per.curve,
      vertical.plot.buffer = vertical.plot.buffer,
      horiz.plot.buffer = horiz.plot.buffer,
      lab.chr = lab.chr,
      lab.chr.1only = lab.chr.1only,
      use.rank = use.rank,
      chr.bg.col = chr.bg.col,
      chr.bg.pch = chr.bg.pch,
      chr.bg.cex = chr.bg.cex,
      chr.id.col = chr.id.col,
      chr.id.cex = chr.id.cex,
      chr.lab.buff = chr.lab.buff,
      chr.abbrev.fun = chr.abbrev.fun,
      gap.prop = gap.prop,
      chr.segm.lwd = chr.segm.lwd,
      chr.segm.col = chr.segm.col,
      chr.buff = chr.buff,
      blk.col = blk.col,
      blk.border = blk.border,
      simplify.poly = simplify.poly)))
}
