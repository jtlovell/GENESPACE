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
#' @param geneID.abbrev.fun function, to strip characters off of gene IDs
#' in riparian plot
#' @param geneid.cex numeric, gene ID character expansion factor in riparian plot
#' @param geneid.offset numeric, amount of space between geneID lines and
#' text in riparian plot
#' @param dodge.geneIDs numeric, amount of distance to dodge geneIDs
#'  in riparian plot
#' @param ortho.col character, to be coerced to a color, specifying the
#' color of lines connecting syntenic orthologs in riparian plot
#' @param scale2dodge numeric, scaling factor for geneID dodging,
#' in riparian plot. Larger values do more dodging.
#' @param genes2plot character, genes / orthology networks to plot
#' in riparian plot
#' @param gene.labels character, labels for geneIDs to plot
#' in riparian plot. Must be of same length as genes2plot.
#' @param gene.colors character, to be coerced to a color.
#' in riparian plot. If specified, must be same length as genes2plot
#' @param lab.chr logical, should chromosomes be labled?
#' @param lab.chr.1only logical, should only the chromosomes of the
#' first genome be labeled?
#' @param chr.segm.col color, what color should the chromosome
#' segments be?
#' @param simplify.poly numeric, how much polygon simplification
#' should happen.
#' @param points.per.curve numeric, specify how complex the curves
#' should be. Higher values yield larger file sizes, but smoother
#' curves
#' @param blks2highlight character, indicating which blocks
#' should be highlighted
#' @param highlight.cols color vector, colors for highlighted blocks
#' @param annotate.genes logical, should genes be annotated?
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
                          genomes,
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
                          simplify.poly = .25,
                          blks2highlight = NULL,
                          highlight.cols = NULL,

                          geneID.abbrev.fun = function(x) gsub("_", "",x),

                          geneid.cex = .5,
                          geneid.offset = 0.2,
                          dodge.geneIDs = .01,
                          ortho.col = "darkorange3",
                          scale2dodge = 2,

                          annotate.genes = NULL,
                          genes2plot = NULL,
                          gene.labels = NULL,
                          gene.colors = NULL){

  omar <- par()$mar
  par(mar = c(1,1,1,1))
  on.exit(par(mar = omar))

  if (length(genomes) != length(chr.list)) {
    stop("genomes must be the same length as chr.list ... \n")
  }

  if(!is.null(blks2highlight)){
    if(is.null(highlight.cols))
      highlight.cols <- pal_genespace(length(blks2highlight))
    if(length(blks2highlight) > 1 & length(highlight.cols) == 1)
      highlight.cols <- rep(highlight.cols, length(blks2highlight))
    trkd <- track_blocks(map = map, blk.ids = blks2highlight,
                         cols = highlight.cols, bg.col = blk.col)
    blk.col <- trkd$blk.cols
    map <- trkd$map
  }


  gff <- format_gffChrlist(gff = gff,
                           chr.list = chr.list,
                           genomes = genomes,
                           use.rank = use.rank,
                           gap.prop = gap.prop)

  combs <- data.table(t(combn(genomes,2)))
  setnames(combs, c("genome1","genome2"))
  map <- merge(map, combs, by = c("genome1","genome2"))

  map <- format_mapChrlist(map = map,
                           gff = gff,
                           chr.list = chr.list,
                           genomes = genomes)
  genomes <- as.character(names(chr.list))

  if (!is.null(genes2plot)) {
    if (is.null(gene.labels)) {
      names(genes2plot) <- genes2plot
    }else{
      if (length(gene.labels) != length(genes2plot)) {
        warning("Length of gene labels != genes2plot, ignoring\n")
        names(genes2plot) <- genes2plot
      }else{
        names(genes2plot) <- gene.labels
      }
    }

    if (is.null(gene.colors)) {
      gene.colors <- rep("black", length(genes2plot))
      names(gene.colors) <- genes2plot
    }else{
      if (length(gene.colors) != length(genes2plot)) {
        warning("Length of gene labels != genes2plot, ignoring\n")
        gene.colors <- rep("black", length(genes2plot))
        names(gene.colors) <- genes2plot
      }else{
        names(gene.colors) <- genes2plot
      }
    }

    if (nrow(gff) == 0) {
      warning("no genes found in specified genomes/chromosomes\n")
      genes2plot <- NULL
    }
  }

  # -- make chromosome coordinates
  genomes <- as.character(names(chr.list))
  fais <- make_fais(genomes = genomes, gff = gff)
  setkey(map, pos1)
  # -- make the block object
  blk <- map[,list(start1 = pos1[1],
                   start2 = pos2[1],
                   end1 = pos1[length(pos1)],
                   end2 = pos2[length(pos2)]),
             by = list(block.id,
                       y1, y2,
                       genome1, genome2,
                       chr1, chr2)]

  blk$blk.col <- blk.col[1]
  if (length(blk.col) > 1 &
      any(names(blk.col) %in% as.character(blk$block.id))) {
    blk.col <- blk.col[names(blk.col) %in% unique(blk$block.id)]
    for (i in names(blk.col)) {
      blk$blk.col[blk$block.id == i] <- blk.col[i]
    }
  }

  # -- make the plot window and y axis
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
  text(labels = genomes,
       x = rep(-.01, length(genomes)),
       y = 0:(length(genomes) - 1),
       adj = c(1, .5))


  # -- for each genome comparison, plot blocks
  draw_blkPolygon(blk = blk,
                  chr.buffer = chr.buff,
                  blk.border = blk.border,
                  simplify.poly = simplify.poly,
                  points.per.curve = points.per.curve)


  # -- add gene paths
  if (!is.null(genes2plot)) {
    ogs2plot <- map$og.id[map$id1 %in% genes2plot |
                            map$id2 %in% genes2plot]
    draw_genePolygon(map = map,
                     genes2plot = genes2plot,
                     chr.buff = chr.buff,
                     simplify.poly = simplify.poly,
                     ortho.col = ortho.col,
                     points.per.curve = points.per.curve)
    if (!is.null(annotate.genes)) {
      genes2plot <- unique(unlist(subset(map, og.id %in% ogs2plot)[,c("id1", "id2")]))
      annotate_riparian(map = map,
                        chr.buff = chr.buff,
                        geneID.abbrev.fun = geneID.abbrev.fun,
                        geneid.offset = geneid.offset,
                        genomes = genomes,
                        dodge.geneIDs = dodge.geneIDs,
                        scale2dodge = scale2dodge,
                        geneid.cex = geneid.cex,
                        gene.colors = gene.colors,
                        genes2plot = genes2plot)
    }
  }

  label_riparian(fais = fais,
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
  return(list(blk = blk, genomes = genomes, fais = fais))
}
