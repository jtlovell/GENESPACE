#' @title plot_riparianGenes
#'
#' @description
#' \code{plot_riparianGenes} plot_riparianGenes
#'
#' @param map data.table, format similar to blast
#' @param rip.out list, output from plot_riparian
#' @param geneID.abbrev.fun function, to strip characters off of gene IDs
#' in riparian plot
#' @param geneid.cex numeric, gene ID character expansion factor in riparian plot
#' @param geneid.offset numeric, amount of space between geneID lines and
#' text in riparian plot
#' @param dodge.geneIDs numeric, amount of distance to dodge geneIDs
#'  in riparian plot
#' @param gene.col character, to be coerced to a color, specifying the
#' color of lines connecting syntenic orthologs in riparian plot
#' @param scale2dodge numeric, scaling factor for geneID dodging,
#' in riparian plot. Larger values do more dodging.
#' @param genes2plot character, genes / orthology networks to plot
#' in riparian plot
#' @param label.col character, to be coerced to a color for gene IDs.
#' @param ... not currently in use
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
plot_riparianGenes <- function(map,
                               rip.out,
                               genes2plot,
                               gene.col = "darkorange3",
                               geneID.abbrev.fun = function(x) gsub("_", "",x),
                               label.col = "black",
                               geneid.cex = .5,
                               geneid.offset = 0.2,
                               dodge.geneIDs = .01,
                               scale2dodge = 2,
                               annotate.genes = FALSE){

  ####################################################################
  # 1. get all unique genes in a block
  m <- data.table(map)
  b <- subset(m, id1 %in% genes2plot | id2 %in% genes2plot)

  genes.in.blk <- unique(unlist(subset(
    m,
    id1 %in% unique(b$id1, b$id2) |
      id2 %in% unique(b$id1, b$id2))[,c("id1","id2")]))
  mb <- subset(m, id1 %in% genes.in.blk & id2 %in% genes.in.blk)
  genes2plot.all <- unique(unlist(mb[,c("id1","id2")]))

  ####################################################################
  # 2. draw the polygon
  with(rip.out$rip.param, draw_genePolygon(
    map = rip.out$map,
    chr.buff = chr.buff,
    genes2plot = genes2plot.all,
    gene.col = gene.col,
    simplify.poly = simplify.poly,
    points.per.curve = points.per.curve))


  if (annotate.genes)
    annotate_riparian(
      map = rip.out$map,
      chr.buff = rip.out$rip.param$chr.buff,
      geneID.abbrev.fun = geneID.abbrev.fun,
      geneid.offset = geneid.offset,
      genomes = genomeIDs,
      dodge.geneIDs = dodge.geneIDs,
      scale2dodge = scale2dodge,
      geneid.cex = geneid.cex,
      gene.colors = label.col,
      genes2plot = genes2plot.all)

  return(list(genes2plot.all = genes2plot.all, map = mb))
}

