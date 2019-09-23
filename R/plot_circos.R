#' @title build plot_circos
#'
#' @description
#' \code{plot_circos} plot_circos
#'
#' @param map data.table, format similar to blast
#' @param gff data.table, format similar to gff
#' @param genomes character, genomes to plot
#' @param chr.list list, named by genome labels (label will
#' replace the genomes element). Must be same length as genomes.
#' Character vector must contain chromosome IDs within the genomes.
#' @param col.list list, with colors that matches chr.list
#' @param genome.cols color vector, to color the genomes.
#' @param palette function, palette function for colors
#' @param use.rank logical, should ranks or bp coordinates be plotted?
#' @param chr.abbrev.fun function, to strip characters off of chromosome IDs
#' @param chr.order numeric, order of chromosomes, passed only to circos
#' @param border.track.height numeric, size of exterior circos track
#' @param link.alpha numeric, transparency of circos links
#' @param gap.ingenome numeric, size of gaps within genomes, overrides gap.deg
#' @param gap.outgenome numeric,  size of gaps among genomes, overrides gap.deg
#' @param gap.deg numeric, degree of gaps, either of length 1, or the number
#' of unique chromosomes
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
plot_circos <- function(map,
                        genomes,
                        gff,
                        chr.list,
                        col.list = NULL,
                        chr.abbrev.fun = function(x) x,
                        chr.order = 1:length(unlist(chr.list)),
                        border.track.height = 0.05,
                        link.alpha = 0.6,
                        gap.ingenome = NULL,
                        gap.outgenome = NULL,
                        use.rank = T,
                        palette = rainbow,
                        genome.cols = greys(length(genomes)),
                        gap.deg = rep(.5, length(unlist(chr.list)))){

  on.exit(circos.clear())

  if (is.null(gap.deg)) {
    gap.deg <- unlist(lapply(chr.list, function(x)
      c(rep(gap.ingenome, (length(x) - 1)),
        gap.outgenome)))
  }

  if (is.null(col.list)) {
    col.list <- chr.list
    for (i in 1:length(col.list))
      col.list[[i]] <- palette(length(col.list[[i]]) + 2)[-c(1, length(col.list[[i]]) + 2)]
  }

  gff <- format_gffChrlist(gff = gff,
                           chr.list = chr.list,
                           genomes = genomes,
                           use.rank = use.rank,
                           gap.prop = gap.prop,
                           do.cumulative = F)
  if (use.rank) {
    m <- map
    m[,start1 := as.integer(frank(start1, ties.method = "dense")),
      by = list(genome1, genome2, chr1)]
    m[,start2 := as.integer(frank(start2, ties.method = "dense")),
      by = list(genome1, genome2, chr2)]
    m[,end1 := as.integer(frank(end1, ties.method = "dense")),
      by = list(genome1, genome2, chr1)]
    m[,end2 := as.integer(frank(end2, ties.method = "dense")),
      by = list(genome1, genome2, chr2)]
  }else{
    m <- map
  }
  map <- format_mapChrlist(map = m,
                           gff = gff,
                           chr.list = chr.list,
                           genomes = genomes,
                           forCircos = T)
  genomes <- names(chr.list)
  # -- make the block object
  blk <- map[,list(start1 = min(pos1),
                   start2 = min(pos2),
                   end1 = max(pos1),
                   end2 = max(pos2)),
             by = list(block.id,
                       genome1, genome2,
                       chr1, chr2)]

  fais <- make_fais(genomes = genomes, gff = gff)
  fais[, chr := chr.abbrev.fun(chr)]
  fais$unique <- with(fais, paste(genome, chr))

  init <- with(fais,
               data.frame(name = unique,
                          start = start,
                          end = end,
                          stringsAsFactors = F))
  if (!is.null(chr.order)) {
    if (length(chr.order) != nrow(init)) {
      warning("chr.order not of same length as chr.list, ignoring\n")
    }else{
      init <- init[chr.order, ]
    }
  }

  circos.par(start.degree = 90,
             clock.wise = TRUE,
             "gap.degree" = gap.deg)
  circos.genomicInitialize(init)
  genome.cols <- NULL
  if (is.null(genome.cols) |
      length(genome.cols) != length(genomes))
    genome.cols <- grey.colors(n = length(genomes))
  fais$color <- NA
  for (i in 1:length(names(chr.list))) {
    fais$color[fais$genome == names(chr.list)[i]] <- genome.cols[i]
  }

  # circos.track(ylim = c(0, 1),
  #              bg.col = "white",
  #              bg.border = NA,
  #              track.height = border.track.height)


  bed <- blk
  col.dt <- data.table(genome = rep(names(col.list),
                                    sapply(col.list, length)),
                       chr = unlist(chr.list),
                       col = sapply(unlist(col.list), add_alpha, alpha = link.alpha))
  col.dt$name <- with(col.dt, paste(genome, chr.abbrev.fun(chr)))
  bed1 <- with(bed,
               data.frame(chr = paste(genome1, chr.abbrev.fun(chr1)),
                          start = start1,
                          end = end1,
                          value = 1,
                          stringsAsFactors = F))
  link.col <- col.dt$col[match(bed1$chr,col.dt$name)]

  bed2 <- with(bed,
               data.frame(chr = paste(genome2, chr.abbrev.fun(chr2)),
                          start = start2,
                          end = end2,
                          value = 1,
                          stringsAsFactors = F))

  circos.genomicLink(bed1[,c(1:3)],
                     bed2[,c(1:3)],
                     col = link.col,
                     border = NA)
}
