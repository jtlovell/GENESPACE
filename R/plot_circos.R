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
#' @param link.col.list list, with colors that matches chr.list
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
                        genomeIDs,
                        gff,
                        chr.list,
                        track.margin = c(0.005, 0.005),
                        link.col.list = NULL,
                        border.col.list = NULL,
                        chr.abbrev.fun = function(x) x,
                        border.track.height = 0.025,
                        link.alpha = 0.6,
                        use.rank = TRUE,
                        genome.cols = greys(length(genomes)),
                        gap.deg = rep(.5, length(unlist(chr.list)))){

  on.exit(circos.clear())

  ####################################################################
  # 1. Get genomic coordinates from gff file
  gff <- format_gffChrlist(
    gff = gff,
    chr.list = chr.list,
    genomes = genomeIDs,
    use.rank = use.rank,
    gap.prop = 0,
    do.cumulative = F)

  ####################################################################
  # 2. Format and (optionally) rerank map
  m <- data.table(map)
  if (use.rank) {
    m[,start1 := as.integer(frank(start1, ties.method = "dense")),
      by = list(genome1, genome2, chr1)]
    m[,start2 := as.integer(frank(start2, ties.method = "dense")),
      by = list(genome1, genome2, chr2)]
    m[,end1 := as.integer(frank(end1, ties.method = "dense")),
      by = list(genome1, genome2, chr1)]
    m[,end2 := as.integer(frank(end2, ties.method = "dense")),
      by = list(genome1, genome2, chr2)]
  }

  map <- format_mapChrlist(
    map = m,
    gff = gff,
    chr.list = chr.list,
    genomes = genomeIDs,
    forCircos = T)

  ####################################################################
  # 3. make the block object
  blk <- map[,list(start1 = min(pos1),
                   start2 = min(pos2),
                   end1 = max(pos1),
                   end2 = max(pos2)),
             by = list(block.id,
                       genome1, genome2,
                       chr1, chr2)]

  ####################################################################
  # 4. Make a fasta index-like object from the gff
  fais <- make_fais(
    genomes = genomeIDs,
    gff = gff)

  fais[, chr := chr.abbrev.fun(chr)]
  fais[, unique := paste(genome, chr)]

  init <- with(fais, data.frame(
    name = unique,
    start = start,
    end = end,
    stringsAsFactors = F))

  ####################################################################
  # 5. Make a data.table with order of chromosomes
  chrl.dt <- rbindlist(lapply(names(chr.list), function(i)
    data.table(name = paste(i, chr.abbrev.fun(chr.list[[i]])),
               col = add_alpha(link.col.list[[i]], alpha = link.alpha))))
  init <- init[match(chrl.dt$name, init$name),]

  if (!is.null(border.col.list)) {
    chrl.bord.dt <- rbindlist(lapply(names(border.col.list), function(i)
      data.table(name = paste(i, chr.abbrev.fun(chr.list[[i]])),
                 col = border.col.list[[i]])))
  }

  ####################################################################
  # 6. Make bed files for each side of the links
  bed <- data.table(blk)
  bed[,name1 := paste(genome1, chr.abbrev.fun(chr1))]
  bed[,name2 := paste(genome2, chr.abbrev.fun(chr2))]
  bed[,ord1 := as.numeric(factor(name1, levels = chrl.dt$name))]
  bed[,ord2 := as.numeric(factor(name2, levels = chrl.dt$name))]
  bed <- subset(bed, ord1 < ord2 | (ord1 == ord2 & start1 < start2))

  bed1 <- with(bed, data.frame(
    chr = name1,
    start = start1,
    end = end1,
    value = 1,
    stringsAsFactors = F))

  bed2 <- with(bed, data.frame(
    chr = name2,
    start = start2,
    end = end2,
    value = 1,
    stringsAsFactors = F))

  ####################################################################
  # 6. Initiate the circos plot
  circos.par(start.degree = 90,
             clock.wise = TRUE,
             track.margin = track.margin,
             "gap.degree" = gap.deg)
  circos.genomicInitialize(init)

  ####################################################################
  # 7. (optionally) plot the border
  if (!is.null(border.col.list)) {
    circos.track(ylim = c(0, .15),
                 bg.col = chrl.bord.dt$col,
                 bg.border = NA,
                 cell.padding = circos.par("cell.padding") - 0.01,
                 track.height = border.track.height)
  }


  ####################################################################
  # 8. Plot the links
  link.col <- chrl.dt$col[match(bed1$chr,chrl.dt$name)]
  circos.genomicLink(bed1[,c(1:3)],
                     bed2[,c(1:3)],
                     col = link.col,
                     border = NA)
}
