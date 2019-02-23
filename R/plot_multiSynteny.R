#' @title build orthology networks
#'
#' @description
#' \code{build_orthoNetworks} Search orthofinder output for genes withou
#' hits in one of the blocks in the network.
#'
#' @param map map results data.table
#' @param chr.order.list produced from match_syntenicChrs
#' @param genes2highlight vector of gene IDs to highlight. Will propagate
#' across entire orthonetwork
#' @param blk.col color of block polygons
#' @param highlight.col color of highlights
#' @param highlight.lwd line thickness of highlights
#' @param chr.col color of chromosome line segments
#' @param chr.lwd line thickness of chromosome line segments
#' @param clean.radius radius of clean blocks to run before construction
#' @param clean.nmappings mappings of clean blocks to run before construction
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
plot_multiSynteny <- function(map,
                              chr.order.list,
                              genes2highlight = NULL,
                              blk.col = rgb(0,0,0,.2),
                              highlight.col = "magenta",
                              highlight.lwd = 1,
                              chr.col = "black",
                              chr.lwd = 2,
                              clean.radius = 5,
                              clean.nmappings = 5,
                              ...){
  ########################################################
  ########################################################
  make_chrDB <- function(map,
                         chr.order.list){

    # -- get chr ids from chr.order.list
    chr.ids = rbindlist(lapply(names(chr.order.list), function(i)
      data.table(genome = i,
                 genome.order = which(names(chr.order.list) == i),
                 chr = chr.order.list[[i]],
                 chr.order = 1:length(chr.order.list[[i]]))))
    setkey(chr.ids, genome, chr)

    # -- convert to long format
    map.long = with(map.in,
                    rbind(
                      data.table(genome = genome1,
                                 id = id1,
                                 chr = chr1,
                                 rank = rank1,
                                 block.id = block.id),
                      data.table(genome = genome2,
                                 id = id2,
                                 chr = chr2,
                                 rank = rank2,
                                 block.id = block.id)))

    # -- Summarize by chromsome
    chr.lens <- map.long[,list(chr.len = max(rank)),
                         by = list(genome, chr)]
    chr.lens[,genome.len := sum(chr.len),
             by = list(genome)]
    setkey(chr.lens, genome, chr)

    out.chr <- merge(chr.ids, chr.lens, by = c("genome","chr"))
    setkey(out.chr, genome.order, chr.order)

    out.chr[, gap.size := c(0,cumsum(rep((genome.len[1] + (genome.len[1]*.01*.N) * .01),.N-1))),
            by = list(genome)]

    chr.md <- rbindlist(lapply(split(out.chr,"genome"), function(x){
      tmp = x[!duplicated(x$chr),c("genome","chr","genome.order","chr.order","chr.len","genome.len")]
      tmp = tmp[order(tmp$chr.order),]
      gap.size = (tmp$genome.len[1]+(tmp$genome.len[1]*.01*nrow(tmp))) * .01
      gap.add = c(0,cumsum(rep(gap.size,nrow(tmp)-1)))
      tmp$chr.start = c(0, cumsum(tmp$chr.len)[-nrow(tmp)])+gap.add
      tmp$chr.end = gap.add + cumsum(tmp$chr.len)
      tmp$max.genome = max(tmp$chr.end)
      tmp$chr.len<-NULL
      tmp$genome.len<-NULL
      return(tmp)
    }))

    chr.md$prop.chr.start = with(chr.md,chr.start/max.genome)
    chr.md$prop.chr.end = with(chr.md, chr.end/max.genome)


    return(chr.md)
  }
  ########################################################
  ########################################################
  make_blockDB <- function(map, chr.db){
    # clean blocks for better presentation
    blk <- clean_blocks(map = map.in,
                        radius = clean.radius,
                        n.mappings = clean.nmappings,
                        rerank = F,
                        verbose = F)$block

    map.tmp <- map.in[,c("genome1","genome2","chr1","chr2")]
    map.tmp <- map.tmp[!duplicated(map.tmp),]
    blk.m <- merge(map.tmp,
                   blk, by = colnames(map.tmp))
    blk.m[,propstart1:=(rankstart1/max(rankend1)),
          by = list(genome1, chr1)]
    blk.m[,propstart2:=(rankstart2/max(rankend2)),
          by = list(genome2, chr2)]
    blk.m[,propend1:=(rankend1/max(rankend1)),
          by = list(genome1, chr1)]
    blk.m[,propend2:=(rankend2/max(rankend2)),
          by = list(genome2, chr2)]

    chr.db1 <- data.table(chr.db)
    setnames(chr.db1, paste0(colnames(chr.db),"1"))
    chr.db2 <- data.table(chr.db)
    setnames(chr.db2, paste0(colnames(chr.db),"2"))

    bm <- merge(chr.db2, blk.m, by = c("genome2","chr2"))
    bm <- merge(chr.db1, bm, by = c("genome1","chr1"))
    bm$proplen.chr1 <- with(bm, prop.chr.end1 - prop.chr.start1)
    bm$proplen.chr2 <- with(bm, prop.chr.end2 - prop.chr.start2)

    bm$tp.start1 <- with(bm, (propstart1*proplen.chr1)+prop.chr.start1)
    bm$tp.start2 <- with(bm, (propstart2*proplen.chr2)+prop.chr.start2)
    bm$tp.end1 <- with(bm, (propend1*proplen.chr1)+prop.chr.start1)
    bm$tp.end2 <- with(bm, (propend2*proplen.chr2)+prop.chr.start2)
    return(bm[,c("genome1","genome2","genome.order1","genome.order2",
                 "chr1","chr2","chr.order1","chr.order2",
                 "block.id","proplen.chr1","proplen.chr2",
                 "prop.chr.start1","prop.chr.start2",
                 "tp.start1","tp.start2","tp.end1","tp.end2")])
  }
  ########################################################
  ########################################################
  make_mapDB <- function(map,
                         chr.db,
                         genes2highlight){
    map <- data.table(map)
    all.genes <- genes2highlight
    for(i in 1:length(genomeIDs)){
      all.genes <- unique(unlist(map[map$id1 %in% all.genes |
                                       map$id2 %in% all.genes,
                                     c("id1","id2")]))
    }
    # clean blocks for better presentation
    map.m <- map
    map.m[,prop1:=(rank1/max(rank1)),
          by = list(genome1, chr1)]
    map.m[,prop2:=(rank2/max(rank2)),
          by = list(genome2, chr2)]

    chr.db1 <- data.table(chr.db)
    setnames(chr.db1, paste0(colnames(chr.db),"1"))
    chr.db2 <- data.table(chr.db)
    setnames(chr.db2, paste0(colnames(chr.db),"2"))

    bm <- merge(chr.db2, map.m, by = c("genome2","chr2"))
    bm <- merge(chr.db1, bm, by = c("genome1","chr1"))
    bm$proplen.chr1 <- with(bm, prop.chr.end1 - prop.chr.start1)
    bm$proplen.chr2 <- with(bm, prop.chr.end2 - prop.chr.start2)

    bm$tp1 <- with(bm, (prop1*proplen.chr1)+prop.chr.start1)
    bm$tp2 <- with(bm, (prop2*proplen.chr2)+prop.chr.start2)

    return(bm[bm$id1 %in% all.genes | bm$id2 %in% all.genes,
              c("genome1","genome2","genome.order1","genome.order2",
                "chr1","chr2","chr.order1","chr.order2",
                "id1","id2",
                "block.id","proplen.chr1","proplen.chr2",
                "prop.chr.start1","prop.chr.start2",
                "tp1","tp2")])
  }
  ########################################################
  ########################################################

  ########################################################
  ########################################################

  ########################################################
  # -- re-rank positions
  map[,rank1 := frank(start1,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr1)]
  map[,rank2 := frank(start2,
                      ties.method = "dense"),
      by = list(genome1, genome2, chr2)]
  ########################################################

  ########################################################
  # -- subset map to comparisons to plot
  map.in = rbindlist(lapply(1:(length(chr.order.list)-1),function(i){
    g1 = names(chr.order.list)[i]
    g2 = names(chr.order.list)[i+1]
    return(map[with(map, genome1 == g1 & genome2 == g2),])
  }))
  ########################################################

  ########################################################
  # -- Make chromsome metadata and plot
  chr.db = make_chrDB(map = map.in,
                      chr.order.list = chr.order.list)

  plot(1:length(genomeIDs),
       seq(from = 0, to = 1, length.out = length(genomeIDs)),
       type = "n",
       axes = F, bty = "n",
       xlab = "",
       ylab = "Relative gene-order position",
       main = "Syntenic blocks by gene-order",
       ... )
  axis(1,
       at = 1:length(genomeIDs),
       labels = genomeIDs,
       las = 2)

  ########################################################

  ########################################################
  # -- Make block metadata and plot
  blk.db = make_blockDB(map = map.in,
                        chr.db = chr.db)
  blk.db = blk.db[with(blk.db, genome.order1 < genome.order2),]
  spl.mat = split(blk.db, "block.id")
  for(x in spl.mat){
    xs = with(x, c(genome.order1, genome.order2, genome.order2, genome.order1))
    ys = with(x, c(tp.start1,tp.start2, tp.end2, tp.end1))
    polygon(x = xs,
            y = ys,
            col = blk.col,
            border = NA)
  }

  ########################################################

  ########################################################
  # -- Make map metadata and plot
  if(!is.null(genes2highlight)){
    map.db = make_mapDB(map = map.in,
                        chr.db = chr.db,
                        genes2highlight = genes2highlight)
    map.db = map.db[with(map.db, genome.order1 < genome.order2),]
    with(map.db, segments(x0 = genome.order1,
                          x1 = genome.order2,
                          y0 = tp1,
                          y1 = tp2,
                          col = highlight.col,
                          lwd = highlight.lwd))
  }
  ########################################################

  ########################################################
  with(chr.db,
       segments(x0 = genome.order,
                x1 = genome.order,
                y0 = prop.chr.start,
                y1 = prop.chr.end,
                col = chr.col,
                lwd = chr.lwd))
  ########################################################

  ########################################################
  return(list(map.db = map.db,
              blk.db = blk.db,
              chr.db = chr.db))
}

