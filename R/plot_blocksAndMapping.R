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
#' @param plotit Logical, Should the plot be drawn?
#' @param returnData Logical, Should the data used to plot be returned?
#' @param ... Not currently in use
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
                                  plotit = T,
                                  returnData = F,
                                  colorSegment = T,
                                  gff.dir = NULL,
                                  blast.results = NULL){

  if(!is.null(gff.dir) & !is.null(blast.results)){
    gff.files <- list.files(gff.dir,
                            full.names = T)
    names(gff.files) <- gsub(".gff3$", "",
                             basename(gff.files))

    parse_gff <- function(gff){
      g <- suppressWarnings(
        data.table::fread(gff,
                          showProgress = F,
                          verbose = F))
      g <- g[g$V3 == "gene", c(9, 1, 4, 5, 7)]
      g$V9 <- sapply(g$V9, function(x) gsub("Name=", "",
                                            strsplit(x, ";")[[1]][2]))
      data.table::setnames(g, c("id", "chr", "start", "end", "strand"))
      return(g)
    }

    gff <- rbindlist(lapply(names(gff.files), function(i){
      tmp <- parse_gff(gff.files[[i]])
      tmp$genome <- i
      tmp$order <- frank(tmp[,c("chr", "start")],
                         ties.method = "random")
      return(tmp)
    }))
    setkey(gff, genome, chr, start)

    gff[,rank := frank(start),
        by = list(genome,chr)]
    gff1 = data.table(gff)
    gff2 = data.table(gff)
    setnames(gff1, paste0(colnames(gff1),1))
    setnames(gff2, paste0(colnames(gff2),2))
    setkey(gff1,genome1,id1)
    setkey(gff2,genome2,id2)

    blast.simp = blast.results[,c("id1","genome1","id2","genome2"), with = F]
    setkey(blast.simp, id2, genome2)
    gffs = merge(gff2, blast.simp)
    setkey(gffs,genome1,id1)
    gffs = merge(gff1, gffs)

    map.simp = map[,c("id1","genome1","id2","genome2","block.id","score"),with = F]
    map.simp$block.id = as.numeric(as.factor(with(map.simp,paste(genome1, genome2, block.id))))

    setkey(map.simp, id1,genome1,id2, genome2)
    setkey(gffs, id1,genome1,id2, genome2)
    ms = merge(gffs, map.simp, all = T)

    tb = make_blocks(ms, rerank = F, drop.NAs = T)
    map <- tb$map
    blk <- tb$block
  }

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
  if(!is.null(chr1toplot)){
    tpb <- tpb[tpb$chr1 %in% chr1toplot,]
    tp <- tp[tp$chr1 %in% chr1toplot,]
  }
  if(!is.null(chr2toplot)){
    tpb <- tpb[tpb$chr2 %in% chr2toplot,]
    tp <- tp[tp$chr2 %in% chr2toplot,]
  }
  sb1 <- split(tpb, tpb$chr1)
  st1 <- split(tp, tp$chr1)
  for(i in names(sb1)){
    sb2 <- split(sb1[[i]], sb1[[i]]$chr2)
    st2 <- split(st1[[i]], st1[[i]]$chr2)
    for(j in names(sb2)){
      t2 <- st2[[j]]
      b2 <- sb2[[j]]
      cols <- rep_len(c("red3", "salmon", "darkorange3", "gold",
                        "grey50", "lightgreen", "forestgreen", "darkgreen",
                        "cyan", "dodgerblue3", "violet", "purple"),
                      length.out = nrow(b2))
      if(colorSegment){
        with(t2, plot(rank1, rank2,
                      col = rgb(0,0,0,.5),
                      pch = 16,
                      cex = .5,
                      ylab = paste(altGenome2plot, j, "gene order"),
                      xlab = paste(ref.id, i, "gene order"),
                      main = main))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 =rankend2,
                          col = cols[as.numeric(as.factor(block.id))],
                          lwd = 2))
      }else{
        with(t2, plot(rank1, rank2,
                      col = cols[as.numeric(as.factor(block.id))],
                      pch = 16,
                      cex = .5,
                      ylab = paste(altGenome2plot, j, "gene order"),
                      xlab = paste(ref.id, i, "gene order"),
                      main = main))
        with(b2, segments(x0 = rankstart1, x1 = rankend1,
                          y0 = rankstart2, y1 =rankend2,
                          col = "black",
                          lwd = 1.5))
      }

      with(b2, text(x = rowMeans(b2[,c("rankstart1","rankend1")]),
                    y = rowMeans(b2[,c("rankstart2","rankend2")]),
                    labels = block.id,
                    col = "black",
                    cex = .5,
                    adj = c(1,-1)))
      if(returnData)
        return(list(t2,
                    b2))
    }
  }
}

