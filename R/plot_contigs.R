#' @title Plot contigs, gaps and telomeres
#' @description
#' \code{plot_contigs} Plots the output from find_contigsGapsTelos
#' @param cgt The output from a single find_contigsGapsTelos call
#' @param cgtList A named list of output from find_contigsGapsTelos
#' @param nColors Integer, the number of colors to cycle through
#' @param nRow Integer, the number of plot facet rows
#' @param nCol Integer, the number of plot facet columns
#' @param palette an R palette function
#' @param ... additional arguments passed to palette
#'
#' @import data.table
#' @import ggplot2
#' @export
plot_contigs <- function(cgt = NULL,
                         cgtList = NULL,
                         nColors = NULL,
                         nRow = NULL,
                         nCol = NULL,
                         palette, ...){

  genome <- type <- x <- end <- start <- index <- y1 <- chr <- y2 <- y <-
    n <- NULL

  if(is.null(cgt) && is.null(cgtList))
    stop("must either provide a single contig/gap/telo object (cgt) or a list of these (cgtList)\n")

  if(!is.null(cgt) && !is.null(cgtList)){
    warning("only provide either cgt or cgtList. Ignoring cgt\n")
    cgt <- NULL
  }
  if(is.null(cgtList))
    cgtList <- list(genome = cgt)

  chk <- function(x) all(c("gaps", "contigs", "telomeres") %in% names(x))
  if(any(!sapply(cgtList, chk)))
    stop("contig/gap/telo object must be a named (gaps, contigs, telomeres) list\n")
  chk <- function(x) all(sapply(x, class) == "GRanges")
  if(any(!sapply(cgtList, chk)))
    stop("all elements of contig/gap/telo object must be GRanges\n")

  if(!is.function(palette))
    stop("palette must be a function\n")
  tst <- palette(1)
  if(!are_colors(tst))
    stop("palette doesn't appear to be an R color palette\n")

  # -- combine the list
  if(is.null(names(cgtList))){
    warning("no names for cgtList provided. assigning arbitrary names\n")
    names(cgtList) <- sprintf("genome%s", 1:length(cgtList))
  }

  tp <- rbindlist(lapply(names(cgtList), function(j){
    x <- cgtList[[j]]
    y <- rbindlist(lapply(names(x), function(i)
      data.table(as.data.frame(x[[i]])[,1:3], type = i)))
    setnames(y, "seqnames", "chr")
    y[,genome := j]
    return(y)
  }))

  # -- get the vector of colors
  if(is.null(nColors)){
    mxn <- with(tp, max(tapply(
      type, chr, FUN = function(x) sum(x == "contigs"))))
    nColors <- floor(ifelse(mxn/3 < 5, 5, mxn / 3))
    if(nColors > 50)
      nColors <- 50
  }

  cols <- palette(nColors, ...)

  # -- get the plotting stuff ready
  tp[,genome := factor(genome, levels = names(cgtList))]
  tpb <- subset(tp, type == "contigs")
  tel <- subset(tp, type == "telomeres")
  tel[,x := ((end + start)/2)/1e6]

  tpb[,index := rep(1:nColors, nrow(tp))[1:.N], by = c("genome", "chr")]
  tpb[,y1 := as.numeric(factor(
    paste(genome, chr), levels = unique(paste(genome, chr))))]
  tpb[,y2 := y1 - .6]
  tpb[,y := (y1 + y2)/2]
  tpb[,genome := factor(genome, levels = unique(genome))]
  ng <- tpb[,list(n = .N, x = max(end)/1e6), by = c("genome", "chr", "y")]
  tmp <- subset(tpb, !duplicated(paste(chr, genome)))
  yv <- tmp$y2; names(yv) <- with(tmp, paste(chr, genome))
  tel[,y2 := yv[paste(chr, genome)]]
  tpl <- subset(tpb, !duplicated(paste(genome, chr)))

  # -- make the plot
  p <- ggplot()+
    geom_rect(data = tpb,
              aes(xmin = start / 1e6, xmax = end / 1e6,
                  ymin = y2, ymax = y1,
                  fill = factor(index)))+
    scale_y_reverse(
      expand = c(0.01, 0.01),
      breaks = tpl$y,
      labels = tpl$chr,
      name = "Chromosome")+
    scale_x_continuous(
      name = "Physical position (Mb); * indicate telomere sequence")+
    geom_text(data = ng,
              aes(x = x, y = y, label = n),
              hjust = -.25, size = 2)+
    geom_point(data = tel, aes(x = x, y = y2),
               shape = "*", col = "red", size = 4)+

    scale_fill_manual(values = cols, guide = "none")+
    theme(panel.background = element_rect(fill = "darkgrey"),
          panel.grid = element_blank())+
    ggtitle(sprintf(
      "Contig positions: each cycle through colors is %s contigs (%s gaps)",
      nColors, nColors - 1))
  if(!is.null(nRow) & !is.null(nCol)){
    p <- p + facet_wrap(~genome, scales = "free", nrow = nRow, ncol = nCol)
  }else{
    if(!is.null(nRow)){
      p <- p + facet_wrap(~genome, scales = "free", nrow = nRow)
    }else{
      if(!is.null(nCol)){
        p <- p + facet_wrap(~genome, scales = "free", ncol = nCol)
      }else{
        p <- p + facet_wrap(~genome, scales = "free")
      }
    }
  }

  print(p)
}
