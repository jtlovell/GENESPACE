#' @title Make a database of block positions
#'
#' @description
#' \code{make_blockDB}
#'
#' @param blk The block object
#' @param map The map object
#' @param genomeIDs Genome identifiers. The first is used as the reference coordinate
#' system
#' @param gff.dir Directory containing the gff3-formatted annotations.
#' @param block.dir Directory to write block data.
#' @param verbose Logical, should updates be printed?
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
make_blockDB = function(blk,
                        map,
                        genomeIDs,
                        gff.dir,
                        block.dir,
                        min.genesInBlock = 10,
                        min.blockBp = 5e4,
                        verbose = T,
                        ...){
  if(verbose)
    cat("Building block database and file structure\n",
        "\t... data will be stored in",block.dir,"\n")

  if(verbose)
    cat("Parsing gff files\n")
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
                       ties.method = "dense")
    return(tmp)
  }))
  setkey(gff, genome, id)

  if(verbose)
    cat("\tSplitting gff files by breakpoints\n")

  gff.spl = sapply(unique(bp$genome), simplify = F, USE.NAMES = T, function(i){
    gg = gff[gff$genome == i,]
    out = sapply(unique(gg$chr), simplify = F, USE.NAMES = T, function(j){
      gg[gg$chr == j,]
    })
    return(out)
  })

  bp.spl = split(bp, bp$breakpoint.id)
  gff.bp = sapply(names(bp.spl), simplify = F, USE.NAMES = T, function(i){
    x = bp.spl[[i]]
    return(rbindlist(lapply(1:nrow(x), function(j){
      y = x[j,]
      go = gff.spl[[y$genome]][[y$chr]]
      go = go[go$end >= y$start & go$start <= y$end,]
    })))
  })


  if(verbose)
    cat("\t... finding physical position of breakpoints\n")

  if(verbose)
    cat("Writing files\n\t... splitting up assemblies into blocks\n")

  if(verbose)
    cat("Splitting up peptide fastas into blocks\n")

  if(verbose)
    cat("Splitting up CDS fastas into blocks\n")

  if(verbose)
    cat("Splitting up full-length transcript fastas into blocks\n")

  if(verbose)
    cat("Splitting up annotations into blocks\n")

  if(verbose)
    cat("Done\n")

}

