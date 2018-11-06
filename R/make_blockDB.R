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
                        verbose = T,
                        ...){
  if(verbose)
    cat("Building block database and file structure\n",
        "\t... data will be stored in",block.dir,"\n")
  if(verbose)
    cat("\t... compiling breakpoints\n")
  rmap = map[map$genome1 == genomeIDs[1],]
  rblk = blk[blk$genome1 == genomeIDs[1],]
  brpts = rbind(data.frame(chr = rblk$chr1, pos = rblk$start1,
                           stringsAsFactors = F),
                data.frame(chr = rblk$chr1, pos = rblk$end1,
                           stringsAsFactors = F))
  brpts = brpts[!duplicated(brpts),]
  brpts = brpts[order(brpts$chr, brpts$pos),]
  brpts = rbindlist(lapply(split(brpts, brpts$chr), function(x){
    x$start = c(x$pos[-nrow(x)],NA)
    x$end = c(x$pos[-1],NA)
    return(x)
  }))
  brpts = data.frame(brpts[complete.cases(brpts),])

  if(verbose)
    cat("\t... splitting map file into breakpoint-informed blocks\n")
  rg.list = lapply(1:nrow(brpts), function(i){
    x = brpts[i,]
    xmap = rmap[with(rmap, chr1 == x$chr &
                       end1 >= x$start &
                       start1 <= x$end),]
    if(nrow(xmap)==0){
      return(NULL)
    }else{
      xblk = make_blocks(xmap)
      xmap = xblk$map
      xblk = xblk$block
      xblk$breakpoint.group = paste0(x$chr,"_",i)
      xmap$breakpoint.group = paste0(x$chr,"_",i)

      xspl = split(xmap$id1, xmap$block.id)
      return(list(block = xblk, map = xmap))
    }
  })
  names(rg.list)<-paste0(brpts$chr,"_",1:nrow(brpts))
  rg.list<-rg.list[!sapply(rg.list, is.null)]


  if(verbose)
    cat("\t... parsing gff files\n")
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
    cat("\t... generating annotation database by breakpoint-blocks\n")

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

