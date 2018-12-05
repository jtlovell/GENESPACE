#' @title Find the coordinates of blocks, based on the original annotations
#'
#' @description
#' \code{find_blkCoords} Reduces the total number of block breakpoints, so that, when
#' blocks are concatenated, there are fewer small orphan blocks.
#'
#' @param blk The block object (data.frame/data.table)
#' @param map The map object (data.frame/data.table)
#' @param genomeIDs Character vector indicating the genome IDs to consider
#' @param gff.dir The directory containing gff3-formatted annotations
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details ...
#' @return A list of block datasets
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
find_blkCoords = function(blk,
                          map,
                          gff.dir,
                          genomeIDs,
                          verbose = T){

  map = data.table(map)
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

  setkey(gff, genome, order)
  rblk = blk[blk$genome1 == genomeIDs[1] & blk$genome2 == genomeIDs[2],]
  gff1 = split(gff[gff$genome == genomeIDs[1],], gff$chr[gff$genome == genomeIDs[1]])
  gff2 = split(gff[gff$genome == genomeIDs[2],], gff$chr[gff$genome == genomeIDs[2]])

  if(verbose)
    cat("Splitting gff by block coordinates\n")
  spl.gff = lapply(1:nrow(rblk), function(i){
    x = rblk[i,]
    rgff = gff1[[x$chr1]]
    agff = gff2[[x$chr2]]
    ro = rgff[with(rgff, end>=x$start1 & start<=x$end1),]
    ro$block.id = x$block.id
    ao = agff[with(agff, end>=x$start2 & start<=x$end2),]
    ao$block.id = x$block.id
    return(list(ref.gff = ro, alt.gff = ao))
  })
  names(spl.gff)<-rblk$block.id

  if(verbose)
    cat("Finding and splitting non-reference/alternative genome blocks\n")
  spl.blk = rbindlist(lapply(1:nrow(rblk), function(i){
    x = rblk[i,]

    rgff = spl.gff[[as.character(x$block.id)]]$ref.gff
    agff = spl.gff[[as.character(x$block.id)]]$alt.gff
    gf = rbind(rgff,agff)
    map.ingff = map[map$id1 %in% gf$id | map$id2 %in% gf$id,]
    map.ingff = map.ingff[with(map.ingff,
                               genome1 == genomeIDs[1] &
                                 genome2 != genomeIDs[2]),]

    out.old = gf[,list(chr = chr[1],
                       start = min(start),
                       end = max(end)),
                 by = list(genome, block.id)]

    if(nrow(map.ingff)>0){
      out.new = map.ingff[,list(chr = chr2[1],
                                start = min(start2),
                                end = max(end2)),
                          by = list(genome2,block.id)]
      setnames(out.new,1,"genome")
      out = rbind(out.old, out.new)
    }else{
      out = out.old
    }
    out$breakpoint.id = i
    return(out)
  }))
  if(verbose)
    cat("Done! Returning dataset with",
        length(unique(spl.blk$breakpoint.id)),"block networks\n")
  return(spl.blk)
}
