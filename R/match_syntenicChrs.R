#' @title match syntenic chromosome IDs
#'
#' @description
#' \code{match_syntenicChrs} Tabulate the number of hits between
#' pairwise genome chromomes and build a database for plotting.
#'
#' @param map map results data.table
#' @param genomeIDs character vector of genomeIDs
#' @param genome1.chrs Chromosome IDs of the first genome IDs to order by.
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
match_syntenicChrs <- function(map,
                               genomeIDs,
                               genome1.chrs,
                               ...){

  m <- map[,c("block.id",
              "genome1","id1","chr1","start1","end1","rank1",
              "genome2","id2","chr2","start2","end2","rank2")]
  m2 <- m[,c(1,8:13,2:7)]
  setnames(m2, colnames(m))
  map <- rbind(m, m2)
  map <- map[!duplicated(map),]

  dspl <- lapply(1:(length(genomeIDs)-1), function(i)
    map[map$genome1 == genomeIDs[i] &
          map$genome2 == genomeIDs[i+1],])

  chr.index <- lapply(dspl, function(x){
    chr.tab <- with(x, table(chr1, chr2))
    wh.row <- apply(chr.tab, 1, which.max)
    wh.col <- apply(chr.tab, 2, which.max)

    g1.df <- rbind(data.frame(genome1 = x$genome1[1],
                              genome2 = x$genome2[1],
                              chr1 = names(wh.row),
                              chr2 = colnames(chr.tab)[wh.row]),
                   data.frame(genome1 = x$genome1[1],
                              genome2 = x$genome2[1],
                              chr1 = rownames(chr.tab)[wh.col],
                              chr2 = names(wh.col)))
    return(data.table(g1.df[!duplicated(g1.df),]))
  })

  nlevs <- genome1.chrs
  for(i in 1:length(chr.index)){
    chr.index[[i]]$chr1 <- factor(chr.index[[i]]$chr1,
                                  levels = nlevs)
    setkey(chr.index[[i]], chr1)
    nlevs <- chr.index[[i]]$chr2[!duplicated(chr.index[[i]]$chr2)]
    chr.index[[i]]$chr2 <- factor(chr.index[[i]]$chr2,
                                  levels = nlevs)
  }

  chr.list <- c(list(as.character(genome1.chrs)),
                lapply(chr.index, function(x) levels(x$chr2)))
  names(chr.list) <- genomeIDs
  return(chr.list)
}
