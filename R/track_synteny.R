#' @title track synteny
#'
#' @description
#' \code{track_synteny} Track most likely position for all genes
#' with at least one syntenic ortholog hit
#'
#' @param map map results data.table
#' @param genomeIDs character vector of genomeIDs
#' @param wind.size Size of window to use to search for the best hit region.
#' Smaller values are prone to larger errors, but larger values are slower.
#' @param quantiles The quanitles of the window locaion to use to infer
#' the position of the gene.
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom compiler cmpfun
#' @export
track_synteny <- function(map,
                          genomeIDs,
                          wind.size = 5,
                          quantiles = c(.4,.6),
                          verbose = T){

  ########################################################
  ########################################################
  make_genomeWindow <- function(map,
                                genome,
                                wind.size = 5){

    # -- Keep only the focal genome
    m = map[map$genome1 == genome & map$genome1 != map$genome2,
            c("block.id","genome1","genome2",
              "id1","id2",
              "chr1","chr2",
              "start1","start2","end1","end2","score")]
    mt = map[map$genome2 == genome & map$genome1 != map$genome2,
             c("block.id","genome2","genome1",
               "id2","id1",
               "chr2","chr1",
               "start2","start1","end2","end1","score")]
    if(nrow(mt)>0){
      setnames(mt, c("block.id","genome1","genome2",
                     "id1","id2",
                     "chr1","chr2",
                     "start1","start2","end1","end2","score"))
      if(nrow(m) > 0){
        m <- data.table(rbind(m, mt))
      }else{
        m <- data.table(mt)
      }
    }

    # -- re-rank by all genes present in map
    spl = split(m, "chr1")
    m <- rbindlist(lapply(spl, function(x){
      x$rank <- frank(x, start1, end1, id1, ties.method = "dense")
      return(x)
    }))

    # -- keep only the single best hits / block and gene ID
    m[ , orderscore := frank(-score,
                             ties.method = "random"),
       by = list(id1,block.id)]

    x = m[m$orderscore == 1, ]
    x$orderscore <- NULL
    x$score <- NULL

    # -- make a set of left offset ranks
    xn = rbindlist(lapply(wind.size:0, function(i){
      y = x
      y$rankn = y$rank-i
      return(y)
    }))

    # -- make a set of right offset ranks
    xp = rbindlist(lapply(1:wind.size, function(i){
      y = x
      y$rankn = y$rank+i
      return(y)
    }))

    # -- combine and return
    xc <- rbind(xn, xp)

    xc[, true.id1 := id1[rank == rankn][1],
       by = list(chr1, rankn)]
    xc <- xc[complete.cases(xc),]
    return(xc)
  }
  ########################################################
  ########################################################
  link_regions <- function(genome.window,
                           quantiles){
    if(length(quantiles) != 2)
      stop("quantiles must be a numeric vector [0,1] of length 2\n")
    xc <- genome.window

    # -- Pull out those with hits
    matched <- xc[xc$rankn == xc$rank,]
    matched$rankn <- NULL
    matched$true.id1 <- NULL

    # -- Pull out those without matches
    no.match <- xc[!with(xc, paste(true.id1, genome2, block.id)) %in%
                     with(matched, paste(id1, genome2, block.id)),
                   c("true.id1", "genome2", "block.id")]
    no.match <- no.match[!duplicated(no.match),]
    nc <- merge(xc,
                no.match,
                by = colnames(no.match))

    # -- Make database of best hits
    xo <- nc[,list(chr2 = chr2[1],
                   rank = rankn[1],
                   start2 = round(quantile(start2, quantiles[1])),
                   end2 = round(quantile(end2, quantiles[2])),
                   n = length(unique(id1))),
             by = list(true.id1, genome2, block.id)]

    setnames(xo,1,"id1")

    mi1 = map[,c("genome1","id1","chr1","start1","end1")]
    mi2 = map[,c("genome2","id2","chr2","start2","end2")]
    setnames(mi2, colnames(mi1))
    mi <- rbind(mi1, mi2)
    mi = mi[!duplicated(mi),]
    setkey(mi, id1)
    setkey(xo, id1)
    out <- merge(mi, xo)
    out$id2 <- with(out, paste(genome2, genome1, 1:nrow(out), sep = "."))
    matched$n <- 0
    return(rbind(matched, out))
  }
  ########################################################
  ########################################################
  make_genomeWindow <- cmpfun(make_genomeWindow)
  link_regions <- cmpfun(link_regions)
  ########################################################
  ########################################################

  ########################################################
  if (verbose)
    cat("Completing pairwise synteny database for every gene in ... \n")
  join2 <- lapply(genomeIDs, function(i){
    if (verbose)
      cat(paste0("\t",i,": "))
    gw = make_genomeWindow(map = map,
                           genome = i,
                           wind.size = wind.size)
    reg <- link_regions(genome.window = gw,
                        quantiles = quantiles)
    if(verbose)
      cat("Found", nrow(reg), "unique links, ",
          sum(reg$n > 0), "have no ortholog\n")
    return(reg)
  })
  join2 <- rbindlist(join2)
  if (verbose)
    cat("\tDone!\n")
  return(join2)
}
