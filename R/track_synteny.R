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
track_synteny <- function(genomeIDs,
                          map,
                          gff,
                          quantiles = c(.4,.6),
                          method = "gff",
                          n.cores = 1,
                          cull2map = TRUE,
                          max.reg.size = 1e5,
                          wind.size = 2,
                          verbose = T){
  if(method == "gff"){
    if(!"array.id" %in% colnames(map))
      stop("If method = gff, must provide map with arrays inferred.\n")
    ########################################################

    ########################################################
    eg <- expand.grid(genomeIDs, genomeIDs)
    eg <- eg[eg[,1] != eg[,2],c(2:1)]
    mapd2 <- rbindlist(lapply(1:nrow(eg), function(i)
      pipe_track(genome1 = eg[i,1],
                 genome2 = eg[i,2],
                 gff = gff,
                 cull2map = F,
                 n.cores = n.cores,
                 map.ta = map.ta,
                 verbose = verbose)))
    ########################################################

    ########################################################
    mapd2$n[!mapd2$no.hit] <- 0
    mapd2$u <- with(mapd2, paste0(genome1, chr1))
    spl = split(mapd2, "u")
    mapd3 <- rbindlist(lapply(spl, function(x) {
      x$rank <- frank(x, start1, end1, id1, ties.method = "dense")
      return(x)
    }))
    mapd3$u <- NULL
    og <- syn.og$map[,c("id1","og1")]
    og <- og[!duplicated(og),]
    mapd4<-merge(mapd3, og, by = "id1", all.x = T)
    mapd4$og<- mapd4$og1
    mapd4$og1<-NULL
    ########################################################

    ########################################################
    map.in <- mapd4[with(mapd4, is.na(id2) & (end2 - start2) < max.reg.size),]
    return(list(nomap = map.in, allmap = mapd4))
  }else{
    ########################################################

    ########################################################
    mog = map[,c("genome1","id1","og1")]
    mog = mog[!duplicated(mog),]
    setnames(mog,3,"og")
    ########################################################

    ########################################################
    if (verbose)
      cat("Completing pairwise synteny database for every gene in ... \n")
    join2 <- lapply(genomeIDs, function(i){
      if (verbose)
        cat(paste0("\t",i,": "))
      gw = make_genomeWindow2(map = map,
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
    mapd4 <- merge(join2, mog, by = c("genome1","id1"), all.x = T)
    map.in <- mapd4[with(mapd4, n > 0 & (end2 - start2) < max.reg.size),]

    if (verbose)
      cat("\tDone!\n")
    return(list(nomap = map.in, allmap = mapd4))
  }
}
