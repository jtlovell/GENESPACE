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

  mog = map[,c("genome1","id1","og1")]
  mog = mog[!duplicated(mog),]
  setnames(mog,3,"og")
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
