#' @title Add ortholog/paralog designation
#'
#' @description
#' \code{add_orthology} Add ortholog/paralog designation from orthofinder
#' with gene trees to the map object
#'
#' @param dir.list The directory list produced by check_environment
#' @param map The map data.frame or data.table
#' @param blast data.table with at least the 10 necessary columns for blast format.
#' @param genomeIDs character vector giving genome IDs to consider.
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param n.cores Number of parallel processes to run, when possible
#' @param rank.buffer The buffer, in gene rank order.
#' @param verbose Logical, should updates be printed
#' @param ... Not currently in use
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom parallel mclapply
#' @export
add_orthology <- function(of.dir,
                          map,
                          verbose = T,
                          block.orthology.threshold = .5){
  if(verbose)
    cat("Pulling orthologs ... \n")
  m <- data.table(map[,c("genome1","genome2","id1","id2","og.id")])
  og.dir <- dirname(list.files(file.path(of.dir, "OrthoFinder"),
                               pattern = "Orthologues_",recursive = T,
                               include.dirs = T, full.names = T))[1]
  ortholog.dirs <- list.files(og.dir, pattern = "^Orthologues_",
                              include.dirs = T, full.names = T)
  og.out <- rbindlist(lapply(ortholog.dirs, function(x){
    fs <- list.files(x, full.names = T)
    if(verbose)
      cat("\tRunning", gsub("Orthologues_","",basename(x)), "... ")

    rd <- rbindlist(lapply(fs, function(y) {

      tmp <- fread(y)

      allgenes <- rbindlist(apply(tmp[,c(2:3),with = F],1,function(z){
        eg <- c(strsplit(z[1],",")[[1]],strsplit(z[2],", ")[[1]])
        return(expand.grid(eg, eg))
      }))
      setnames(allgenes,c("id1","id2"))
      allgenes[,genome1 := colnames(tmp)[2]]
      allgenes[,genome2 := colnames(tmp)[3]]

      return(allgenes)
    }))

    if(verbose)
      cat("Done!\n")
    return(rd)
  }))

  og.out[,is.ortholog := TRUE]
  map.out <- merge(og.out, map, by = colnames(og.out)[1:4], all.y = T)
  map.out$is.ortholog[map.out$id1 == map.out$id2] <- TRUE
  map.out$is.ortholog[is.na(map.out$is.ortholog)] <- FALSE

  # map.out[,prop.og := sum(is.ortholog) / .N,by = "block.id"]
  # map.out[,block.is.ortholog := prop.og > block.orthology.threshold]

  return(map.out)
}
