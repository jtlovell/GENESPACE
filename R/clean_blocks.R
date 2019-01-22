#' @title Clean form_syntenicBlocks output
#'
#' @description
#' \code{clean_blocks} Clusters hits and drops low-confidence blocks.
#'
#' @param map the map data.table or data.frame
#' @param blk the block data.table or data.frame
#' @param rerank logical, should the ranks be re-calculated prior to cleaning?
#' @param radius numeric, what should the radius of 2d density clustering be?
#' @param n.mappings numeric, how many mappings are required for a cluster?
#' @param n.cores The number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @import parallel
#' @import dbscan
#' @export
clean_blocks <- function(map,
                         rerank = TRUE,
                         radius = 10,
                         n.mappings = 5,
                         n.cores = 1,
                         clean.by.unique.genes = T,
                         min.unique.genes = ceiling(n.mappings/2),
                         min.unique.og = 2,
                         clean.by.og = F,
                         verbose = T){

  map <- data.table(map)
  setkey(map, chr1, chr2, start1, start2)
  if(rerank){
    map[,rank1 := frank(start1,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr2)]
  }

  if (verbose)
    cat("Cleaning / merging via dbscan\n")

  map$unique = with(map, paste(genome1, genome2, chr1, chr2))
  map$unique.genome = with(map, paste(genome1, genome2))

  setkey(map, chr1, chr2, start1, start2)


  if (verbose)
    cat("Calculating 2d densities for ... \n\t")
  spl.gen = split(map, "unique.genome")
  merged_map<-rbindlist(lapply(spl.gen, function(x){
    g1 = x$genome1[1]
    g2 = x$genome2[1]
    if (verbose)
      cat(g1,"-->",g2,paste0("(initial hits = ",nrow(x),") ... "))
    spl.map <- split(x, "unique")
    chr.map <- rbindlist(mclapply(spl.map, mc.cores = n.cores, mc.preschedule = F, function(tmp){
      x <- run_dbs(y = tmp[, c("rank1","rank2"), with = F],
                   eps.radius = radius,
                   mappings = n.mappings)
      tmp$block.id <- x$cluster
      return(tmp)
    }))
    chr.map <- chr.map[chr.map$block.id != 0, ]
    if (verbose)
      cat(nrow(chr.map), "hits in",
          length(unique(paste(chr.map$unique, chr.map$block.id))),"blocks\n\t")
    return(chr.map)
  }))



  merged_map$block.id <- with(merged_map,
                              as.numeric(as.factor(paste(unique, block.id))))

  if(clean.by.unique.genes){
    merged_map[,n.genes1 := length(unique(id1)), by = list(block.id)]
    merged_map[,n.genes2 := length(unique(id2)), by = list(block.id)]
    merged_map <- merged_map[with(merged_map,
                                    n.genes1 >= min.unique.genes &
                                    n.genes2 >= min.unique.genes),]
  }
  if(clean.by.og){
    merged_map[,n.orthogroups := length(unique(orthogroup)), by = list(block.id)]
    merged_map <- merged_map[with(merged_map,
                                  n.orthogroups >= min.unique.og),]
  }


  merged_blk <- make_blocks(map = merged_map,
                            rename.blocks = T,
                            rerank = T,
                            clean.columns = T,
                            ties.method = "dense")

  if (verbose)
    cat("Cleaned n blocks / mappings =",
        nrow(merged_blk$block),
        "/",
        nrow(merged_blk$map),"\n\tDone!\n")

  return(merged_blk)
}
