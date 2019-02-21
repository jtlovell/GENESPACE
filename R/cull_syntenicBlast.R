#' @title Cull blast to syntenic regions
#'
#' @description
#' \code{cull_syntenicBlast} An internal function, designed to be called
#' by extend_blocks and find_syntenicOrthogs.
#'
#' @param map the map data.table or data.frame. This is used to infer
#' syntenic regions.
#' @param blast the blast dataset to screen for syntenic hits
#' @param gff gff data.table
#' @param rank.buffer The buffer, in rank order of gff annotation
#' data.table.
#' @param plotit Logical, should plots be made? If so, and n.cores
#' is 1, then plotting is done.
#' @param n.cores The number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom dbscan frNN
#' @export
cull_syntenicBlast <- function(map,
                               blast,
                               gff,
                               rank.buffer = 250,
                               verbose = T,
                               plotit = F,
                               n.cores = 1){
  #######################################################
  #######################################################
  cull_blast2MapChr <- function(map, blast){
    map$unique.genome <- with(map,
                              paste(genome1, genome2))
    blast$unique.genome <- with(blast,
                                paste(genome1, genome2))
    map$unique.chr <- with(map,
                           paste(unique.genome,
                                 chr1, chr2))
    blast$unique.chr <- with(blast,
                             paste(unique.genome,
                                   chr1, chr2))

    ugm <- map$unique.chr
    ugb <-  blast$unique.chr
    gug <- intersect(unique(ugb), unique(ugm))
    wh.blast <- which(ugb %in% gug)

    blast <- data.table(blast[wh.blast,])
    return(list(map = map,
                blast = blast))
  }
  #######################################################
  #######################################################
  cull_blast2NewIds <- function(blast, map){
    setkey(map, id1, id2)
    setkey(blast, id1, id2)
    t.map <- map[,c("id1","id2")]
    t.map$in.map <- TRUE

    blast <- merge(blast, t.map, all.x = T)
    blast <- blast[is.na(blast$in.map),colnames(blast), with = F]
    return(blast)
  }
  #######################################################
  #######################################################
  find_whichInBuffer <- function(x,
                                 y,
                                 which.in.blk,
                                 rank.buffer,
                                 ...){

    nn <- frNN(x = data.frame(x,
                              y),
               eps = rank.buffer,
               ...)

    all.near.blk <- unique(unlist(nn$id[which.in.blk]))

    return(all.near.blk[order(all.near.blk)])
  }
  #######################################################
  #######################################################
  find_hitsInBuffer <- function(map,
                                blast,
                                rank.buffer,
                                verbose,
                                plotit,
                                n.cores = 1,
                                ...){

    spl.map <- split(map, "unique.genome")
    spl.blast <- split(blast, "unique.genome")
    ns <- unique(names(spl.blast))
    ns <- ns[ns %in% unique(names(spl.map))]

    res.by.genome <- lapply(ns, function(i){
      i.map = spl.map[[i]]
      i.blast = spl.blast[[i]]
      if(verbose)
        cat(paste0("\t",i," ... (new.hits = ",
                   nrow(i.blast),", ","map.size = ",nrow(i.map),")"))


      spl.i.map <- split(i.map, "unique.chr")
      spl.i.blast <- split(i.blast, "unique.chr")
      nis <- unique(names(spl.i.blast))
      nis <- nis[nis %in% unique(names(spl.i.map))]

      ids2keep <- rbindlist(mclapply(nis, mc.cores = n.cores, function(j){
        j.map = spl.i.map[[j]]
        j.blast = spl.i.blast[[j]]
        j.blast <- j.blast[with(j.blast,
                                rank1 >= (min(j.map$rank1)-(rank.buffer*2)) &
                                  rank2 >= (min(j.map$rank2)-(rank.buffer*2)) &
                                  rank1 <= (max(j.map$rank1)+(rank.buffer*2)) &
                                  rank2 <= (max(j.map$rank2)+(rank.buffer*2))), ]
        j.out <- rbind(j.map[,c("id1","id2")],
                       j.blast[,c("id1","id2")])
        x = c(j.map$rank1, j.blast$rank1)
        y = c(j.map$rank2, j.blast$rank2)
        wh = 1:nrow(j.map)

        tokeep <- find_whichInBuffer(x = x,
                                     y = y,
                                     which.in.blk = wh,
                                     rank.buffer = rank.buffer,
                                     ...)
        if(plotit){
          map.tp = cbind(j.map$rank1, j.map$rank2)
          int = intersect(tokeep, (nrow(j.map)+1):length(x))
          plot(x, y, pch = ".",
               xlab = paste(j.map$genome1[1],
                            j.map$chr1[1],"gene order"),
               ylab = paste(j.map$genome2[1],
                            j.map$chr2[1],"gene order"))
          points(map.tp, cex = .5, col = "dodgerblue")
          points(x[int], y[int], cex = .3, col = "darkred")
        }

        return(j.out[tokeep,])
      }))

      ids2keep <- ids2keep[!duplicated(ids2keep),]
      if(verbose)
        cat(" returning", nrow(ids2keep),"\n")
      return(ids2keep)
    })
    out <- rbindlist(res.by.genome)
    setkey(out, id1, id2)
    return(out)
  }
  #######################################################

  #######################################################
  #######################################################
  cull_blast2MapChr <- cmpfun(cull_blast2MapChr)
  cull_blast2NewIds <- cmpfun(cull_blast2NewIds)
  find_whichInBuffer <- cmpfun(find_whichInBuffer)
  find_hitsInBuffer <- cmpfun(find_hitsInBuffer)
  #######################################################
  #######################################################

  #######################################################
  #######################################################
  if (verbose)
    cat("Dropping chromosome combinations in blast not found in map... ")
  tmp <- cull_blast2MapChr(blast = blast, map = map)
  m.blast <- tmp$blast
  map <- tmp$map
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Dropping blast hits in the map already ... ")
  t.blast <- cull_blast2NewIds(blast = m.blast, map = map)
  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Making new blast and map with gff-based ranks ... ")
  id.names <- c("genome1","genome2","id1","id2","chr1","chr2","start1","start2","end1","end2")
  all.ids <- rbind(map[, id.names, with = F],
                   blast[, id.names, with = F])
  all.ids <- all.ids[!duplicated(all.ids),]
  all.ids[,rank1 := frank(start1, ties.method = "dense"),
          by = list(genome1, chr1)]
  all.ids[,rank2 := frank(start2, ties.method = "dense"),
                by = list(genome2, chr2)]
  setkey(all.ids, id1, id2)
  setkey(map, id1, id2)
  setkey(blast, id1, id2)
  r.map <- merge(all.ids, map[,c("id1","id2")])
  r.blast <- merge(all.ids, blast[,c("id1","id2")])

  if (verbose)
    cat("Done!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Finding neighbors to mappings in blast hits ... completed:\n")
  r.map$unique.genome <- with(r.map,
                            paste(genome1, genome2))
  r.blast$unique.genome <- with(r.blast,
                              paste(genome1, genome2))
  r.map$unique.chr <- with(r.map,
                         paste(unique.genome,
                               chr1, chr2))
  r.blast$unique.chr <- with(r.blast,
                           paste(unique.genome,
                                 chr1, chr2))
  ids2keep <- find_hitsInBuffer(map = r.map,
                                blast = r.blast,
                                rank.buffer = rank.buffer,
                                verbose = verbose,
                                plotit = plotit)

  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reformatting blast output ... ")

  setkey(blast, id1, id2)

  out.blast <- merge(ids2keep, blast)
  re.out <- with(out.blast,
                 rerank_fromIDs(id1 = id1,
                                id2 = id2,
                                gff = gff))

  if (verbose)
    #######################################################
  cat("Done!\n")
  return(re.out)
}
