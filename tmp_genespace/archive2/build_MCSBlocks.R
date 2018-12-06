#' @title Build out the MCScanX collinear blocks
#'
#' @description
#' \code{build_MCSBlocks} Runs the underlying `run_MCScanX` function both across
#' and within pairwise genome comparisons.
#'
#' @param blast.results R object contain the blast results
#' @param MCScanX.params Final arameters to pass to MCScanX
#' @param MCScanX.params.init Initial parameters to pass to MCScanX. These should be less stringent
#' than MCScanX.params
#' @param mcscanx.input.dir Directory containing the MCScanX-formatted mapping files
#' @param abbrevs Genome abbreviations
#' @param onlyPairwise Logical, should the genomes only be considered in a pairwise fashion?
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
build_MCSBlocks <- function(blast.results,
                            MCScanX.param.init = "-a -s 5 -m 25 -w 5",
                            MCScanX.params = "-a -s 20 -m 10 -w 5",
                            abbrevs,
                            mcscanx.input.dir,
                            verbose = TRUE,
                            onlyPairwise = TRUE,
                            eps.radius = 5e4,
                            ...){


  find_unmappedBlast <- function(blk,
                                 map,
                                 blast){

    blk$uniq <- with(blk, paste(genome1, genome2, chr1, chr2))
    map$uniq <- with(map, paste(genome1, genome2, chr1, chr2))
    blast$uniq <- with(blast, paste(genome1, genome2, chr1, chr2))
    map$uniq.map <- with(map, paste(genome1, genome2, id1, id2))
    blast$uniq.map <- with(blast, paste(genome1, genome2, id1, id2))

    blast.noblk <- blast[!blast$uniq.map %in% map$uniq.map,]

    spl.blk = split(blk, blk$uniq)
    blast.noblk <- blast.noblk[blast.noblk$uniq %in% unique(map$uniq),]

    spl.blast <- split(blast.noblk, factor(blast.noblk$uniq))
    blast.noblk <- rbindlist(lapply(names(spl.blast), function(i){
      x.blk <- spl.blk[[i]]
      x.blast <- spl.blast[[i]]
      if(nrow(x.blast) > 1){
        inblk <- sapply(1:nrow(x.blk), function(j){
          x.blast$start1 <= x.blk$end1[j] & x.blast$end1 >= x.blk$start1[j] &
            x.blast$start2 <= x.blk$end2[j] & x.blast$end2 >= x.blk$start2[j]
        })
        outblk <- x.blast[rowSums(inblk)==0,]
      }else{
        outblk <- x.blast
      }
    }))

    allu <- unique(c(map$uniq.map, blast.noblk$uniq.map))
    return(blast[blast$uniq.map %in% allu,])
  }

  add_gapHits2Block = function(blast,
                               map,
                               blk,
                               eps.radius){
    test = find_unmappedBlast(blk = blk,
                              map = map,
                              blast = all.blast.results)
    b = blk
    b$id1 = "block"
    b$id2 = "block"
    test = map
    test$rankstart1 = test$rank1
    test$rankstart2 = test$rank2
    test$rankend1 = test$rank1
    test$rankend2 = test$rank2
    tmp = rbind(b, test, fill = T)[,1:17,with = F]
    tmp$unique = with(tmp, paste(genome1, genome2, chr1, chr2))
    spl = split(tmp, tmp$unique)
    maps2keep = rbindlist(lapply(spl, function(x){
      y1 = x[,c("unique","id1","id2","start1","start2")]
      y2 = x[,c("unique","id1","id2","end1","end2")]
      setnames(y1,4:5,c("pos1","pos2"))
      setnames(y2,4:5,c("pos1","pos2"))
      y = rbind(y1,y2)

      nn = frNN(y[,c("pos1","pos2"), with  =F], eps = eps.radius)
      dbs = dbscan(nn, minPts = 2)$cluster
      y$init.clus = dbs
      grp.clus = unique(y$init.clus[y$id1 == "block"])
      grp.clus = grp.clus[grp.clus!=0]
      yo = y[y$init.clus %in% grp.clus,]
      return(yo[,c("id1","id2"),with = F])
    }))
    ret = merge(test,maps2keep)
    return(ret)
  }


  if(verbose)
    cat("Initial MCScanX run on",
        nrow(blast.results),
        "culled BLAST hits... \n\t")
  if(!onlyPairwise){
    mcscan.cull <- run_MCScanX(
      blast.results = blast.results,
      abbrevs = abbrevs,
      mcscanx.input.dir = mcscanx.input.dir,
      MCScanX.params = MCScanX.param.init,
      verbose = FALSE)
    if(verbose)
      cat("Clustered",
          nrow(mcscan.cull),
          "BLAST hits into",
          length(unique(mcscan.cull$block.id)),
          "colinear blocks\n")

    blast.results.mcs <- merge(blast.results,
                               mcscan.cull[,1:2,with = F])
    if(verbose)
      cat("Secondary MCScanX run ...\n\t")
    mcscan.cull <- run_MCScanX(
      blast.results = blast.results.mcs,
      abbrevs = abbrevs,
      mcscanx.input.dir = mcscanx.input.dir,
      MCScanX.params = MCScanX.param.init,
      verbose = FALSE)
    if(verbose)
      cat("Clustered",
          nrow(mcscan.cull),
          "BLAST hits into",
          length(unique(mcscan.cull$block.id)),
          "colinear blocks\n")

    blast.results.mcs <- merge(blast.results,
                               mcscan.cull[,1:2,with = F])
  }else{
    mcscan.cull <- rbindlist(
      lapply(split(blast.results,
                   with(blast.results,
                        paste(genome1, genome2))), function(x){
                          run_MCScanX(
                            blast.results = x,
                            abbrevs = abbrevs,
                            mcscanx.input.dir = mcscanx.input.dir,
                            MCScanX.params = MCScanX.param.init,
                            verbose = FALSE)
                        }))

    blast.results.mcs <- merge(blast.results,
                               mcscan.cull[,1:2,with = F])
  }

  if(verbose)
    cat("MCScanX run within pairwise genome comparisons...\n\t")

  mcscan.cull <- rbindlist(
    lapply(split(blast.results.mcs,
                 with(blast.results.mcs,
                      paste(genome1, genome2))), function(x){
                        run_MCScanX(
                          blast.results = x,
                          abbrevs = abbrevs,
                          mcscanx.input.dir = mcscanx.input.dir,
                          MCScanX.params = MCScanX.params,
                          verbose = FALSE)
                      }))

  blast.results.mcs <- merge(blast.results,
                             mcscan.cull[,1:2,with = F])

  mcscan.list <- lapply(
    split(blast.results.mcs,
          with(blast.results.mcs,
               paste(genome1, genome2))), function(x){
                 mctmp = run_MCScanX(
                   blast.results = x,
                   abbrevs = abbrevs,
                   mcscanx.input.dir = mcscanx.input.dir,
                   MCScanX.params = MCScanX.params,
                   verbose = FALSE)
                 return(make_blocks(mctmp))
               })

  blk <- rbindlist(lapply(mcscan.list, function(x) x$block))
  map <- rbindlist(lapply(mcscan.list, function(x) x$map))

  if(verbose)
    cat("Clustered",
        nrow(map),
        "BLAST hits into",
        nrow(blk),
        "colinear blocks\n")
  if(verbose)
    cat("Finding BLAST hits in gaps between blocks...\n\t")

  blast.unmap <- add_gapHits2Block(map = map,
                                   blk = blk,
                                   blast = blast.results,
                                   eps.radius = eps.radius)
  blast.unmap <- blast.unmap[,colnames(blast.results.mcs), with = F]
  blast.unmap <- rbind(blast.results.mcs, blast.unmap)
  if(verbose)
    cat("Added",
        nrow(blast.unmap)-nrow(map),
        "BLAST hits\n")

  if(verbose)
    cat("Final MCScanX run within pairwise genome comparisons...\n\n")
  mcscan.cull <- rbindlist(
    lapply(split(blast.unmap,
                 with(blast.unmap,
                      paste(genome1, genome2))), function(x){
                        run_MCScanX(
                          blast.results = x,
                          abbrevs = abbrevs,
                          mcscanx.input.dir = mcscanx.input.dir,
                          MCScanX.params = MCScanX.params,
                          verbose = FALSE)
                      }))
  blast.results.mcs <- merge(blast.results,
                             mcscan.cull[, 1:2, with = F])

  mcscan.list <- lapply(split(blast.results.mcs,
                              with(blast.results.mcs,
                                   paste(genome1, genome2))), function(x){
                                     mctmp = run_MCScanX(
                                       blast.results = x,
                                       abbrevs = abbrevs,
                                       mcscanx.input.dir = mcscanx.input.dir,
                                       MCScanX.params = MCScanX.params,
                                       verbose = FALSE)
                                     return(make_blocks(mctmp))
                                   })

  blk <- make_blocks(rbindlist(lapply(mcscan.list, function(x) x$map)))
  map <- blk$map
  blk <- blk$block

  if(verbose)
    cat("Done! ... Clustered",
        nrow(map),
        "BLAST hits into",
        nrow(blk),
        "colinear blocks\n")
  return(list(block = blk,
              map = map))
}
