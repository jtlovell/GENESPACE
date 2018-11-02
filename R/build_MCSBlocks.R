#' @title Iteratively run MCScanX program
#'
#' @description
#' \code{build_MCSBlocks} Runs the underlying `run_MCScanX` function both across
#' and within pairwise genome comparisons.
#'
#' @param blast.results R object contain the blast results
#' @param MCScanX.params Parameters to pass to MCScanX
#' @param mcscanx.input.dir Directory containing the MCScanX-formatted mapping files
#' @param abbrevs Genome abbreviations
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
                            MCScanX.params = "-a -s 10 -m 5 -w 2",
                            abbrevs,
                            mcscanx.input.dir,
                            verbose = TRUE,
                            ...){

  find_unmappedBlast <- function(blk,
                                map,
                                blast){

    blk$uniq = with(blk, paste(genome1, genome2, chr1, chr2))
    map$uniq = with(map, paste(genome1, genome2, chr1, chr2))
    blast$uniq = with(blast, paste(genome1, genome2, chr1, chr2))
    map$uniq.map = with(map, paste(genome1, genome2, id1, id2))
    blast$uniq.map = with(blast, paste(genome1, genome2, id1, id2))

    blast.noblk = blast[!blast$uniq.map %in% map$uniq.map,]

    spl.blk = split(blk, blk$uniq)
    blast.noblk = blast.noblk[blast.noblk$uniq %in% unique(map$uniq),]

    spl.blast = split(blast.noblk, factor(blast.noblk$uniq))
    blast.noblk = rbindlist(lapply(names(spl.blast), function(i){
      x.blk = spl.blk[[i]]
      x.blast = spl.blast[[i]]
      if(nrow(x.blast)>1){
        inblk = sapply(1:nrow(x.blk), function(j){
          x.blast$start1 <= x.blk$end1[j] & x.blast$end1 >= x.blk$start1[j] &
            x.blast$start2 <= x.blk$end2[j] & x.blast$end2 >= x.blk$start2[j]
        })
        outblk = x.blast[rowSums(inblk)==0,]
      }else{
        outblk = x.blast
      }
    }))

    allu = unique(c(map$uniq.map, blast.noblk$uniq.map))
    return(blast[blast$uniq.map %in% allu,])
  }

  if(verbose)
    cat("Initial MCScanX run on",nrow(blast.results), "culled BLAST hits... \n\t")
  mcscan.cull = run_MCScanX(
    blast.results = blast.results,
    abbrevs = abbrevs,
    mcscanx.input.dir = mcscanx.input.dir,
    MCScanX.params = MCScanX.params,
    verbose = FALSE)
  if(verbose)
    cat("Clustered",nrow(mcscan.cull),"BLAST hits into",
        length(unique(mcscan.cull$block.id)), "colinear blocks\n")

  blast.results.mcs = merge(blast.results, mcscan.cull[,1:2,with = F])
  if(verbose)
    cat("Secondary MCScanX run ...\n\t")
  mcscan.cull = run_MCScanX(
    blast.results = blast.results.mcs,
    abbrevs = abbrevs,
    mcscanx.input.dir = mcscanx.input.dir,
    MCScanX.params = MCScanX.params,
    verbose = FALSE)
  if(verbose)
    cat("Clustered",nrow(mcscan.cull),"BLAST hits into",
        length(unique(mcscan.cull$block.id)), "colinear blocks\n")

  blast.results.mcs = merge(blast.results, mcscan.cull[,1:2,with = F])
  if(verbose)
    cat("MCScanX run within pairwise genome comparisons...\n\t")

  mcscan.cull = rbindlist(lapply(split(blast.results.mcs, with(blast.results.mcs, paste(genome1, genome2))),
                       function(x){
                         run_MCScanX(
                           blast.results = x,
                           abbrevs = abbrevs,
                           mcscanx.input.dir = mcscanx.input.dir,
                           MCScanX.params = MCScanX.params,
                           verbose = FALSE)
                       }))
  blast.results.mcs = merge(blast.results, mcscan.cull[,1:2,with = F])

  mcscan.list = lapply(split(blast.results.mcs, with(blast.results.mcs, paste(genome1, genome2))),
                       function(x){
                         mctmp = run_MCScanX(
                           blast.results = x,
                           abbrevs = abbrevs,
                           mcscanx.input.dir = mcscanx.input.dir,
                           MCScanX.params = MCScanX.params,
                           verbose = FALSE)
                         return(make_blocks(mctmp))
                       })

  blk = rbindlist(lapply(mcscan.list, function(x) x$block))
  map = rbindlist(lapply(mcscan.list, function(x) x$map))

  if(verbose)
    cat("Clustered",nrow(map),"BLAST hits into",
        nrow(blk), "colinear blocks\n")
  if(verbose)
    cat("Finding BLAST hits in gaps between blocks...\n\t")

  blast.unmap = find_unmappedBlast(map = map, blk = blk, blast = blast.results)
  blast.unmap = blast.unmap[,colnames(blast.results.mcs), with = F]
  if(verbose)
    cat("Added",nrow(blast.unmap)-nrow(map),"BLAST hits\n")

  if(verbose)
    cat("Final MCScanX run within pairwise genome comparisons...\n\n")
  mcscan.cull = rbindlist(lapply(split(blast.unmap, with(blast.unmap, paste(genome1, genome2))),
                                 function(x){
                                   run_MCScanX(
                                     blast.results = x,
                                     abbrevs = abbrevs,
                                     mcscanx.input.dir = mcscanx.input.dir,
                                     MCScanX.params = MCScanX.params,
                                     verbose = FALSE)
                                 }))
  blast.results.mcs = merge(blast.results, mcscan.cull[,1:2,with = F])

  mcscan.list = lapply(split(blast.results.mcs, with(blast.results.mcs, paste(genome1, genome2))),
                       function(x){
                         mctmp = run_MCScanX(
                           blast.results = x,
                           abbrevs = abbrevs,
                           mcscanx.input.dir = mcscanx.input.dir,
                           MCScanX.params = MCScanX.params,
                           verbose = FALSE)
                         return(make_blocks(mctmp))
                       })

  blk = rbindlist(lapply(mcscan.list, function(x) x$block))
  map = rbindlist(lapply(mcscan.list, function(x) x$map))

  if(verbose)
    cat("Done! ... Clustered",nrow(map),"BLAST hits into",
        nrow(blk), "colinear blocks\n")
  return(list(block = blk, map = map))
}
