#' @title Extend and complete syntenic mappings
#'
#' @description
#' \code{extend_blocks} Extend and complete syntenic mappings.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL default, taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param use.score.cull.blast logical, should the blasts be from the
#' score.cull directory, or just the genome.cull directory?
#' @param ... Not currently in use
#' @param rank.buffer numeric, the radius to search within for syntenic
#' mappings.
#'
#' @details ...
#'
#' @return A blast data.table, without block.id information
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
extend_blocks <- function(gff,
                          genomeIDs,
                          dir.list,
                          map,
                          verbose = T,
                          MCScanX.s.param = 5,
                          MCScanX.m.param = 50,
                          MCScanX.path,
                          use.score.cull.blast = T,
                          return.reg.only = F,
                          final.search = F,
                          rank.buffer){

  if (use.score.cull.blast) {
    of.dir <- dir.list$cull.score.blast
  }else{
    of.dir <- dir.list$cull.blast
  }

  gff <- data.table(gff[,c("id","chr","start","end","strand","genome","order")])

  if(verbose)
    cat("Importing orthofinder database and blast results ...")
  sids <- read_speciesIDs(
    of.dir = of.dir,
    genomeIDs = genomeIDs)
  gids <- read_geneIDs(
    of.dir = of.dir,
    species.num.id = sids,
    gff = gff)


  blast <- import_allBlasts(
    orthofinder.dir = of.dir,
    gff = gff,
    genomeIDs = genomeIDs)

  bcols <- c("genome1","genome2","id1","id2")
  allblastsyn <- blast[!duplicated(blast[,bcols,with =F]),]
  map <- merge(
    allblastsyn,
    map[!duplicated(map[,bcols,with  =F]),bcols,with =F],
    by = bcols)

  if(verbose)
    cat("Done!\nSubsetting blast hits to chromosomes found in syntenic blocks ... ")
  chrs <- map[,c("genome1","genome2","chr1","chr2")]
  chrs <- chrs[!duplicated(chrs),]
  blast.more <- merge(chrs, blast, by = colnames(chrs))

  all.map <- merge(map[,c("genome1","genome2","id1","id2")],
                   blast.more,
                   by = c("genome1","genome2","id1","id2"))

  all.blast <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = mirror_map(blast.more),
    intergenome.only = F)

  all.map <- reduce_recipBlast(
    genomeIDs = genomeIDs,
    blast = mirror_map(all.map),
    intergenome.only = F)

  all.blast[,what := "blast"]
  all.map[,what := "map"]

  cn <- colnames(all.blast)[colnames(all.blast) %in% colnames(all.map)]

  z <- rbind(all.map[, cn, with = F],
             all.blast[, cn, with = F])
  z <- z[!duplicated(z[,c("genome1","genome2","id1","id2")]),]

  if(verbose)
    cat("Done!\nReranking blast hits by position in gff ... ")
  zid <- merge(gff,
               with(z,  data.table(
                 genome = c(genome1, genome2),
                 id = c(id1, id2))),
               by = c("genome","id"))
  zid <- zid[!duplicated(zid),]
  setkey(zid, genome, chr, start, end)
  zid[,rank := frank(start, ties.method = "dense"),
      by = list(genome, chr)]
  zid1 <- with(zid,
               data.table(genome1 = genome,
                          id1 = id,
                          rank1 = rank))
  zid2 <- with(zid,
               data.table(genome2 = genome,
                          id2 = id,
                          rank2 = rank))
  zo <- merge(zid1,
              merge(zid2,
                    z,
                    by = c("genome2","id2")),
              by = c("genome1","id1"))
  spl.map <- split(zo, by = c("genome1","genome2","chr1","chr2"))

  if(verbose)
    cat("Done!\nSearching for syntenic hits within", length(spl.map),"pairwise chromosome combinations ...\n")

  out <- rbindlist(lapply(names(spl.map), function(i){
    if(which(names(spl.map) == i) %% 100 == 0 & verbose)
      cat("\tCompleted",which(names(spl.map) == i),"/",length(spl.map),"\n")
    z <- spl.map[[i]]
    wh <- find_whichInBuffer(x = z$rank1,
                             y = z$rank2,
                             which.in.blk = which(z$what == "map"),
                             rank.buffer = rank.buffer)
    return(data.table(z[wh,]))
  }))
  if(verbose)
    cat("\tDone!\nPulling collinear hits ... ")
  out <- merge(blast.more,
               out[,c("genome1","genome2","id1","id2")],
               by = c("genome1","genome2","id1","id2"))

  if(!return.reg.only){
    syntenic.blast <- run_MCScanX(
      blast = mirror_map(out),
      gff = gff,
      mcscan.dir = dir.list$mcscanx,
      overwrite.output.dir = T,
      genomeIDs = genomeIDs,
      MCScanX.path = MCScanX.path,
      MCScanX.s.param = MCScanX.s.param,
      MCScanX.m.param = MCScanX.m.param,
      verbose = F)

    if(verbose)
      cat("Done!\nCompleting subgraphs ... ")

    comp2 <- complete_graph(
      map = syntenic.blast,
      gff = gff,
      ignore.self = T,
      verbose = F)

    if(final.search){
      if(verbose)
        cat("Done!\nRe-Searching for syntenic hits near collinear blocks  ... ")
      zo <- merge(comp2[,c("genome1","genome2","id1","id2")],
                  zo,
                  by = c("genome1","genome2","id1","id2"))
      spl.map <- split(zo, by = c("genome1","genome2","chr1","chr2"))
      out <- rbindlist(lapply(names(spl.map), function(i){
        z <- spl.map[[i]]
        wh <- find_whichInBuffer(x = z$rank1,
                                 y = z$rank2,
                                 which.in.blk = which(z$what == "map"),
                                 rank.buffer = rank.buffer)
        return(data.table(z[wh,]))
      }))

      bl.iden <- subset(blast.more[,c("genome1","genome2","id1","id2")],
                        genome1 == genome2 & id1 == id2)
      out2 <- merge(blast.more,
                    rbind(out[,c("genome1","genome2","id1","id2")],
                          bl.iden),
                    by = c("genome1","genome2","id1","id2"))
      if(verbose)
        cat("Done!\n")
      out <- data.table(out2)
    }else{
      out <- data.table(comp2)
    }

  }

  return(out)
}
