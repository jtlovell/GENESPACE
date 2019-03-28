#' @title Pseudogene utility functions
#' @description
#' \code{psg_utils} Five utilities functions meant for internal calls in compareGeneSpace
#' @name syn_utils
#'
#' @param map map results data.table
#' @param assembly.dir path to assembly fastas
#' @param tmp.dir path to temp directory
#' @param peptide.dir path to peptide directory
#' @param buffer numeric, the number of basepairs outside of the range to look at
#' @param genomeIDs character, indicating genomeIDs to consider.
#' @param clean logical, should the intermediate files be removed?
#' @param max.window.bp Numeric, the maximum size of region to blast
#' @param gff.spl a list of gff objects, split by genome and chr
#' @param diamond.sensitive logical, should diamond blastx be run in the
#' sensitive mode?
#'
#' @param verbose logical, should updates be printed?
#' @param ... not currently in use
#'
#' @note \code{syn_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{syn_utils} returns its own arguments.
#'
#' @title cull intermediate map file
#' @description
#' \code{cull_gpMap} cull intermediate map file
#' @rdname syn_utils
#' @import data.table
#' @export
cull_gpMap <- function(gpmap,
                       alt.map){
  gpmap$id.pass<- with(gpmap, score.pass & !is.na(score.pass) & is.na(id2))
  inmap <- alt.map[,c("id1","id2","og","block.id",
                      "array.id","array.key","score")]
  inmap$in.map <- TRUE
  mapout2 <- merge(gpmap, inmap, by = c("id1","id2"),
                   all = T)
  mapout2$in.map[is.na(mapout2$in.map)]<-F
  mapout2$no.hit <- with(mapout2, id.pass & score.pass & n.pass)

  mapout3 <- mapout2[with(mapout2, no.hit | !is.na(id2)),
                     c("genome1","genome2",
                       "og","block.id","array.id","array.key",
                       "id1","id2",
                       "chr1","chr2",
                       "start1","end1",
                       "start2","end2",
                       "score","n","no.hit")]
  return(mapout3)
}

#' @title run dbs
#' @description
#' \code{rdbs} run dbs clustering without dropping any hits
#' @rdname syn_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
rdbs <- function(x, eps){
  nn <- frNN(data.frame(x,x),
             eps = eps)
  return(dbscan(nn, minPts = 0)$cluster)
}

#' @title find genomic positions
#' @description
#' \code{find_genomicPos} find genomic positions of unmapped blast hits
#' @rdname syn_utils
#' @import data.table
#' @export
find_genomicPos <- function(gff1,
                            gff2,
                            alt.map,
                            n.cores = 1,
                            wind.size = 2,
                            rel.sup.thresh = 1,
                            abs.sup.thresh = 0.75,
                            n.thresh = 1){
  splg <- split(gff1, "chr")
  map2gen <- rbindlist(lapply(splg, function(x){
    xo = x[,c("genome","id","chr","start","end","map.rank")]
    setnames(xo, paste0(colnames(xo),"1"))

    xw <- make_genomeWindow(x,
                            wind.size = wind.size)

    spl = split(xw$id,
                xw$id2)
    sup.pos <- rbindlist(mclapply(names(spl), mc.cores = n.cores, function(i)
      pull_altGenomePos(id1 = i,
                        id2 = spl[[i]],
                        alt.gff = gff2,
                        alt.map = alt.map)))
    sup.out <- merge(xo,
                     sup.pos,
                     by = "id1")
    return(sup.out)
  }))

  map2gen$score.pass <- with(map2gen,
                             (abs.sup >= abs.sup.thresh &
                                rel.sup >= abs.sup.thresh) |
                               rel.sup >= rel.sup.thresh)
  map2gen$n.pass <- with(map2gen, n > n.thresh)

  map2gen$pos2.near10Mb <- round(map2gen$start2,-7)
  map2gen$pos2.near100rank <- with(map2gen,
                                   round((rankend2+rankstart2)/2,-2))

  return(map2gen)
}

#' @title convert a multi-genome gff into a list of two
#' @description
#' \code{make_2gffs} convert a multi-genome gff into a list of two
#' @rdname syn_utils
#' @import data.table
#' @export
make_2gffs <- function(gff,
                       genome1,
                       genome2,
                       genes1.inarray,
                       genes1.noarray,
                       genes2.inarray,
                       genes2.noarray){

  # -- Subset to focal genomes
  g1 <- gff[gff$genome == genome1,]
  g2 <- gff[gff$genome == genome2,]

  # -- add in annotation about location of arrays
  g1$in.map <- g1$id %in% genes1.noarray
  g2$in.map <- g2$id %in% genes2.noarray
  g1$in.ta <- g1$id %in% genes1.inarray
  g2$in.ta <- g2$id %in% genes2.inarray

  # -- count number of genes by chromsome, cull ones without any hits
  g1[,n.inmap := sum(in.map), by = list(chr)]
  g2[,n.inmap := sum(in.map), by = list(chr)]
  g1 <- g1[g1$n.inmap > 0,]
  g2 <- g2[g2$n.inmap > 0,]
  g1$n.inmap <- NULL
  g2$n.inmap <- NULL

  # -- get mean position of each gene
  g1$mean.pos <- with(g1, (start+end)/2)
  g2$mean.pos <- with(g2, (start+end)/2)

  # -- re-order within chromosome
  g1$order <- NULL
  g2$order <- NULL
  g1[,order := frank(mean.pos, ties.method = "random"), by = list(chr)]
  g2[,order := frank(mean.pos, ties.method = "random"), by = list(chr)]

  # -- make ranks by cumsum of non-array gene hits
  g1$torank = g1$in.map
  g2$torank = g2$in.map
  g1$torank[with(g1, !torank & order == 1)]<-TRUE
  g2$torank[with(g2, !torank & order == 1)]<-TRUE
  setkey(g1, chr, mean.pos)
  setkey(g2, chr, mean.pos)
  g1[,map.rank := cumsum(torank), by = list(chr)]
  g2[,map.rank := cumsum(torank), by = list(chr)]
  g1$torank = NULL
  g2$torank = NULL

  return(list(gff1 = g1, gff2 = g2))
}

#' @title pull alternative genome position
#' @description
#' \code{pull_altGenomePos} pull alternative genome position of a given set of genes
#' @rdname syn_utils
#' @import data.table
#' @export
pull_altGenomePos <- function(id1,
                              id2,
                              alt.gff,
                              alt.map,
                              clus.rank.dist =  100){
  idt = id2[!is.na(id2)]
  zid = alt.map$id2[alt.map$id1 %in% idt]
  z = alt.gff[alt.gff$id %in% zid,]
  if(length(zid) == 0){
    return(data.table(clus = NA,
                      genome2 = alt.gff$genome[1],
                      chr2 = NA,
                      start2 = NA,
                      end2 = NA,
                      rankstart1 = NA,
                      rankend2 = NA,
                      n = NA,
                      rel.sup = NA,
                      abs.sup = NA,
                      id1 = id1))
  }else{
    if(length(zid) == 1){
      z$clus <- 1
    }else{
      z[,clus := rdbs(x = map.rank, eps = clus.rank.dist),
        by = list(chr)]
      z$clus <- with(z, as.numeric(as.factor(paste(chr, clus))))
    }
    z[, n := .N, by = list(clus)]
    z[,rel.sup := n/max(n)]
    z[,abs.sup := n/(length(unique(idt)))]
    o <- z[,list(genome2 = genome[1],
                 chr2 = chr[1],
                 start2 = min(start),
                 end2 = max(end),
                 rankstart2 = min(map.rank),
                 rankend2 = max(map.rank),
                 n = n[1],
                 rel.sup = rel.sup[1],
                 abs.sup = abs.sup[1]),
           by = list(clus)]
    o$id1 <- id1
    return(o)
  }
}

#' @title add region around focal genes
#' @description
#' \code{make_genomeWindow} add region around focal genes
#' @rdname syn_utils
#' @import data.table
#' @export
make_genomeWindow <- function(x,
                              wind.size = 5){

  idv <- c(
    sapply(wind.size:1, function(i)
      c(rep(NA, i), x$id[1:(length(x$id)-i)])),
    x$id,
    sapply(1:wind.size, function(i)
      c(x$id[(i+1):(length(x$id))],rep(NA, i)))
  )

  # -- make a set of left offset ranks
  xn = rbindlist(lapply(wind.size:0, function(i){
    y = x
    y$rankn = y$map.rank-i
    return(y)
  }))

  # -- make a set of right offset ranks
  xp = rbindlist(lapply(1:wind.size, function(i){
    y = x
    y$rankn = y$map.rank+i
    return(y)
  }))

  # -- combine and return
  xc <- rbind(xn, xp)
  xc$id2 <- idv
  xc$true.id <- xc$id == idv

  # xc <- xc[!is.na(xc$true.id),]
  return(xc)
}

#' @title pull map for two genomes
#' @description
#' \code{pull_map} pull map for two genomes
#' @rdname syn_utils
#' @import data.table
#' @export
pull_map <- function(map,
                     genome1,
                     genome2){

  m <- map[,c("genome1","genome2",
              "id1","id2",
              "og1","block.id","array.id","score",
              "start1","end1","start2","end2")]
  setnames(m, "og1","og")
  m2 <- m[,c(2,1,4,3,5:8,11:12,9:10)]
  setnames(m2, colnames(m))
  m <- rbind(m, m2)
  m <- m[!duplicated(m),]
  wh = which(m$genome1 == genome1 &
               m$genome2 == genome2)
  mt <- data.table(m[wh,])
  mt[,array.key := frank(-score, ties.method = "random") == 1,
     list(array.id)]
  return(mt)
}

#' @title pipeline to track hits by gff method
#' @description
#' \code{pipe_track} pipeline to track hits by gff method
#' @rdname syn_utils
#' @import data.table
#' @export
pipe_track <- function(map.ta,
                       gff,
                       cull2map = F,
                       genome1,
                       genome2,
                       n.cores = 1,
                       wind.size = 2,
                       verbose = T){
  if(verbose)
    cat(paste0("\t", genome1," --> ",genome2, ": "))
  # -- build out map files.
  mt <- pull_map(map = map.ta,
                 genome1 = genome1,
                 genome2 = genome2)

  mo = mt[with(mt, is.na(array.id) | array.key),]
  ma = mt[with(mt, !(is.na(array.id) | array.key)),]
  ma1 <- unique(ma$id1)
  ma2 <- unique(ma$id2)
  mo1 <- unique(mo$id1)
  mo2 <- unique(mo$id2)

  if(cull2map){
    gff <- gff[gff$id %in% unique(c(map.ta$id1, map.ta$id2))]
  }
  # -- make gffs for each genome with non-array ranks
  gff.list <- make_2gffs(gff = gff,
                         genome1 = genome1,
                         genome2 = genome2,
                         genes1.inarray = ma1,
                         genes2.inarray = ma2,
                         genes1.noarray = mo1,
                         genes2.noarray = mo2)
  if(verbose)
    cat("n.genes =", nrow(gff.list$gff1), " ... ")

  # -- Pull the inferred best hit position for each gene
  map2gen <- find_genomicPos(gff1 = gff.list$gff1,
                             gff2 = gff.list$gff2,
                             alt.map = mo,
                             n.cores = n.cores,
                             wind.size = wind.size)
  # -- Merge initial hits with new g2 ranks
  mapout <- merge_gpMap(gff2 = gff.list$gff2,
                        alt.map = mt,
                        gpmap = map2gen)

  gout <- cull_gpMap(alt.map = mt,
                     gpmap = mapout)
  if(verbose)
    cat("n. regions / n. blast =",
        sum(gout$no.hit), "/",
        sum(!gout$no.hit), "\n")
  return(gout)
}

#' @title merge_gpMap
#' @description
#' \code{merge_gpMap} merge_gpMap
#' @rdname syn_utils
#' @import data.table
#' @export
merge_gpMap <- function(gff2,
                        alt.map,
                        gpmap){
  g2m <- merge(alt.map[,c("id1","id2")],
               with(gff2,
                    data.table(id2 = id,
                               chr2 = chr,
                               pos2.near10Mb = round(start,-7),
                               pos2.near100rank = round(map.rank,-2),
                               map.rank2 = map.rank)),
               by = "id2")
  g2m <- g2m[!duplicated(g2m),]
  map3gen1 <- merge(gpmap[complete.cases(gpmap),],
                    g2m[,-5,with = F],
                    by = c("id1", "chr2", "pos2.near10Mb"),
                    all = T)
  map3gen2 <- merge(gpmap[complete.cases(gpmap),],
                    g2m[,-4,with = F],
                    by = c("id1", "chr2", "pos2.near100rank"),
                    all = T)
  mns <- intersect(colnames(map3gen1), colnames(map3gen2))
  mapout <- rbind(map3gen1[,mns, with = F],
                  map3gen2[,mns, with = F])
  mapout <- mapout[!duplicated(mapout),]
  return(mapout)
}

#' @title make genome window
#' @description
#' \code{make_genomeWindow} make genome window for map method
#' @rdname syn_utils
#' @import data.table
#' @export
make_genomeWindow2 <- function(map,
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
     by = list(id1, block.id)]

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

#' @title link regions for map method
#' @description
#' \code{link_regions} link regions for map method
#' @rdname syn_utils
#' @import data.table
#' @export
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

  mi1 = genome.window[,c("genome1","id1","chr1","start1","end1")]
  mi2 = genome.window[,c("genome2","id2","chr2","start2","end2")]
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

#' @title summarize_mapByArray
#' @description
#' \code{summarize_mapByArray} summarize_mapByArray
#' @rdname syn_utils
#' @import data.table
#' @export
summarize_mapByArray <- function(map, verbose = TRUE){

  m <- map[!is.na(map$array.id),]
  mo <- map[is.na(map$array.id),]
  adb <- with(m, data.table(
    genome1 = c(genome1, genome1),
    genome2 = c(genome2, genome2),
    id = c(id1, id2),
    block.id = c(block.id, block.id),
    og = c(og1, og2),
    chr = c(chr1, chr2),
    start = c(start1, start2),
    end = c(end1, end2),
    score = c(score, score)))
  if(verbose)
    cat("n. initial hits =", nrow(adb))
  adb <- data.table(adb[!duplicated(adb),])
  setkey(adb, score)
  adb <- adb[,tail(.SD, 1),
             by = list(genome1, genome2, id, block.id, og, chr, start, end)]
  adb[, pos.rank := frank((start+end)/2, ties.method = "dense"),
      by = list(genome1, genome2, block.id, og)]
  adb[, array.n := .N, by = list(genome1, genome2, block.id, og)]
  adb[, dist2mid := abs(ceiling(.N/2)-pos.rank),
      by = list(genome1, genome2, block.id, og)]
  adb[, score.rank := frank(-score, ties.method = "dense"),
      by = list(genome1, genome2, block.id, og)]
  adb[, score.dist.rank := frank(((score.rank*2) + dist2mid), ties.method = "random"),
      by = list(genome1, genome2, block.id, og)]
  setkey(adb, genome1, genome2, chr, block.id, score.dist.rank)
  ai <- adb[,head(.SD, 1),
            by = list(genome1, genome2, block.id, og)]
  if(verbose)
    cat(" -->", nrow(ai), "representative hits")

  ai1 <- with(ai,
              data.table(genome1 = genome1,
                         genome2 = genome2,
                         block.id = block.id,
                         og1 = og,
                         id1 = id))
  ai2 <- with(ai,
              data.table(genome1 = genome1,
                         genome2 = genome2,
                         block.id = block.id,
                         og2 = og,
                         id2 = id))


  out1 <- merge(ai1, m, by = colnames(ai1))
  setkey(out1, score)
  out1 <- out1[,tail(.SD, 1),
               by = list(genome1, genome2, id1, block.id, og1)]

  out2 <- merge(ai2, m, by = colnames(ai2))
  setkey(out2, score)
  out2 <- out2[,tail(.SD, 1),
               by = list(genome1, genome2, id2, block.id, og2)]

  return(list(condensed = rbind(out1, out2, mo),
              array.db = adb))
}

#' @title track_hits
#' @description
#' \code{track_hits} track_hits
#' @rdname syn_utils
#' @import data.table
#' @export
track_hits <- function(map, gff.spl, max.window.bp = 2e5, verbose){
  spl <- split(map, by = "block.id")
  evy <- round(length(spl),-2)/10
  evy <- min(evy, 10)
  if(verbose)
    cat("Tracking unmapped hits across", length(spl), "blocks ...\n")
  out.all <- rbindlist(lapply(1:length(spl), function(k){
    if(verbose)
      if(k %% evy == 0)
        cat("\tCompleted",k,"/",length(spl),"\n")
    mblk = spl[[k]]
    gff1 <- gff.spl[[mblk$genome1[1]]][[mblk$chr1[[1]]]]
    gff1 <- gff1[gff1$end >= min(mblk$start1) & gff1$start <= max(mblk$end1),]

    gff2 <- gff.spl[[mblk$genome2[1]]][[mblk$chr2[[1]]]]
    gff2 <- gff2[gff2$end >= min(mblk$start2) & gff2$start <= max(mblk$end2),]

    gff1$inblk = gff1$id %in% mblk$id1
    gff2$inblk = gff2$id %in% mblk$id2

    setkey(gff1, chr, start, end, id)
    setkey(gff2, chr, start, end, id)

    gff1$rankblk.start <- cumsum(gff1$inblk)
    gff2$rankblk.start <- cumsum(gff2$inblk)
    gff1$rankblk.end <- gff1$rankblk.start+1
    gff2$rankblk.end <- gff2$rankblk.start+1

    find1 <- gff1[!gff1$inblk,]
    look1 <- gff1[gff1$inblk,]
    find2 <- gff2[!gff2$inblk,]
    look2 <- gff2[gff2$inblk,]

    if(nrow(look1) <= 1 | nrow(find1) < 1 ){
      out1 <- NULL
    }else{
      findlist1 <- rbindlist(apply(find1[,c("rankblk.start","rankblk.end")], 1, function(x){
        ids <- look1$id[with(look1, rankblk.start %in% x)]
        ids2 <- mblk$id2[mblk$id1 %in% ids]
        l2 <- look2[look2$id %in% ids2,]
        return(with(l2, data.table(look.genome = genome[1],
                                   look.chr = chr[1],
                                   look.start = min(start),
                                   look.end = max(end))))
      }))
      out1 <- cbind(find1, findlist1)
    }

    if(nrow(look2) <= 1 | nrow(find2) < 1 ){
      out2 <- NULL
    }else{
      findlist2 <- rbindlist(apply(find2[,c("rankblk.start","rankblk.end")], 1, function(x){
        ids <- look2$id[with(look2, rankblk.start %in% x)]
        ids1 <- mblk$id1[mblk$id2 %in% ids]
        l1 <- look1[look1$id %in% ids1,]
        return(with(l1, data.table(look.genome = genome[1],
                                   look.chr = chr[1],
                                   look.start = min(start),
                                   look.end = max(end))))
      }))
      out2 <- cbind(find2, findlist2)
    }

    if(is.null(out1) & is.null(out2)){
      out <- NULL
    }else{
      out <- rbind(out1, out2)

      out$block.id <- mblk$block.id[1]
      out$look.width = with(out, abs(look.start - look.end))
    }
    return(out)
  }))

  out.small <- out.all[out.all$look.width <= max.window.bp,]
  if(verbose)
    cat("Done! ... Found",nrow(out.all),"total hits,",
        nrow(out.small),"are in windows smaller than",max.window.bp,"bp\n")
  return(out.small)
}

#' @title track_synHits
#' @description
#' \code{track_synHits} track_synHits
#' @rdname syn_utils
#' @import data.table
#' @export
track_synHits <- function(map,
                          gff,
                          verbose = TRUE,
                          max.window.bp = 1e5){

  map[ , is.sameReg := any(id1 %in% id2) & all(genome1 %in% genome2),
       by = list(block.id)]
  map <- map[!map$is.sameReg,]
  if(verbose)
    cat("Condensing arrays ... ")
  sum.array <- summarize_mapByArray(map = map, verbose = verbose)
  map <- sum.array$condensed
  if(verbose)
    cat(" ... Done!\n")

  if(verbose)
    cat("Subsetting and indexing gff to hits in map ... ")
  all.genes <- unique(unlist(map[,c("id1","id2")]))

  gff.inblk <- gff[gff$id %in% all.genes,]
  gff.inblk[, rank := frank((start+end)/2, ties.method = "random"),
            by = list(genome, chr)]
  gff.spl <- lapply(split(gff.inblk, by = "genome"), function(x) split(x, by = "chr"))
  if(verbose)
    cat("Done!\n")

  trk <- track_hits(map = map,
                    max.window.bp = max.window.bp,
                    gff.spl = gff.spl,
                    verbose = T)

  allmap <- map[,c("genome1","genome2","block.id","id1","id2")]
  allmap$og <- map$og1
  gff1 <- data.table(gff.inblk[,c("genome","id","chr","start","end","strand","rank")])
  gff2 <- data.table(gff1)
  setnames(gff1, paste0(colnames(gff1),"1"))
  setnames(gff2, paste0(colnames(gff2),"2"))

  allmap <- merge(gff1, merge(gff2, allmap,
                              by = c("genome2","id2")),
                  by = c("genome1","id1"))

  nomap <- with(trk,
                data.table(genome1 = genome,
                           id1 = id,
                           chr1 = chr,
                           start1 = start,
                           end1 = end,
                           strand1 = strand,
                           rank = rank,
                           genome2 = look.genome,
                           chr2 = look.chr,
                           start2 = look.start,
                           end2 = look.end,
                           block.id = block.id))
  am.og <- with(allmap,
                data.table(id1 = c(id1, id2),
                           og = c(og, og)))
  am.og <- am.og[!duplicated(am.og$id1),]
  nomap <- merge(nomap,
                 am.og,
                 by = "id1")
  nomap[,id2:=paste(genome1, genome2, block.id, frank(start1, ties.method = "random"), sep = "."),
        by = list(genome1, genome2, block.id)]


  return(list(allmap = allmap,
              nomap = nomap,
              array.db = sum.array$array.db))
}
