#' @title selection utility functions
#' @description
#' \code{genespace_utils} Functions that allow selection stat calculation in GENESPACE
#' @name genespace_utils
#'
#' @param add.gff logical, should
#' @param add.metadata logical, should
#' @param alpha numeric, length
#' @param bg.col character, length
#' @param blast map-formatted data.table, without blockID info.
#' @param blast.dir file path to
#' @param blast.file.in file path to
#' @param blast.file.out file path to
#' @param blast.ids character, length
#' @param blk block-formatted data.table
#' @param blk.border numeric, length
#' @param blk.ids character, length
#' @param border.color character, length
#' @param check.ogs logical, should
#' @param chr character, length
#' @param chr.abbrev.fun function to
#' @param chr.bg.cex numeric, length
#' @param chr.bg.col character, length
#' @param chr.bg.pch numeric, length
#' @param chr.buff numeric, length
#' @param chr.buffer numeric, length
#' @param chr.id.cex numeric, length
#' @param chr.id.col character, length
#' @param chr.lab.buff numeric, length
#' @param chr.list list of chromosomes to plot
#' @param chr.segm.col character, length
#' @param chr.segm.lwd numeric, length
#' @param clean.columns logical, should
#' @param col character, length
#' @param cols character, length
#' @param comb genome combinations
#' @param cull.blast.dir file path to
#' @param dir.list list of file paths to
#' @param do.cumulative logical, should
#' @param dodge.geneIDs numeric, length
#' @param dodge.x numeric, length
#' @param drop.NAs logical, should
#' @param e1 numeric, length
#' @param e2 numeric, length
#' @param end numeric, length
#' @param eps.radius numeric, length
#' @param fais fai-like data.table
#' @param fasta.dir file path to
#' @param fill.color character, length
#' @param forCircos logical, should
#' @param gap.prop numeric, length
#' @param gene.colors character, length
#' @param gene.dict1 gene dictionary
#' @param gene.dict2 gene dictionary
#' @param geneID.abbrev.fun function to
#' @param geneid.cex numeric, length
#' @param geneid.offset numeric, length
#' @param genes2plot character, length
#' @param genomeIDs character, length
#' @param genomes character, length
#' @param gff gff-like data.table
#' @param gff.dir file path to
#' @param gff.file file path to
#' @param id.db database of gene IDs / numbers
#' @param is.peptide logical, should
#' @param keep.best.id1.hit logical, should
#' @param keep.geneNum logical, should
#' @param lab.chr logical, should
#' @param lab.chr.1only logical, should
#' @param m.param numeric, length
#' @param map map-formatted data.table
#' @param map.bychr map-formatted data.table split by chr
#' @param mappings numeric, length
#' @param maxn numeric, length
#' @param mcs.file file path to
#' @param mcscan.dir file path to
#' @param mcscan.param character, length
#' @param MCScanX.path file path to
#' @param min.blockSize numeric, length
#' @param min.dist2end numeric, length
#' @param min.perc.iden numeric, length
#' @param min.score numeric, length
#' @param n.cores numeric, length
#' @param n.mapping numeric, length
#' @param n.mappings numeric, length
#' @param n.of.cores numeric, length
#' @param n.out numeric, length
#' @param n.reps numeric, length
#' @param n.sample numeric, length
#' @param num numeric, length
#' @param of.blast blast file from orthofinder
#' @param of.dir file path to
#' @param of.ids character, length
#' @param only.orthogroups logical, should
#' @param ortho.col character, length
#' @param parse_fastaHeader.FUN function to
#' @param pattern character, length
#' @param peptide.dir file path to
#' @param points.per.curve numeric, length
#' @param prop.of.best numeric, length
#' @param pw.of logical, should
#' @param radius numeric, length
#' @param rank.buffer numeric, length
#' @param rename.blocks logical, should
#' @param rerank logical, should
#' @param return.start logical, should
#' @param s1 numeric, length
#' @param s2 numeric, length
#' @param scale.it logical, should
#' @param scale2dodge numeric, length
#' @param silent.mcs logical, should
#' @param simplify.poly logical, should
#' @param start numeric, length
#' @param syn.blast syntenic blast, map format data.table
#' @param syn.ortho.map syntenic orthologous blast, map format data.table
#' @param ties.method character, length
#' @param tmp.dir file path to
#' @param use.rank logical, should
#' @param use.recip logical, should
#' @param use.topn logical, should
#' @param verbose logical, should
#' @param w2ki numeric vector specifying which genome comparisons to keep
#' @param which.in.blk numeric, length
#' @param whichAttr numeric, length
#' @param x numeric or character, length
#' @param y numeric, length
#' @param y.end numeric, length
#' @param y.start numeric, length
#'
#' @note \code{genespace_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespace_utils} returns its own arguments.
#'
#'

#' @title track_blocks
#' @description
#' \code{track_blocks} track blocks
#' @rdname genespace_utils
#' @export
track_blocks <- function(map,
                         blk.ids,
                         cols,
                         bg.col){
  ulist <- unlist(lapply(blk.ids, function(x)
    unique(subset(map, block.id == x)$og.id)))
  if (any(duplicated(ulist)))
    ulist <- ulist[!duplicated(ulist)]

  if (length(cols) != length(blk.ids))
    stop("cols and blk.ids must be of the same length")

  blk.cols <- bg.col
  for (i in 1:length(blk.ids)) {
    bi <- blk.ids[i]
    ci <- cols[i]

    mb <- subset(map, block.id == bi)
    ob <- unique(mb$og.id)

    map$block.id <- with(map, ifelse(og.id %in% ob, paste0(block.id, "_inblk",i,"col"), block.id))
    bo <- rep(ci, length(grep(paste0("_inblk",i,"col"), map$block.id)))
    names(bo) <- map$block.id[grep(paste0("_inblk",i,"col"), map$block.id)]
    bo <- bo[!duplicated(names(bo))]
    blk.cols <- c(blk.cols, bo)
  }
  return(list(map = map,
              blk.cols = blk.cols))
}









#' @title drop_mapDup
#' @description
#' \code{drop_mapDup} drop_mapDup
#' @rdname genespace_utils
#' @export
drop_mapDup <- function(map){
  map[,ns := -score]
  setkey(map, ns)
  out <- map[, head(.SD, 1), by = list(block.id, id1)]
  out$ns <- NULL
  return(out)
}

#' @title extend_blk2chrend
#' @description
#' \code{extend_blk2chrend} extend_blk2chrend
#' @rdname genespace_utils
#' @export
extend_blk2chrend <- function(blk,
                              gff,
                              map,
                              genomeIDs,
                              min.dist2end){

  gl <- gff[,list(min = min(start),
                  max = max(end)),
            by = list(genome, chr)]

  spl.gff <- split(gff, by = "genome")
  spl.map <- split(map, by = c("genome1","genome2"))

  for (x in spl.map) {
    g1 <- subset(spl.gff[[x$genome1[1]]], id %in% x$id1)
    g2 <- subset(spl.gff[[x$genome2[1]]], id %in% x$id2)
    gm1 <- g1[,list(g.min = id[frank(start, ties.method = "random") <= min.dist2end],
                    g.max = id[frank(-end, ties.method = "random") <= min.dist2end]),
              by = list(genome, chr)]
    gm2 <- g2[,list(g.min = id[frank(start, ties.method = "random") <= min.dist2end],
                    g.max = id[frank(-end, ties.method = "random") <= min.dist2end]),
              by = list(genome, chr)]

    gm1 <- merge(gl, gm1, by = c("genome", "chr"))
    gm2 <- merge(gl, gm2, by = c("genome", "chr"))

    for (j in 1:nrow(gm1)) {
      wh.blks <- unique(x$block.id[x$id1 %in% gm1$g.min[j]])
      blk$start1[with(blk, block.id %in% wh.blks)] <- gm1$min[j]

      wh.blks <- unique(x$block.id[x$id1 %in% gm1$g.max[j]])
      blk$end1[with(blk, block.id %in% wh.blks)] <- gm1$max[j]
    }

    for (j in 1:nrow(gm2)) {
      wh.blks <- unique(x$block.id[x$id2 %in% gm2$g.min[j]])
      blk$start2[with(blk, block.id %in% wh.blks)] <- gm2$min[j]

      wh.blks <- unique(x$block.id[x$id2 %in% gm2$g.max[j]])
      blk$end2[with(blk, block.id %in% wh.blks)] <- gm2$max[j]
    }
  }

  return(blk)
}

#' @title pull_duplicatesInChr
#' @description
#' \code{pull_duplicatesInChr} pull_duplicatesInChr
#' @rdname genespace_utils
#' @export
pull_duplicatesInChr <- function(map.bychr,
                                 gff){
  dup.blks <- character()
  b <- map.bychr[,list(start1 = min(start1),
                       end1 = max(end1),
                       start2 = min(start2),
                       end2 = max(end2)),
                 by = list(genome1,genome2,chr1,chr2,block.id)]
  if (nrow(b) > 1) {
    b$block.id <- as.character(b$block.id)

    ovl <- 1
    nb <- nrow(b)
    while (ovl > .1 & nb > 1) {
      genes.byBlk <- lapply(1:nrow(b), function(i){
        subset(gff, genome == b$genome1[i] & chr == b$chr1[i] &
                 start >= b$start1[i] & end <= b$end1[i])$id
      })
      names(genes.byBlk) <- as.character(b$block.id)

      eg <- data.table(expand.grid(as.character(b$block.id),
                                   as.character(b$block.id),
                                   stringsAsFactors = F)[,c(2,1)])
      eg <- subset(eg, Var1 != Var2)

      eg$prop.overlap <- apply(eg, 1, function(x)
        length(intersect(genes.byBlk[[x[1]]],
                         genes.byBlk[[x[2]]]))/length(genes.byBlk[[x[1]]]))
      eg <- eg[order(-eg$prop.overlap),]
      ovl <- eg$prop.overlap[1]
      if (eg$prop.overlap[1] > 0.1) {
        dup.blks <- c(dup.blks, eg$Var2[1])
        b <- b[!b$block.id %in% dup.blks]
        nb <- nrow(b)
      }
    }
  }
  return(dup.blks)
}

#' @title find_blocksToExtend
#' @description
#' \code{find_blocksToExtend} find_blocksToExtend
#' @rdname genespace_utils
#' @export
find_blocksToExtend <- function(blk,
                                chr.end.buffer = 5,
                                ovl.buffer = 3,
                                same.chr.buffer = 500,
                                verbose = T){

  if(verbose)
    cat("Finding adjacent blocks ...")
  adj.blk1 <- blk[,find_adjBlks1(.SD,
                                 ovl.buffer = ovl.buffer,
                                 same.chr.buffer = same.chr.buffer),
                  by = list(genome1, genome2, chr1)]
  adj.blk2 <- blk[,find_adjBlks2(.SD,
                                 ovl.buffer = ovl.buffer,
                                 same.chr.buffer = same.chr.buffer),
                  by = list(genome1, genome2, chr2)]
  adjb <- merge(adj.blk1, adj.blk2, by = c("genome1","genome2","block.id"))
  adjb <- subset(adjb, block.id != blk.left | block.id != blk.right |
                   block.id != blk.down | block.id != blk.up)
  rs <- rowSums(is.na(adjb[,c("blk.left","blk.right","blk.up","blk.down")]))
  blks.adj <- adjb[rs <= 3,]
  if(verbose)
    cat("Done\n")

  m <- melt(blks.adj, id.vars = c("genome1","genome2","chr1","chr2"),
            measure.vars = c("block.id","blk.left","blk.right","blk.up","blk.down"))
  m[,variable:= NULL]
  m <- m[!duplicated(m),]
  spl.adj <- split(m, by = c("genome1","genome2","chr1","chr2"))
  blk[,chr.rankend1 := max(rankend1),
      by = list(genome1, genome2, chr1)]
  blk[,chr.rankend2 := max(rankend2),
      by = list(genome1, genome2, chr2)]
  blk[,n.inchr := .N,
      by = list(genome1, genome2, chr1, chr2)]
  blks.solo <- subset(blk, n.inchr == 1 &
                        rankstart1 > chr.end.buffer &
                        (rankend1 - chr.end.buffer) < chr.rankend1 &
                        rankstart2 > chr.end.buffer &
                        (rankend2 - chr.end.buffer) < chr.rankend2)

  blks2test <- unique(c(blks.solo$block.id, blks.adj$block.id))
  return(blks2test)
}

#' @title find_leftBlk
#' @description
#' \code{find_leftBlk} find_leftBlk
#' @rdname genespace_utils
#' @export
find_leftBlk <- function(st.rank, end.ranks, chrs, chr, ovl.buffer, same.chr.buffer){

  is.left.onchr <- st.rank > (end.ranks - ovl.buffer) & chrs == chr
  is.left.offchr <- st.rank > (end.ranks - ovl.buffer) & chrs != chr
  if(sum(is.left.onchr | is.left.offchr) == 0){
    out <- NA
  }else{
    if(sum(is.left.offchr) == 0 | sum(is.left.onchr) == 0){
      ens <- end.ranks[is.left.onchr | is.left.offchr]
      out <- which(end.ranks == ens[which.min(st.rank - ens)])[1]
    }else{
      min.dist.onchr <- min(st.rank - end.ranks[is.left.onchr])
      min.dist.offchr <- min(st.rank - end.ranks[is.left.offchr])
      if(min.dist.onchr < min.dist.offchr + same.chr.buffer){
        ens <- end.ranks[is.left.onchr]
        out <- which(end.ranks == ens[which.min(st.rank - ens)] &
                       is.left.onchr)[1]
      }else{
        ens <- end.ranks[is.left.offchr]
        out <- which(end.ranks == ens[which.min(st.rank - ens)] &
                       is.left.offchr)[1]
      }
    }
  }
  return(out)
}

#' @title find_rightBlk
#' @description
#' \code{find_rightBlk} find_rightBlk
#' @rdname genespace_utils
#' @export
find_rightBlk <- function(st.ranks, end.rank, chrs, chr, ovl.buffer, same.chr.buffer){

  is.right.onchr <- st.ranks > (end.rank - ovl.buffer) & chrs == chr
  is.right.offchr <- st.ranks > (end.rank - ovl.buffer) & chrs != chr
  if(sum(is.right.onchr | is.right.offchr) == 0){
    out <- NA
  }else{
    if(sum(is.right.offchr) == 0 | sum(is.right.onchr) == 0){
      sts <- st.ranks[is.right.onchr | is.right.offchr]
      out <- which(st.ranks == sts[which.min(sts - end.rank)])[1]
    }else{
      min.dist.onchr <- min(st.ranks[is.right.onchr] - end.rank)
      min.dist.offchr <- min(st.ranks[is.right.offchr] - end.rank)
      if(min.dist.onchr < min.dist.offchr + same.chr.buffer){
        sts <- st.ranks[is.right.onchr]
        out <-  which(st.ranks == sts[which.min(sts - end.rank)] &
                        is.right.onchr)[1]
      }else{
        sts <- st.ranks[is.right.offchr]
        out <-  which(st.ranks == sts[which.min(sts - end.rank)] &
                        is.right.offchr)[1]
      }
    }
  }
  return(out)
}

#' @title find_adjBlks1
#' @description
#' \code{find_adjBlks1} find_adjBlks1
#' @rdname genespace_utils
#' @export
find_adjBlks1 <- function(b, ovl.buffer, same.chr.buffer){
  which.left <- sapply(1:nrow(b), function(i)
    find_leftBlk(st.rank = b$rankstart1[i], end.ranks = b$rankend1,
                 chr = b$chr2[i], chrs = b$chr2,
                 ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  which.right <- sapply(1:nrow(b), function(i)
    find_rightBlk(st.ranks = b$rankstart1, end.rank = b$rankend1[i],
                  chr = b$chr2[i], chrs = b$chr2,
                  ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))

  return(data.table(block.id = b$block.id,
                    start1 = b$start1,
                    end1 = b$end1,
                    blk.left = b$block.id[which.left],
                    blk.right = b$block.id[which.right],
                    blk.left.end = b$end1[which.left],
                    blk.right.start = b$start1[which.right]))
}

#' @title find_adjBlks2
#' @description
#' \code{find_adjBlks2} find_adjBlks2
#' @rdname genespace_utils
#' @export
find_adjBlks2 <- function(b, ovl.buffer, same.chr.buffer){
  which.left <- sapply(1:nrow(b), function(i)
    find_leftBlk(st.rank = b$rankstart2[i], end.ranks = b$rankend2,
                 chr = b$chr1[i], chrs = b$chr1,
                 ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  which.right <- sapply(1:nrow(b), function(i)
    find_rightBlk(st.ranks = b$rankstart2, end.rank = b$rankend2[i],
                  chr = b$chr1[i], chrs = b$chr1,
                  ovl.buffer = ovl.buffer, same.chr.buffer = same.chr.buffer))
  return(data.table(block.id = b$block.id,
                    start2 = b$start2,
                    end2 = b$end2,
                    blk.down = b$block.id[which.left],
                    blk.up = b$block.id[which.right],
                    blk.down.end = b$end1[which.left],
                    blk.up.start = b$start1[which.right]))
}

#' @title add_orthology
#' @description
#' \code{add_orthology} add_orthology
#' @rdname genespace_utils
#' @export
add_orthology <- function(of.dir,
                          map,
                          verbose,
                          block.orthology.threshold){
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

  map.out[,prop.og := sum(is.ortholog) / .N,by = "block.id"]
  map.out[,block.is.ortholog := prop.og > block.orthology.threshold]

  return(map.out)
}
