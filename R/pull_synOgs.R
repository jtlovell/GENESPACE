#' @title Post-synteny orthofinder runs
#' @description
#' \code{pull_synOgs} Functions to run synteny-constrained orthofinder
#' @name pull_synOgs
#'
#' @param path2orthofinder single file path, pointing to the orthofinder
#' executable. If orthofinder is in the path (i.e. if R was launched within a
#' conda env.) then the default 'orthofinder' is sufficient.
#' @param gsAnnot named character vector of file.paths to gff-like annotation
#' files. The names in this vector must match genome1/genome2 columns in the
#' syntenyParams data.table.
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run) and ploidy
#' (integer vector of genome ploidies named by genomeIDs).
#' @param synParam file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param hits data.table containing annotated blast-format pairwise hits
#' @param blkSize integer of length 1 specifying the minimum size for a syntenic
#' block and the -s 'size' MCScanX parameter
#' @param nGaps integer, number of gaps allowed within an MCScanX block
#' @param path2mcscanx character string coercible to a file.path.
#' @param nCores integer of length 1 specifying the number of parallel processes
#' to run
#'
#' @details ...


#' @title Extract and save syntenic orthogroups from pairwise blast hits
#' @description
#' \code{pull_synOgs} Main pipeline to go from pairwise blast hits to an
#' integrated database of syntenic orthology.
#' @rdname pull_synOgs
#' @import data.table
#' @importFrom parallel mclapply
#' @export
pull_synOgs <- function(gsParam,
                        gsAnnot,
                        synParam,
                        dropInterleavesSmallerThan = 2,
                        overwrite = FALSE,
                        inferMissingPos = TRUE){

  blkFiles <- with(synParam, file.path(
    gsParam$synteny,
    sprintf("%s_vs_%s_synog.txt.gz", genome1, genome2)))
  if(all(file.exists(blkFiles)) & !overwrite){
    warning("allSynOgBlks.txt.gz and synOgClus.txt.gz exist in /synteny\n\tSet overwrite = TRUE to run anyways.")
    return(blkFiles)
  }else{
    if(gsParam$verbose)
      cat("Building orthogroup and gff annotation objects ... ")

    ##############################################################################
    # 1. Split global orthogroups by synteny
    ##############################################################################
    # read in the orthogroup vector
    genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]
    gffAll <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), blastDir = gsParam$blast)
    ov <- gffAll$ord; names(ov) <- gffAll$ofID
    cv <- gffAll$chr; names(cv) <- gffAll$ofID
    sv <- gffAll$start; names(sv) <- gffAll$ofID
    ev <- gffAll$end; names(ev) <- gffAll$ofID
    ogv <- pull_orthogroups(gsParam$blast, genomeIDs = genomeIDs)$ogv

    # read in the syntenic hits, subset to only the same orthogroup
    synOg <- rbindlist(mclapply(synParam$synHitFiles, mc.cores = gsParam$nCores, function(i)
      subset(fread(i, select = c("ofID1","ofID2","og1","og2")), og1 == og2)[,c("ofID1","ofID2")]))

    # make graph
    synOg[,og := clus_igraph(ofID1, ofID2)]
    synOg <- synOg[,list(u = unique(c(ofID1, ofID2))), by = "og"]
    sogv <- synOg$og; names(sogv) <- synOg$u
    rm(list = "synOg")

    ##############################################################################
    # 2. Calculate RBH hits that do not have an ortholog within each group
    ##############################################################################
    # read in the orthogroup vector
    if(gsParam$verbose)
      cat("Done!\nPulling RBHs that are missing from orthogroups ... ")
    rbhs <- rbindlist(mclapply(1:nrow(synParam), mc.cores = gsParam$nCores, function(i){
      sp <- synParam[i,]
      hits <- fread(sp$synHitFiles, select = c("ofID1", "ofID2", "score", "blkID"))
      hits[,isSelf := any(ofID1 == ofID2), by = "blkID"]
      hits[,`:=`(og1 = sogv[ofID1], og2 = sogv[ofID2])]
      rbh <- pull_rbh(hits)[,c("ofID1", "ofID2", "score", "blkID")]
      return(rbh)
    }))


    ##############################################################################
    # 3. Make final blocks and coordinates from orthogroups and RBHs
    ##############################################################################

    if(gsParam$verbose)
      cat("Done!\nForming blocks on RBH and OG hits",
          "\n\tGenome Comp.\tOG hits (RBHs): reps; nBlks (drop, split, alt.):final; n missing\n")
    blkCoords <- rbindlist(lapply(1:nrow(synParam), function(i){
      sp <- synParam[i,]
      if(gsParam$verbose){
        s1 <- substr(sp$genome1, 1, 5)
        s2 <- substr(sp$genome2, 1, 5)
        if(nchar(s1) < 5)
          s1 <- sprintf("%s%s", s1, paste(rep(" ", 5-nchar(s1)), collapse = ""))
        if(nchar(s2) < 5)
          s2 <- sprintf("%s%s", s2, paste(rep(" ", 5-nchar(s2)), collapse = ""))
        cat(sprintf("\t%s vs %s:  ",s1,s2))
      }

      # read in blasts, keeping only same orthogroups or RBHs
      hits <- fread(sp$synHitFiles, select = c("ofID1", "ofID2", "score", "blkID"))
      hits[,`:=`(og1 = sogv[ofID1], og2 = sogv[ofID2],
                 isRBH = paste(ofID1, ofID2) %in% paste(rbhs$ofID1, rbhs$ofID2))]
      hits <- subset(hits, og1 == og2 | isRBH)

      if(gsParam$verbose)
        with(hits, cat(sprintf("%s (%s)", length(og1), sum(isRBH))))
      hits[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2],
                 chr1 = cv[ofID1],  chr2 = cv[ofID2],
                 genome1 = sp$genome1, genome2 = sp$genome2)]

      # get representatives
      hitRep <- prune_hits2reps(
        hits = hits,
        blkSize = sp$blkSize,
        nGaps = sp$nGaps,
        path2mcscanx = gsParam$path2mcscanx,
        nCores = gsParam$nCores)
      if(gsParam$verbose)
        cat(sprintf(": %s;", nrow(hitRep)))

      # make blocks
      hitBlk <- make_blks(
        hits = hitRep,
        nCores = gsParam$nCores,
        path2mcscanx = gsParam$path2mcscanx,
        blkSize = sp$blkSize,
        nGaps = sp$nGaps)
      hitBlk[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2], blkID = paste0("blk_", blkID))]
      ub <- unique(hitBlk$blkID)

      # split blocks that were interleaved
      hitBlk <- split_intlvBlks(
        hits = hitBlk,
        blkSize = dropInterleavesSmallerThan)
      us <- unique(hitBlk$blkID)
      usl <- sapply(us, function(x) strsplit(x, " ")[[1]][1])
      if(gsParam$verbose)
        cat(sprintf(
          " %s (%s, %s, %s): %s",
          length(ub), sum(!ub %in% usl), sum(duplicated(usl)),
          sum(!us %in% ub), length(us)))

      setkey(hitBlk, chr1, chr2, ord1, ord2)
      hitBlk[,blkID := sprintf("%s_%s_%s", genome1[1], genome2[1], as.numeric(as.factor(blkID)))]
      hitBlk[,`:=`(ord1 = ov[ofID1], ord2 = ov[ofID2])]

      # plot the hits
      pdf(file.path(gsParam$results, sprintf("%s_vs_%s_synBlkHits.pdf", sp$genome1, sp$genome2)),
        height = 8, width = 8)
      pltdat <- plot_hits(
        hits = hitBlk,
        gsParam = gsParam,
        gsAnnot = gsAnnot,
        reorderChrs = FALSE)
      dev.off()
      ##############################################################################
      # 4. Get blk coordinates and place genes therein.
      ##############################################################################

      # pull synog hits in coordinates

      blkCoord <- hitBlk[,list(s1 = min(ord1, na.rm = T), e1 = max(ord1, na.rm = T),
                               s2 = min(ord2, na.rm = T), e2 = max(ord2, na.rm = T),
                               u = paste(chr1[1], chr2[1])),
                         by = c("chr1","chr2","blkID")]

      hitBlk[,`:=`(sv1 = sv[ofID1], ev1 = ev[ofID1],
                   sv2 = sv[ofID2], ev2 = ev[ofID2])]
      bcout <- hitBlk[,list(
        ordStart1 = min(ord1, na.rm = T),
        ordEnd1 = max(ord1, na.rm = T),
        ordStart2 = min(ord2, na.rm = T),
        ordEnd2 = max(ord2, na.rm = T),
        bpStart1 = min(sv1, rm = T),
        bpEnd1 = max(ev1, rm = T),
        bpStart2 = min(sv2, rm = T),
        bpEnd2 = max(ev2, rm = T),
        nHits = uniqueN(c(ofID1, ofID2)),
        orient = ifelse(length(ord1) <= 1, "+",
                        ifelse(cor(jitter(ord1),
                                   jitter(ord2)) > 0,"+", "-")),
        firstGene1 = ofID1[which.min(ord1)],
        firstGene2 = ofID2[which.min(ord2)],
        lastGene1 = ofID1[which.max(ord1)],
        lastGene2 = ofID2[which.max(ord2)]),
        by = c("chr1","chr2","blkID")]

      bcout[,`:=`(genome1 = sp$genome1, genome2 = sp$genome2)]
      blkCoord[,u := paste(chr1, chr2)]
      hits[,u := paste(chr1, chr2)]
      setkey(hits, ord1, ord2)
      spl <- split(hits, by = "u")
      gspl1 <- split(subset(gffAll, genome == sp$genome1), by = "chr")
      gspl2 <- split(subset(gffAll, genome == sp$genome2), by = "chr")

      ogo <- rbindlist(mclapply(1:nrow(blkCoord), mc.cores = gsParam$nCores, function(j){
        x <- blkCoord[j,]

        y <- spl[[x$u]]
        y <- subset(y, ord1 >= x$s1 & ord1 <= x$e1 & ord2 >= x$s2 & ord2 <= x$e2)
        if(!inferMissingPos){
          return(y)
        }else{
          if(nrow(y) == 0){
            return(NULL)
          }else{
            y$ofID2[y$isRBH] <- paste0(y$ofID2[y$isRBH], "RBH")
            o <- y[,c("ofID1","ofID2","ord1","ord2")]

            if(any(o$ofID1 == o$ofID2)){
              o <- y[,c("ofID1","ofID2")]
              o[,blkID := x$blkID]

              return(o)
            }else{
              # infer positions of missing genes on genome1
              g1 <- subset(gspl1[[x$chr1]], ord >= x$s1 & ord <= x$e1)
              g1 <- with(subset(g1, !ofID %in% y$ofID1), data.table(ofID1 = ofID, ord1 = ord))
              g2 <- subset(gspl2[[x$chr2]], ord >= x$s2 & ord <= x$e2)
              g2 <- with(subset(g2, !ofID %in% y$ofID2), data.table(ofID2 = ofID, ord2 = ord))

              if(nrow(g1) > 0){
                tmp1 <- rbind(o, g1, fill = T)
                tmp1[,ord2 := round(infer_pos(ord1 = as.numeric(ord2), ord2 = ord1),1)]
                v <- is.na(tmp1$ofID2)
                tmp1$ofID2[v] <- sprintf(
                  "%s:%s:%s:%s:%s",  sp$genome2, x$chr2, x$blkID, tmp1$ord2[v], 1:sum(v))
              }else{
                tmp1 <- NULL
              }

              if(nrow(g2) > 0){
                tmp2 <- rbind(o, g2, fill = T)
                tmp2[,ord1 := round(infer_pos(ord1 = as.numeric(ord1), ord2 = ord2),1)]
                v <- is.na(tmp2$ofID1)
                tmp2$ofID1[v] <- sprintf(
                  "%s:%s:%s:%s:%s",  sp$genome1,x$chr1, x$blkID, tmp2$ord1[v], 1:sum(v))
              }else{
                tmp2 <- NULL
              }
              out <- rbind(o, tmp1, tmp2)[,c("ofID1","ofID2")]
              out <- subset(out, !duplicated(out))
              out[,blkID := x$blkID]

              return(out)
            }
          }
        }
      }))

      # pull missing genes in coordinates
      nb <- with(ogo, sum(grepl(":",paste(ofID1, ofID2), fixed = T)))
      if(gsParam$verbose)
        cat(sprintf("; %s", nb))
      cat("\n")
      fwrite(
        ogo,
        file = blkFiles[i],
        sep = "\t")
      return(bcout)
    }))

    ##############################################################################
    # 5. Produce low-memory output vectors and pangenome
    ##############################################################################

    # parse and write blks and orthogroups
    fwrite(blkCoords, sep = "\t", quote = F,
           file = file.path(gsParam$results,
                            "blockCoordinates_geneOrder.txt.gz"))
    return(synogFiles = blkFiles)
  }
}

#' @title build high-confidence blocks
#' @description
#' \code{make_blks} re-form blocks within pre-existing blocks, specified with
#' a blkID column in the hits object. Uses MCScanX to prune, then dbscan to
#' make the actual blocks.
#' @rdname pull_synOgs
#' @import data.table
#' @export
make_blks <- function(hits, nCores, blkSize, nGaps, path2mcscanx){
  splHits <- split(hits, by = "blkID")
  out <- rbindlist(mclapply(splHits, mc.cores = nCores, function(y){
    ib <- y$blkID[1]
    y[,`:=` (ord1 = frank(ord1, ties.method = "dense"),
             ord2 = frank(ord2, ties.method = "dense"))]
    if(any(y$ofID1 %in% y$ofID2)){
      return(y)
    }else{
      y[,blkID := run_mcscanx(
        hits = y, blkSize = blkSize,
        nGaps = nGaps, path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)]]

      if("blkID" %in% colnames(y)){
        y <- subset(y, !is.na(blkID))
        if(nrow(y) == 0){
          y <- NULL
        }else{
          y[,`:=` (ord1 = frank(ord1, ties.method = "dense"),
                   ord2 = frank(ord2, ties.method = "dense"))]
          y[,blkID := dbscan(frNN(cbind(ord1, ord2), eps = sqrt(blkSize^2+blkSize^2)/2,),
                             minPts = 2)$cluster]
        }
      }
      if("blkID" %in% colnames(y)){
        y <- subset(y, !is.na(blkID))
        if(nrow(y) == 0){
          y <- NULL
        }else{
          y[,blkID := paste(ib, blkID)]
        }
        return(y)
      }else{
        return(NULL)
      }
    }
  }))
  out[,blkID := as.numeric(factor(blkID, levels = unique(blkID)))]
  return(out)
}

#' @title prune hits to representatives by orthogroup
#' @description
#' \code{prune_hits2reps} Condense arrays to hits that are a) are nearest the
#' middle of the array and b) in mcscanx blocks.
#' @rdname pull_synOgs
#' @import data.table
#' @export
prune_hits2reps <- function(hits, blkSize, nGaps, path2mcscanx, nCores){
  spl <- split(hits, by = "blkID")
  out <- rbindlist(mclapply(spl, mc.cores = nCores, function(y){
    y[,clus := clus_igraph(ofID1, ofID2)]
    ib <- y$blkID[1]
    if(any(y$ofID1 %in% y$ofID2)){
     return(subset(y, ofID1 == ofID2))
    }else{
      y[,`:=`(tmp1 = ord1, tmp2 = ord2)]
      y[,`:=` (ord1 = frank(ord1, ties.method = "dense"),
               ord2 = frank(ord2, ties.method = "dense"))]
      y[,`:=` (mord1 = as.integer(median(ord1, na.rm = T)),
               mord2 = as.integer(median(ord2, na.rm = T))),
        by = "clus"]
      y[,`:=`(sumDi = sqrt(abs(mord1 - ord1)^2 + abs(mord2 - ord2)^2))]
      setkey(y, sumDi)
      y <- subset(y, !duplicated(clus))
      y[,`:=` (ord1 = frank(ord1, ties.method = "dense"),
               ord2 = frank(ord2, ties.method = "dense"),
               mord1 = NULL, mord2 = NULL, sumDi = NULL)]
      if(nrow(y) < 2){
        return(NULL)
      }else{
        if(nrow(y) <= blkSize){
          suppressWarnings(y[,tmp := run_mcscanx(
            hits = y, blkSize = nrow(y)-1,
            nGaps = nGaps, path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)]])
        }else{
          suppressWarnings(y[,tmp := run_mcscanx(
            hits = y, blkSize = blkSize,
            nGaps = nGaps, path2mcscanx = path2mcscanx)[paste(ofID1, ofID2)]])
        }
        if("tmp" %in% colnames(y)){
          y <- subset(y, !is.na(tmp))
          y[,`:=`(tmp = NULL, ord1 = NULL, ord2 = NULL)]
          y[,`:=`(ord1 = tmp1, ord2 = tmp2, tmp1 = NULL, tmp2 = NULL)]
          return(y)
        }
      }
    }
  }), fill = T)
  return(out)
}

#' @title pull reciprocal best scoring hits
#' @description
#' \code{pull_rbh} Given a hits object with a 'blockID' column and og1/og2,
#' parse out the reciprocal best hits that are not duplicated and do not have
#' identity in an intergenomic/interhomeolog to a gene in the same orthogroup.
#' @rdname pull_synOgs
#' @import data.table
#' @export
pull_rbh <- function(hits){
  hits <- subset(hits, !isSelf)
  hcols <- colnames(hits)
  hits[,anyOg1 := any(og1 == og2), by = c("blkID", "ofID1")]
  hits[,anyOg2 := any(og1 == og2), by = c("blkID", "ofID2")]
  hits[,maxScr1 := max(score), by = c("blkID", "ofID1")]
  hits[,maxScr2 := max(score), by = c("blkID", "ofID2")]
  ogb1 <- with(subset(hits, anyOg1), paste(blkID, ofID1))
  ogb2 <- with(subset(hits, anyOg2), paste(blkID, ofID2))
  rbh <- subset(hits, !paste(blkID, ofID1) %in% ogb1 & !paste(blkID, ofID2) %in% ogb2 & maxScr1 == score & maxScr2 == score)
  rbhd1 <- with(rbh, paste(ofID1, blkID)[duplicated(paste(ofID1, blkID))])
  rbhd2 <- with(rbh, paste(ofID2, blkID)[duplicated(paste(ofID2, blkID))])
  rbh <- subset(rbh, !paste(ofID1, blkID) %in% rbhd1 & !paste(ofID2, blkID) %in% rbhd2)
  return(rbh[,hcols, with = F])
}

#' @title split interleaved blocks
#' @description
#' \code{split_intlvBlks} Given hits and a blkID column, check overlapping
#' coordinates that do not have duplicated hits. Split these into their own
#' blocks, retaining blocks with > blkSize number of unique hits.
#' @rdname pull_synOgs
#' @import data.table
#' @export
split_intlvBlks <- function(hits, blkSize = 2){

  i.start <- i.end <- x <- ndup <- blkID <- i.blkID <- rlID <- ord1 <- NULL
  end <- start <- ord1 <- chr1 <- ord2 <- chr2 <- i.x <- rl <- NULL

  blks <- data.table(hits)
  blkCols <- colnames(blks)
  setkey(blks, blkID)
  ub <- unique(blks$blkID)

  # check overlapping on the first genome
  b <- blks[,list(start = min(ord1), end = max(ord1)), by = c("chr1","blkID")]
  setkey(b, chr1, start, end)
  bo <- subset(foverlaps(b, b), blkID != i.blkID)

  if(nrow(bo) > 0){
    # get list of genes in each block
    bu <- with(blks, tapply(ofID1, blkID, function(x) list(unique(x))))

    # for each overlapping pair of blocks, count the number of duplicates
    bo[,ndup := sum(duplicated(c(bu[[blkID]], bu[[i.blkID]]))),
       by = c("blkID","i.blkID")]

    # subset to overlaps with no dups with the larger block as the first block
    bo[,`:=` (i.x = jitter(i.end), x = jitter(end),
              sz = (end - start) + (i.end - i.start))]
    setorder(bo, -sz)
    n0 <- subset(bo, ndup == 0 & (x-start) >= (i.x-i.start))

    while(nrow(n0) > 0){
      # subset to at most one unique entry for each block
      tmp <- subset(n0, !duplicated(blkID))
      tmp <- subset(tmp, !i.blkID %in% blkID)
      tmp <- subset(tmp, !duplicated(i.blkID))

      tblkID <- unique(c(tmp$blkID, tmp$i.blkID))
      tblk <- subset(blks, blkID %in% tblkID)
      bspl <- split(tblk, by = "blkID")
      tmp <- subset(tmp, complete.cases(tmp))
      tmp <- subset(tmp, blkID %in% names(bspl) & i.blkID %in% names(bspl))

      # for each unique entry of overlapping blocks, pull the hits
      splt <- rbindlist(lapply(1:nrow(tmp), function(i){
        x <- NULL
        ji <- as.character(tmp$blkID[i])
        if(length(ji) == 1 && ji %in% names(bspl)){
          xj <- bspl[[ji]]
        }else{
          xj <- NULL
        }
        ki <- as.character(tmp$i.blkID[i])
        if(length(ki) == 1 && ki %in% names(bspl)){
          xk <- bspl[[ki]]
        }else{
          xk <- NULL
        }
        x <- rbind(xj, xk)
        if(nrow(x) > 0 && !is.null(x)){
          setkey(x, ord1)
          mrl <- (-1)

          # split the run of genes, tossing those up to the min block size
          while(mrl < blkSize){
            x[,rl := add_rle(blkID)]
            mrl <- ifelse(nrow(x) == 0,blkSize, min(x$rl))
            if(mrl < blkSize){
              x[,rlID := add_rle(blkID, which = "id")]
              smb <- with(x, tapply(rl, rlID, min))
              smb <- smb[which.min(smb)]
              x <- subset(x, rlID != names(smb))
            }
          }
          if(nrow(x) >= blkSize){
            x[,rlID := add_rle(blkID, which = "id")]
            x[,blkID := paste(blkID, rlID)]
            return(x)
          }
        }
      }))
      if(!is.null(splt) && nrow(splt) > 0){
        blks <- rbind(subset(blks, !blkID %in% names(bspl)),
                      splt[,colnames(blks), with = F])
      }else{
        blks <- subset(blks, !blkID %in% names(bspl))
      }



      b <- blks[,list(start = min(ord1), end = max(ord1)), by = c("chr1","blkID")]
      setkey(b, chr1, start, end)
      bo <- subset(foverlaps(b, b), blkID != i.blkID)

      bu <- with(blks, tapply(ofID1, blkID, function(x) list(unique(x))))
      if(nrow(bo) == 0){
        n0 <- data.table(bo)
      }else{
        # for each overlapping pair of blocks, count the number of duplicates
        bo[,ndup := sum(duplicated(c(bu[[blkID]], bu[[i.blkID]]))),
           by = c("blkID","i.blkID")]

        # subset to overlaps with no dups with the larger block as the first block
        bo[,`:=` (i.x = jitter(i.end), x = jitter(end),
                  sz = end - start, i.sz = i.end - i.start)]
        setorder(bo, -sz, -i.sz)
        n0 <- subset(bo, ndup == 0 & (x-start) >= (i.x-i.start))
      }
    }
  }

  b <- blks[,list(start = min(ord2), end = max(ord2)), by = c("chr2","blkID")]
  setkey(b, chr2, start, end)
  bo <- subset(foverlaps(b, b), blkID != i.blkID)

  if(nrow(bo) > 0){
    # get list of genes in each block
    bu <- with(blks, tapply(ofID2, blkID, function(x) list(unique(x))))

    # for each overlapping pair of blocks, count the number of duplicates
    bo[,ndup := sum(duplicated(c(bu[[blkID]], bu[[i.blkID]]))),
       by = c("blkID","i.blkID")]

    # subset to overlaps with no dups with the larger block as the first block
    bo[,`:=` (i.x = jitter(i.end), x = jitter(end),
              sz = (end - start) + (i.end - i.start))]
    setorder(bo, -sz)
    n0 <- subset(bo, ndup == 0 & (x-start) >= (i.x-i.start))

    while(nrow(n0) > 0){
      # subset to at most one unique entry for each block
      tmp <- subset(n0, !duplicated(blkID))
      tmp <- subset(tmp, !i.blkID %in% blkID)
      tmp <- subset(tmp, !duplicated(i.blkID))

      tblkID <- unique(c(tmp$blkID, tmp$i.blkID))
      tblk <- subset(blks, blkID %in% tblkID)
      bspl <- split(tblk, by = "blkID")
      tmp <- subset(tmp, complete.cases(tmp))
      tmp <- subset(tmp, blkID %in% names(bspl) & i.blkID %in% names(bspl))

      # for each unique entry of overlapping blocks, pull the hits
      splt <- rbindlist(lapply(1:nrow(tmp), function(i){
        x <- NULL
        ji <- as.character(tmp$blkID[i])
        if(length(ji) == 1 && ji %in% names(bspl)){
          xj <- bspl[[ji]]
        }else{
          xj <- NULL
        }
        ki <- as.character(tmp$i.blkID[i])
        if(length(ki) == 1 && ki %in% names(bspl)){
          xk <- bspl[[ki]]
        }else{
          xk <- NULL
        }
        x <- rbind(xj, xk)
        if(nrow(x) > 0 && !is.null(x)){
          setkey(x, ord2)
          mrl <- (-1)

          # split the run of genes, tossing those up to the min block size
          while(mrl < blkSize){
            x[,rl := add_rle(blkID)]
            mrl <- ifelse(nrow(x) == 0,blkSize, min(x$rl))
            if(mrl < blkSize){
              x[,rlID := add_rle(blkID, which = "id")]
              smb <- with(x, tapply(rl, rlID, min))
              smb <- smb[which.min(smb)]
              x <- subset(x, rlID != names(smb))
            }
          }
          if(nrow(x) >= blkSize){
            x[,rlID := add_rle(blkID, which = "id")]
            x[,blkID := paste(blkID, rlID)]
            return(x)
          }
        }
      }))

      if(!is.null(splt) && nrow(splt) > 0){
        blks <- rbind(subset(blks, !blkID %in% names(bspl)),
                      splt[,colnames(blks), with = F])
      }else{
        blks <- subset(blks, !blkID %in% names(bspl))
      }

      b <- blks[,list(start = min(ord2), end = max(ord2)), by = c("chr2","blkID")]
      setkey(b, chr2, start, end)
      bo <- subset(foverlaps(b, b), blkID != i.blkID)

      bu <- with(blks, tapply(ofID2, blkID, function(x) list(unique(x))))
      if(nrow(bo) == 0){
        n0 <- data.table(bo)
      }else{
        # for each overlapping pair of blocks, count the number of duplicates
        bo[,ndup := sum(duplicated(c(bu[[blkID]], bu[[i.blkID]]))),
           by = c("blkID","i.blkID")]

        # subset to overlaps with no dups with the larger block as the first block
        bo[,`:=` (i.x = jitter(i.end), x = jitter(end),
                  sz = end - start, i.sz = i.end - i.start)]
        setorder(bo, -sz, -i.sz)
        n0 <- subset(bo, ndup == 0 & (x-start) >= (i.x-i.start))
      }
    }
  }
  return(blks)
}

#' @title split interleaved blocks
#' @description
#' \code{infer_pos} Given hits and a blkID column, check overlapping
#' coordinates that do not have duplicated hits. Split these into their own
#' blocks, retaining blocks with > blkSize number of unique hits.
#' @rdname pull_synOgs
#' @import data.table
#' @export
infer_pos <- function(ord1, # named order vector from gff
                      ord2){ # all alt ofIDs within blk coordinates
  o12 <- subset(
    data.table(o1 = ord1, o2 = ord2, index = 1:length(ord1), key = "o2"),
    !is.na(o2))
  if(sum(complete.cases(o12)) < 2){
    return(NA)
  }else{
    whNaSt <- which(is.na(o12$o1) & o12$o2 == min(o12$o2, na.rm = T))
    whNaEn <- which(is.na(o12$o1) & o12$o2 == max(o12$o2, na.rm = T))
    o12$o2[whNaSt] <-  o12$o2[whNaSt] + .001
    o12$o2[whNaEn] <-  o12$o2[whNaEn] - .001
    setkey(o12, o2)
    o12[,runID := add_rle(is.na(o1), which = "id")]
    oi <- subset(o12, !is.na(o1))
    om <- subset(o12, is.na(o1))
    si <- oi[,list(r = o1[1], l = o1[length(o1)]), by = "runID"]
    rv <- si$r; names(rv) <- as.character(si$runID-1)
    lv <- si$l; names(lv) <- as.character(si$runID+1)
    om[,l := lv[as.character(runID)]]
    om[,r := rv[as.character(runID)]]
    om[,o1 := seq(from = l[1], to = r[1], length.out = .N+2)[2:(.N+1)],
          by = "runID"]
    ord1[om$index] <- om$o1
    return(ord1)
  }
}
