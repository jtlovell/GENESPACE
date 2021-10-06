#' @title Orthofinder-calling routines from within R
#' @description
#' \code{rerun_orthofinder} Sets of functions to run orthofinder from within
#' the R environment. Used for blockwise and self-array routines.
#' @name rerun_orthofinder
#'
#' @param gsParam list of genespace parameters
#' @param writeDir file path character string

#' \cr
#' If called, \code{rerun_orthofinder} returns its own arguments.
#'
#'

#' @title blkwise_orthofinder
#' @description
#' \code{blkwise_orthofinder} blkwise_orthofinder
#' @rdname rerun_orthofinder
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom Biostrings readAAStringSet
#' @export
blkwise_orthofinder <- function(gsParam,
                                genomeIDs = NULL,
                                minGenes4of = 40,
                                verbose = NULL){

  #############################################################################
  # 1.Checking
  #############################################################################
  # -- get various parameters
  if(is.null(genomeIDs))
    genomeIDs <- gsParam$genomes$genomeIDs
  if(is.null(verbose))
    verbose <- gsParam$params$verbose

  # -- make sure orthofinder has been run and blast results exist
  if(is.na(gsParam$paths$blastDir))
    gsParam <- find_orthofinderResults(gsParam)

  # -- read in gff
  gffFile <- file.path(gsParam$paths$results, "gffWithOgs.txt.gz")
  gff <- subset(fread(gffFile, showProgress = F), genome %in% genomeIDs)
  ov <- gff$ord; sv <- gff$start; ev <- gff$end
  names(ov) <- names(sv) <- names(ev) <- gff$ofID

  # -- get the synParams in
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, gnum1 <= gnum2)
  setkey(synp, gnum1, gnum2)

  # -- check and report the number of blocks that will be excluded.
  if(verbose)
    cat("Checking blocks\n")
  blks <- count_nBlks(
    gsParam = gsParam, genomeIDs = genomeIDs, minGenes4of = minGenes4of)
  selfBlks <- subset(blks, grepl("self", blkID))

  if(verbose)
    cat(sprintf(
      "Running synteny-constrained orthofinder for %s genome combinations\n",
      nrow(synp)))
  cat("\tGen. comb. (type): private / single region / multiple regions\n")
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  pepspl <- sapply(genomeIDs, USE.NAMES = T, simplify = F, function(i)
    readAAStringSet(file.path(gsParam$paths$blastDir,
                              sprintf("Species%s.fa", ofSpId[i]))))
  blnk <- paste(rep(" ", 16), collapse = "")
  test <- lapply(1:nrow(synp), function(i){
    fs <- with(synp[i,], file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",genome1, genome2)))
    hits <- fread(
      fs, na.strings = c("NA", ""), showProgress = F,
      select = c("regID", "regBuffer","regAnchor","ofID1","ofID2","gen1","gen2",
                 "chr1","chr2","start1","start2","end1","end2","ord1","ord2","og"))
    setnames(hits, c("regID", "regBuffer", "regAnchor"),
             c("blkID", "blkBuffer", "blkAnchor"))
    blks <- calc_blkCoords(hits)
    blk1 <- with(blks, data.table(
      chr = chr1, start = startBp1, end = endBp1, blkID = blkID))
    blk2 <- with(blks, data.table(
      genome = gen2, chr = chr2, start = minBp2, end = maxBp2, blkID = blkID))
    g1 <- with(subset(gff, genome == hits$gen1[1]), data.table(
      chr = chr, start = (end + start)/2, end = (end + start)/2, ofID = ofID))
    g2 <- with(subset(gff, genome == hits$gen2[1]), data.table(
      chr = chr, start = (end + start)/2, end = (end + start)/2, ofID = ofID))
    setkey(blk1, chr, start, end)
    setkey(blk2, chr, start, end)
    setkey(g1, chr, start, end)
    setkey(g2, chr, start, end)

    fo1 <- foverlaps(blk1, g1)[,c("ofID", "blkID")]
    fo2 <- foverlaps(blk2, g2)[,c("ofID", "blkID")]
    foExp1 <- fo1[,list(n = uniqueN(blkID)), by = "ofID"]
    foExp2 <- fo2[,list(n = uniqueN(blkID)), by = "ofID"]

    gcnt <- merge(gff, rbind(fo1, fo2), by = "ofID", allow.cartesian = T)[,c("ofID","og","synOg","blkID")]
    oghits <- subset(hits, !is.na(og))
    fo1[,hasOg := paste(ofID, blkID) %in% paste(oghits$ofID1, oghits$blkID)]
    fo2[,hasOg := paste(ofID, blkID) %in% paste(oghits$ofID2, oghits$blkID)]

    ngenomes <- with(synp[i,], uniqueN(c(genome1, genome2)))
    setnames(fo1, 1:2, c("ofID1", "blkID1"))
    setnames(fo2, 1:2, c("ofID2", "blkID2"))
    fo1x <- data.table(fo1)
    setnames(fo1x, 1:2, c("ofID2", "blkID2"))
    fo2x <- data.table(fo2)
    setnames(fo2x, 1:2, c("ofID1", "blkID1"))

    if(ngenomes == 1){
      h <- read_blast(
        path = gsParam$paths$blastDir, onlyIDScore = F,
        ofID1 = ofSpId[synp$genome1[i]], ofID2 = ofSpId[synp$genome1[i]])
    }else{
      h <- read_blast(
        path = gsParam$paths$blastDir, onlyIDScore = F,
        ofID1 = ofSpId[synp$genome1[i]], ofID2 = ofSpId[synp$genome1[i]])
      h1 <- h[,c(2,1,3:6,8,7,10,9,11,12)]
      setnames(h1, colnames(h))
      h <- rbind(h, h1)
      setorder(h, -V12)
      h <- subset(h, !duplicated(paste(V1, V2)))
      setnames(h, 1:2, c("ofID1", "ofID2"))
      h1 <- merge(h, fo1x[,1:2], by = "ofID2", allow.cartesian = T)
      h1 <- merge(h1, fo1[,1:2], by = "ofID1", allow.cartesian = T)
      h00 <- subset(h1, blkID1 == blkID2)

      h <- rbind(
        read_blast(
          path = gsParam$paths$blastDir, onlyIDScore = F,
          ofID1 = ofSpId[synp$genome1[i]], ofID2 = ofSpId[synp$genome2[i]]),
        read_blast(
          path = gsParam$paths$blastDir, onlyIDScore = F,
          ofID1 = ofSpId[synp$genome2[i]], ofID2 = ofSpId[synp$genome1[i]]))
      h1 <- h[,c(2,1,3:6,8,7,10,9,11,12)]
      setnames(h1, colnames(h))
      h <- rbind(h, h1)
      setorder(h, -V12)
      h <- subset(h, !duplicated(paste(V1, V2)))
      setnames(h, 1:2, c("ofID1", "ofID2"))
      h1 <- merge(h, fo2[,1:2], by = "ofID2", allow.cartesian = T)
      h1 <- merge(h1, fo1[,1:2], by = "ofID1", allow.cartesian = T)
      h01 <- subset(h1, blkID1 == blkID2)
      h10 <- h01[,c(2,1,3:6,8,7,10,9,11:14)]
      setnames(h10, colnames(h01))

      h <- read_blast(
        path = gsParam$paths$blastDir, onlyIDScore = F,
        ofID1 = ofSpId[synp$genome2[i]], ofID2 = ofSpId[synp$genome2[i]])
      h1 <- h[,c(2,1,3:6,8,7,10,9,11,12)]
      setnames(h1, colnames(h))
      h <- rbind(h, h1)
      setorder(h, -V12)
      h <- subset(h, !duplicated(paste(V1, V2)))
      setnames(h, 1:2, c("ofID1", "ofID2"))
      h1 <- merge(h, fo2[,1:2], by = "ofID2", allow.cartesian = T)
      h1 <- merge(h1, fo2x[,1:2], by = "ofID1", allow.cartesian = T)
      h11 <- subset(h1, blkID1 == blkID2)

      setnames(h00, 1:2, c("V1", "V2"))
      setnames(h01, 1:2, c("V1", "V2"))
      setnames(h10, 1:2, c("V1", "V2"))
      setnames(h11, 1:2, c("V1", "V2"))
      spl00 <- lapply(split(h00, by = "blkID1"), function(x) x[,1:12])
      spl01 <- lapply(split(h01, by = "blkID1"), function(x) x[,1:12])
      spl10 <- lapply(split(h10, by = "blkID1"), function(x) x[,1:12])
      spl11 <- lapply(split(h11, by = "blkID1"), function(x) x[,1:12])
      ns <- intersect(intersect(intersect(names(spl00), names(spl01)), names(spl10)), names(spl11))
    }

    inBlkOgs <- rbindlist(mclapply(ns, mc.cores = gsParam$params$nCores, mc.preschedule = F, function(j){
      tmpDir <- file.path(gsParam$params$wd, sprintf("%s_og4inBlkTMPdir", j))
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      og <- run_ofFromObj(
        blast00 = spl00[[j]],
        blast01 = spl01[[j]],
        blast10 = spl10[[j]],
        blast11 = spl11[[j]],
        pep0 = pepspl[[synp$genome1[i]]],
        pep1 = pepspl[[synp$genome2[i]]],
        writeDir = tmpDir)
      if(dir.exists(tmpDir))
        unlink(tmpDir, recursive = T)
      og[,blkID := j]
      return(og)
    }))
    ogv <- with(inBlkOgs, paste(og, blkID)); names(ogv) <- inBlkOgs$ofID
    h01[,`:=`(og1 = ogv[V1], og2 = ogv[V2])]
    hog <- subset(h01, og1 == og2)
    hog <- with(hog, data.table(ofID = c(V1, V2), blkID = blkID1, og = og1))
    fog <- data.table(ofID = c(fo1$ofID1, fo2$ofID2), blkID = c(fo1$blkID1, fo2$blkID2))
    mog <- merge(hog, fog, by = c("ofID", "blkID"), all = T)
    setorder(mog, ofID, blkID, og, na.last = T)
    mog <- subset(mog, !duplicated(mog))
    inBlkOgs[,n := uniqueN(ofID), by = c("blkID","og")]
    ogOut <- subset(inBlkOgs, complete.cases(inBlkOgs))
    og1 <- subset(ogOut, n == 1)
    tmp1 <- with(og1, data.table(
      blkID = blkID, og = og, ofID1 = ofID, ofID2 = ofID))
    og2 <- subset(ogOut, n == 2)
    tmp2 <- og2[,list(ofID1 = ofID[1], ofID2 = ofID[2]), by = c("blkID", "og")]
    tmp2 <- rbind(tmp2, with(tmp2, data.table(
      blkID = blkID, og = og, ofID1 = ofID2, ofID2 = ofID1)))
    og3 <- subset(ogOut, n > 2)
    tmp3 <- og3[,CJ(ofID1 = ofID, ofID2 = ofID, unique = T), by = c("blkID","og")]
    out <- subset(rbind(tmp1, tmp2, tmp3), ofID1 != ofID2)
    out <- subset(out, !duplicated(out))

    grph <- with(out, clus_igraph(ofID1, ofID2))
    g[,blkOg := grph[ofID]]
    g$blkOg[is.na(g$blkOg)] <- max(grph + 1):(sum(is.na(g$blkOg))+ max(grph))
    inblk <- g[,list(n = uniqueN(tmp)), by = "blkOg"]

    hits[,`:=`(ogInblk = grph[ofID1], tmp = grph[ofID2])]
    hits$ogInblk[with(hits, ogInblk != tmp)] <- NA

    if(verbose)
      cat(sprintf(
        "\t%s(in blk): %s / %s / %s / %s / %s\n",
        blnk,
        sum(inblk$n == 1), sum(inblk$n == 2), sum(inblk$n == 3), sum(inblk$n == 4), sum(inblk$n > 4)))
    h <- subset(hits, !is.na(ogInblk))[,c("ofID1", "ofID2", "ogInblk")]
    return(h)
  })
}


#' @title read_synBlastInReg
#' @description
#' \code{read_synBlastInReg} read_synBlastInReg
#' @rdname rerun_orthofinder
#' @import data.table
#' @importFrom dbscan dbscan frNN
#' @importFrom Biostrings readAAStringSet
#' @export
read_synBlastInReg <- function(gff, pwSynParam, gsParam, hits, minGenes4of){
  nCores <- gsParam$params$nCores
  ofSpId <- read_orthofinderSpeciesIDs(gsParam$paths$blastDir)
  syn <- data.table(pwSynParam)
  synBuff <- syn$synBuff

  # -- subset gff
  g1 <- subset(gff, genome == syn$genome1)
  g2 <- subset(gff, genome == syn$genome2)
  ov1 <- g1$ord; names(ov1) <- g1$ofID
  ov2 <- g2$ord; names(ov2) <- g2$ofID
  g1l <- split(g1, by = "chr")
  g2l <- split(g2, by = "chr")

  # -- read in hits
  fs <- with(syn, file.path(
    gsParam$paths$results,
    sprintf("%s_%s_synHits.txt.gz",genome1, genome2)))
  h <- subset(hits, complete.cases(hits))
  h <- subset(h, blkBuffer)

  # -- calculate the block coordinates
  b <- calc_blkCoords(subset(h, blkAnchor))
  b <- subset(b, nHits1 > minGenes4of & nHits2 > minGenes4of)

  # -- read in raw blast
  bl <- read_blast(
    path = gsParam$paths$blastDir,
    ofID1 = ofSpId[syn$genome1],
    ofID2 = ofSpId[syn$genome2],
    onlyIDScore = F)

  blg1 <- read_blast(
    path = gsParam$paths$blastDir,
    ofID1 = ofSpId[syn$genome1],
    ofID2 = ofSpId[syn$genome1],
    onlyIDScore = F)

  blg2 <- read_blast(
    path = gsParam$paths$blastDir,
    ofID1 = ofSpId[syn$genome2],
    ofID2 = ofSpId[syn$genome2],
    onlyIDScore = F)

  # -- subset the gff to blocks
  gl <- mclapply(1:nrow(b), mc.cores = gsParam$params$nCores, function(j){
    u <- with(subset(h, blkID == b$blkID[j] & blkAnchor), paste(ofID1, ofID2))
    g1t <- subset(
      g1l[[as.character(b$chr1[j])]],
      ord >= b$startOrd1[j] & ord <= b$endOrd1[j])
    setkey(g1t, start)
    g1t[,`:=`(o = 1:.N, blkOfID = sprintf("0_%s", (1:.N)-1))]

    g2t <- subset(
      g2l[[as.character(b$chr2[j])]],
      ord >= b$minOrd2[j] & ord <= b$maxOrd2[j])
    setkey(g2t, start)
    g2t[,o := 1:.N]
    g2t[,`:=`(o = 1:.N, blkOfID = sprintf("1_%s", (1:.N)-1))]

    ov1 <- g1t$o; gv1 <- g1t$blkOfID; names(ov1) <- names(gv1) <- g1t$ofID
    ov2 <- g2t$o; gv2 <- g2t$blkOfID; names(ov2) <- names(gv2) <- g2t$ofID

    bl12 <- subset(bl, V1 %in% g1t$ofID)
    bl12 <- subset(bl12, V2 %in% g2t$ofID)
    bl12[,`:=`(ord1 = ov1[V1], ord2 = ov2[V2], isAnchor = paste(V1, V2) %in% u)]

    nn <- with(bl12, frNN(x = data.frame(ord1, ord2), eps = synBuff))
    wh <- unique(c(which(bl12$isAnchor), unlist(nn$id[bl12$isAnchor])))
    bl12[,inBuffer := 1:.N %in% wh]
    bl12 <- subset(bl12, inBuffer)[,1:12]
    bl21 <- bl12[,c(2,1,3:6,8,7,10,9,11,12)]
    setnames(bl21, colnames(bl12))

    buff <- sqrt(synBuff^2 * 2) + 1

    bl11 <- subset(blg1, V1 %in% unique(bl12$V1))
    bl11 <- subset(bl11, V2 %in% unique(bl12$V1))
    bl11[,`:=`(ord1 = ov1[V1], ord2 = ov1[V2])]
    bl11 <- subset(bl11, abs(ord1 - ord2) <= synBuff*2)
    bl11[,`:=`(ord1 = NULL, ord2 = NULL)]

    bl22 <- subset(blg2, V1 %in% unique(bl12$V2))
    bl22 <- subset(bl22, V2 %in% unique(bl12$V2))
    bl22[,`:=`(ord1 = ov2[V1], ord2 = ov2[V2])]
    bl22 <- subset(bl22, abs(ord1 - ord2) <= synBuff*2)
    bl22[,`:=`(ord1 = NULL, ord2 = NULL)]

    return(list(
      bl11 = bl11, bl12 = bl12, bl21 = bl21, bl22 = bl22,
      blkCoord = b[j,]))
  })
  return(gl)
}

#' @title Run orthofinder from R blast data.table
#' @description
#' \code{run_ofFromObj} General engine for blockwise and syntenic array runs
#' @rdname rerun_orthofinder
#' @import data.table
#' @importFrom Biostrings writeXStringSet
#' @export
run_ofFromObj <- function(blast00,
                          blast01,
                          blast10,
                          blast11,
                          pep0,
                          pep1,
                          writeDir){

  if(dir.exists(writeDir))
    stop(sprintf("%s exists. Specify non-existing directory\n",
                 writeDir))
  dir.create(writeDir)

  if(colnames(blast00)[1] != "ofID1"){
    setnames(blast00, 1:2, c("ofID1", "ofID2"))
    setnames(blast01, 1:2, c("ofID1", "ofID2"))
    setnames(blast10, 1:2, c("ofID1", "ofID2"))
    setnames(blast11, 1:2, c("ofID1", "ofID2"))
  }

  # -- make gene ID dictionaries
  id0 <- unique(c(blast00$ofID1, blast00$ofID2, blast01$ofID1, blast10$ofID2))
  id1 <- unique(c(blast11$ofID1, blast11$ofID2, blast01$ofID2, blast10$ofID1))
  id0 <- id0[order(id0)]
  id1 <- id1[order(id1)]
  names(id0) <- sprintf("0_%s", (1:length(id0))-1)
  names(id1) <- sprintf("1_%s", (1:length(id1))-1)

  # -- ensure that all ids are in the peptide files
  id0 <- id0[id0 %in% names(pep0)]
  id1 <- id1[id1 %in% names(pep1)]

  # -- rename peptides and invert dictionary
  p0 <- pep0[id0]; names(p0) <- names(id0)
  p1 <- pep1[id1]; names(p1) <- names(id1)
  di0 <- names(id0); di1 <- names(id1); names(di0) <- id0; names(di1) <- id1

  # -- write the peptide files / fake diamond dbs
  writeXStringSet(p0, filepath = file.path(writeDir, "Species0.fa"))
  writeXStringSet(p1, filepath = file.path(writeDir, "Species1.fa"))
  cat("NA", file = file.path(writeDir, "diamondDBSpecies0.dmnd"))
  cat("NA", file = file.path(writeDir, "diamondDBSpecies1.dmnd"))

  # -- write the species and sequence IDs
  sid <- data.table(
    of = paste0(c(names(id0), names(id1)), ":"),
    id = c(id0, id1))
  fwrite(
    sid, file = file.path(writeDir, "SequenceIDs.txt"),
    sep = " ", quote = F, row.names = F, col.names = F)
  cat(
    c("0: species1.fa", "1: species2.fa"),
    sep = "\n", file = file.path(writeDir, "SpeciesIDs.txt"))

  # -- rename the blast files
  bl00 <- subset(blast00, ofID1 %in% names(di0) & ofID2 %in% names(di0))
  bl01 <- subset(blast01, ofID1 %in% names(di0) & ofID2 %in% names(di1))
  bl10 <- subset(blast10, ofID1 %in% names(di1) & ofID2 %in% names(di0))
  bl11 <- subset(blast11, ofID1 %in% names(di1) & ofID2 %in% names(di1))
  bl00[,`:=`(ofID1 = di0[ofID1], ofID2 = di0[ofID2])]
  bl01[,`:=`(ofID1 = di0[ofID1], ofID2 = di1[ofID2])]
  bl10[,`:=`(ofID1 = di1[ofID1], ofID2 = di0[ofID2])]
  bl11[,`:=`(ofID1 = di1[ofID1], ofID2 = di1[ofID2])]
  bl00 <- subset(bl00, complete.cases(bl00[,1:12]))[,1:12]
  bl01 <- subset(bl01, complete.cases(bl01[,1:12]))[,1:12]
  bl10 <- subset(bl10, complete.cases(bl10[,1:12]))[,1:12]
  bl11 <- subset(bl11, complete.cases(bl11[,1:12]))[,1:12]

  # -- write the blasts
  fwrite(
    bl00, file = file.path(writeDir, "Blast0_0.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl01, file = file.path(writeDir, "Blast0_1.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl10, file = file.path(writeDir, "Blast1_0.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)
  fwrite(
    bl11, file = file.path(writeDir, "Blast1_1.txt.gz"),
    sep = "\t", quote = F, row.names = F, col.names = F)

  # -- run orthofinder
  com <- sprintf("orthofinder -b %s -og -a 1 -t 1  1>/dev/null 2>&1", writeDir)
  system(com)

  # -- find the files
  ogf <- order_filesByMtime(
    path = writeDir,
    recursive = T,
    pattern = "Orthogroups.tsv")[1]

  # -- read the orthogroups.tsv file and process
  ogdt <- fread(ogf, showProgress = F, verbose = F)
  ogdt <- melt(
    ogdt, id.vars = "Orthogroup", variable.name = "genome", value.name = "id")
  ogdt <- ogdt[,list(id = strsplit(id, ",")[[1]]), by = c("Orthogroup", "genome")]
  ogdt[,`:=`(genome = NULL, ofID = trimws(id), id = NULL,
           og = trimws(Orthogroup), Orthogroup = NULL)]
  ogdt <- subset(ogdt, !duplicated(ogdt))
  hasDup <- subset(ogdt, ofID %in% subset(ogdt, duplicated(ofID))$ofID)
  if(nrow(hasDup) > 1){
    m <- merge(hasDup, hasDup, by = "ofID", all = T, allow.cartesian = T)
    ci <- clus_igraph(m$og.x, m$og.y)
    ci <- ci[!duplicated(names(ci))]
    ogdt[,og := ifelse(og %in% names(ci), ci[og], og)]
    ogdt <- subset(ogdt, !duplicated(ogdt))
  }
  ogdt[,og := as.integer(factor(og, unique(og)))]

  nog <- c(id0, id0)
  nog <- nog[!nog %in% ogdt$ofID]

  # -- return data.table of ogs
  if(length(nog) > 0){
    ogdt <- rbind(ogdt, data.table(
      ofID = nog, og = sprintf("NOG_%s",1:length(nog))))
  }

  return(ogdt)
}

#' @title count_nBlks
#' @description
#' \code{count_nBlks} Gcount_nBlks
#' @rdname rerun_orthofinder
#' @import data.table
#' @export
count_nBlks <- function(gsParam, genomeIDs, minGenes4of){
  synp <- data.table(gsParam$params$synteny)
  synp[,`:=`(gnum1 = match(genome1, genomeIDs),
             gnum2 = match(genome2, genomeIDs))]
  synp <- subset(synp, gnum1 <= gnum2)
  setkey(synp, gnum1, gnum2)
  blks <- rbindlist(mclapply(1:nrow(synp), mc.cores = gsParam$params$nCores, function(i){
    fs <- with(synp[i,], file.path(
      gsParam$paths$results,
      sprintf("%s_%s_synHits.txt.gz",genome1, genome2)))
    h <- fread(
      fs, na.strings = c("NA", ""), showProgress = F,
      select = c("regID", "regBuffer","regAnchor","ofID1","ofID2","gen1","gen2","chr1","chr2","start1","start2","end1","end2","ord1","ord2"))
    h <- subset(h, complete.cases(h))
    h <- subset(h, regBuffer)
    h[,`:=`(blkID = regID, regID = NULL)]
    b <- calc_blkCoords(subset(h, regAnchor))
    return(b)
  }))
  if(all(blks$nHits1 >= minGenes4of) && all(blks$nHits2 >= minGenes4of)){
    cat(sprintf("\tAll %s regions have > %s genes\n", nrow(blks), minGenes4of))
  }else{
    bd <- subset(blks, nHits1 < minGenes4of | nHits2 < minGenes4of)
    bk <- subset(blks, nHits1 >= minGenes4of | nHits2 >= minGenes4of)
    cat(sprintf("\tUsing synteny-constrained global OGs for %s/%s regions/hits, with < %s hits\n",
                nrow(bd), sum(bd$nHits1), minGenes4of))
    cat(sprintf("\tRe-building syntenic orthogroups for %s regions with %s hits\n",
                nrow(bk), sum(bk$nHits1)))
  }
  return(blks)
}

