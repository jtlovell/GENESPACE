#' @title Build datasets for genespace pangenome annotation
#' @description
#' \code{build_pgDb} Inferred reference positions of all syntenic
#' orthogroups and decodes these into a pan-genome annotation representation
#' of gene copy number and presence absence variation.
#' @name build_pgDb
#'
#' @param synParam file.path to the directory storing the input orthofinder
#' blast files and orthogroups.tsv. The orthogroups file can be in its
#' original subdirectory. Genesppace will only use the most recently modified
#' occurance of orthogroups.tsv in all subdirectories of blastDir.
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#' @param gsAnnot named character vector of file.paths to gff-like annotation
#' files. The names in this vector must match genome1/genome2 columns in the
#' syntenyParams data.table.
#' @param synOgFiles syntenic orthogroup files generated from pull_synOgs
#' @details ...
#'
#' @note \code{build_pgDb} is a generic name for the functions documented.
#' \cr
#' If called, \code{build_pgDb} returns its own arguments.
#'
#' @title Pipeline to build a pangenome
#' @description
#' \code{build_pgDb} Combines syntenic orthogroups and inferred positions of
#' missing genes into a dataset containing all information needed to print
#' a pan-genome annotationo
#' @rdname build_pgDb
#' @import data.table
#' @export
build_pgDb <- function(maxDistBtwHits,
                       gsParam,
                       gsAnnot,
                       synParam,
                       synOgFiles){

  ##############################################################################
  # 1. load annotations and other information
  ##############################################################################
  if(gsParam$verbose)
    cat("Building syntenic pangenome annotation graph ...\n\tLoading annotations and orthogroups ... ")
  genomeIDs <- gsParam$genomeIDs[!gsParam$genomeIDs %in% gsParam$outgroup]
  gffAll <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), blastDir = gsParam$blast)

  ##############################################################################
  # 2. read in the syntenic hits, subset to only the same orthogroup
  ##############################################################################
  if(gsParam$verbose)
    cat("Done!\n\tBuilding syntenic orthogroup subgraphs ... ")
  sogv <- extract_synOGs(synParam = synParam, gsParam = gsParam)

  ##############################################################################
  # 3. Add in reciprocal best hits into the network
  ##############################################################################
  if(gsParam$verbose)
    cat("Done!\n\tAdding reciprocal best hits to subgraphs ... ")
  rbhs <- extract_rbhFromSynOgFiles(
    gsParam = gsParam,
    synParam = synParam,
    synOgv = sogv)
  ov <- gffAll$ord; cv <- gffAll$chr; gv <- gffAll$genome; idv <- gffAll$id
  bpv <- with(gffAll, (start + end)/2)
  names(ov) <- names(cv) <- names(gv) <- names(idv) <- names(bpv) <- gffAll$ofID
  rbhs[,`:=`(genome = gv[gsub("*","",ofID, fixed = T)],
             chr = cv[gsub("*","",ofID, fixed = T)],
             ord = ov[gsub("*","",ofID, fixed = T)],
             bp = bpv[gsub("*","",ofID, fixed = T)],
             id = sprintf("%s*",idv[gsub("*","",ofID, fixed = T)]))]
  ov <- cv <- gv <- idv <- bpv <- NULL

  ##############################################################################
  # 4. summarize inferred positions of missing regions
  ##############################################################################
  if(gsParam$verbose)
    cat("Done!\n\tAdding inferred positions of missing orthogroup members ... ")
  missNoRBH <- combine_inferredPos(
    synOgFiles = synOgFiles,
    rbhs = rbhs,
    gsParam = gsParam,
    synOgv = sogv,
    maxDistBtwHits = maxDistBtwHits)

  ##############################################################################
  # 5. get physical position of absences
  ##############################################################################
  if(gsParam$verbose)
    cat("Done!\n\tConverting gene-rank order to physical positions ... ")
  spl <- split(missNoRBH, by = c("genome","chr"))
  splg <- split(gffAll, by = c("genome","chr"))
  missNoRBH <- rbindlist(lapply(names(spl), function(i){
    x <- spl[[i]]
    y <- splg[[i]]
    z <- rbind(with(x, data.table(ofID = ofID, ord = ord, bp = NA)),
               with(y, data.table(ofID = ofID, ord = ord, bp = (start + end)/2)))
    setkey(z, ord)
    z[,bp := round(infer_pos(ord1 = bp, ord2 = ord),0)]
    x <- merge(x, z[,c("ofID","bp")], all.x = T, by = "ofID", allow.cartesian = T)
    x[,id := sprintf("%s:%s", chr, bp)]
    x <- subset(x, !duplicated(x[,c("ofID","og","genome")]))
    return(x)
  }))
  spl <- splg <- NULL

  ##############################################################################
  # 6. Combine rbh, ogs and absence positions into pangenome
  ##############################################################################
  pangenomeFile <- file.path(gsParam$results, "genespace_pangenomeLong.txt.gz")
  if(gsParam$verbose)
    cat("Done!\n\tDone! .. Pangenome written to ... ")
  a <- gffAll[,c("ofID","genome","chr","ord","id")]
  a[,`:=`(og =  sogv[ofID], bp = (gffAll$start + gffAll$end)/2,
          genome = factor(genome, levels = genomeIDs))]
  a <- subset(a, !is.na(og)); setkey(a, genome, ord)
  oglev <- unique(a$og)
  a[,`:=`(genome = as.character(genome))]
  pg <- rbind(
    a[,c("og","id","genome","chr","bp")],
    rbhs[,c("og","id","genome","chr","bp")],
    missNoRBH[,c("og","id","genome","chr","bp")])
  a <- missNoRBH <- NULL
  pg[,og := factor(og, levels = oglev)]
  setkey(pg, og, genome, bp)
  pg[,og := as.character(og)]
  if(gsParam$verbose)
    cat(sprintf("/results/%s\n", basename(pangenomeFile)))
  fwrite(pg, sep = "\t", quote = F, file = pangenomeFile)

  ##############################################################################
  # 7. Pull and write out non-syntenic orthogroups
  ##############################################################################
  if(gsParam$verbose)
    cat("Pulling non-syntenic orthogroup members by genome \n")
  gv <- gffAll$genome; idv <- gffAll$id; names(idv) <- names(gv) <- gffAll$ofID
  nonSynOgFiles <- sapply(genomeIDs, function(i){
    if(gsParam$verbose)
      cat(sprintf("\t%s: ", i))
    ns <- pull_nonSynOGs(
      gsParam = gsParam,
      rbhs = rbhs,
      synParam = synParam,
      refGenome = i)
    ns <- ns[,list(nsog = unlist(nsOgGenes)), by = "ofID"]
    ns[,`:=`(refGenome = gv[ofID], id = idv[ofID], ofID = NULL,
             genome = gv[nsog], nonSynOgID = idv[nsog], nsog = NULL)]
    ns[,id := factor(id, levels = idv)]
    setkey(ns, id)
    ns[,id := as.character(id)]
    fileOut <- file.path(gsParam$results, sprintf("%s_nonSynOgs.txt.gz", i))
    fwrite(ns, sep = "\t", quote = F, file = fileOut)
    if(gsParam$verbose)
      with(ns, cat(sprintf("%s non-syntenic og hits --> /results/%s\n",
                           uniqueN(paste(genome, nonSynOgID)), basename(fileOut))))
    return(fileOut)
  })

  # setnames(ns, "ofID", "pgID")
  # pgw <- merge(pgw, ns, by = "pgID", all.x = T)

  if(gsParam$verbose)
    cat("\tDone!")
  return(c(pangenomeFile = pangenomeFile, nonSynOgFiles = nonSynOgFiles))
}

#' @title Pipeline to build a pangenome
#' @description
#' \code{build_pgDb} Combines syntenic orthogroups and inferred positions of
#' missing genes into a dataset containing all information needed to print
#' a pan-genome annotationo
#' @rdname build_pgDb
#' @import data.table
#' @importFrom parallel mclapply
#' @export
extract_synOGs <- function(synParam, gsParam){
  sf <- synParam$synHitFiles
  synOg <- rbindlist(mclapply(sf, mc.cores = gsParam$nCores, function(i)
    subset(fread(
      i,
      select = c("ofID1","ofID2","og1","og2")),
      og1 == og2)[,c("ofID1","ofID2")]))

  synOg[,og := clus_igraph(ofID1, ofID2)]
  synOg <- synOg[,list(u = unique(c(ofID1, ofID2))), by = "og"]
  sogv <- synOg$og; names(sogv) <- synOg$u
  return(sogv)
}

#' @title Pipeline to build a pangenome
#' @description
#' \code{build_pgDb} Combines syntenic orthogroups and inferred positions of
#' missing genes into a dataset containing all information needed to print
#' a pan-genome annotationo
#' @rdname build_pgDb
#' @import data.table
#' @importFrom parallel mclapply
#' @export
extract_rbhFromSynOgFiles <- function(synParam, gsParam, synOgv){
  if(class(synOgFiles) == "character"){
    sf <- synOgFiles
  }else{
    sf <- synParam$synHitFiles
  }

  rbhs <- rbindlist(mclapply(sf, mc.cores = gsParam$nCores, function(i){
    x <- fread(i, select = 1:2)
    x <- subset(x, grepl("RBH$", ofID2))
    x[,`:=`(og = synOgv[ofID1], ofID1 = NULL,
            ofID = gsub("RBH","",ofID2), ofID2 = NULL)]
    return(x)
  }))
  rbhs[,`:=`(ofID = paste0(ofID,"*"))]
  return(rbhs)
}

#' @title Pipeline to build a pangenome
#' @description
#' \code{build_pgDb} Combines syntenic orthogroups and inferred positions of
#' missing genes into a dataset containing all information needed to print
#' a pan-genome annotationo
#' @rdname build_pgDb
#' @import data.table
#' @importFrom parallel mclapply
#' @importFrom dbscan dbscan frNN
#' @export
combine_inferredPos <- function(synOgFiles, rbhs, maxDistBtwHits, gsParam, synOgv){
  d <- rbindlist(mclapply(synOgFiles, mc.cores = gsParam$nCores, function(i){
    x <- fread(i, select = 1:2)
    x1 <- subset(x, grepl(":", ofID1, fixed = T))
    x2 <- subset(x, grepl(":", ofID2, fixed = T))
    if(nrow(x1) == 0 & nrow(x2) == 0){
      return(NULL)
    }else{
      x <- rbind(x1, with(x2, data.table(ofID1 = ofID2, ofID2 = ofID1)))
      x[,og := synOgv[ofID2]]
      x <- subset(x, !is.na(og))
      x <- subset(with(subset(x, grepl(":", paste(ofID1, ofID2))),
                       data.table(og = c(og, og), loc = c(ofID1, ofID2))),
                  grepl(":",loc))
      x[,c("genome","chr","blk","ord","index"):=tstrsplit(loc,":")]
      x[,`:=`(index = NULL, loc = NULL, blk = NULL, ord = as.numeric(ord))]
      return(x)
    }
  }))
  dm <- merge(d, with(rbhs, data.table(
    og = og, genome = genome, chr = chr, rbhOrd = ord)),
    by = c("og","genome","chr"), allow.cartesian = T)
  suppressWarnings(dm[,minDiff := min(abs(rbhOrd - ord)), by = c("og","genome","chr")])

  closRbh <- with(subset(dm, minDiff <= maxDistBtwHits), paste(og, genome, chr, ord))
  missNoRBH <- subset(d, !paste(og, genome, chr, ord) %in% closRbh)
  missNoRBH[,rng := diff(range(ord)), by = c("og", "genome","chr")]
  missNoRBH[,dbs := 1]
  xs <- subset(missNoRBH, rng > maxDistBtwHits)
  xr <- subset(missNoRBH, rng <= maxDistBtwHits)
  if(nrow(xs) > 1){
    xs[,dbs := dbscan(frNN(cbind(
      ord, ord), eps = maxDistBtwHits), minPts = 1)$cluster,
      by = c("og", "genome","chr")]
  }

  missNoRBH <- rbind(xs, xr)
  missNoRBH <- missNoRBH[,list(ord = median(ord)), by = c("og", "genome","chr","dbs")]
  missNoRBH[,`:=`(ord = round(ord, 1), dbs = NULL, ofID = sprintf("%s:%s",chr, ord))]

  return(missNoRBH)
}

#' @title Pull non-syntenic orthogroup members
#' @description
#' \code{pull_nonSynOGs} Parses non-syntenic orthogroup members against a
#' reference genome and places these in the context of the pan-genome
#' annotation.
#' @rdname build_pgDb
#' @import data.table
#' @importFrom parallel mclapply
#' @export
pull_nonSynOGs <- function(gsParam, synParam, refGenome, rbhs){
  spid <- read_orthofinderSpeciesIDs(gsParam$blast)
  ogv <- pull_orthogroups(gsParam$blast)$ogv
  rbhsv <- with(rbhs, data.table(og = og, rbhofID = gsub("*", "", ofID, fixed = T)))
  wh <- which(synParam$genome1 == refGenome | synParam$genome2 == refGenome)
  nsOg <- rbindlist(mclapply(wh, mc.cores = gsParam$nCores, function(i){
    x <- read_blast(
      path = gsParam$blast,
      ofID1 = spid[synParam$genome1[i]],
      ofID2 = spid[synParam$genome2[i]])[,`:=`(og1 = ogv[ofID1], og2 = ogv[ofID2])]
    x <- subset(x, og1 == og2)[,c("ofID1","ofID2")]
    y <- fread(file.path(gsParam$synteny, sprintf(
      "%s_vs_%s_synHits.txt.gz", synParam$genome1[i], synParam$genome2[i])),
      select = 1:2)
    u <- c(paste(y$ofID1, y$ofID2), paste(y$ofID2, y$ofID1))
    u <- u[!duplicated(u)]
    x <- subset(x, !paste(ofID1, ofID2) %in% u)
    if(synParam$genome2[i] == refGenome)
      x <- with(x, data.table(ofID1 = ofID2, ofID2 = ofID1))
    setnames(x, "ofID1", "ofID")
    return(x)
  }))
  nsOg <- nsOg[,list(nsOgGenes = list(unique(ofID2))), by = "ofID"]
  return(nsOg)
}
