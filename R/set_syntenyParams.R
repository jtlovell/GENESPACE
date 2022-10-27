#' @title Set synteny parameters and parse required data
#' @description
#' \code{set_syntenyParams} Generate all data needed to run synteny. This
#' includes the synteny parameters, combined/annotated bed file and annotated
#' blast files.
#' @name set_syntenyParams
#'
#' @param gsParam A list of genespace parameters. This should be created
#' by init_genespace.
#' @param overwrite logical, should the combBed.txt and synHits be re-generated
#' if they exist. Default is FALSE.
#'
#' \cr
#' If called, \code{set_syntenyParams} returns its own arguments.
#'
#' @details set_synParam can be called directly to generate all data needed to
#' run synteny. This MUST be run before synteny, to ensure that all data are
#' in the right format and ready to go. The first step follows the functionality
#' in < v1, generating a matrix of query v. target genomes and populates that
#' with all needed parameters. The second step follows what was `annotate_gff`
#' in < v1, and is now called `annotate_bed`. In each case, this function
#' (detailed below) is designed to be used internally, but can be called by the
#' user if desired. The third step follows the first step of `synteny` in < v1,
#' where the blast files are read in and populated with the required positional
#' etc. info from the combined bed file. This functionality is also offered as
#' an exported function `annotate_blast`, but should just be called within
#' `set_syntenyParam`.
#'
#' The combined bed file is annotated with a variety of additional
#' information. These data are required for synteny and pangenome functions. The
#' biggest steps are 1) finding physically dispersed orthogroups (those which
#' hit lots of places across the genome) and flagging these so they won't be
#' used in synteny detection, 2) determining whether you should use phylo.
#' hierarchical orthogroups (HOGs) or the traditional orthogroup.tsv (OGs)
#' methods, and 3) determining membership and representatives of tandem arrays,
#' which here are defined as physically proximate orthogroup members. The
#' function is verbose and prints several sets of data that may be useful to
#' contextualize your run and ensure that you are using sufficiently high
#' quality genomes and annotations.
#'
#' The combined bed file contains the first four bed columns (chr, start, end,
#' name), and the following (columns 4:16):
#' \enumerate{
#' \item genome unique genomeID, taken from the individual bed files
#' \item ord integer with the order of genes
#' \item unique orthofinder ID, taken from 'SequenceIDs.txt'
#' \item pepLen: amino acid length of the gene CDS
#' \item arrayID: tandem array ID, unique for each genome, og and chr.
#' \item isArayRep: logical flag whether a gene is the representative for array
#' \item globOG: global orthogroup id (from 'Orthogroups.tsv')
#' \item globHOG: global phylogenetically hierarchical orthogroup id (N0.tsv)
#' \item synOG: og, split by synteny
#' \item inblkOG: synOG, merged across within-block orthogroups
#' \item noAnchor: logical flag whether this gene can ever be a syntenic anchor
#' \item og: either globOG or globHOG, depending on if useHOG = FALSE or TRUE
#' }
#'
#'
#'

#' @title set_syntenyParams
#' @description
#' \code{set_syntenyParams} set_syntenyParams
#' @rdname set_syntenyParams
#' @import data.table
#' @export
set_syntenyParams <- function(gsParam, overwrite){

  # -- if there is not combined bed file, make it
  bedFile <- file.path(gsParam$paths$results, "combBed.txt")
  if(file.exists(bedFile) & !overwrite)
    cat("**NOTE** annotated/combined bed file exists, not overwriting\n")
  if(file.exists(bedFile) & overwrite)
    cat(strwrap(
      "**NOTE** annotated/combined bed file exists but overwrite = TRUE,
      re-making it", indent = 0, exdent = 8), sep = "\n")
  if(!file.exists(bedFile))
    cat("combining and annotating the bed files\n")
  if((file.exists(bedFile) & overwrite) || !file.exists(bedFile)){
    tmp <- annotate_bed(gsParam = gsParam)
    tmp <- NULL
  }

  # if there are no blast md in the gsParam
  if(!"annotBlastMd" %in% names(gsParam))
    cat("**NOTE** could not find annotated blast files, making them.\n")
  if("annotBlastMd" %in% names(gsParam) && overwrite)
    cat("**NOTE** annotated blast files exist, but overwrite = TRUE,
      re-making them")
  if(("annotBlastMd" %in% names(gsParam) && overwrite) ||
    !"annotBlastMd" %in% names(gsParam))
    gsParam <- annotate_blast(gsParam = gsParam)

  # -- add ploidy
  ploidy1 <- ploidy2 <- query <- target <- synBuff <- NULL
  blMd <- data.table(gsParam$annotBlastMd)
  blMd[,`:=`(ploidy1 = gsParam$ploidy[query],
             ploidy2 = gsParam$ploidy[target])]

  # -- add all parameters
  pnames <- c(
    "orthofinderInBlk", "blkSize", "nGaps", "synBuff", "onlyOgAnchors",
    "nSecondaryHits", "blkSizeSecond", "nGapsSecond", "onlyOgAnchorsSecond")
  for(i in pnames)
    blMd[[i]] <- gsParam$params[[i]]

  blMd <- data.table(blMd)
  blMd[, `:=`(synRad = sqrt(2) * synBuff,
              inBufferRadius = synBuff/2)]

  if(!is.na(gsParam$outgroup))
    blMd <- subset(blMd, !(query %in% gsParam$outgroup | target %in% gsParam$outgroup))

  gsParam$annotBlastMd <- blMd
  return(gsParam)
}

#' @title annotate_bed
#' @description
#' \code{annotate_bed} annotate_bed
#' @rdname set_syntenyParams
#' @import data.table
#' @importFrom stats median
#' @export
annotate_bed <- function(gsParam){

  add_ofID2bed <- function(bed, gsParam){
    genomeNum <- genome <- id <- ofID <- NULL
    seqid <- read_orthofinderSequenceIDs(
      file.path(gsParam$paths$results, "SequenceIDs.txt"))
    speid <- read_orthofinderSpeciesIDs(
      file.path(gsParam$paths$results, "SpeciesIDs.txt"))
    seqid[,genome := names(speid)[match(genomeNum, speid)]]
    di <- seqid$ofID; names(di) <- with(seqid, paste(genome, id))
    bed[,ofID := di[paste(genome, id)]]
    return(bed)
  }

  add_pepLen2bed <- function(bed,
                             gsParam){
    id <- pepLen <- NULL
    genomeIDs <- unique(bed$genome)
    pepFiles <- file.path(gsParam$paths$peptide, sprintf("%s.fa", genomeIDs))
    pepFiles <- pepFiles[file.exists(pepFiles)]
    if(length(pepFiles) > 1){
      naa <- character()
      for(i in pepFiles){
        gid <- gsub(".fa$", "", basename(i))
        tmp <- get_nAA(i)
        names(tmp) <- paste(gid, names(tmp))
        naa <- c(naa, tmp)
      }
      bed[,pepLen := naa[paste(genome, id)]]
    }else{
      bed[,pepLen := NA]
    }
    return(bed)
  }

  add_ogs2gff <- function(gsParam, gff){
    ogs <- parse_ogs(gsParam)
    gff <- merge(
      gff, ogs,
      by = c("genome","id"), all.x = T)
    gff$ogID[is.na(gff$ogID)] <- paste0("NOG",1:sum(is.na(gff$ogID)))
    setnames(gff, "ogID", "globOG")
    return(gff)
  }

  add_og2bed <- function(bed, gsParam){
    genome <- id <- NULL
    ogPath <- file.path(gsParam$paths$results, "Orthogroups.tsv")
    ogs <- parse_ogs(path = ogPath, genomeIDs = genomeIDs)
    di <- ogs$ogID; names(di) <- with(ogs, paste(genome, id))
    bed[,globOG := di[paste(genome, id)]]
    wh <- which(is.na(bed$globOG))
    lab <- gsub(" ", "0", align_charRight(c("00000", 1:length(wh))))[-1]
    bed$globOG[wh] <- sprintf("NoOG%s", lab)
    return(bed)
  }

  add_hog2bed <- function(bed, gsParam){
    genome <- id <- NULL
    ogPath <- file.path(gsParam$paths$results, "N0.tsv")
    ogs <- parse_hogs(path = ogPath, genomeIDs = genomeIDs)
    di <- ogs$hogID; names(di) <- with(ogs, paste(genome, id))
    bed[,globHOG := di[paste(genome, id)]]
    wh <- which(is.na(bed$globHOG))
    lab <- gsub(" ", "0", align_charRight(c("00000", 1:length(wh))))[-1]
    bed$globHOG[wh] <- sprintf("N0.NoHOG%s", lab)
    return(bed)
  }

  add_nOGplaces2bed <- function(bed, arrayJump){
    genome <- chr <- og <- ord <- jumpLeft <- NULL
    setkey(bed, genome, chr, og, ord)
    bed[,jumpLeft := c(arrayJump + 1, diff(ord)), by = c("genome", "chr", "og")]
    bed[,nOGPlaces := sum(jumpLeft > arrayJump), by = c("genome","og")]
    bed[,`:=`(jumpLeft = NULL)]
    setkey(bed, genome, ord)
    return(bed)
  }

  add_array2bed <- function(bed, maxPlaces, synBuff, maxIter = 10){

    # -- we can make arrays for everything except that we want to exclude the
    # arrays that are huge and problematic
    genome <- chr <- id <- arrayID <- nOGPlaces <- tord <- NULL
    tmp <- subset(bed, nOGPlaces <= maxPlaces)
    tmp[,tord := as.numeric(ord)]

    # -- set up the iteration
    cnt <- 1
    diffn <- 1
    tmp[,arrayID := sprintf(
      "tmp%s", as.integer(as.factor(paste(genome, chr, id))))]
    while(cnt <= maxIter && diffn > 0){
      # -- for each iteration, calculate clusters by the size of jump between
      # genes larger than the synBuffer
      cnt <- cnt + 1
      initn <- uniqueN(tmp$arrayID)
      genome <- chr <- og <- ord <- jumpLeft <- clus <- n <- NULL
      setkey(tmp, genome, chr, og, tord)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      tmp[,jumpLeft := c(synBuff + 1, diff(tord)), by = c("genome", "chr", "og")]
      tmp[,clus := as.integer(jumpLeft > synBuff), by = c("genome", "chr")]
      tmp[,clus := cumsum(clus), by = c("genome", "chr", "og")]
      tmp[,`:=`(
        arrayID = sprintf("tmp%s",
                          as.integer(as.factor(paste(genome, chr, og, clus)))),
        jumpLeft = NULL, clus = NULL)]

      # -- get the new order of genes based on array ID
      tmp[,n := .N, by = "arrayID"]
      tmp1 <- subset(tmp, n == 1)
      tmp2 <- subset(tmp, n > 1)
      tmp2[,tord := as.numeric(tord)]
      tmp2[,tord := mean(as.numeric(tord)), by = "arrayID"]
      tmp <- rbind(tmp1, tmp2)
      tmp[,tord := frank(tord, ties.method = "dense"), by = "genome"]
      newn <- uniqueN(tmp$arrayID)
      diffn <- initn - newn
    }

    # -- relabel
    ogID <- NULL
    lab <- gsub(" ", "0",
                align_charRight(
                  as.numeric(factor(tmp$arrayID,
                                    levels = unique(tmp$arrayID)))))
    arrv <- sprintf("Arr%s", lab); names(arrv) <- tmp$ofID
    bed[,arrayID := arrv[ofID]]

    tmp <- subset(bed, is.na(arrayID))
    wh <- with(tmp,
               as.numeric(as.factor(paste(genome, chr, og))))
    lab <- gsub(" ", "0", align_charRight(wh))
    wh <- which(is.na(bed$arrayID))
    bed$arrayID[wh] <- sprintf("NoArr%s", lab)
    return(bed)
  }

  add_arrayReps2bed <- function(bed){
    n <- ord <- medOrd <- medDiff <- pepLen <- arrayID <- isRep <-
      isArrayRep <- ofID <- NULL
    bed[,n := .N, by = "arrayID"]
    tmp <- subset(bed, n > 1)
    tmp[,medOrd := as.numeric(median(ord)), by = "arrayID"]
    tmp[,medDiff := as.numeric(abs(medOrd - ord))]
    setorder(tmp, arrayID, medDiff, -pepLen)
    tmp$noAnchor[duplicated(tmp$arrayID)] <- TRUE
    tmp[,isRep := !duplicated(arrayID)]
    bed[,isArrayRep := ofID %in% tmp$ofID[tmp$isRep] | n == 1]
    tmp <- subset(bed, isArrayRep)
    tmp[,ord := frank(ord, ties.method = "dense"), by = "genome"]
    di <- tmp$ord; names(di) <- tmp$ofID
    return(bed)
  }

  check_ploidyOgMatch <- function(ploidy, bed){
    genome <- nGenome <- nGenes <- ofID <- n <- nOGs <- smPloid <- NULL
    cat("Checking which type of orthogroup most closely matches ploidy\n")
    tmp <- data.table(genome = names(ploidy), ploidy = ploidy)
    cat("Checking global orthogroups (OGs) ... \n")
    bed[,nGenome := uniqueN(genome), by = "globOG"]
    cat(sprintf("\t%s OGs (of %s) with %s genes (of %s) hit all genomes\n",
                uniqueN(bed$globOG[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$globOG),
                uniqueN(bed$ofID[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$ofID)))
    cnts <- subset(bed, nGenome == nrow(tmp))
    cnts <- cnts[,list(nGenes = uniqueN(ofID)), by = c("globOG", "genome")]
    cnts <- cnts[,list(n = .N), by = c("genome", "nGenes")]
    ogCnt <- merge(tmp, cnts, by = "genome", all = T)
    ogCnt[,smPloid := ploidy == nGenes]
    ogCnt <- ogCnt[,list(nOGs = sum(n), type = "OG"), by = c("genome", "smPloid")]

    cat("Checking phylogenetic hierarchical orthogroups (HOGs) ... \n")
    bed[,nGenome := uniqueN(genome), by = "globHOG"]
    cat(sprintf("\t%s HOGs (of %s) with %s genes (of %s) hit all genomes\n",
                uniqueN(bed$globHOG[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$globHOG),
                uniqueN(bed$ofID[bed$nGenome == nrow(tmp)]),
                uniqueN(bed$ofID)))
    cnts <- subset(bed, nGenome == nrow(tmp))
    cnts <- cnts[,list(nGenes = uniqueN(ofID)), by = c("globHOG", "genome")]
    cnts <- cnts[,list(n = .N), by = c("genome", "nGenes")]
    hogCnt <- merge(tmp, cnts, by = "genome", all = T)
    hogCnt[,smPloid := ploidy == nGenes]
    hogCnt <- hogCnt[,list(nOGs = sum(n), type = "HOG"), by = c("genome", "smPloid")]
    cnts <- rbind(ogCnt, hogCnt)
    tp <- dcast(cnts, genome ~ type + smPloid, value.var = "nOGs")
    tp <- rbind(tp,
                matrix(c(genome = "total",
                         colSums(tp[,-1,with = F])),nrow =  1), use.names = F)

    pOG <- pHOG <- genome <- HOG <- OG <- HOG_TRUE <- HOG_FALSE <- OG_TRUE <-
      OG_FALSE <- NULL
    tp[,`:=`(pHOG = as.numeric(HOG_TRUE)/(as.numeric(HOG_FALSE) + as.numeric(HOG_TRUE)),
             pOG = as.numeric(OG_TRUE)/(as.numeric(OG_FALSE) + as.numeric(OG_TRUE)))]
    tp[,`:=`(pHOG = round(pHOG*100, 2), pOG = round(pOG*100,2))]

    md <- tp[,list(genome = align_charLeft(c("", genome)),
                   HOG = align_charLeft(c("match, not, % (HOGs)",
                                          sprintf("%s, %s, %s", HOG_TRUE, HOG_FALSE, pHOG))),
                   OG = align_charLeft(c("match, not, % (OGs)",
                                         sprintf("%s, %s, %s", OG_TRUE, OG_FALSE, pOG))))]
    cat("##############################################################################\n")
    nu <- apply(md, 1, function(x){
      x <- unlist(x)
      cat(sprintf("%s: %s || %s\n", x[1], x[2], x[3]))
    })
    cat("##############################################################################\n")
    tmp <- subset(tp, genome == "total")
    if(tmp$pHOG > tmp$pOG){
      useHOGs <- TRUE
      cat("More HOGs match patterns of ploidy than OGs. Setting useHOGs = TRUE")
    }else{
      useHOGs <- FALSE
      cat("More OGs match patterns of ploidy than HOGs. Setting useHOGs = FALSE")
    }
    return(useHOGs)
  }

  ##############################################################################
  # 0. get the parameters in order
  resDir <- gsParam$paths$results

  genomeIDs <- gsParam$genomeIDs
  if(!is.na(gsParam$outgroup))
    genomeIDs <- genomeIDs[genomeIDs %in% gsParam$outgroup]
  useHOGs <- gsParam$params$useHOGs
  nGaps <- gsParam$params$nGaps
  synBuff <- gsParam$params$synBuff
  blkSize <- gsParam$params$blkSize
  arrayJump <- gsParam$params$arrayJump
  maxOgPlaces <- gsParam$params$maxOgPlaces
  maxOgPlaces <- gsParam$ploidy * 8
  arrayJump <- gsParam$params$synBuff / 2
  bed <- NULL

  ##############################################################################
  # 2. combine bed files
  # -- 2.1 get the file paths
  bedDir <- gsParam$paths$bed
  bedFiles <- file.path(bedDir, sprintf("%s.bed", genomeIDs))
  if(!all(file.exists(bedFiles)))
    stop(sprintf("could not find all bed files in %s\n", bedDir))

  # -- 2.2 read in the bed files
  bed <- rbindlist(lapply(bedFiles, function(i){
    tmp <- fread(
      i, col.names = c("chr", "start", "end", "id"),
      colClasses = c("character", "integer", "integer", "character"),
      showProgress = FALSE, verbose = FALSE)
    tmp[,genome := gsub(".bed$", "", basename(i))]
    return(tmp)
  }))

  # -- 2.3 add in gene orders
  chrn <- n <- chr <- ord <- NULL
  bed[,chrn := as.numeric(gsub('\\D+','', chr))]
  bed[,n := .N, by = c("genome", "chr")]
  setorder(bed, genome, chrn, chr, -n, start)
  bed[,ord := 1:.N, by = "genome"]

  # -- 2.4 add in all other columns
  noAnchor <- arrayID <- chr <- n <- ofID <- pepLen <- arrayID <- isArrayRep <-
    globOG <- globHOG <- synOG <- inblkOG <- NULL
  bed[,`:=`(
    chrn = NULL, n = NULL, ofID = NA, pepLen = NA, arrayID = NA,
    isArrayRep = NA,  globOG = NA, globHOG = NA, synOG = NA,
    inblkOG = NA, noAnchor = NA)]

  ##############################################################################
  # 3. Add in necesary data columns if they do not exist
  # -- 3.1 orthofinder IDs
  bed <- add_ofID2bed(bed = bed, gsParam = gsParam)

  # -- 3.2 peptide lengths
  bed <- add_pepLen2bed(bed = bed, gsParam = gsParam)

  # -- 3.3 orthogroup IDs
  bed <- add_og2bed(bed = bed, gsParam = gsParam)

  # -- 3.4 hier orthogroup IDs
  bed <- add_hog2bed(bed = bed, gsParam = gsParam)

  ##############################################################################
  # 4. Arrays and multi-placed orthogroups
  if(is.null(useHOGs) || is.na(useHOGs))
    useHOGs <- check_ploidyOgMatch(bed = bed, ploidy = gsParam$ploidy)

  og <- globHOG <- globOG <- genome <- nOGPlaces <- isArrayRep <- tmp <-
    n <- maxPlaces <- NULL
  if(useHOGs){
    bed[,og := globHOG]
  }else{
    bed[,og := globOG]
  }

  # -- 4.1 n. orthogroup placements
  bed <- add_nOGplaces2bed(bed = bed, arrayJump = arrayJump)
  bed[,maxPlaces := maxOgPlaces[genome]]

  # -- 4.2 find the arrays
  bed <- add_array2bed(
    bed = bed, maxPlaces = maxPlaces, synBuff = synBuff, maxIter = 10)

  # -- 4.3 choose the array representative genes
  bed <- add_arrayReps2bed(bed)

  # -- 4.4 determine un-anchorable genes
  bed[,tmp := uniqueN(og), by = c("genome", "chr")]
  bed[,noAnchor := tmp < blkSize | nOGPlaces > maxPlaces | !isArrayRep]
  bed[,`:=`(tmp = NULL, n = NULL)]

  ##############################################################################
  # 5. Clean up and summarize
  n <- nchr <- nGenes <- nHOG <- og <- ofID <- nGenesInArray <- nOGs <-
    nGenesSmallScaff <- nSmallScaff <- nDispArray <- genome <- narr <-
    nGenome <- NULL
  cat(strwrap(
    "\n\nBelow, printing a table with QC information for each genome. Abbreviations:
  'nArr' = n. tandem arrays, 'smChrs' = scaffolds smaller than `2 * blkSize`,
  'dispOG' = over-dispersed OG, by default, hitting > `8 * ploidy` unique
  places in the genome. If there are more than ~5-10% of genes in either of the
  last two columns, you may have an issue with your genome or annotation.",
    indent = 0, exdent = 0), sep = "\n")
  bed[,n := uniqueN(ofID), by = c("genome", "arrayID")]
  bed[,nchr := uniqueN(og), by = c("genome", "chr")]
  md <- bed[,list(nGenes = uniqueN(ofID), nOGs = uniqueN(globOG),
                  nHOG = uniqueN(globHOG), narr = uniqueN(arrayID[n > 1]),
                  nGenesInArray = uniqueN(ofID[n > 1]),
                  nGenesSmallScaff = sum(nchr <= blkSize * 2),
                  nSmallScaff = uniqueN(chr[nchr <= blkSize * 2]),
                  nDispArray = sum(nOGPlaces > maxPlaces)), by = "genome"]
  md <- md[,list(genome = align_charLeft(c("", genome)),
                 nGenes = align_charLeft(c("nGenes",nGenes)),
                 nOGs = align_charLeft(c("nOG", nOGs)),
                 nHOG = align_charLeft(c("nHOG", nHOG)),
                 narr = align_charLeft(c("nArr", narr)),
                 nGenesInArray = align_charLeft(c("genes", nGenesInArray)),
                 nGenesSmallScaff = align_charLeft(c("genes", nGenesSmallScaff)),
                 nSmallScaff = align_charLeft(c("smChrs",nSmallScaff)),
                 nDispArray = align_charLeft(c("dispOG", nDispArray)))]
  cat("##############################################################################\n")
  nu <- apply(md, 1, function(x){
    x <- unlist(x)
    cat(sprintf("%s: %s (%s / %s) || %s (%s) || %s (%s) || %s\n",
                x[1], x[2], x[3], x[4], x[5], x[6], x[8], x[7], x[9]))
  })

  bed[,`:=`(
    nGenome = NULL, nOGPlaces = NULL, maxPlaces = NULL, n = NULL, nchr = NULL)]
  write_combBed(
    x = bed, filepath = file.path(gsParam$paths$results, "combBed.txt"))
  return(bed)
}

#' @title annotate_blast
#' @description
#' \code{annotate_blast} annotate_blast
#' @rdname set_syntenyParams
#' @import data.table
#' @importFrom grDevices pdf dev.off rgb
#' @export
annotate_blast <- function(gsParam){
  ##############################################################################
  # 1. Build the blast file metadata

  # -- 1.1 read in the combined bed
  bed <- read_combBed(file.path(gsParam$paths$results, "combBed.txt"))


  # -- 1.2 get the blast files together
  genome1 <- genome2 <- n1 <- n2 <- blastFile <- query <- target <- hasBlast <-
    type <- NULL
  tmp <- gsParam$ofFiles$blast
  if(is.null(tmp)){
    gsParam$ofFiles <- find_gsResults(
      resultsDir = gsParam$paths$results,
      genomeIDs = gsParam$genomeIDs)
    tmp <- gsParam$ofFiles$blast
  }

  nseqs <- table(bed$genome)
  tmp[,`:=`(n1 = nseqs[genome1], n2 = nseqs[genome2])]
  tmp[,query := ifelse(n1 >= n2, genome1, genome2)]
  tmp[,target := ifelse(query == genome2, genome1, genome2)]
  tmp[,hasBlast := !is.na(blastFile)]
  setorder(tmp, -hasBlast, query, target)
  tmp[,type := c("queryBlast", "targetBlast")[1:.N], by = c("query", "target")]
  synMd <- dcast(tmp, query + target ~ type, value.var = "blastFile")

  ##############################################################################
  # 2. Process each set
  # -- for each unique set in synMd
  blNames <- c(
    "ofID1", "ofID2", "pid", "length", "mismatches","gapopenings",
    "queryStart", "queryEnd", "subjectStart", "subjectEnd", "Evalue", "bitScore")
  blNamesR <- c(
    "ofID2", "ofID1", "pid", "length", "mismatches","gapopenings",
    "subjectStart", "subjectEnd", "queryStart", "queryEnd", "Evalue", "bitScore")
  bdNames <- c(
    "chr", "start", "end", "id", "ofID", "ord", "genome", "og", "noAnchor", "isArrayRep")
  bed1 <- data.table(bed[,bdNames, with = F])
  setnames(bed1, paste0(names(bed1), "1"))
  bed2 <- data.table(bed[,bdNames, with = F])
  setnames(bed2, paste0(names(bed2), "2"))
  annotBlastFile <- NULL
  synMd[,annotBlastFile := NA]
  for(i in 1:nrow(synMd)){

    # -- 2.1 read in the blast file(s)
    bl <- fread(
      synMd$queryBlast[i],
      showProgress = F,
      na.strings = c("", "NA"),
      col.names = blNamesR)
    if(!is.na(synMd$targetBlast[i])){
      blr <- fread(
        synMd$targetBlast[i],
        showProgress = F,
        na.strings = c("", "NA"),
        col.names = blNames)
      bl <- rbind(blr, use.names = T)
      bl <- bl[,blNames, with = F]
    }

    # -- 2.2 drop duplicated hits, keeping the highest score
    bitScore <- ofID1 <- ofID2 <- og1 <- og2 <- og <- noAnchor1 <- ns1 <- ns2 <-
      noAnchor2 <- noAnchor <- rnd1 <- rnd2 <- ngene1 <- ngene2 <- sameOg <-
      ord1 <- ord2 <- ancOrd1 <- ancOrd2 <- n <- NULL
    setorder(bl, -bitScore)
    bl <- subset(bl, !duplicated(paste(ofID1, ofID2)))
    cat(sprintf("\t%s vs. %s: total hits = %s",
                synMd$query[i], synMd$target[i], nrow(bl)))

    # -- 2.3 merge with bed information
    bl <- merge(bed1, merge(bed2, bl, by = "ofID2"), by = "ofID1")

    # -- 2.4 get anchor and og information
    bl[,sameOg := og1 == og2 & !is.na(og1) & !is.na(og2)]
    bl[,noAnchor := noAnchor1 | noAnchor2]
    bl[,`:=`(og1 = NULL, og2 = NULL, noAnchor1 = NULL, noAnchor2 = NULL)]
    cat(sprintf(", same og = %s\n", sum(bl$sameOg)))

    bl[,`:=`(isAnchor = NA, inBuffer = NA, regID = NA, blkID = NA, lgBlkID = NA)]
    blf <- file.path(
      gsParam$paths$syntenicHits,
      sprintf("%s_vs_%s.synBlast.txt.gz", synMd$query[i], synMd$target[i]))

    write_synBlast(bl, filepath = blf)
    synMd$annotBlastFile[i] <- blf
    # -- 2.5 summarize the hits, make plots
    hits <- data.table(bl)
    hits[,`:=`(rnd1 = round_toInteger(ord1, 50), rnd2 = round_toInteger(ord2, 50))]
    un1 <- uniqueN(hits$ofID1)
    un2 <- uniqueN(hits$ofID2)
    if(un1 > un2){
      ht <- 8
      wd <- ht * (un1/un2)
    }else{
      wd <- 8
      ht <- wd * (un2/un1)
    }

    dpFile <- file.path(gsParam$paths$dotplots,
                        sprintf("%s_vs_%s.rawHits.pdf",
                                synMd$query[i], synMd$target[i]))
    pdf(dpFile, height = ht, width = wd*1.07)
    hits[,ngene1 := uniqueN(ofID1), by = "chr1"]
    hits[,ngene2 := uniqueN(ofID2), by = "chr2"]
    hc <- subset(hits, sameOg & ngene1 > 50 & ngene2 > 50)
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    hc[,ns2 := sum(n), by = c("chr2")]
    hc$n[hc$n > 20] <- 20
    p1 <- ggplot(subset(hc), aes(rnd1, rnd2, col = n)) +
      geom_point(pch = ".") +
      scale_color_viridis_c(begin = .1, trans = "log10") +
      scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
      scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = rgb(1,1,1,.2), size = .2, linetype = 2),
            panel.spacing = unit(.1, "mm"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = 90, family = "Helvetica", size = 5),
            strip.text.y.left = element_text(angle = 0, family = "Helvetica", size = 5),
            axis.title = element_text( family = "Helvetica", size = 6),
            plot.title = element_text( family = "Helvetica", size = 7)) +
      facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
      labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome1[1], uniqueN(hits$ofID1)),
           y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome2[1], uniqueN(hits$ofID2)),
           title = "Blast hits where query and target are in the same OG, 50-gene window thresholded at 20 hits")

    hc <- subset(hits, bitScore > 100 & ngene1 > 50 & ngene2 > 50)
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    hc[,ns2 := sum(n), by = c("chr2")]
    hc[,ns2 := sum(n), by = c("chr2")]
    hc$n[hc$n > 20] <- 20
    p0 <- ggplot(subset(hc), aes(rnd1, rnd2, col = n)) +
      geom_point(pch = ".") +
      scale_color_viridis_c(begin = .1, trans = "log10") +
      scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
      scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = rgb(1,1,1,.2), size = .2, linetype = 2),
            panel.spacing = unit(.1, "mm"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = 90, family = "Helvetica", size = 5),
            strip.text.y.left = element_text(angle = 0, family = "Helvetica", size = 5),
            axis.title = element_text( family = "Helvetica", size = 6),
            plot.title = element_text( family = "Helvetica", size = 7)) +
      facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
      labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome1[1], uniqueN(hits$ofID1)),
           y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome2[1], uniqueN(hits$ofID2)),
           title = "All blast hits, with bitscore > 100, 50-gene window thresholded at 20 hits")

    hc <- subset(hits, !noAnchor & sameOg & ngene1 > 50 & ngene2 > 50)
    hc <- hc[,list(n = .N), by = c("chr1", "chr2", "rnd1", "rnd2")]
    hc[,ns2 := sum(n), by = c("chr2")]
    hc[,ns2 := sum(n), by = c("chr2")]
    hc$n[hc$n > 20] <- 20
    p2 <- ggplot(subset(hc), aes(rnd1, rnd2, col = n)) +
      geom_point(pch = ".") +
      scale_color_viridis_c(begin = .1, trans = "log10") +
      scale_x_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd1), by = 1e3))+
      scale_y_continuous(expand = c(0,0), breaks = seq(from = 1e3, to = max(hc$rnd2), by = 1e3))+
      theme(panel.background = element_rect(fill = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(color = rgb(1,1,1,.2), size = .2, linetype = 2),
            panel.spacing = unit(.1, "mm"),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(angle = 90, family = "Helvetica", size = 5),
            strip.text.y.left = element_text(angle = 0, family = "Helvetica", size = 5),
            axis.title = element_text( family = "Helvetica", size = 6),
            plot.title = element_text( family = "Helvetica", size = 7)) +
      facet_grid(chr2 ~ chr1, scales = "free", space = "free", as.table = F, switch = "both")+
      labs(x = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome1[1], uniqueN(hits$ofID1)),
           y = sprintf("%s: gene rank order position (%s genes with blast hits), gridlines every 1000 genes",
                       hits$genome2[1], uniqueN(hits$ofID2)),
           title = "Only blast hits to be considered for initial syntenic anchors, 50-gene window thresholded at 20 hits")

    print(p0)
    print(p1)
    print(p2)
    dev.off()
  }
  gsParam$annotBlastMd <- synMd
  return(gsParam)
}
