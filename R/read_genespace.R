#' @title Genespace plotting routines
#' @description
#' \code{read_genespace} Genespace human-readable data import routines
#' @name read_genespace
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
#'
#' @details ...
#'
#' @note \code{read_genespace} is a generic name for the functions documented.
#' \cr
#' If called, \code{read_genespace} returns its own arguments.
#'
#' @title Read in hits between two genomes
#' @description
#' \code{read_syntenicHits} Read in hits between two genomes
#' @rdname read_genespace
#' @import data.table
#' @export


#' @title Read in the pangenome database
#' @description
#' \code{read_pangenome} Read in the pangenome database
#' @rdname read_genespace
#' @import data.table
#' @export
read_pangenome <- function(gsParam,
                           genomeIDs,
                           includePrivateGenes = TRUE,
                           arraySep = ",",
                           RBHFlag = "*",
                           nonSynOGFlag = "+"){

}

#' @title import a gff-like file
#' @description
#' \code{read_gff} reads a genespace-formatted gff-like annotation file into
#' memory
#' @rdname utils
#' @export
read_gff <- function(gffFiles){
  genome <- NULL
  gff <- rbindlist(lapply(names(gffFiles), function(i){
    x <- fread(
      gffFiles[[i]],
      key = c("chr", "start","end","strand"))
    x[,genome := i]
    return(x)
  }))
  gffCols <- c("ord","chr","genome","start", "end", "strand","id")
  if(!all(gffCols %in% colnames(gff)))
    stop(paste(gffCols, collapse = ", "),
         " must all be column names in gff\n")
  return(gff)
}

#' @title Convert orthofinder metadata to vector
#' @description
#' \code{pull_orthogroups} Parse orthogroup file, pull orthofinder gene
#' IDs and return a named vector of orthogroups
#' @rdname utils
#' @import data.table
#' @export
pull_orthogroups <- function(path, genomeIDs = NULL){
  id <- ofID <- genomeNum <- genome <- NULL

  # find orthogroup.tsv file
  tsvFile <- file.path(path, "Orthogroups.tsv")[1]
  if(!file.exists(tsvFile))
    stop("Can't find the orthogroups file:", tsvFile,"\n")

  # read and parse species and seq ids
  specIDs <- read_orthofinderSpeciesIDs(path)
  if(is.null(genomeIDs))
    genomeIDs <- names(specIDs)
  specIDs <- specIDs[genomeIDs]
  seqIDs <- subset(read_orthofinderSequenceIDs(path), genomeNum %in% unique(specIDs))

  gv <- names(specIDs); names(gv) <- as.character(specIDs)
  seqIDs[,genome :=  gv[as.character(genomeNum)]]
  idv <- seqIDs$ofID; names(idv) <- with(seqIDs, paste(genome, id))

  # convert orthogroup.tsv into data.table
  og <- fread(
    tsvFile,
    verbose = F,
    showProgress = F)[,c("Orthogroup", genomeIDs),with = F]
  og <- melt(
    og,
    id.vars = "Orthogroup",
    measure.vars = names(og)[-1],
    value.name = "id",
    variable.name = "genome")
  og <- subset(og, id != "")
  og[,id := list(strsplit(id, ","))[[1]]]
  og <- og[,list(id = unlist(id)),
           by = c("Orthogroup","genome")]
  og[,id := trimws(id)]

  # convert to vector and return
  og[,ofID:= idv[paste(genome, id)]]
  out <- og$Orthogroup; names(out) <- og$ofID
  seqIDs[,og := out[ofID]]
  seqIDs[,`:=`(genomeNum = NULL, geneNum = NULL)]
  seqIDs$og[is.na(seqIDs$og)] <- paste0("NoOG_", 1:sum(is.na(seqIDs$og)))

  return(list(ogv = out, ogdt = seqIDs))
}

#' @title Read orthofinder species IDs
#' @description
#' \code{read_orthofinderSpeciesIDs} Parses the SpeciesIDs.txt file into a
#' data.table and returns to R.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSpeciesIDs <- function(path){
  genome <- NULL
  si <- fread(
    file.path(path, "SpeciesIDs.txt"),
    sep = ":",
    header = F,
    col.names = c("genomeNum", "genome"),
    colClasses = c("numeric", "character"))
  si[,genome := gsub(".fa", "", genome, fixed = T)]
  sio <- si$genomeNum; names(sio) <- si$genome
  return(sio)
}

#' @title Read orthofinder sequence IDs
#' @description
#' \code{read_orthofinderSequenceIDs} Reads the sequence
#' IDs:gene name dictionary into memory.
#' @rdname utils
#' @import data.table
#' @export
read_orthofinderSequenceIDs <- function(path){
  ofID <- NULL
  gi <- fread(
    file.path(path, "SequenceIDs.txt"),
    header = F,
    sep = ":",
    col.names = c("ofID","id"),
    colClasses = c("character","character"))
  gi[,c("genomeNum","geneNum") := tstrsplit(ofID, "_", type.convert = T)]
  return(gi)
}

#' @title Read orthofinder blast file
#' @description
#' \code{read_blast} Reads in a single pairwise orthofinder-formaatted blast
#' file
#' @rdname utils
#' @import data.table
#' @export
read_blast <- function(blFile = NULL,
                       ofID1 = NULL,
                       ofID2 = NULL,
                       path = NULL,
                       onlyIDScore = TRUE){
  V12 <- score <- NULL
  if(is.null(blFile)){
    blFile <- file.path(path, paste0("Blast", ofID1, "_", ofID2,".txt.gz"))
  }
  if(!file.exists(blFile))
    stop("cannot find ", blFile, "\n")

  if(!onlyIDScore){
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE)
    g1 <- strsplit(bl$V1[1], "_")[[1]][1]
    g2 <- strsplit(bl$V2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3:6,8,7,10,9,11,12)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -V12)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }else{
    bl <-  fread(
      blFile,
      showProgress = FALSE,
      verbose = FALSE,
      select = c(1,2,12),
      col.names = c("ofID1","ofID2","score"))
    g1 <- strsplit(bl$ofID1[1], "_")[[1]][1]
    g2 <- strsplit(bl$ofID2[1], "_")[[1]][1]

    if(g1 == g2){
      tmp <- data.table(bl[,c(2,1,3)])
      setnames(tmp, colnames(bl))
      tmp <- tmp[,colnames(bl),with = F]
      bl <- rbind(bl, tmp)
      setorder(bl, -score)
      bl <- subset(bl, !duplicated(bl[,c(1:2)]))
    }
  }

  return(bl)
}


write_phytozome <- function(gsParam, gsAnnot){
  outdir <- file.path(gsParam$results, "phytozome")
  dir.create(outdir)
  blksFile <- file.path(gsParam$results, "blockCoordinates_geneOrder.txt.gz")
  blk <- fread(blksFile)

  genomeIDs <- unique(c(blk$genome1, blk$genome2))
  gff <- add_ofID2gff(read_gff(gsAnnot$gff[genomeIDs]), blastDir = gsParam$blast)

  sogFile <- list.files(gsParam$synteny, pattern = "_synog.txt.gz", full.names = T)
  so <- rbindlist(lapply(sogFile, fread))
  so <- subset(so, ofID1 %in% gff$ofID & ofID2 %in% gff$ofID)
  ov <- gff$ord; names(ov) <- gff$ofID
  cv <- gff$chr; names(cv) <- gff$ofID
  sv <- gff$start; names(sv) <- gff$ofID
  ev <- gff$end; names(ev) <- gff$ofID
  gv <- gff$genome; names(gv) <- gff$ofID
  iv <- gff$id; names(iv) <- gff$ofID
  so[,`:=`(genome1 = gv[ofID1], genome2 = gv[ofID2],
           chr1 = cv[ofID1], chr2 = cv[ofID2],
           start1 = sv[ofID1], start2 = sv[ofID2],
           end1 = ev[ofID1], end2 = ev[ofID2])]
  so <- rbind(so, with(so, data.table(
    ofID1 = ofID2, ofID2 = ofID1, blkID = blkID,
    genome1 = genome2, genome2 = genome1, chr1 = chr2, chr2 = chr1,
    start1 = start2, start2 = start1, end1 = end2, end2 = end1)))
  so <- subset(so, !duplicated(so))
  bo <- rbind(
    with(subset(blk, orient == "+"), data.table(
      proteomeA = c(genome1, genome2),
      chromA = c(chr1, chr2),
      startA = c(bpStart1, bpStart2),
      endA = c(bpEnd1, bpEnd2),
      proteomeB = c(genome2, genome1),
      chromB = c(chr2, chr1),
      startB = c(bpStart2, bpStart1),
      endB = c(bpEnd2, bpEnd1),
      strand = c(orient, orient),
      blkID = c(blkID, blkID))),
    with(subset(blk, orient == "-"), data.table(
      proteomeA = c(genome1, genome2),
      chromA = c(chr1, chr2),
      startA = c(bpStart1, bpStart2),
      endA = c(bpEnd1, bpEnd2),
      proteomeB = c(genome2, genome1),
      chromB = c(chr2, chr1),
      startB = c(bpEnd2, bpEnd1),
      endB = c(bpStart2, bpStart1),
      strand = c(orient, orient),
      blkID = c(blkID, blkID))))
  bo <- subset(bo, !duplicated(bo))
  splb <- split(bo, by = "proteomeA")
  splo <- split(so, by = "genome1")
  nu <- lapply(names(splb), function(i){
    x <- splb[[i]]
    y <- splo[[i]]
    splx <- split(x, by = "proteomeB")
    sply <- split(y, by = "genome2")
    dirx <- file.path(outdir, i)
    dir.create(dirx)
    dir.create(file.path(dirx, "alignments"))
    for(j in names(splx)){
      bcoords <- splx[[j]]
      oc <- sply[[j]]
      splk <- split(oc, by = "blkID")
      fwrite(bcoords[,1:9, with = F],
             file = file.path(dirx, sprintf("%s_%s.csv", i, j)))
      for(k in 1:nrow(bcoords)){
        w <- bcoords[k,]
        print(w)
        z <- with(splk[[w$blkID]], data.table(
          geneA = iv[ofID1], startA = start1, endA = end1,
          geneB = iv[ofID2], startB = start2, endB = end2))
        fn <- with(w, sprintf("%s__%s:%s-%s_%s__%s:%s-%s.csv",
                              proteomeA, chromA, startA, endA,
                              proteomeB, chromB, startB, endB))
        fwrite(z, file = file.path(dirx, "alignments", fn))
      }
    }
  })
}
