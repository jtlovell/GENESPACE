#' @title Build orthofinder database
#'
#' @description
#' \code{build_ofDb} Simplified blast database construction for and orthogroup
#' construction in orthofinder.
#'
#' @param gsAnnot list of length 2, containing the genespace annotation paths
#' 'gff' and 'peptide' -- file path character vector with the locations of
#' the peptide and gff-like annotation files. Each element is named by the
#' associated genomeID in gsParam
#' @param gsParam A list of genespace parameters. This should be created
#' by setup_genespace, but can be built manually. Must have the following
#' elements: blast (file.path to the original orthofinder run), synteny (
#' file.path to the directory where syntenic results are stored), genomeIDs (
#' character vector of genomeIDs).
#'
#' @details Build pairwise blast database without running both sides of
#' pairwise blasts. The genome with more genes is set as the query, then to
#' complete the blast database, a mirrored blast file is written with
#' the query genome as the target. Self blast runs are done as usual.
#'
#' @return A data.table that is needed for many downstream analyses. This has
#' metadata for each pair of genomes with the following columns:
#' \enumerate{
#' \item genome1/2: the genome IDs for each pairwise run
#' \item uGenes1/2: the number of unique genes for each genome
#' \item query/target: the identity of query/target genomes in each blast run
#' \item run/mirrorBlast: logical whether the blast was run or mirrored
#' \item gn1/2: orthofinder genomeID numbers for each genome
#' \item db/fa/1/2/blFile: diamond database, fasta files and blast files.
#' }
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#' @import R.utils
#' @import data.table
#' @export
build_OFDb <- function(gsAnnot,
                       gsParam){

  gn1 <- gn2 <- nGenes2 <- nGenes1 <- u <- genome2 <- genome1 <- NULL
  fa1 <- fa2 <- db2 <- db1 <- mirrorBlast <- runBlast <- target <- query <- NULL
  tiebreak <- invertFile <- blFile <- NULL

  pepFiles <- gsAnnot$peptide
  gffFiles <- gsAnnot$gff
  nCores <- gsParam$nCores
  path2diamond <- gsParam$path2diamond
  path2orthofinder <- gsParam$path2orthofinder
  verbose <- gsParam$verbose
  blastDir <- gsParam$blast

  ##############################################################################
  # -- make directories
  if (verbose)
    cat("\tSetting up orthofinder environment ... ")
  if (dir.exists(blastDir))
    unlink(blastDir, recursive = T)
  dir.create(blastDir)

  genomeIDs <- names(pepFiles)
  ##############################################################################
  # -- decompress peps to a tempdir, format for orthofinder, copy to blastdir

  if(gsParam$orthofinderMethod == "fast"){
    tmpDir <- file.path(
      getwd(),
      paste0("tmp_", gsub("[^A-Za-z0-9]", "", Sys.time())))
    if (dir.exists(tmpDir))
      unlink(tmpDir, recursive = T)
    dir.create(tmpDir)
    on.exit(expr = unlink(tmpDir, recursive = T))

    for(i in pepFiles)
      gunzip(
        filename = i,
        destname = file.path(tmpDir, gsub(".gz$","",basename(i))),
        remove = F)

    com <- paste(path2orthofinder,
                 "-f", tmpDir,
                 "-t", nCores,
                 "-S diamond -op 1>/dev/null 2>&1")
    system(com)

    ofInputLoc <- dirname(list.files(
      path = tmpDir,
      pattern = "SequenceIDs.txt",
      recursive = T,
      full.names = T))
    ofInputFiles = list.files(path = ofInputLoc, full.names = T)
    nu <- file.copy(ofInputFiles, blastDir)

    ##############################################################################
    # -- choose the files to blast
    # -- the query genome should be the larger of the two.
    cmb <- data.table(genome1 = genomeIDs, genome2 = genomeIDs)
    cmb <- cmb[,CJ(genome1, genome2)]
    nGenes <- sapply(pepFiles, function(x) length(readAAStringSet(x)))
    cmb[,u := paste(
      genomeIDs[genomeIDs %in% c(genome1, genome2)], collapse =" vs. "),
      by = c("genome1","genome2")]
    cmb[,nGenes1 := nGenes[genome1]]
    cmb[,nGenes2 := nGenes[genome2]]
    cmb[,tiebreak :=
          as.numeric(factor(genome1, levels = genomeIDs)) -
          as.numeric(factor(genome2, levels = genomeIDs))]
    cmb[,query := ifelse(nGenes1 > nGenes2,
                         genome1,
                         ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                genome1,
                                genome2)),
        by = "u"]
    cmb[,target := ifelse(nGenes1 > nGenes2,
                          genome2,
                          ifelse(nGenes1 == nGenes2 & tiebreak > 0,
                                 genome2,
                                 genome1)),
        by = "u"]
    cmb[,runBlast := genome1 == query]
    cmb[,mirrorBlast := genome1 == target & genome1 != query]

    si <- read_orthofinderSpeciesIDs(blastDir)
    cmb[,gn1 := si[genome1]]
    cmb[,gn2 := si[genome2]]
    cmb[,db1 := file.path(blastDir, paste0("diamondDBSpecies",gn1,".dmnd"))]
    cmb[,db2 := file.path(blastDir, paste0("diamondDBSpecies",gn2,".dmnd"))]
    cmb[,fa1 := file.path(blastDir, paste0("Species",gn1,".fa"))]
    cmb[,fa2 := file.path(blastDir, paste0("Species",gn2,".fa"))]
    cmb[,blFile := file.path(blastDir, paste0("Blast",gn1,"_",gn2,".txt.gz"))]
    cmb[,invertFile := file.path(blastDir, paste0("Blast",gn2,"_",gn1,".txt.gz"))]

    ##############################################################################
    # -- run blasts for each of the primary runs
    # -- mirror blasts over where necessary
    if (verbose)
      cat("Done!\n\tDiamond blastp searches ...\n")
    cmb$invertFile[with(cmb, genome1 == genome2)] <- NA
    cmbbl <- subset(cmb, runBlast)

    for (i in 1:nrow(cmbbl)) {
      if (verbose)
        cat("\t\tRunning", i,"/", nrow(cmbbl), "... ")
      x <- cmbbl[i,]
      com <-  with(x, sprintf(
        "%s blastp --quiet -e %s -p %s --compress 1 -d %s -q %s -o %s",
        path2diamond, .1, nCores, db2, fa1, blFile))
      system(com)

      if (!is.na(x$invertFile)) {
        tmp <- fread(x$blFile, verbose = F, showProgress = F)
        tmp <- tmp[,c(2,1,3:6,8,7,10,9,11,12)]
        fwrite(tmp, sep = "\t",
               quote = FALSE,
               col.names = FALSE,
               row.names = FALSE,
               file = x$invertFile,
               showProgress = FALSE,
               verbose = FALSE)
      }
      if (verbose)
        cat("Done!\n")
    }

    ##############################################################################
    # -- run orthofinder
    if (verbose)
      cat("\tRunning orthofinder ... ")
    com <- paste(path2orthofinder,
                 "-b", blastDir,
                 "-a", nCores,
                 "-t", nCores,
                 "-og")
    if (!verbose)
      com <- paste(com, "1>/dev/null 2>&1")
    system(com)

    if(verbose)
      cat("\nCleaning up results ... ")
    tsvFile <- order_filesByMtime(
      path = blastDir,
      pattern = "Orthogroups.tsv",
      recursive = T)
    if(length(tsvFile) > 0){
      tsvFile <- tsvFile[1]
    }
    file.copy(tsvFile, blastDir)
    unlink(file.path(blastDir, "Orthofinder"), recursive = T)
  }else{
    for(i in pepFiles)
      gunzip(
        filename = i,
        destname = file.path(blastDir, gsub(".gz$","",basename(i))),
        remove = F)

    com <- paste(path2orthofinder,
                 "-f", blastDir,
                 "-a", nCores,
                 "-t", nCores,
                 "-og")
    if (!verbose)
      com <- paste(com, "1>/dev/null 2>&1")
    system(com)

    if(verbose)
      cat("\nCleaning up results ... ")
    tsvFile <- order_filesByMtime(
      path = blastDir,
      pattern = "Orthogroups.tsv",
      recursive = T)[1]
    blFiles <- order_filesByMtime(
      path = blastDir,
      pattern = "Blast",
      recursive = T)
    blFiles <- blFiles[!duplicated(basename(tsvFile))]
    dmndFiles <- order_filesByMtime(
      path = blastDir,
      pattern = "diamondDBSpecies",
      recursive = T)
    dmndFiles <- dmndFiles[!duplicated(basename(dmndFiles))]
    seqIDFile <- file.path(dirname(blFiles[1]),"SequenceIDs.txt")
    speIDFile <- file.path(dirname(blFiles[1]),"SpeciesIDs.txt")
    faFiles <- file.path(dirname(blFiles[1]),
                         sprintf("Species%s.fa", 0:(length(gsParam$genomeIDs)-1)))
    nu <- sapply(c(tsvFile, blFiles, seqIDFile, speIDFile, faFiles, dmndFiles), function(x)
      file.copy(x, blastDir))
    for(i in gsParam$genomeIDs)
      unlink(file.path(blastDir, sprintf("%s.fa", i)))
    unlink(file.path(blastDir, "Orthofinder"), recursive = T)
  }
  if(verbose)
    cat("Done!\n")
  return(com)
}
