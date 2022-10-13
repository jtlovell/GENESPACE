#' @title run_mcscanx
#'
#' @description
#' \code{run_mcscanx} run_mcscanx
#'
#' @param hits hits
#' @param nGaps integer of length 1, specifying the -m param to mcscanx
#' for the primary MCScanX run. This acts on the results from the initial
#' MCScanX run.
#' @param blkSize integer of length 1, specifying the -s param to mcscanx
#' @param tmpDir tmpDir
#' @param MCScanX_hCall MCScanX_hCall
#'
#' @details xxx
#'
#' @return xxx
#' \enumerate{
#' \item xxx
#' }
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @import data.table
#' @export
run_mcscanx <- function(hits,
                        blkSize,
                        nGaps,
                        tmpDir,
                        MCScanX_hCall){
  setDTthreads(1)
  mcsID1 <- mcsID2 <- NULL
  ##############################################################################
  # parameter argument checking
  hitCols <- c("ofID1","ofID2","chr1","chr2","ord1","ord2","score")
  if(!all(hitCols %in% colnames(hits)))
    stop("Looks like the hits data.table is misformed. It must be a data.table",
         "comparing two genomes with the columns:", hitCols,"\n")
  if(length(blkSize) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(length(nGaps) != 1)
    stop("blkSize must be a integer of length 1\n")
  if(!is.integer(blkSize)){
    blkSize <- as.integer(blkSize)
    if(is.na(blkSize))
      stop("Cannot coerce blkSize to an integer\n")
  }
  if(!is.integer(nGaps)){
    nGaps <- as.integer(nGaps)
    if(is.na(nGaps))
      stop("Cannot coerce nGaps to an integer\n")
  }
  if(!file.exists(MCScanX_hCall))
    stop("Cannot find MCScanX_h executable in", MCScanX_hCall,"\n")

  ##############################################################################
  # set tmp directory
  tmpd <- file.path(
    tmpDir,
    paste0("tmp_", paste(
      sample(c(letters,LETTERS), 20, replace = T),
      collapse = "")))
  if(dir.exists(tmpd))
    unlink(tmpd, recursive = T)
  dir.create(tmpd)
  on.exit(expr = unlink(tmpd, recursive = T))

  ##############################################################################
  # convert gene locations to 'gff' like mcscanx file
  ord1 <- ord2 <- chr1 <- chr2 <- ofID1 <- ofID2 <- NULL
  u1 <- subset(hits, !duplicated(ofID1))
  u2 <- subset(hits, !duplicated(ofID2))
  setkey(u1, chr1, ord1)
  setkey(u2, chr2, ord2)
  u1[,mcsID1 := paste0("aa", as.numeric(factor(chr1, levels = unique(chr1))))]
  u2[,mcsID2 := paste0("bb", as.numeric(factor(chr2, levels = unique(chr2))))]

  mcs1 <- u1[,c("mcsID1", "ofID1","ord1","ord1")]
  setnames(mcs1, c("chr", "id","start","end"))

  mcs2 <- u2[,c("mcsID2", "ofID2","ord2","ord2")]
  setnames(mcs2, c("chr", "id","start","end"))
  mcsGffIn <- rbind(mcs1, mcs2)

  ##############################################################################
  # convert hits to mcscanx_h blast input
  mcsBlsIn <- hits[,c("ofID1", "ofID2", "score")]
  # mcsBlsIn[,score := 1]

  mcsBlsIn[,ofID2 := paste0(ofID2, "xxxx")]
  mcsGffIn$id[grepl("bb", mcsGffIn$chr)] <- paste0(mcsGffIn$id[grepl("bb", mcsGffIn$chr)],"xxxx")

  blFile <- file.path(tmpd, "mcs.homology")
  gfFile <- file.path(tmpd, "mcs.gff")
  colFile <- file.path(tmpd, "mcs.collinearity")

  fwrite(
    mcsGffIn, file = gfFile, sep = "\t", quote = FALSE, col.names = FALSE,
    showProgress = FALSE, verbose = FALSE)
  fwrite(
    mcsBlsIn, file = blFile, sep = "\t", quote = FALSE, col.names = FALSE,
    showProgress = FALSE, verbose = FALSE)

  ##############################################################################
  # run mcscanx_h
  mcsCom <- sprintf(
    "-a -b 2 -c 2 -m %s -s %s %s",
    nGaps, blkSize, file.path(tmpd, "mcs"))
  comout <- system2(MCScanX_hCall, mcsCom, stdout = TRUE, stderr = TRUE)

  ##############################################################################
  # parse collinearity file
  idg <- strsplit(as.character(hits$ofID1[1]), "_")[[1]][1]
  suppressWarnings(collin <- fread(
    cmd = sprintf("cat %s | grep %s_ | grep :", colFile, idg),
    col.names = c("blkID","gn1","gn2"),
    select = 1:3,
    showProgress = FALSE,
    header = FALSE))

  if(nrow(collin) > 1){
    blkID <- gn1 <- gn2 <- NULL
    collin[,blkID := as.numeric(sapply(blkID, function(x)
      strsplit(x, "-")[[1]][1])) + 1]
    collin[,gn2 := gsub("xxxx$", "", gn2)]
    mcsb <- collin$blkID
    names(mcsb) <- with(collin, paste(gn1, gn2))
    return(mcsb)
  }
}
