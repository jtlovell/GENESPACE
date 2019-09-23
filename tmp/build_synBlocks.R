#' @title build_synBlocks
#'
#' @description
#' \code{build_synBlocks} build_synBlocks
#'
#' @param dir.list dir.list
#' @param genomeIDs genomeIDs
#' @param mcscan.dir path to the directory where MCScanX will be run
#' @param gap.multiplier parameters to supply MCScanX
#' @param MCScanX.path the path to the the MCScanX program. If in the path,
#' just use "MCScanX".
#' @param min.blockSize min.blockSize
#' @param radius radius
#' @param clean.before.mcscanx clean.before.mcscanx
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
build_synBlocks <- function(dir.list,
                            genomeIDs,
                            n.cores = 6,
                            gff,
                            min.blockSize = 5,
                            m.param = 25,
                            MCScanX.path,
                            verbose = T,
                            only.orthogroups = T,
                            use.topn = T,
                            method = "pairwise"){

  MCScanX.tool <- file.path(MCScanX.path,"MCScanX")

  #######################################################

  #######################################################
  if(method == "pairwise"){
    if (verbose)
      cat("Running pairwise orthofinder calls ... \n")
    comb <- combn(genomeIDs, 2, simplify = F)
    wh.2keep <- sapply(genomeIDs, USE.NAMES = T, function(x)
      min(which(sapply(comb, function(y) x %in% y))))

    pw.map <- rbindlist(lapply(1:length(comb), function(i){
      x <- comb[[i]]
      w2ki <- wh.2keep[wh.2keep == i]

      if (verbose)
        cat(paste0("\t",x[1]), "<-->", x[2],"... ")
      blast <- rerun_pairwiseOF(dir.list = dir.list,
                                genomeIDs = x,
                                n.cores = n.cores,
                                gff = gff,
                                verbose = F)
      out <- mirror_blast(blast = blast,
                          w2ki = w2ki,
                          genomes = x)
      if (verbose)
        cat("found", nrow(out), "hits", with(out, sum(og1 == og2)), "in orthogroups\n")

      out <- subset(out, og1 == og2)
      out[, og.id := og1]
      out$og1 <- NULL
      out$og2 <- NULL

      return(out)
    }))
    if (verbose)
      cat("\tDone!\n")
    #######################################################
  }else{
    stop("Only pairwise blocks currently implemented")
  }
  #######################################################
  if (verbose)
    cat("Forming syntenic blocks via MCScanX ... \n")
  mcsp <- paste("-a -s", min.blockSize,
                "-m", m.param)

  ug12 <- data.table(pw.map[,c("genome1","genome2")])
  ug12 <- ug12[!duplicated(ug12),]
  wh1 <- grep("1$", colnames(pw.map))
  wh2 <- grep("2$", colnames(pw.map))
  who <- which(!grepl("1$|2$", colnames(pw.map)))

  n1 <- colnames(pw.map)[wh1]
  n2 <- colnames(pw.map)[wh2]
  no <- colnames(pw.map)[who]
  m1 <- pw.map[,c(n1,n2,no), with = F]
  m2 <- pw.map[,c(n2,n1,no), with = F]
  setnames(m2, c(n1,n2,no))
  pw.map <- rbind(m1, m2)
  pw.map <- pw.map[!duplicated(pw.map[,c("id1","id2")]),]

  pw.map <- merge(ug12, pw.map, by = c("genome1","genome2"))
  syn.blks <- pipe_mcscanx(blast = pw.map,
                           gff = gff,
                           dir.list = dir.list,
                           genomeIDs = genomeIDs,
                           MCScanX.path = MCScanX.tool,
                           mcscan.dir = dir.list$mcscanx,
                           mcscan.param = mcsp,
                           verbose = T)
  if (verbose)
    cat("\tDone!")
  out <- list(blast = pw.map,
              map = syn.blks$map,
              blk = syn.blks$block)

  if(verbose)
    cat("\nFormatting and returning block and map objects ... ")
  m <- data.table(out$map)
  m[,block.id := paste0("blk_",as.numeric(as.factor(paste(genome1, genome2, chr1, chr2, block.id))))]
  tab <- table(m$block.id)
  good.blks <- names(tab)[tab >= min.blockSize]

  wh1 <- grep("1$", colnames(m))
  wh2 <- grep("2$", colnames(m))
  who <- which(!grepl("1$|2$", colnames(m)))

  m <- subset(m, block.id %in% good.blks)
  n1 <- colnames(m)[wh1]
  n2 <- colnames(m)[wh2]
  no <- colnames(m)[who]
  m1 <- m[,c(n1,n2,no), with = F]
  m2 <- m[,c(n2,n1,no), with = F]
  setnames(m2, c(n1,n2,no))
  mt <- rbind(m1, m2)

  mo <- make_blocks(mt[!duplicated(mt[,c("id1","id2")]),], clean.columns = F)
  out$map <- mo$map
  out$blk <- mo$block

  if(verbose)
    cat("Done!\n")
  return(out)
}
