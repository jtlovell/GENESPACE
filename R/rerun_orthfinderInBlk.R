#' @title Run the orthofinder program
#'
#' @description
#' \code{rerun_orthofinderInBlk} A simple wrapper to run orthofinder from R.
#'
#' @param init.results Results from initial orthofinder run (via parse_orthofinder)
#' @param blk block file from make blocks or whatever
#' @param genomeIDs the genome ids.
#' @param dir.list List of directories, created by check_environment
#' @param n.cores numeric, number of parallel processes to spawn.
#' @param ... Additional arguments passed on to run_orthofinder
#' @details Nothing yet

#' @return Nothing, writes results to the cull.blast.dir directory
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
rerun_orthofinderInBlk = function(blk,
                                  init.results,
                                  dir.list,
                                  genomeIDs,
                                  verbose = T,
                                  n.cores = 1,
                                  ...){

  #######################################################
  #######################################################
  if (verbose)
    cat("Importing annotations ... ")
  gff <- init.results$gff
  gff.wNum <- mclapply(1:nrow(blk), mc.cores = n.cores, function(i)
    pull_gff(gff = gff,
             blk.line = blk[i,],
             gene.index = init.results$ortho.info$gene.index))
  if (verbose)
    cat("Done!\n")


  gffn <- gff.wNum
  genenum.list <- lapply(gff.wNum, function(x)
    unique(c(x[[1]]$gene.num, x[[2]]$gene.num)))

  #######################################################
  #######################################################
  if (verbose)
    cat("Parsing blast files ... ")

  remake_blast(blast.dir = dir.list$blast,
               cull.blast.dir = dir.list$tmp,
               genenum.list = genenum.list,
               n.cores = n.cores,
               verbose = verbose)
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Running orthofinder from culled blast files ...")

  run_orthofinder(
    peptide.dir = NULL,
    tmp.dir = dir.list$tmp,
    blast.dir = dir.list$cull.blast,
    og.threads = n.cores,
    og.silent = F,
    verbose = T)
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Reading blast results ...")
  all.blasts <- read_allBlast(blast.dir = dir.list$cull.blast)
  all.blkgff <- rbindlist(unlist(gff.wNum, recursive = F))
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Building orthogroup networks into gff annotations ...")
  gffog <- make_mapFromOGs(gff.wNum = gff.wNum,
                           cull.blast.dir = dir.list$cull.blast)
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Reading in cds fastas and building sequence length metadata ...")
  cds.fastas <- do.call(c,lapply(genomeIDs, function(i)
    readDNAStringSet(file.path(dir.list$cds,
                               paste0(i,".fa")))))
  cds.md <- data.table(id = names(cds.fastas),
                       length = width(cds.fastas),
                       stringsAsFactors = F)
  setkey(cds.md, id)
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Combining block information with cds and gff metadata ...")
  md.list <- make_blockMetadata(cds.md = cds.md,
                                gffog = gffog,
                                blk = blk)
  spl.gff <- md.list$spl.gff
  blk.md <- md.list$blk.md
  blk.md2 <- md.list$blk.md2
  if (verbose)
    cat("Done!\n")

  #######################################################
  #######################################################
  if (verbose)
    cat("Finding orphan genes in blocks ...")

  orphan.md <- with(md.list,
                    find_orphans(spl.gff = spl.incompleteGff,
                                 blk.md2 = simple.blk.metadata,
                                 n.cores = n.cores))
  if (verbose)
    cat("Done!\n")


  #######################################################
  #######################################################
  if (verbose)
    cat("Loading annotations into memory ...")
  annot <- load.annotations(genomeIDs = genomeIDs,
                            cds.dir = dir.list$cds,
                            peptide.dir = dir.list$peptide,
                            assembly.dir = dir.list$assembly)
  if (verbose)
    cat("Done!\n")
  return(list(annot = annot,
              orhan.md = orphan.md,
              md.list = md.list,
              cds.md = cds.md,
              gffog = gffog,
              cds.fastas = cds.fastas))

}
