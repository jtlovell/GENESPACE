#' @title Synteny-constrained orthology pipeline
#'
#' @description
#' \code{find_synOrthos} Subset blast hits to syntenic regions and
#' re-run orthofinder.
#'
#' @param dir.list The directory list produced by check_environment
#' @param map The map data.frame or data.table
#' @param blast data.table with at least the 10 necessary columns for blast format.
#' @param genomeIDs character vector giving genome IDs to consider.
#' @param gff The gff-like data.table or data.frame produced by
#' form_syntenicBlocks. Can also be made by hand - just a parsed gff
#' file with the following columns: 'id' (gene identifier), 'chr',
#' 'start', 'end', 'strand', 'genome' (matching an element in genomeIDs),
#' 'order' (gene order within that genome).
#' @param n.cores Number of parallel processes to run, when possible
#' @param rank.buffer The buffer, in gene rank order.
#' @param verbose Logical, should updates be printed
#' @param ... Not currently in use
#' @details None yet

#' @return A 4-element list of block, map, blast output and
#' orthofinder output.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom parallel mclapply
#' @export
find_synOrthos <- function(dir.list,
                           gff,
                           map,
                           genomeIDs,
                           n.cores = 1,
                           verbose = T,
                           run.gene.trees = F,
                           min.homology.score = 50,
                           min.homology.prop.of.best = 0.8,
                           min.perc.iden = NULL,
                           enforce.brkpts = TRUE,
                           silent.orthofinder = FALSE,
                           block.orthology.threshold = 0.5,
                           do.msa = FALSE,
                           blk.rank.buffer = 5,
                           rank.buffer = 50){

  if (verbose)
    cat("Preparing environment ...\n\tRe-making orthofinder input for culled data ... ")
  of.dir <- dir.list$syn.blast

  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = of.dir,
               peptide.dir = dir.list$peptide,
               genomeIDs = genomeIDs,
               verbose = F)

  if (verbose)
    cat("Done!\n\tLoading all blast files into memory ... ")
  blast <- read_allBlasts(gff = gff,
                          keep.geneNum = T,
                          add.gff = T,
                          check.ogs = F,
                          blast.dir = dir.list$cull.blast,
                          genomeIDs = genomeIDs,
                          verbose = F)
  blast.in <- data.table(blast)

  if(verbose)
    cat("Done!\n\tRe-ranking blast and map gene-order ... ")
  map <- map[!duplicated(map[,c("id1","id2")]),]
  ci <- c(intersect(colnames(map), colnames(blast)),"block.id")
  tmp <- rbind(data.table(map[,ci[ci %in% colnames(map)], with = F], type = "map"),
               data.table(blast[,ci[ci %in% colnames(blast)], with = F], type = "blast"),
               fill = T)
  tmp <- rerank_fromIDs(map = tmp, gff = gff, cull2genes.inmap = T)

  tmp2 <- make_blocks(subset(tmp, type == "map"))
  blast <- subset(tmp, type == "blast")
  map <- data.table(tmp2$map)
  blk <- data.table(tmp2$block)

  tp <- subset(map, block.id %in% b$block.id)
  with(tp, plot(rank1, rank2, col = as.numeric(as.factor(block.id))))

  if(verbose)
    cat("Done!\nExtending block coordinates by",rank.buffer,"... ")
  blast <- cull_blast2blk(blast = blast,
                          blk = blk,
                          blk.rank.buffer = blk.rank.buffer)
  tp <- subset(blast, block.id %in% b$block.id)
  with(tp, plot(rank1, rank2, col = as.numeric(as.factor(block.id))))

  if(verbose)
    cat("Done!\nPulling blast results hits within", rank.buffer,
        "hits of syntenic hits by block\n")
  blast.out <- cull_synBlast(map = map,
                             blast = blast,
                             blast.in = blast.in,
                             rank.buffer = rank.buffer,
                             verbose = verbose)
  tp <- subset(blast.out, block.id %in% b$block.id)
  with(tp, plot(rank1, rank2, col = as.numeric(as.factor(block.id))))

  if(verbose)
    cat("\tDone!\nPrepping blast for orthofinder ... ")
  nothing <- io_blast4of(of.dir = of.dir,
                         gff = gff,
                         blast = blast.out,
                         genomeIDs = genomeIDs)

  #######################################################
  if (verbose)
    cat("Done!\nRe-running orthofinder on culled blast hits ...\n")

  if (run.gene.trees) {
    if (do.msa) {
      com <- paste("orthofinder", "-b", of.dir,
                   "-a", n.cores, "-M msa")
    }else{
      com <- paste("orthofinder", "-b", of.dir,
                   "-a", n.cores)
    }
  }else{
    com <- paste("orthofinder", "-b", of.dir,
                 "-a", n.cores, "-og")
  }

  if (silent.orthofinder)
    com <- paste(com, "1>/dev/null 2>&1")

  system(com)

  if (verbose)
    cat("Done!\nCompiling results ... ")
  og <- read_ogs(of.dir, gff = gff)
  og1 <- data.table(id1 = og$id,
                    og1 = og$og)
  og2 <- data.table(id2 = og$id,
                    og2 = og$og)

  if ("og1" %in% colnames(blast.out))
    blast.out$og1 <- NULL
  if ("og2" %in% colnames(blast.out))
    blast.out$og2 <- NULL

  yo <- merge(og1, merge(og2, blast.out, by = "id2"), by = "id1")

  y.ortho <- subset(yo, og1 == og2)
  y.ortho$og.id <-  gsub(":", "", y.ortho$og1, fixed = T)
  y.ortho$og1 <- y.ortho$og2 <-  NULL
  if (verbose)
    cat("Done!\n")

  syn.homos <- pull_synHomos(
    syn.blast = blast.out,
    syn.ortho.map = y.ortho,
    min.score = min.homology.score,
    min.perc.iden = min.perc.iden,
    prop.of.best = min.homology.prop.of.best)

  if (run.gene.trees)
    y.ortho <- add_orthology(map = y.ortho,
                             of.dir = of.dir,
                             verbose = T,
                             block.orthology.threshold = block.orthology.threshold)

  y.ortho <- rerank_fromIDs(map = y.ortho,
                            gff = gff,
                            cull2genes.inmap = T)
  tp <- subset(y.ortho, block.id %in% b$block.id)
  with(tp, plot(rank1, rank2, col = as.numeric(as.factor(block.id))))
  out <- make_blocks(y.ortho,
                     rerank = F,
                     rename.blocks = F,
                     add.metadata = F,
                     clean.columns = F)

  return(list(map = out$map,
              blk = out$block,
              blast = blast.out,
              syn.homologs = syn.homos))
}
