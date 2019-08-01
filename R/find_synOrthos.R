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
                           enforce.brkpts = TRUE,
                           silent.orthofinder = FALSE,
                           block.orthology.threshold = 0.5,
                           rank.buffer = 100){

  if(verbose)
    cat("\nFormatting and returning block and map objects ... ")
  m <- data.table(map)
  m[,block.id := paste0("blk_",as.numeric(as.factor(paste(genome1, genome2, chr1, chr2, block.id))))]

  wh1 <- grep("1$", colnames(m))
  wh2 <- grep("2$", colnames(m))
  who <- which(!grepl("1$|2$", colnames(m)))

  n1 <- colnames(m)[wh1]
  n2 <- colnames(m)[wh2]
  no <- colnames(m)[who]
  m1 <- m[,c(n1,n2,no), with = F]
  m2 <- m[,c(n2,n1,no), with = F]
  setnames(m2, c(n1,n2,no))
  mt <- rbind(m1, m2)

  mo <- make_blocks(mt[!duplicated(mt[,c("id1","id2")]),], clean.columns = F)
  map <- mo$map
  blk <- mo$block

  if(verbose)
    cat("Done!\nPreparing environment ... \n")

  of.dir <- dir.list$syn.blast
  mb <- make_blocks(map, rename.blocks = F)
  map <- data.table(mb$map)
  blk <- data.table(mb$block)

  of.dir <- dir.list$syn.blast

  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = of.dir,
               peptide.dir = dir.list$peptide,
               genomeIDs = genomeIDs)

  blast <- read_allBlasts(gff = gff,
                          keep.geneNum = T,
                          add.gff = T,
                          check.ogs = F,
                          blast.dir = dir.list$cull.blast,
                          genomeIDs = genomeIDs,
                          verbose = T)


  spl.blast <- split(blast, by = c("genome1","genome2","chr1","chr2"))
  spl.blk <- split(blk, by = c("genome1","genome2","chr1","chr2"))
  blast <- rbindlist(lapply(names(spl.blk), function(i){
    y <- spl.blk[[i]]
    x <- spl.blast[[i]]

    xc <- rbindlist(lapply(1:nrow(y), function(j){
      int1s <- findInterval(x$start1, range(y[j,c("start1","end1")])) == 1
      int1e <- findInterval(x$end1, range(y[j,c("start1","end1")])) == 1
      int2s <- findInterval(x$start2, range(y[j,c("start2","end2")])) == 1
      int2e <- findInterval(x$end2, range(y[j,c("start2","end2")])) == 1
      wh <- which((int1s | int1e) & (int2s | int2e))
      if(length(wh)<1){
        return(NULL)
      }else{
        return(data.table(x[wh,], block.id = y$block.id[j]))
      }
    }))
    return(xc)
  }))


  map$rank1 <- NULL
  map$rank2 <- NULL
  map$what <- "map"
  blast$what <- "blast"

  # combn <- rbind(data.table(t(combn(genomeIDs,2))),
  #                rbindlist(lapply(genomeIDs, function(x)
  #                  data.table(V1 = x, V2 = x))))
  # setnames(combn, c("genome1","genome2"))
  # map.sim <- merge(combn, map, by = c("genome1","genome2"))
  # self.blk <- map.sim[,list(is.self = any(id1 %in% id2)),
  #                     by = list(block.id)]
  # self.blk <- self.blk$block.id[self.blk$is.self]
  # map.sim <- subset(map.sim, !block.id %in% self.blk)
  spl.map <- split(map, by = "block.id")
  spl.blast <- split(blast, by = "block.id")
  cn <- colnames(map)

  if(verbose)
    cat("Pulling blast results hits within", rank.buffer, "hits of syntenic hits by block\n")
  idl <- rbindlist(lapply(1:length(names(spl.map)), function(j){
    i = names(spl.map)[j]
    if(verbose)
      if(j %% 100 == 0)
        cat("\tCompleted",j,"/",length(spl.map),"blocks\n")

    y <- spl.map[[i]]
    x <- spl.blast[[i]][,cn, with = F]
    z <- rbind(y,x)
    z <- z[!duplicated(z[,c("id1","id2")]),]

    wh <- find_whichInBuffer(x = frank(z$start1, ties.method = "dense"),
                             y = frank(z$start2, ties.method = "dense"),
                             which.in.blk = which(z$what == "map"),
                             rank.buffer = rank.buffer)
    return(data.table(id1 = z$id1[wh],
                      id2 = z$id2[wh]))
  }))
  if(verbose)
    cat("\tDone!\n")

  blast.out <- merge(idl, blast, by = c("id1","id2"))
  blast.out$what <- NULL

  gi1 <- read_geneIDs(of.dir = of.dir, gff = gff)[,c("id","gene.num")]
  gi2 <- data.table(gi1)
  setnames(gi1, c("id1","gn1"))
  setnames(gi2, c("id2","gn2"))
  blast.out$gn2 <- blast.out$gn1 <- NULL
  blast.out <- merge(gi1, merge(gi2, blast.out, by = "id2"), by = "id1")

  #######################################################
  si <- read_speciesIDs(of.dir = of.dir, genomeIDs = genomeIDs)
  combs <- expand.grid(si$genome, si$genome)
  combs.n <- expand.grid(si$genome.num, si$genome.num)
  combs.file <- file.path(of.dir, paste0("Blast",combs.n[,1],"_", combs.n[,2],".txt"))
  cols2write <- c("gn1","gn2","perc.iden","align.length","n.mismatch","n.gapOpen",
                  "q.start","q.end","s.start","s.end","eval","score")

  for (i in 1:nrow(combs)) {
    tmp <- subset(blast.out,
                  genome1 == combs[i,1] &
                    genome2 == combs[i,2])[,cols2write, with = F]
    write.table(tmp, sep = "\t",
                file = combs.file[i],
                quote = F,
                col.names = F,
                row.names = F)

  }

  #######################################################
  if (verbose)
    cat("Re-running orthofinder on culled blast hits ...\n")
  if (verbose)
    cat("\tRunning orthofinder ... ")

  if (run.gene.trees) {
    com <- paste("orthofinder", "-b", of.dir,
                 "-a", n.cores)
  }else{
    com <- paste("orthofinder", "-b", of.dir,
                 "-a", n.cores,
                 "-og")
  }
  if(silent.orthofinder){
    com <- paste(com, "1>/dev/null 2>&1")
  }

  system(com)
  if (verbose)
    cat("Done!\n\tCompiling results ... ")
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
  y.ortho = subset(yo, og1 == og2)
  y.ortho$og.id <-  gsub(":","",y.ortho$og1, fixed = T)
  y.ortho$og1 <- y.ortho$og2 <-  NULL
  if (verbose)
    cat("Done!\n")


  syn.homos <- pull_synHomos(
    syn.blast = blast.out,
    syn.ortho.map = y.ortho,
    min.score = min.homology.score,
    prop.of.best = min.homology.prop.of.best)

  if(run.gene.trees){
    y.ortho <- add_orthology(map = y.ortho,
                             of.dir = of.dir,
                             verbose = T,
                             block.orthology.threshold = block.orthology.threshold)
  }
  out <- make_blocks(y.ortho,rerank = T,rename.blocks = F,add.metadata = F,clean.columns = F)
  return(list(map = out$map,
              blast = blast.out,
              blk = out$block,
              syn.homologs = syn.homos))
}
