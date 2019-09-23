#' @title Assign syntenic homologs
#'
#' @description
#' \code{assign_synHomologs} Assign syntenic homologs as paralogs, orthologs,
#' or just homologs.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param syn.blast data.table, optional, containing the
#' syntenic blast hits to feed to orthofinder.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL (default), taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param rank.buffer numeric, the radius to search within for syntenic
#' mappings.
#' @param min.homolog.score numeric, minimum blast bit score to be
#' included
#' @param n.cores integer length 1, the number of parallel processes
#' to run.
#' @param min.perc.iden numeric, minimum blast proprotion of sequence
#' identity to be included
#' @param min.prop.of.best numeric, minimum proportion of top score for
#' each gene to be included.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use

#' @details ...
#'
#' @return A blast data.table, with an additional column indicating
#' whether each is is homologous, paralogous or orthologous
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom igraph clusters graph_from_data_frame
#' @export
assign_synHomologs <- function(map,
                               syn.blast = NULL,
                               genomeIDs,
                               gff,
                               dir.list,
                               rank.buffer,
                               verbose,
                               min.homolog.score = 50,
                               min.perc.iden = NULL,
                               min.prop.of.best = .8,
                               n.cores = 1){

  add_orthology <- function(of.dir,
                            map,
                            verbose = T){
    if(verbose)
      cat("Pulling orthologs ... \n")
    m <- data.table(map[,c("genome1","genome2","id1","id2","og.id")])
    og.dir <- dirname(list.files(file.path(of.dir, "OrthoFinder"),
                                 pattern = "Orthologues_",recursive = T,
                                 include.dirs = T, full.names = T))[1]
    ortholog.dirs <- list.files(og.dir, pattern = "^Orthologues_",
                                include.dirs = T, full.names = T)
    og.out <- rbindlist(lapply(ortholog.dirs, function(x){
      fs <- list.files(x, full.names = T)
      if(verbose)
        cat("\tRunning", gsub("Orthologues_","",basename(x)), "... ")

      rd <- rbindlist(lapply(fs, function(y) {

        tmp <- fread(y)

        allgenes <- rbindlist(apply(tmp[,c(2:3),with = F],1,function(z){
          eg <- c(strsplit(z[1],",")[[1]],strsplit(z[2],", ")[[1]])
          return(expand.grid(eg, eg))
        }))
        setnames(allgenes,c("id1","id2"))
        allgenes[,genome1 := colnames(tmp)[2]]
        allgenes[,genome2 := colnames(tmp)[3]]

        return(allgenes)
      }))

      if(verbose)
        cat("Done!\n")
      return(rd)
    }))

    og.out[,is.ortholog := TRUE]
    map.out <- merge(og.out, map, by = colnames(og.out)[1:4], all.y = T)
    map.out$is.ortholog[map.out$id1 == map.out$id2] <- TRUE
    map.out$is.ortholog[is.na(map.out$is.ortholog)] <- FALSE

    return(map.out)
  }

  pull_synHomos <- function(syn.blast,
                            syn.ortho.map,
                            min.score,
                            prop.of.best,
                            min.perc.iden){

    if (is.na(min.score)) {
      sb <- subset(syn.blast, perc.iden >= min.perc.iden)

    }else{
      sb <- subset(syn.blast, score >= min.score)
    }

    sb[,propscore1 := score/max(score),
       by = list(genome2, chr2, id1)]
    sb[,propscore2 := score/max(score),
       by = list(genome1, chr1, id2)]
    sb <- subset(sb,
                 propscore1 > prop.of.best |
                   propscore2 > prop.of.best)
    mu <- with(syn.ortho.map, unique(paste(id1, id2)))
    so <- subset(sb, !paste(id1, id2) %in% mu)
    return(so)
  }


  io_blast4of <- function(of.dir, gff, blast, genomeIDs){
    gi1 <- read_geneIDs(of.dir = of.dir, gff = gff)[,c("id","gene.num")]
    gi2 <- data.table(gi1)
    setnames(gi1, c("id1","gn1"))
    setnames(gi2, c("id2","gn2"))
    blast$gn2 <- blast$gn1 <- NULL
    blast <- merge(gi1, merge(gi2, blast, by = "id2"), by = "id1")

    #######################################################
    si <- read_speciesIDs(of.dir = of.dir, genomeIDs = genomeIDs)
    combs <- expand.grid(si$genome, si$genome)
    combs.n <- expand.grid(si$genome.num, si$genome.num)
    combs.file <- file.path(of.dir, paste0("Blast",combs.n[,1],"_", combs.n[,2],".txt"))
    cols2write <- c("gn1","gn2","perc.iden","align.length","n.mismatch","n.gapOpen",
                    "q.start","q.end","s.start","s.end","eval","score")
    for (i in 1:nrow(combs)) {
      tmp <- subset(blast,
                    genome1 == combs[i,1] &
                      genome2 == combs[i,2])[,cols2write, with = F]
      write.table(tmp, sep = "\t",
                  file = combs.file[i],
                  quote = F,
                  col.names = F,
                  row.names = F)
    }
  }

  if (is.null(syn.blast)) {
    if(verbose)
      cat("Extracting syntenic blast hits for each unique chromosome combination ... \n")
    syn.blast <- extend_blocks(
      gff = gff,
      genomeIDs = genomeIDs,
      map = map,
      dir.list = dir.locs,
      use.score.cull.blast = F,
      rank.buffer = rank.buffer,
      verbose = F)
  }else{
    cat("Using previously calculated syntenic blast results ... ")
  }
  if(verbose)
    cat("\tDone!\nLoading all blast hits into memory ... ")
  all.blast <- read_allBlasts(
    gff = gff,
    genomeIDs = genomeIDs,
    of.dir = dir.locs$cull.blast,
    check.ogs = F,
    add.gff = T,
    verbose = F,
    keep.geneNum = T)

  if(verbose)
    cat("Done!\nBuilding new orthofinder database ... ")
  sb <- with(syn.blast,
             rbind(data.table(id1 = id1, id2 = id2),
                   data.table(id1 = id2, id2 = id1)))
  sb <- sb[!duplicated(sb),]
  syn.all.blast <- merge(all.blast, sb, by = c("id1","id2"))

  of.dir <- dir.list$syn.blast
  make_newOFdb(tmp.dir = dir.list$tmp,
               cull.blast.dir = of.dir,
               peptide.dir = dir.list$peptide,
               genomeIDs = genomeIDs,
               verbose = F)
  if(verbose)
    cat("Done!\nPrepping blast for orthofinder ... ")
  nothing <- io_blast4of(of.dir = of.dir,
                         gff = gff,
                         blast = syn.all.blast,
                         genomeIDs = genomeIDs)
  if(verbose)
    cat("Done!\nRunning orthofinder ... \n#######################\n")
  com <- paste("orthofinder", "-b", of.dir,
               "-a", n.cores)
  system(com)

  if (verbose)
    cat("\n#######################\nDone!\nCompiling orthogroups ... ")
  og <- read_ogs(of.dir, gff = gff)
  og1 <- data.table(id1 = og$id,
                    og1 = og$og)
  og2 <- data.table(id2 = og$id,
                    og2 = og$og)

  yo <- merge(og1, merge(og2, syn.all.blast, by = "id2"), by = "id1")

  y.ortho <- subset(yo, og1 == og2)
  y.ortho$og.id <-  gsub(":", "", y.ortho$og1, fixed = T)
  y.ortho$og1 <- y.ortho$og2 <-  NULL
  if(verbose)
    cat("Done!\nPulling syntenic homologs ... ")
  syn.homos <- pull_synHomos(
    syn.blast = syn.blast,
    syn.ortho.map = y.ortho,
    min.score = min.homolog.score,
    min.perc.iden = min.perc.iden,
    prop.of.best = min.prop.of.best)
  if(verbose)
    cat("Done!\nSplitting orthologs and paralogs ... ")
  syn.homos[,hit.type := "syntenic.homolog"]
  y.out <- add_orthology(of.dir = of.dir, map = y.ortho)
  y.out[,hit.type := ifelse(is.ortholog, "ortholog","paralog")]
  y.out <- y.out[,colnames(y.out) %in% colnames(syn.homos), with = F]
  syn.homos <- syn.homos[,colnames(syn.homos) %in% colnames(y.out), with = F]

  out <- mirror_map(map = rbind(y.out, syn.homos),
                    keycol = "hit.type")

  if(verbose)
    cat("Done!\n")
  return(out)
}
