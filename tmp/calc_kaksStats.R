#' @title Extend and complete syntenic mappings
#'
#' @description
#' \code{calc_kaksStats} Extend and complete syntenic mappings.
#'
#' @param map data.table, containing the merged gff and blast results
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param cds.align.files vector of file paths, pointing to the CDS
#' alignment files
#' @param infer.consensus logical, should consensus sequences be
#' calculated?
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' If NULL default, taken as all unique elements in the 'genome' column
#' of the gff data.table.
#' @param ogs2test vector of orthogroup ids to consider
#' @param delimer charcter string, length 1, to paste, then split
#' genomeIDs and geneIDs together for alignment IDs
#' @param batch.size numeric, length 1, how many ogs should be analyzed
#' simultaneously
#' @param mafft.params character vector, length 1, specifying parameters
#' to paass to mafft. Highly recommended to no use --threads > 1.
#' @param n.cores numeric, length 1, number of parallel processes to run
#' @param min.cdsMSA.length numeric, length 1, minimum length of CDS
#' alignments to calculate stats.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#'
#' @details ...
#'
#' @return A blast data.table, without block.id information
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom seqinr kaks read.alignment reverse.align consensus
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet
#' @importFrom parallel mclapply
#' @export
calc_kaksStats <- function(map,
                           dir.list,
                           genomeIDs,
                           cds.align.files = NULL,
                           min.cdsMSA.length = 100,
                           n.cores = 1,
                           just.do.alignments = F,
                           ogs2test = NULL,
                           delimer = NULL,
                           batch.size = NULL,
                           infer.consensus = TRUE,
                           mafft.params = "--retree 2 --maxiterate 0",
                           verbose = TRUE){

  ##########################################################
  ##########################################################
  run_seqinrStats <- function(alg){
    kaks.out <- kaks(alg,  verbose = T, rmgap = F)

    cols <- c("ka","ks","vka","vks",
              "a0","a2","a4",
              "b0","b2","b4",
              "l0","l2","l4")
    coln <- c("ka","ks","var.ka","var.ks",
              "tr.0f","tr.2f","tr.4f",
              "tv.0f","tv.2f","tv.4f",
              "freq.0f","freq.2f","freq.4f")

    kaks.comb <- rbindlist(lapply(1:length(cols), function(i)
      data.table(
        melt(
          data.matrix(
            kaks.out[[cols[i]]])),
        test = coln[i])))

    kaks.comb <- dcast(kaks.comb,
                       Var1 + Var2 ~ test,
                       value.var = "value")
    kaks.comb <- subset(kaks.comb, Var1 != Var2)
    setnames(kaks.comb, 1:2, c("genome_id1","genome_id2"))
    kaks.comb[,kaks := ka / ks]
    kaks.comb[,var.kaks := abs(kaks) * sqrt((var.ka/ka)^2 + (var.ks/ks)^2)]
    return(kaks.comb)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  make_ugenesByOG <- function(map){
    og.dt <- with(subset(map, !is.na(og)), rbind(
      data.table(genome = as.character(genome1),
                 id = id1,
                 og = og),
      data.table(genome = as.character(genome2),
                 id = id2,
                 og = og)))
    og.dt <- og.dt[!duplicated(og.dt),]
    og.dt[,genome_id := paste0(genome, "_", id)]
    spl <- split(og.dt$genome_id, og.dt$og)
    return(spl)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  align_cds2pep <- function(pep.unal.fa,
                            pep.align.fa,
                            cds.unal.fa,
                            cds.align.fa,
                            mafft.params){

    if (grepl("--thread",mafft.params, fixed = T) &
        grepl("--quiet", mafft.params, fixed = T)) {
      mc <- "mafft"
    }else{
      if (grepl("--thread",mafft.params, fixed = T)) {
        mc <- "mafft --quiet"
      }else{
        if (grepl("--quiet", mafft.params, fixed = T)) {
          mc <- "mafft --thread 1"
        }else{
          mc <- "mafft --thread 1 --quiet"
        }
      }
    }

    mafft.com <- paste(
      mc, mafft.params,
      pep.unal.fa, ">", pep.align.fa)
    system(mafft.com)

    nothin <- reverse.align(
      nucl.file = cds.unal.fa,
      protaln.file = pep.align.fa,
      out.file = cds.align.fa,
      forceDNAtolower = F)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  make_alignFilePaths <- function(align.dir,
                                  og.ids,
                                  batch.size){
    file.loc <- data.table(og = og.ids)
    file.loc[,pep.unal.file := file.path(align.dir, paste0(og,".pep.unal.fa"))]
    file.loc[,pep.align.file := file.path(align.dir, paste0(og,".pep.align.fa"))]
    file.loc[,cds.unal.file := file.path(align.dir, paste0(og,".cds.unal.fa"))]
    file.loc[,cds.align.file := file.path(align.dir, paste0(og,".cds.align.fa"))]
    file.loc[,batch.id := ceiling(1:nrow(file.loc)/batch.size)]
    return(split(file.loc, by = "batch.id"))
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  pipe_kaks <- function(cds.algns,
                        batch.size,
                        n.cores){
    batches <- ceiling(1:length(cds.algns) / batch.size)
    cds.algns.batch <- split(cds.algns, batches)
    selstats <- lapply(1:length(cds.algns.batch), function(i){
      out <- mclapply(
        cds.algns.batch[[i]],
        mc.cores = n.cores,
        run_seqinrStats)
      if (verbose)
        cat("\tCompleted batch", i,
            "/", length(cds.algns.batch), "\n")
      return(out)
    })
    outss <- unlist(selstats, recursive = F)
    return(outss)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  add_consensus <- function(algn){
    con <- paste(
      consensus(algn,
                method = "majority"),
      collapse = "")
    algn$nb <- algn$nb + 1
    algn$nam <- c(algn$nam, "consensus")
    algn$seq <- c(algn$seq, con)
    return(algn)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  import_revAlignCDS <- function(align.paths,
                                 cds.align.files = NULL,
                                 n.cores,
                                 batch.size,
                                 min.cdsMSA.length){
    if (is.null(cds.align.files)) {
      align.paths[,align.worked := file.exists(cds.align.file)]
      align.paths <- subset(align.paths, align.worked)
      ap <- align.paths$cds.align.file
    }else{
      align.worked <- file.exists(cds.align.files)
      ap <- cds.align.files[align.worked]
    }

    cds.algns <- mclapply(
      ap,
      mc.cores = n.cores,
      function(x)
        read.alignment(x,
                       format = "fasta"))

    algns.len <- sapply(cds.algns, function(x) nchar(x$seq[[1]]))
    algns.is <- sapply(cds.algns, function(x) substr(x$seq[[1]],1,1) == ">")
    algns.len[algns.is] <- 0
    cds.algns2test <- cds.algns[algns.len >= min.cdsMSA.length]
    algns.len <- sapply(cds.algns2test, function(x) nchar(x$seq[[1]]))

    return(cds.algns2test[order(-algns.len)])
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  convert_gid2map <- function(g_id.map, delimer){
    g_id.map[,genome1 := sapply(as.character(genome_id1), function(x)
      strsplit(x, delimer, fixed = T)[[1]][1])]
    g_id.map[,genome2 := sapply(as.character(genome_id2), function(x)
      strsplit(x, delimer, fixed = T)[[1]][1])]
    g_id.map[,id1 := sapply(as.character(genome_id1), function(x)
      strsplit(x, delimer, fixed = T)[[1]][2])]
    g_id.map[,id2 := sapply(as.character(genome_id2), function(x)
      strsplit(x, delimer, fixed = T)[[1]][2])]
    g_id.map[,genome_id1 := NULL]
    g_id.map[,genome_id2 := NULL]
    return(g_id.map)
  }
  ##########################################################
  ##########################################################

  ##########################################################
  ##########################################################
  align_cdsPep <- function(align.dir,
                           cds.dir,
                           pep.dir,
                           map,
                           genomeIDs,
                           genome.id.delimer,
                           n.cores,
                           mafft.params,
                           batch.size){
    if (verbose)
      cat("Building file paths for alignments ... ")
    ug.list <- make_ugenesByOG(map = map)
    file.paths <- make_alignFilePaths(
      og.ids = names(ug.list),
      align.dir = align.dir,
      batch.size = batch.size)

    if (verbose)
      cat("Done!\nReading in peptide and cds sequences ... ")

    peps <- do.call(c, lapply(genomeIDs, function(x){
      y <- readAAStringSet(file.path(pep.dir, paste0(x,".fa")))
      names(y) <- paste0(x,genome.id.delimer,names(y))
      return(y)
    }))
    cdss <- do.call(c, lapply(genomeIDs, function(x){
      y <- readAAStringSet(file.path(cds.dir, paste0(x,".fa")))
      names(y) <- paste0(x,genome.id.delimer,names(y))
      return(y)
    }))

    if (verbose)
      cat("Done!\nWriting sequences for",
          length(file.paths), "groups of",
          batch.size,"orthogroups\n")
    nothing <- lapply(1:length(file.paths), function(i){
      x <- file.paths[[i]]
      sl <- split(x, by = "og")
      og.ids <- x$og
      ugl <- ug.list[og.ids]
      uniq.genes <- unique(unlist(ugl))
      cds.tmp <- cdss[uniq.genes]
      pep.tmp <- peps[uniq.genes]
      for (j in 1:nrow(x)) {
        o <- x$og[j]
        u <- ugl[[o]]
        l <- sl[[o]]
        co <- cds.tmp[u]
        po <- pep.tmp[u]
        writeXStringSet(co, filepath = l$cds.unal.file)
        writeXStringSet(po, filepath = l$pep.unal.file)
      }
      if (verbose)
        cat("\tCompleted", i,"/",length(file.paths),"batches\n")
    })

    if (verbose)
      cat("\tDone!\nAligning peptide sequences ...\n")

    dont <- lapply(1:length(file.paths), function(i){
      if (verbose)
        cat("\tRunning alignments on", i, "/",
            length(file.paths)," batches\n")
      x <- file.paths[[i]]
      don <- mclapply(
        1:nrow(x),
        mc.cores = n.cores,
        mc.preschedule = F,
        function(j){
          align_cds2pep(
            pep.unal.fa = x$pep.unal.file[j],
            pep.align.fa = x$pep.align.file[j],
            cds.unal.fa = x$cds.unal.file[j],
            cds.align.fa = x$cds.align.file[j],
            mafft.params = mafft.params)
        })
    })
    align.paths <- rbindlist(file.paths)
    align.paths[,batch.id := NULL]
    return(align.paths)
  }
  ##########################################################
  ##########################################################

  map <- subset(map, genome1 %in% genomeIDs & genome2 %in% genomeIDs)

  ug <- unique(unlist(map[,c("genome1","genome2","id1","id2")]))
  if (is.null(delimer))
    delimer <- ifelse(
      !any(grepl("_",ug, fixed = T)),
      "_",
      ifelse(
        !any(grepl("|",ug, fixed = T)),
        "|",
        ifelse(
          !any(grepl("_xxxx_",ug, fixed = T)),
          "_xxxx_",
          NA)))

  if(is.null(cds.align.files)){
    ##########################################################
    pep.dir <- dir.list$peptide
    cds.dir <- dir.list$cds
    results.dir <- dir.list$results
    align.dir <- file.path(results.dir, "alignments")
    if (!dir.exists(align.dir))
      dir.create(align.dir)

    if (!"og" %in% colnames(map))
      stop("Can't find the orthogroup (og) column in map, has assign_homologs been run?\n")

    if (!is.null(ogs2test)) {
      if (!any(ogs2test %in% map$og))
        stop("ogs2test must represent a vector of og ids, which are stored in the map og column\n")
      map <- subset(map, og %in% ogs2test)
    }

    if (is.na(delimer))
      stop("choose a string for delimer that is not found in genomeIDs or geneIDs\n")

    nog <- length(unique(map$og))
    if (is.null(batch.size))
      batch.size <- nog / 50

    ##########################################################
    align.paths <- align_cdsPep(
      align.dir = align.dir,
      pep.dir = pep.dir,
      cds.dir = cds.dir,
      genomeIDs = genomeIDs,
      genome.id.delimer = delimer,
      map = map,
      mafft.params = mafft.params,
      n.cores = n.cores,
      batch.size = batch.size)
    cds.align.files <- align.paths$cds.align.file
    if(verbose)
      cat("\tDone!\n")
  }else{
    align.paths <- NULL
  }

  if(just.do.alignments){
    out.ss <- NULL
  }else{
    ##########################################################
    if (verbose)
      cat("Reading CDS alignments into memory ... ")
    cds.algns2test <- import_revAlignCDS(
      cds.align.files = cds.align.files,
      align.paths = NULL,
      n.cores = n.cores,
      batch.size = batch.size,
      min.cdsMSA.length = min.cdsMSA.length)

    if (infer.consensus) {
      if(verbose)
        cat("Done!\nInferring consensus ancestral sequence ...\n")
      batchv <- ceiling(1:length(cds.algns2test)/batch.size)
      spl <- split(cds.algns2test, batchv)
      con.out <- unlist(lapply(1:length(spl), function(i){
        if(verbose)
          cat("\tRunnning batch",i,"/",length(spl),"\n")
        con.tmp <- mclapply(
          spl[[i]],
          mc.cores = n.cores,
          add_consensus)
        return(con.tmp)
      }), recursive = F)
      cds.algns2test <- con.out
      conseq <- lapply(con.out, function(x) x$seq[length(x$seq)])
    }else{
      conseq <- NULL
    }

    ##########################################################
    if (verbose)
      cat("\tDone!\nCalculating selection stats ... \n")
    kaks.out <- rbindlist(
      pipe_kaks(
        cds.algns = cds.algns2test,
        batch.size = batch.size,
        n.cores = 4))
    kaks.out <- subset(kaks.out, genome_id1 != "consensus")
    kaks.out[,genome_id2 := gsub("consensus",paste0(delimer,"consensus"), genome_id2)]

    ##########################################################
    if (verbose)
      cat("\tDone!\nProcessing output ... ")
    all.selstats <- convert_gid2map(
      kaks.out,
      delimer = delimer)

    if (infer.consensus) {
      all.constats <- subset(all.selstats, id2 == "consensus")
      all.selstats <- subset(all.selstats, id2 != "consensus")
    }

    out.ss <- merge(
      map,
      all.selstats,
      by = c("genome1", "genome2", "id1", "id2"))


    if (infer.consensus) {
      cno <- colnames(out.ss)[!colnames(out.ss) %in% colnames(all.constats)]
      out.cs1 <- merge(
        out.ss[,c("genome1","id1",cno), with = F],
        all.constats,
        by = c("genome1","id1"))

      cnnull <- colnames(out.cs1)[-grep("1$", colnames(out.cs1))]
      cnnull <- cnnull[!cnnull %in% c("id2","og",colnames(all.constats))]

      for(i in cnnull)
        out.cs1[[i]] <- ""
      out.cs1 <- out.cs1[!duplicated(out.cs1),]
      out.ss <- rbind(out.ss, out.cs1)
    }
  }


  if (verbose)
    cat("Done!\n")
  return(list(out.ss, align.paths))
}
