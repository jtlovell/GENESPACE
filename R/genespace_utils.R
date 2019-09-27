#' @title GENESPACE utility functions
#' @description
#' \code{genespace_utils} Functions for GENESPACE
#' @name genespace_utils

#' @param map data.table, containing the merged gff and blast results
#' @param n.mappings numeric, the number of mappings required within a
#' given radius. Length must match that of radius
#' @param radius numeric, the radius to search within for syntenic
#' mappings. Length must match that of n.mappings.
#' @param blast data.table, containing the merged gff and blast results.
#' Unlike the 'map' object, which really just needs merged gff coordinates,
#' this must have all the blast8 columns. See details.
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param dir.list list, containing the paths to various genespace
#' subdirectories.
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#' @param id.db genome ID dictionary
#' @param cull.blast.dir file.path, to the subdirectory
#' where the culled blast results should be written
#' @param blast.file.in file.path, to the blast file to read in
#' @param blast.file.out file.path, to the filename to which the
#' blast results should be written
#' @param maxn integer length 1, indicating the maximum number of
#' hits to be retained per query.
#' @param of.dir file.path, to the subdirectory containing
#' blast results
#' @param gene.dict1  data.table with column names V1, new1, keyed on
#' V1. V1 contains the original geneIDs, new1 contains the new geneIDs
#' @param gene.dict2  data.table with column names V2, new2, keyed on
#' V2. V2 contains the original geneIDs, new2 contains the new geneIDs
#' @param min.score numeric length 1, specifying the minimum blast
#' bit score for a hit to be retained.
#' @param n.cores integer length 1, the number of parallel processes
#' to run.
#' @param tmp.dir file.path, to the subdirectory
#' where the temporary results should be written
#' @param peptide.dir file.path, to the subdirectory containing
#' the parsed peptide files
#'
#'
#' @note \code{genespace_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{genespace_utils} returns its own arguments.
#'
#'


#' @title reduce_recipBlast
#' @description
#' \code{reduce_recipBlast} reduce_recipBlast
#' @rdname genespace_utils
#' @import data.table
#' @export
read_ofBlast <- function(gids,
                         blast.files){
  all.blast <- rbindlist(lapply(blast.files, function(x)
    fread(x,
          col.names = c("gn1", "gn2", "perc.iden", "align.length",
                        "n.mismatch", "n.gapOpen", "q.start",
                        "q.end", "s.start",
                        "s.end", "eval", "score"),
          key = "gn2")))
  if("gene.num" %in% colnames(gids))
    setnames(gids, "gene.num","gn")
  gids1 <- data.table(gids)
  setnames(gids1, paste0(colnames(gids),"1"))
  gids2 <- data.table(gids)
  setnames(gids2, paste0(colnames(gids),"2"))
  all.blast <- merge(gids1,
                     merge(gids2,
                           all.blast,
                           by = "gn2"),
                     by = "gn1")
  return(all.blast)
}

#' @title reduce_recipBlast
#' @description
#' \code{reduce_recipBlast} reduce_recipBlast
#' @rdname genespace_utils
#' @import data.table
#' @export
reduce_recipBlast <- function(genomeIDs,
                              blast,
                              intergenome.only = F){
  ug <- data.table(t(combn(genomeIDs,2)))
  setnames(ug, c("genome1","genome2"))

  if(!intergenome.only)
    ug <- rbind(ug,
              data.table(genome1 = genomeIDs,
                         genome2 = genomeIDs))
  out <- merge(blast,
               ug,
               by = c("genome1","genome2"))
  return(out)
}

#' @title pull_orthologs
#' @description
#' \code{pull_orthologs} pull_orthologs
#' @rdname genespace_utils
#' @import data.table
#' @export
pull_orthologs <- function(of.dir){
  ortholog.dirs <- list.files(
    file.path(of.dir, "OrthoFinder"),
    pattern = "Orthologues_",
    recursive = T,
    include.dirs = T,
    full.names = T)
  og.files <- unlist(
    lapply(ortholog.dirs, list.files, full.names = T))

  if(length(og.files) == 0)
    stop("Could not find orthologue files ... did you do a full orthofinder run?\n")

  og.ids <- rbindlist(lapply(og.files, function(x){
    og2 <- fread(x)
    gn <- colnames(og2)
    mel <- melt(og2, id.vars = "Orthogroup")
    setnames(mel, 1:2, c("og","genome"))
    lon <- mel[,list(id = unlist(sapply(value, function(x) strsplit(x, ", ")))),
               by = list(og, genome)]
    l1 <- with(subset(lon, genome == gn[2]),
               data.table(og = og,
                          genome1 = genome,
                          id1 = id))
    l2 <- with(subset(lon, genome == gn[3]),
               data.table(og = og,
                          genome2 = genome,
                          id2 = id))
    mer <- merge(l1,
                 l2,
                 by = "og",
                 allow.cartesian = T)
    mer[,is.ortholog := TRUE]
    return(mer)
  }))

  og.ids <- og.ids[!duplicated(og.ids[,c("genome1","genome2","id1","id2")]),]
  return(og.ids)
}


#' @title write_ofBlast2file
#' @description
#' \code{write_ofBlast2file} write_ofBlast2file
#' @rdname genespace_utils
#' @import data.table
#' @export
write_ofBlast2file <- function(blast,
                               output.dir){
  cols2write <- c("gn1","gn2","perc.iden","align.length","n.mismatch","n.gapOpen",
                  "q.start","q.end","s.start","s.end","eval","score")
  if (!all(c("genome.num1", "genome.num2", cols2write) %in% colnames(blast)))
    stop("blast must have genome.num1 and genome.num2 columns\n")

  spl <- split(blast,
               by = c("genome.num1",
                      "genome.num2"))

  ns <- names(spl)[order(names(spl))]

  names(ns) <- paste0("Blast",
                      gsub(".", "_",ns, fixed = T),
                      ".txt")

  for (i in 1:length(ns)) {
    n <- ns[i]
    fn <- names(ns)[i]
    write.table(
      spl[[n]][,cols2write, with = F],
      sep = "\t",
      file = file.path(output.dir, fn),
      quote = F,
      col.names = F,
      row.names = F)
  }
}

#' @title merge_gffWithBlast
#' @description
#' \code{merge_gffWithBlast} merge_gffWithBlast
#' @rdname genespace_utils
#' @import data.table
#' @export
merge_gffWithBlast <- function(gff,
                               blast,
                               include.gene.num,
                               include.genome.num,
                               mirror){
  if(mirror)
    blast <- mirror_map(blast)

  if(include.gene.num & include.genome.num){
    g1 <- with(gff, data.table(
      genome1 = genome,
      id1 = id,
      genome.num1 = genome.num,
      gn1 = gene.num))

    g2 <- with(gff, data.table(
      genome2 = genome,
      id2 = id,
      genome.num2 = genome.num,
      gn2 = gene.num))

    if("genome.num1" %in% colnames(blast))
      blast[,genome.num1 := NULL]
    if("genome.num2" %in% colnames(blast))
      blast[,genome.num2 := NULL]
    if("gn1" %in% colnames(blast))
      blast[,gn1 := NULL]
    if("gn2" %in% colnames(blast))
      blast[,gn2 := NULL]

  }else{
    if(include.gene.num){
      g1 <- with(gff, data.table(
        genome1 = genome,
        id1 = id,
        gn1 = gene.num))

      g2 <- with(gff, data.table(
        genome2 = genome,
        id2 = id,
        gn2 = gene.num))

      if("gn1" %in% colnames(blast))
        blast[,gn1 := NULL]
      if("gn2" %in% colnames(blast))
        blast[,gn2 := NULL]

    }else{
      if(include.genome.num){
        g1 <- with(gff, data.table(
          genome1 = genome,
          id1 = id,
          genome.num1 = genome.num))

        g2 <- with(gff, data.table(
          genome2 = genome,
          id2 = id,
          genome.num2 = genome.num))

        if("genome.num1" %in% colnames(blast))
          blast[,genome.num1 := NULL]
        if("genome.num2" %in% colnames(blast))
          blast[,genome.num2 := NULL]

      }else{
        g1 <- with(gff, data.table(
          genome1 = genome,
          id1 = id))

        g2 <- with(gff, data.table(
          genome2 = genome,
          id2 = id))
      }
    }
  }


  blast <- merge(
    g1,
    merge(
      g2,
      blast,
      by = c("genome2","id2")),
    by = c("genome1","id1"))

  blast <- mirror_map(blast)
  return(blast)
}

#' @title prep_ofDB
#' @description
#' \code{prep_ofDB} prep_ofDB
#' @rdname genespace_utils
#' @import data.table
#' @export
prep_ofDB <- function(tmp.dir,
                      orig.dir,
                      peptide.dir,
                      output.dir,
                      n.cores,
                      genomeIDs,
                      verbose,
                      min.score,
                      copy2dir = NULL){
  if (verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")
  make_newOFdb(
    tmp.dir = tmp.dir,
    output.dir = output.dir,
    peptide.dir = peptide.dir,
    verbose = verbose,
    n.cores = n.cores,
    genomeIDs = genomeIDs)

  fs <- list.files(path = output.dir,
                   full.names = T)
  if (!is.null(copy2dir))
    if (dir.exists(copy2dir))
    cpd <- file.copy(from = fs,
                     to = copy2dir)
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = orig.dir,
                             genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = output.dir,
                             genomeIDs = genomeIDs)

  id.db <- merge(old.ids, new.ids, by = "genome")
  setnames(id.db,2:3,c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       of.dir = orig.dir,
                       cull.blast.dir = output.dir)
  #
  old.genes <- read_geneIDs(of.dir = orig.dir,
                            gff = gff,
                            species.num.id = old.ids)
  new.genes <- read_geneIDs(of.dir = output.dir,
                            gff = gff,
                            species.num.id = new.ids)
  genes <- merge(old.genes[,c("genome","id","gene.num")],
                 new.genes[,c("genome","id","gene.num")],
                 by = c("genome","id"))
  g1 <- with(genes,
             data.table(V1 = gene.num.x,
                        new1 = gene.num.y,
                        key = "V1"))
  g2 <- with(genes,
             data.table(V2 = gene.num.x,
                        new2 = gene.num.y,
                        key = "V2"))
  if (verbose)
    cat("\tDone!\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")

  d <- lapply(1:nrow(map.db), function(i){
    if (verbose)
      cat(paste0("\t", map.db$genome1[i]),
          "-->", map.db$genome2[i], "... ")
    bl <- readRename_blastGenes(
      gene.dict1 = g1,
      gene.dict2 = g2,
      verbose = verbose,
      min.score = min.score,
      blast.file.in = map.db$filename[i],
      blast.file.out = map.db$new.filename[i])
    if (verbose)
      cat("Done!\n")
  })
  return(map.db)
}

#' @title find_ofFiles
#' @description
#' \code{find_ofFiles} find_ofFiles
#' @rdname genespace_utils
#' @import data.table
#' @export
find_ofFiles <- function(of.dir){
  blast.loc <- dirname(
    list.files(
      of.dir,
      pattern = "SequenceIDs",
      recursive = T,
      full.names = T)[1])

  ortho.loc <- dirname(
    list.files(
      of.dir,
      pattern = "Orthogroups.txt",
      recursive = T,
      full.names = T)[1])

  blast.files <- list.files(
    blast.loc,
    pattern = "Blast*",
    full.names = T)

  fa.files <- list.files(
    blast.loc,
    pattern = "Species*",
    full.names = T)

  fa.files <- fa.files[grep(".fa$", fa.files)]

  dmnd.files <- list.files(
    blast.loc,
    pattern = "diamondDBSpecies*",
    full.names = T)

  og.files <- file.path(
    ortho.loc,
    "Orthogroups.txt")
  sp.id.files <- file.path(
    blast.loc,
    "SpeciesIDs.txt")
  seq.id.files <- file.path(
    blast.loc,
    "SequenceIDs.txt")

  files <- c(blast.files,
             fa.files,
             dmnd.files,
             og.files,
             sp.id.files,
             seq.id.files)
  return(files)
}

#' @title check_gsDir
#' @description
#' \code{check_gsDir} check_gsDir
#' @rdname genespace_utils
#' @import data.table
#' @export
check_gsDir <- function(dir2check, overwrite.output.dir){
  if (!dir.exists(dir2check)) {
    dir.create(dir2check)
  }else{
    if (overwrite.output.dir) {
      unlink(dir2check, recursive = T)
      dir.create(dir2check)
    }else{
      if (length(dir(dir2check)) > 0)
        stop(dir2check, "is not empty and overwrite.output.dir = F\n\t",
             "Change the output directory if you wish to retain these files\n\t",
             "Set overwrite.output.dir = TRUE to overwrite this directory\n")
    }
  }
}


#' @title read_ogs
#' @description
#' \code{read_ogs} fread_ogs
#' @rdname genespace_utils
#' @import data.table
#' @export
read_ogs <- function(of.dir,
                     gff){

  ortho.loc <- dirname(list.files(of.dir,
                                  pattern = "Orthogroups.tsv",
                                  recursive = T,
                                  full.names = T)[1])
  og2 <- fread(file.path(ortho.loc,
                         "Orthogroups.tsv"))
  gn <- colnames(og2)
  mel <- melt(og2, id.vars = "Orthogroup")
  setnames(mel, 1:2, c("og","genome"))
  lon <- mel[,list(id = unlist(sapply(value, function(x) strsplit(x, ", ")))),
             by = list(og, genome)]
  l1 <- with(lon, data.table(og = og, genome1 = genome, id1 = id))
  l2 <- with(lon, data.table(og = og, genome2 = genome, id2 = id))
  mer <- merge(l1, l2, by = "og", allow.cartesian = T)
  return(mer)
}

#' @title make_newOFdb
#' @description
#' \code{make_newOFdb} make_newOFdb
#' @rdname genespace_utils
#' @import data.table
#' @export
make_newOFdb <- function(tmp.dir,
                         output.dir,
                         peptide.dir,
                         genomeIDs,
                         verbose = T,
                         n.cores = 1){
  if (verbose)
    cat("\tCleaning tmp dir ... ")
  unlink(tmp.dir, recursive = T)
  dir.create(tmp.dir)

  unlink(output.dir, recursive = T)
  dir.create(output.dir)

  if (verbose)
    cat("Done!\n\tCopying peptide fastas ... ")
  cpd <- file.copy(from = file.path(peptide.dir,
                                    paste0(genomeIDs,".fa")),
                   to = tmp.dir)

  if (verbose)
    cat("Done!\n\tRe-making orthofinder input format ... ")
  system(paste("orthofinder", "-f", tmp.dir,
               "-a", n.cores, "-S diamond",
               "-op >/dev/null 2>&1"))

  if (verbose)
    cat("Done!\n\tMoving results to cull.blast directory ... ")
  blast.loc <- dirname(
    list.files(tmp.dir,
               pattern = "SequenceIDs",
               recursive = T,
               full.names = T)[1])
  fa.files <- list.files(path = blast.loc,
                         pattern = "Species*",
                         full.names = T)
  fa.files <- fa.files[grep(".fa$", fa.files)]
  dmnd.files <- list.files(path = blast.loc,
                           pattern = "diamondDBSpecies*",
                           full.names = T)

  sp.id.files <- file.path(blast.loc, "SpeciesIDs.txt")
  seq.id.files <- file.path(blast.loc, "SequenceIDs.txt")
  files <- c(fa.files,
             dmnd.files,
             sp.id.files,
             seq.id.files)
  nu <- file.copy(files,
                  output.dir)
  if (verbose)
    cat("Done!\n")
}


#' @title readRename_blastGenes
#' @description
#' \code{readRename_blastGenes} readRename_blastGenes
#' @rdname genespace_utils
#' @import data.table
#' @export
readRename_blastGenes <- function(gene.dict1,
                                  gene.dict2,
                                  blast.file.in,
                                  blast.file.out,
                                  min.score = NULL,
                                  verbose = T){
  d <- fread(blast.file.in, key = "V2")
  if (verbose)
    cat("Processing", nrow(d), "raw hits ... ")
  if (!is.null(min.score)) {
    d <- subset(d, V12 >= min.score)
    if (verbose)
      cat(nrow(d),"with score >",min.score,"... ")
  }
  d <- merge(gene.dict2, d)
  setkey(d, V1)
  d <- merge(gene.dict1, d)
  d$V1 <- NULL
  d$V2 <- NULL
  setnames(d, 1:2, c("V1", "V2"))

  d$V13 <- d$V12 * (-1)
  setkey(d, V13)
  d <- d[!duplicated(d[, c("V1", "V2")])]
  d$V13 <- NULL
  if (verbose)
    cat(paste0("(", nrow(d), " unique) "))

  write.table(d,
              file = blast.file.out,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}


######################################################################
######################################################################
# read_allBlasts functions
######################################################################

#' @title read_speciesIDs
#' @description
#' \code{read_speciesIDs} read_speciesIDs
#' @rdname genespace_utils
#' @import data.table
#' @export
read_speciesIDs <- function(of.dir,
                            genomeIDs){
  si <- read.delim(file.path(of.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  si <- si[match(genomeIDs, si$genome),]
  return(data.table(si))
}

#' @title make_mapDB
#' @description
#' \code{make_mapDB} make_mapDB
#' @rdname genespace_utils
#' @import data.table
#' @export
make_mapDB <- function(id.db,
                       of.dir,
                       cull.blast.dir){
  sm <- with(id.db, cbind(expand.grid(n.old, n.old,
                                      stringsAsFactors = F),
                          expand.grid(genome, genome,
                                      stringsAsFactors = F)))
  sm2 <- with(id.db, cbind(expand.grid(n.new, n.new,
                                       stringsAsFactors = F),
                           expand.grid(genome, genome,
                                       stringsAsFactors = F)))
  colnames(sm) <- c("n.old1", "n.old2", "genome1", "genome2")
  colnames(sm2) <- c("n.new1", "n.new2", "genome1", "genome2")
  sm <- merge(sm, sm2, by = c("genome1", "genome2"))
  sm$filename <- file.path(of.dir,
                           paste0("Blast",
                                  sm$n.old1, "_", sm$n.old2,
                                  ".txt"))
  sm$new.filename <- file.path(cull.blast.dir,
                               paste0("Blast",
                                      sm$n.new1, "_", sm$n.new2,
                                      ".txt"))
  return(data.table(sm))
}

#' @title read_geneIDs
#' @description
#' \code{read_geneIDs} read_geneIDs
#' @rdname genespace_utils
#' @import data.table
#' @export
read_geneIDs <- function(of.dir,
                         gff,
                         species.num.id){
  gi <- fread(file.path(of.dir,
                        "SequenceIDs.txt"),
              sep = ":",
              stringsAsFactors = F,
              header = F,
              strip.white = T,
              col.names = c("gene.num", "id"))
  gi[,genome.num := as.integer(sapply(gene.num, function(x)
    strsplit(x,"_")[[1]][1]))]
  gi <- merge(species.num.id, gi, by = "genome.num")

  setkey(gi, genome, id)
  setkey(gff, genome, id)

  gi <- merge(gff, gi)
  return(gi)
}

#' @title import_ofResults
#' @description
#' \code{import_ofResults} import_ofResults
#' @rdname genespace_utils
#' @import data.table
#' @export
import_ofResults <- function(gff,
                             of.dir,
                             verbose = T,
                             genomeIDs){
  gff <- data.table(gff)
  #######################################################
  gz <- list.files(of.dir,
                   pattern = ".gz$")
  if (length(gz) > 0) {
    if (verbose)
      cat("\tDecompressing blast results\n")
    system(paste("gunzip -f",
                 file.path(of.dir,
                           "*.gz")))
  }
  #######################################################
  if (verbose)
    cat("\tReading Species IDs\n")
  si <- read.delim(file.path(of.dir,
                             "SpeciesIDs.txt"),
                   sep = ":",
                   stringsAsFactors = F,
                   header = F,
                   strip.white = T,
                   col.names = c("genome.num",
                                 "genome"))
  si$genome <- gsub(".fa$", "", si$genome)
  rownames(si) <- si$genome
  #######################################################
  sm <- cbind(expand.grid(si$genome.num,
                          si$genome.num,
                          stringsAsFactors = F),
              expand.grid(si$genome,
                          si$genome,
                          stringsAsFactors = F))
  colnames(sm) <- c("n1","n2","genome1","genome2")
  for (i in rev(genomeIDs))
    sm$ref[sm$genome1 == i | sm$genome2 == i] <- i
  sm$alt <- ifelse(sm$genome1 == sm$ref, sm$genome2, sm$genome1)
  sm$filename <- file.path(of.dir,
                           paste0("Blast",
                                  sm$n1, "_", sm$n2,
                                  ".txt"))
  sm$unique <- paste(sm$ref, sm$alt)
  sm$map.rank <- as.numeric(factor(sm$genome1,
                                   levels = genomeIDs))
  #######################################################
  if (verbose)
    cat("\tReading gene IDs\n")
  sequence.index <- fread(file.path(of.dir,
                                    "SequenceIDs.txt"),
                          sep = ":",
                          stringsAsFactors = F,
                          header = F,
                          strip.white = T,
                          col.names = c("gene.num", "id"))
  #######################################################
  if (verbose)
    cat("\tReading orthogroup networks\n")
  og <- readLines(file.path(of.dir,
                            "Orthogroups.txt"))
  og <- lapply(og, function(x) strsplit(x, " ")[[1]])
  ons <- sapply(og, function(x) x[1])
  names(og) <- ons
  og <- lapply(og, function(x) x[-1])

  if (verbose)
    cat("\tBuilding orthogroup data.table\n")
  og2 <- readLines(file.path(of.dir,
                             "Orthogroups.txt"))
  og2 <- lapply(og2, function(x) strsplit(x, " ")[[1]])
  og.name <- sapply(og2, function(x) x[1])
  og.length <- sapply(og2, length) - 1
  og.ids <- sapply(og2, function(x) x[-1])
  og2 <- data.table(block.id = NA,
                    og = rep(og.name, og.length),
                    id = unlist(og.ids),
                    stringsAsFactors = F)
  #######################################################
  if (verbose)
    cat("\tCompiling metadata\n")

  gffi <- gff[, c("id","genome")]
  setkey(gffi, "id")
  setkey(og2, "id")
  ogo <- merge(gffi, og2)
  setkey(ogo, "block.id", "genome", "og", "id")
  ogo[, og.n.genes := length(unique(id)), by = list(block.id, og)]
  ogo[, og.n.genomes := length(unique(genome)), by = list(block.id, og)]

  ogo.meta = ogo[!duplicated(ogo[, -c(1:2), with = F]),-c(1:2), with = F]
  #######################################################
  return(list(orthogroups = og,
              species.mappings = sm,
              species.index = si,
              orthogroup.metadata = ogo.meta,
              orthogroup.data.table = ogo,
              gene.index = sequence.index))
}

######################################################################
######################################################################
# remake_ofInput functions
######################################################################

#' @title cull_blastByScore
#' @description
#' \code{cull_blastByScore} cull_blastByScore
#' @rdname genespace_utils
#' @import data.table
#' @export
cull_blastByScore <- function(blast.file.in,
                              blast.file.out,
                              maxn,
                              verbose = T){
  d <- fread(blast.file.in, key = "V12")
  if (verbose)
    cat("Read in", nrow(d), "hits ... ")
  do <- d[,tail(.SD, maxn), by = list(V1)]
  if (verbose)
    cat(nrow(do),
        "top", maxn, "hits ... ")
  write.table(do,
              file = blast.file.out,
              sep = "\t",
              row.names = F,
              col.names = F,
              quote = F)
}

#' @title find_whichInBuffer
#' @description
#' \code{find_whichInBuffer} find_whichInBuffer
#' @rdname genespace_utils
#' @import data.table
#' @importFrom dbscan frNN
#' @export
find_whichInBuffer <- function(x,
                               y,
                               which.in.blk,
                               rank.buffer){

  nn <- frNN(x = data.frame(x, y),
             eps = rank.buffer)

  all.near.blk <- unique(unlist(nn$id[which.in.blk]))

  return(all.near.blk[order(all.near.blk)])
}

#' @title run_dbs
#' @description
#' \code{run_dbs} run_dbs
#' @rdname genespace_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @export
run_dbs <- function(map,
                    radius,
                    n.mappings){
  nn <- frNN(
    data.frame(map[, c("rank1", "rank2"), with = F]),
    eps = radius)
  dbs <- dbscan(nn,
                minPts = n.mappings)
  map[,cluster := dbs$cluster]
  return(map)
}

#' @title mirror_map
#' @description
#' \code{mirror_map} mirror_map
#' @rdname genespace_utils
#' @import data.table
#' @export
mirror_map <- function(map, keycol = NULL){
  cn <- colnames(map)
  cn1 <- cn[c(grep("1$",cn),grep("2$", cn),which(!grepl("1$|2$", cn)))]
  cn2 <- cn[c(grep("2$",cn),grep("1$", cn),which(!grepl("1$|2$", cn)))]
  out1 <- data.table(map[, cn1, with = F])
  out2 <- data.table(map[, cn2, with = F])
  setnames(out2, cn1)
  out3 <- rbind(out1, out2)
  if (!is.null(keycol))
    if (keycol %in% colnames(out3))
      setkeyv(out3, keycol)

  out3 <- out3[!duplicated(out3[, c("id1", "id2")]), ]
  return(out3)
}

#' @title clean_it
#' @description
#' \code{clean_it} clean_it
#' @rdname genespace_utils
#' @import data.table
#' @export
clean_it <- function(map,
                     genomeIDs,
                     rerank,
                     radius,
                     n.mappings,
                     verbose){

  map <- data.table(map)
  setkey(map, chr1, chr2, start1, start2)

  if (rerank) {
    map[,rank1 := frank(start1,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr1)]
    map[,rank2 := frank(start2,
                        ties.method = "dense"),
        by = list(genome1, genome2, chr2)]
  }

  if (!"unique.genome" %in% colnames(map)) {
    map$unique.genome <- with(map, paste(genome1, genome2))
  }
  if (!"unique" %in% colnames(map)) {
    map$unique <- with(map, paste(genome1, genome2, chr1, chr2))
  }

  if (is.null(genomeIDs))
    genomeIDs <- unique(gff$genome)

  map <- subset(map,
                genome1 %in% genomeIDs &
                  genome2 %in% genomeIDs)
  map[, genome1 := factor(genome1, levels = genomeIDs)]
  map[, genome2 := factor(genome2, levels = genomeIDs)]
  setkey(map, genome1, genome2)
  map[, genome1 := as.character(genome1)]
  map[, genome2 := as.character(genome2)]

  spl.gen <- split(map, by = c("genome1", "genome2"))
  merged_map <- rbindlist(lapply(spl.gen, function(x){
    g1 = x$genome1[1]
    g2 = x$genome2[1]

    if (verbose)
      cat(paste0("\t", g1), "-->", g2,
          paste0("(initial hits = ", nrow(x) , ") ... "))
    spl.map <- split(x, by = c("chr1", "chr2"))

    chr.map <- rbindlist(lapply(spl.map, function(tmp){
      x <- run_dbs(
        map = tmp[, c("rank1", "rank2"), with = F],
        radius = radius,
        n.mappings = n.mappings)

      tmp[, block.id := x$cluster]
      return(tmp)
    }))

    chr.map <- chr.map[chr.map$block.id != 0, ]
    if (verbose)
      cat(nrow(chr.map), "hits in",
          length(unique(paste(chr.map$unique,
                              chr.map$block.id))),
          "blocks\n")
    return(chr.map)
  }))

  merged_map[,block.id := as.numeric(
    as.factor(
      paste(genome1, genome2,
            chr1, chr2, block.id)))]

  merged_blk <- make_blocks(
    map = merged_map,
    rename.blocks = T,
    rerank = rerank,
    clean.columns = F,
    ties.method = "dense")
  return(merged_blk)
}
