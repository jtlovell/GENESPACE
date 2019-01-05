#' @title Import blast results
#' @description
#'  \code{import_ofBlast} Import blast results from orthofinder
#' @param species.mappings data.frame of species mappings and file locations
#' produced by import_ofResults
#' @param orthogroups data.table of orthogroups
#' produced by import_ofResults
#' @param genomeIDs character, genome identifiers
#' @param gff data.table containing the parsed gff annotation data
#' @param gene.index gene index from import_ofResults
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details Nothing yet
#' @return new list of block and map data.tables
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
import_ofBlast <- function(species.mappings,
                           genomeIDs,
                           orthogroups,
                           gene.index,
                           gff,
                           verbose = T,
                           ...){
  #######################################################
  #######################################################
  read_ofBlast <- function(genomeID1,
                           genomeID2,
                           species.mappings){
    sm <- species.mappings[with(species.mappings,
                                (genome1 == genomeID1 &
                                   genome2 == genomeID2) |
                                  (genome2 == genomeID1 &
                                     genome1 == genomeID2)),]
    sm <- sm[order(sm$n1), ]
    filenames<- sm$filename
    filenames <- filenames[]
    if(length(filenames) == 1){
      suppressWarnings(
        blast.in <- fread(filenames,
                          showProgress = F))
    }else{
      suppressWarnings(
        b1 <- fread(filenames[1],
                    showProgress = F))

      suppressWarnings(
        b2 <- fread(filenames[2],
                    showProgress = F))
      b2 <- data.table(b2[, c(2, 1, 3:6, 9:10, 7:8, 11:12)])
      setnames(b2, colnames(b1))

      blast.in <- rbind(b1, b2)
    }
    setnames(blast.in, c("gn1", "gn2", "perc.iden",
                         "align.length", "n.mismatch",
                         "n.gapOpen", "q.start", "q.end",
                         "s.start", "s.end",
                         "eval", "score"))
    return(blast.in)
  }
  #######################################################
  #######################################################
  parse_ofBlast <- function(blast,
                            parsed.ogs,
                            gff1,
                            gff2,
                            verbose = T){
    # - Drop duplicates
    blast <- data.table(blast)
    g1 <- parsed.ogs$g1
    g2 <- parsed.ogs$g2

    blast$neg.score <- blast$score * (-1)
    setkey(blast, "gn1", "gn2", "neg.score")

    blast <- blast[!duplicated(blast[, c("gn1","gn2"), with = F])]

    blast$neg.score <- NULL

    if (verbose)
      cat(paste0("(unique hits = ", nrow(blast),")"))

    # - Merge with gene IDs
    setkey(blast, "gn2")
    bl1 <- merge(g2, blast)
    setkey(bl1, "gn1")
    bl2 <- merge(g1, bl1)

    blast <- data.table(bl2[with(bl2, og1 == og2),])

    if (verbose)
      cat("\t",nrow(blast), "hits in orthogroups ... ")

    # - Merge with gffs
    gf1 <- data.table(gff1)
    gf2 <- data.table(gff2)
    setkey(blast, "id2")
    setkey(gf2, "id2")
    setkey(gf1, "id1")
    blast <- merge(gf2, blast)
    setkey(blast, "id1")
    blast <- merge(gf1, blast)

    if (verbose)
      cat("Done!\n")
    return(blast)
  }
  #######################################################
  #######################################################
  parse_orthogroups <- function(orthogroups,
                                gff,
                                gene.index){
    wh2 <- which(sapply(orthogroups, length) > 1)
    g2 <- orthogroups[wh2]
    l2 <- sapply(g2, length)
    n2 <- rep(names(l2), l2)

    g1 <- data.table(og = n2,
                     id = unlist(g2),
                     stringsAsFactors = F)
    setkey(g1, "id")

    setkey(gff, "id")
    g1a = merge(gff, g1)
    g1a[, n.genomes := length(unique(genome)),
        by = list(og)]
    g1a <- g1a[g1a$n.genomes > 1,]
    g1 <- g1[g1$id %in% g1a$id, ]

    setkey(gene.index, "id")
    g1 <- merge(g1, gene.index)
    g2 <- data.table(g1)

    setnames(g1, c("id1","og1", "gn1"))
    setnames(g2, c("id2", "og2", "gn2"))
    setkey(g1, "gn1")
    setkey(g2, "gn2")
    return(list(g1 = g1,
                g2 = g2))
  }
  #######################################################
  #######################################################

  gene.index <- data.table(gene.index)
  #######################################################

  #######################################################
  if (verbose)
    cat("Parsing gff annotations ... ")
  gff1 <- data.table(gff)
  gff2 <- data.table(gff)
  setnames(gff1,c("id1", "chr1", "start1", "end1",
                  "strand1", "genome1", "order1"))
  setnames(gff2,c("id2", "chr2", "start2", "end2",
                  "strand2", "genome2", "order2"))
  setkey(gff1, "id1")
  setkey(gff2, "id2")

  spl.gff1 <- split(gff1, "genome1")
  spl.gff2 <- split(gff2, "genome2")

  if (verbose)
    cat("Done\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Parsing orthogroups\n")
  gs <- parse_orthogroups(orthogroups = orthogroups,
                          gff = gff,
                          gene.index = gene.index)

  if (verbose)
    cat("Done\n")
  #######################################################

  #######################################################
  if (verbose)
    cat("Importing BLAST results ... \n")

  comb <- combn(genomeIDs, 2, simplify = F)
  for(i in genomeIDs){
    comb[[length(comb)+1]]<-rep(i,2)
  }

  blast <- rbindlist(lapply(comb, function(x){
    if(verbose)
      cat(paste0("\t", x[1]),"-->", x[2],"")
    blast <- read_ofBlast(genomeID1 = x[1],
                          genomeID2 = x[2],
                          species.mappings = species.mappings)
    blast <- parse_ofBlast(blast = blast,
                           parsed.ogs = gs,
                           gff1 = spl.gff1[[x[1]]],
                           gff2 = spl.gff2[[x[2]]])
    return(blast)
  }))
  if(verbose)
    cat("\tDone!\n")
  return(blast)
}
