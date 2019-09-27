#' @title Read all blast files
#'
#' @description
#' \code{read_allBlasts} Read all blast files in an orthofinder output directory
#'
#' @param gff data.table, containing the parsed gff annotations,
#' as produced by import_gff.
#' @param genomeIDs character vector, specifying the genomeIDs to include.
#' @param of.dir file.path, to the subdirectory containing
#' blast results
#' @param check.ogs logical, should the orthogroups be considered at all?
#' @param add.gff logical, should annotation data be added?
#' @param keep.geneNum logical, should the orthofinder gene ID be retained?
#' @param verbose logical, should updates be printed to the console?
#' @param ... Not currently in use
#'
#' @details Small and dispersed blocks are dropped using 2-dimensional
#' clustering. Essentially, any hits that are not near n.mappings hits
#' within a specified radius, are dropped. The remaining hits are clustered
#' following standard DBScan methods.
#'
#' @return A list of length 2, block and map, as output by make_blocks.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
read_allBlasts <- function(gff,
                           genomeIDs,
                           of.dir,
                           check.ogs = T,
                           add.gff = F,
                           keep.geneNum = F,
                           simple = F,
                           blast.files = NULL,
                           verbose = T){

  if(simple){
    if(is.null(blast.files))
      stop("must specify blast files if simple\n")
    all.blast <- rbindlist(lapply(blast.files, fread,
                                  col.names = c("gn1", "gn2", "perc.iden", "align.length",
                                                "n.mismatch", "n.gapOpen", "q.start",
                                                "q.end", "s.start",
                                                "s.end", "eval", "score"),
                                  key = "gn2"))
    sp.ids <- read_speciesIDs(of.dir = of.dir,
                               genomeIDs = genomeIDs)
    of.geneIndex <- read_geneIDs(of.dir = of.dir,
                                 species.num.id = sp.ids,
                                 gff = gff)
    setnames(of.geneIndex, "gene.num","gn")
    of1 <- data.table(of.geneIndex)
    setnames(of1, paste0(colnames(of1), "1"))
    of2 <- data.table(of.geneIndex)
    setnames(of2, paste0(colnames(of2), "2"))
    out <- merge(of2, all.blast, by = "gn2")
    out <- merge(of1, out, by = "gn1")
    return(out)
  }else{
    sp.ids <- read_speciesIDs(of.dir = of.dir,
                              genomeIDs = genomeIDs)
    if (check.ogs) {
      if (verbose)
        cat("Importing orthofinder database ... ")
      ofdat <- import_ofResults(
        gff = gff,
        genomeIDs = genomeIDs,
        of.dir = of.dir,
        verbose = F)
      if (verbose)
        cat("Done!\n")
      of.geneIndex <- ofdat$gene.index
      of.blastFiles <- ofdat$species.mappings
      of.blastFiles <- subset(of.blastFiles, genome1 %in% genomeIDs & genome2 %in% genomeIDs)
    }else{
      of.blastFiles <- list(filename = list.files(of.dir,
                                                  pattern = "^Blast",
                                                  full.names = T))
    }

    if (verbose)
      cat("Reading all blast files into memory ... ")
    all.blast <- rbindlist(lapply(of.blastFiles$filename, fread,
                                  col.names = c("gn1", "gn2", "perc.iden", "align.length",
                                                "n.mismatch", "n.gapOpen", "q.start",
                                                "q.end", "s.start",
                                                "s.end", "eval", "score"),
                                  key = "gn2"))
    of.geneIndex <- read_geneIDs(of.dir = of.dir,
                                 species.num.id = sp.ids,
                                 gff = gff)
    setnames(of.geneIndex, "gene.num","gn")
    if (verbose)
      cat("Done!\nParsing results and merging with geneIDs ... ")
    if (!add.gff) {
      of.geneIndex1 <- data.table(of.geneIndex)
      of.geneIndex2 <- data.table(of.geneIndex)

      setnames(of.geneIndex1, paste0(colnames(of.geneIndex1),"1"))
      setnames(of.geneIndex2, paste0(colnames(of.geneIndex2),"2"))

      m <- merge(of.geneIndex2, all.blast, by = "gn2")
      m <- merge(of.geneIndex1, m, by = "gn1")

      if (!keep.geneNum) {
        m$gn1 <- NULL
        m$gn2 <- NULL
      }
    }else{
      if (add.gff) {
        gff1 <- data.table(of.geneIndex)
        gff2 <- data.table(of.geneIndex)
        setnames(gff1, paste0(colnames(gff1), "1"))
        setnames(gff2, paste0(colnames(gff2), "2"))
        setkey(gff2, gn2)
        setkey(all.blast, gn2)
        m <- merge(gff2, all.blast)
        setkey(gff1, gn1)
        setkey(m, gn1)
        m <- merge(gff1, m)
        if (!keep.geneNum) {
          m$gn1 <- NULL
          m$gn2 <- NULL
        }
      }
    }

    if (verbose)
      cat("Done!\n")
    return(m)
  }
}
