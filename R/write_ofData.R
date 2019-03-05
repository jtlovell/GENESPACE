#' @title Prep data for orthofinder
#'
#' @description
#' \code{write_ofData} Write blast, move peptide fastas, make diamond databases,
#' write species index and sequence index files.
#'
#' @param blast A data.table with at least the 10 necessary blast columns
#' @param peptide.dir The directory containing the peptide fasta annotations
#' @param of.dir The directory to write results
#' @param genomeIDs Character vector giving the genome IDs.
#' @param verbose Logical, should updates be printed?
#' @param ... Not currently in use
#' @details None yet

#' @return Nothing, just writes to file
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom compiler cmpfun
#' @importFrom utils write.table
#' @export
  write_ofData <- function(blast,
                           genomeIDs,
                           of.dir,
                           peptide.dir,
                           verbose = T,
                           ...){

    ########################################################
    ########################################################
    extract_seqIDs <- function(blast, genomeIDs){
      ug <- with(blast,
                 data.table(genome = c(genome1, genome2),
                            chr = c(chr1,chr2),
                            pos = c(start1, start2),
                            id = c(id1, id2),
                            stringsAsFactors = F))
      ug$genome <- factor(ug$genome, levels = genomeIDs)
      setkey(ug, genome, chr, pos)
      ug <- ug[!duplicated(ug$id),]
      ug <- ug[ug$genome %in% genomeIDs,]
      ug <- ug[with(ug, order(genome, chr, pos)),]
      ugenes <- data.table(rbindlist(lapply(split(ug, "genome"), function(x){
        x$num <- 0:(nrow(x) - 1)
        return(x)
      })))

      ugenes$genome.num <- as.numeric(factor(ugenes$genome, levels = genomeIDs)) - 1
      ugenome <- ugenes[!duplicated(ugenes[,c("genome","genome.num")]),c("genome","genome.num")]

      seqs <- with(ugenes, paste0(genome.num, "_", num, ": ", id))
      spes <- with(ugenome, paste0(genome.num, ": ", genome, ".fa"))

      dmnd.ids <- with(ugenome, paste0("diamondDBSpecies", genome.num))
      names(dmnd.ids) <- ugenome$genome

      fa.ids <- with(ugenome, paste0("Species", genome.num, ".fa"))
      names(fa.ids) <- ugenome$genome

      ugenome <- data.frame(ugenome)
      rownames(ugenome) <- ugenome$genome
      eg <- as.matrix(expand.grid(genomeIDs, genomeIDs))
      eg <- data.table(genome1 = eg[,1],
                       genome2 = eg[,2],
                       genome.num1 = ugenome[eg[,1], "genome.num"],
                       genome.num2 = ugenome[eg[,2], "genome.num"])
      return(list(species.ids = spes,
                  sequence.ids = seqs,
                  dmnd.ids = dmnd.ids,
                  fa.ids = fa.ids,
                  blast.ids = data.frame(eg),
                  gene.dict = with(ugenes,
                                   data.table(id = id,
                                              num = paste0(genome.num, "_", num))),
                  id.list = with(ugenes, split(id, genome)),
                  num.list = with(ugenes, split(paste0(genome.num, "_", num), genome))))
    }
    ########################################################
    ########################################################
    write_ofBlast <- function(blast,
                              blast.ids,
                              of.dir){
      bl <- blast[,c("id1", "id2",
                     "perc.iden", "align.length",
                     "n.mismatch","n.gapOpen",
                     "q.start", "q.end",
                     "s.start", "s.end",
                     "eval", "score",
                     "genome1","genome2")]
      bl1 <- data.table(bl)
      bl2 <- data.table(bl[,c(2,1,3:6,8,7,10,9,11:12,14,13)])
      setnames(bl2, colnames(bl1))
      bl <- rbind(bl1, bl2)
      bl$negscore <- bl$score * (-1)
      setkey(bl, negscore)
      bl <- bl[!duplicated(bl[,c("id1", "id2")]),]
      bl$negscore <- NULL

      d1 <- with(of.ids$gene.dict, data.table(id1 = id, num1 = num))
      d2 <- with(of.ids$gene.dict, data.table(id2 = id, num2 = num))
      setkey(bl, id2)
      setkey(d2, id2)
      setkey(d1, id1)
      blt <-  merge(d2, bl)
      setkey(blt, id1)
      bl <- merge(d1, blt)[,c(2, 4, 5:16)]

      eg <- data.frame(of.ids$blast.ids)
      for (i in 1:nrow(eg)) {
        x1 = eg[i, 1]
        x2 = eg[i, 2]
        n1 = eg[i, 3]
        n2 = eg[i, 4]
        tmp.bl <- bl[with(bl, genome1 == x1 & genome2 == x2), -c(13:14)]
        write.table(tmp.bl,
                    file = file.path(of.dir,
                                     paste0("Blast", n1, "_", n2, ".txt")),
                    quote = F, col.names = F,
                    row.names = F, sep = "\t")
      }
    }

    #######################################################
    #######################################################
    write_ofBlast <- cmpfun(write_ofBlast)
    extract_seqIDs <- cmpfun(extract_seqIDs)
    #######################################################
    #######################################################

    ########################################################
    if (verbose)
      cat("Extracting sequence and species IDs ... ")
    of.ids <- extract_seqIDs(blast = blast, genomeIDs = genomeIDs)
    cat(of.ids$species.ids,
        sep = "\n",
        file = file.path(of.dir, "SpeciesIDs.txt"))
    cat(of.ids$sequence.ids,
        sep = "\n",
        file = file.path(of.dir, "SequenceIDs.txt"))
    if (verbose)
      cat("Done!\n")
    ########################################################

    ########################################################
    if (verbose)
      cat("Parsing and databasing peptide sequences ... ")
    pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
      readAAStringSet(file.path(peptide.dir,
                                paste0(i, ".fa")))))
    pep.spl <- sapply(names(of.ids$id.list), USE.NAMES = T, simplify = F, function(i){
      x <- pep.fastas[of.ids$id.list[[i]]]
      names(x) <- of.ids$num.list[[i]]
      return(x)
    })

    pep.files <- sapply(names(pep.spl), function(i){
      outf <- file.path(of.dir, of.ids$fa.ids[i])
      outdb <- file.path(of.dir, of.ids$dmnd.ids[i])
      writeXStringSet(pep.spl[[i]], filepath = outf)
      system(paste("diamond makedb --quiet",
                   "--in", outf,
                   "-d", outdb))
      return(list(fa = outf,
                  db = outdb))
    })
    if (verbose)
      cat("Done!\n")
    ########################################################

    ########################################################
    if (verbose)
      cat("Converting blast results to orthofinder-formatted text files ... ")
    write_ofBlast(blast = blast,
                  blast.ids = of.ids$blast.ids,
                  of.dir = of.dir)
    if (verbose)
      cat("Done!\n")
    ########################################################
  }
