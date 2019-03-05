#' @title Run the exonerate program
#'
#' @description
#' \code{run_exonerate} A simple wrapper to run orthofinder from R.
#'
#' @param ass.fasta.file The path to the assembly fasta to search
#' @param pep.fasta.file The path to the directory where the blast results should be stored
#' @param bed The path to the directory where temporary files will be stored
#' then deleted
#' @param min.score The number of threads for blast run
#' @param tmp.dir The number of threads used for orthogroup construction
#' @param ... Not currently in use
#' @details  ...

#' @return Nothing, writes results to the blast.dir directory
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom compiler cmpfun
#' @export
run_exonerate <- function(locs,
                          verbose = T,
                          clean = T,
                          min.score = 20){

  ########################################################
  ########################################################
  call_exonerate <- function(ass.fasta.file,
                             pep.fasta.file,
                             min.score = 20){
    system(paste("exonerate",
                 "--model protein2genome",
                 "--score", min.score,
                 "--softmasktarget --alignmentwidth 1000000",
                 "--bestn 1",
                 "--dpmemory 2000",
                 "--maxintron 10000",
                 "--showalignment no ",
                 "--showtargetgff yes",
                 "--showvulgar no",
                 "--showcigar no",
                 '--ryo ">\n%tas\n"',
                 "--query", pep.fasta.file,
                 "--target", ass.fasta.file),
           intern = T)
  }
  ########################################################
  ########################################################
  proc_exonerate <- function(exonerate_output,
                             bed){
    exo <- exonerate_output

    if (length(grep("^>", exo)) == 0) {
      exo <- DNAStringSet("")
      names(exo) <- bed$id
      gff1 <- data.frame(exon.start = NA,
                         exon.end = NA,
                         strand = NA,
                         align.info = NA)
      gffo <- data.frame(bed,
                         gff1)
    }else{
      end <- which(exo == "")[1] + 1
      exo <- exo[1:end]
      gff1 <- do.call(rbind,
                      lapply(exo[grepl(bed$id, exo) &
                                   grepl("exon\t", exo)], function(x)
                                     strsplit(x, "\t")[[1]][c(4:5, 7, 9)]))
      gff1 <- data.frame(gff1,
                         stringsAsFactors = F)
      colnames(gff1) <- c("exon.start",
                          "exon.end",
                          "strand",
                          "align.info")
      gffo <- data.frame(bed, gff1)

      wh <- grep("^>", exo) + 1
      exo <- DNAStringSet(paste(exo[wh:(length(exo) - 2)],
                                collapse = ""))
      names(exo) <- bed$id
    }
    return(list(cds.seq = exo,
                gff = gffo))
  }
  #######################################################
  #######################################################
  proc_exonerate <- cmpfun(proc_exonerate)
  call_exonerate <- cmpfun(call_exonerate)
  #######################################################
  #######################################################

  if (verbose) {
    n <- nrow(locs)
    nb <- ifelse(n < 5000, 100,
                 ifelse(n < 10000, 500,1000))
  }

  exon.out <- lapply(1:nrow(locs), function(i){

    if (verbose)
      if (i %% nb == 0)
        cat(i,"/",nrow(locs),"\n\t")
    x <- locs[i]
    bedt <- data.frame(x[,c("chr","start","end","reg.name")])
    colnames(bedt)[4] <- "id"

    exonerate.out <- call_exonerate(ass.fasta.file = x$fa.file,
                                    pep.fasta.file = x$pep.file,
                                    min.score = min.score)
    exon <- proc_exonerate(exonerate_output = exonerate.out,
                           bed = bedt)
    return(exon)
  })
  exon.cds <- do.call(c, lapply(exon.out, function(x) x$cds.seq))
  gff.gene <- rbindlist(lapply(exon.out, function(x) x$gff))
  gff.gene$ex.start <- with(gff.gene,
                            as.numeric(start) + as.numeric(exon.start))
  gff.gene$ex.end <- with(gff.gene,
                          as.numeric(start) + as.numeric(exon.end))
  gff.out <- gff.gene[,list(start = min(ex.start),
                            end = max(ex.end)),
                      by = list(id, chr, strand)]

  if(clean)
    unlink(c(locs$fa.file, locs$pep.file))

  return(list(gff.region = gff.out,
              gff.cds = gff.gene,
              cds = exon.cds))
}
