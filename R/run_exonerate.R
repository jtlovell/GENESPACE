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
#' @export
run_exonerate <- function(ass.fasta.file,
                          pep.fasta.file,
                          bed,
                          min.score = 20,
                          tmp.dir){


  out1 = system(paste("exonerate --model protein2genome",
                      "--score",min.score,
                      "--softmasktarget --alignmentwidth 1000000",
                      "--bestn 1",
                      "--dpmemory 2000 --maxintron 10000",
                      "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                      '--ryo ">\n%tas\n"',
                      "--query",pep.fasta.file,
                      "--target", ass.fasta.file),
                intern = T)

  if(length(grep("^>", out1))==0){
    out1 = DNAStringSet("")
    names(out1) <- bed$id
    gff1 = data.frame(exon.start = NA,
                      exon.end = NA,
                      strand = NA,
                      align.info = NA)
    gffo = data.frame(bed,gff1)
  }else{
    end = which(out1 == "")[1]+1
    out1 = out1[1:end]
    gff1 = do.call(rbind,lapply(out1[grepl(bed$id, out1) & grepl("exon\t", out1)],
                                function(x) strsplit(x,"\t")[[1]][c(4:5,7,9)]))
    gff1 = data.frame(gff1,
                      stringsAsFactors = F)
    colnames(gff1)<-c("exon.start","exon.end","strand","align.info")
    gffo = data.frame(bed,gff1)

    wh = grep("^>",out1)+1
    out1 = DNAStringSet(paste(out1[wh:(length(out1)-2)], collapse = ""))
    names(out1)<-bed$id
  }
  return(list(cds.seq = out1, gff = gffo))
}
