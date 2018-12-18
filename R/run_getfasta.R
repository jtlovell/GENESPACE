#' @title Run the exonerate program
#'
#' @description
#' \code{run_getfasta} A simple wrapper to run orthofinder from R.
#'
#' @param ass.file The path to the assembly fasta to search
#' @param bed.file The path to the directory where the blast results should be stored
#' @param bed The path to the directory where temporary files will be stored
#' then deleted
#' @param fa.file The number of threads for blast run
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
run_getfasta <- function(bed,
                        bed.file,
                        ass.file,
                        fa.file){
  write.table(bed,
              file=bed.file,
              quote=F,
              sep="\t", row.names=F, col.names=F)
  system(paste("bedtools getfasta -fi",
               ass.file, "-bed", bed.file,"-name -fo",fa.file))
}
