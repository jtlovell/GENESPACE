#' @title Run the exonerate program
#'
#' @description
#' \code{run_getfasta} A simple wrapper to run orthofinder from R.
#'
#' @param bed.dt data.table or data.frame containing bed-formatted coordinates
#' to pull. Column names must include "chr", "start", "end", and "reg.name".
#' @param tmp.dir directory path to the tempory directory
#' @param assembly.dir directory path to the assembly fasta files
#' @param n.cores numeric, the number of parallel processes to run.
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#' @details  ...

#' @return A data.table with the same format as bed.dt, but with an
#' additional column with the file path to the fasta file.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @importFrom compiler cmpfun
#' @export
run_getfasta <- function(bed.dt,
                         write.dir,
                         tmp.dir,
                         assembly.dir,
                         n.cores = 1,
                         verbose = T){
  setnames(bed.dt,
           c("chr","start","end","genome","name"))

  bed.dt$start <- bed.dt$start - 1
  bed.dt$end <- bed.dt$end - 1

  ########################################################
  ########################################################
  get_fasta <- function(bed,
                        bed.file,
                        ass.file,
                        fa.file){

    write.table(bed,
                file = bed.file,
                quote = F,
                sep = "\t",
                row.names = F,
                col.names = F)

    system(paste("bedtools getfasta",
                 "-fi", ass.file,
                 "-bed", bed.file,
                 "-name+ -fo", fa.file))
  }
  ########################################################
  ########################################################
  get_fasta <- cmpfun(get_fasta)
  ########################################################
  ########################################################
  fa.locs <- unlist(lapply(1:nrow(bed.dt), function(i){

    bed <- bed.dt[i,c("chr","start","end","genome")]
    bed.id <- bed.dt$name[i]

    ass.file <- file.path(assembly.dir,
                          paste0(bed$genome, ".fa"))
    fa.file <- file.path(write.dir,
                         paste0(bed.id, ".fa"))
    bed.file <- file.path(tmp.dir,
                          paste0(bed.id, ".bed"))

    get_fasta(bed = bed,
              ass.file = ass.file,
              fa.file = fa.file,
              bed.file = bed.file)

    return(fa.file)
  }))

  bed.dt$fa.file <- fa.locs
  return(bed.dt)

}
