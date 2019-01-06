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
#' @export
run_getfasta <- function(bed.dt,
                         tmp.dir,
                         assembly.dir,
                         n.cores = 1,
                         verbose = T){

  bed.dt <- data.table(bed.dt)

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
                 "-name -fo", fa.file))
  }
  ########################################################
  ########################################################

  if (dir.exists(tmp.dir))
    unlink(tmp.dir,
           recursive = T)

  dir.create(tmp.dir)

  nb <- ifelse(nrow(bed.dt) < 1000, 100,
               ifelse(nrow(bed.dt) < 5000, 500,
                      ifelse(nrow(bed.dt) < 10000, 1000, 5000)))

  fa.locs <- unlist(mclapply(1:nrow(bed.dt), mc.cores = n.cores, function(i){
    if (verbose)
      if (i %% nb == 0)
        cat(i, "/",
            nrow(bed.dt), "\n\t")

    x <- bed.dt[i,]
    bed <- x[, c("chr", "start", "end", "reg.name")]

    ass.file <- file.path(assembly.dir,
                          paste0(x$genome, ".fa"))
    fa.file <- file.path(tmp.dir,
                         paste0(x$reg.name, ".fa"))
    bed.file <- file.path(tmp.dir,
                          paste0(x$reg.name, ".bed"))

    get_fasta(bed = bed,
                 ass.file = ass.file,
                 fa.file = fa.file,
                 bed.file = bed.file)

    return(fa.file)
  }))

  bed.dt$fa.file <- fa.locs
  return(bed.dt)

}
