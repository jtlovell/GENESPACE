#' @title Convert input data formats
#'
#' @description
#' \code{convert_input2v1} Takes an existing GENESPACE <= v0.9.4 run, converts
#' file formats and copies to a new directory which can be used for GENESPACE >=
#' v1.0.0.
#'
#' @param existingDir character string coercible to a file.path, pointing to an
#' existing GENESPACE <= v0.9.4 run.
#' @param v1Dir character string coercible to a file.path, pointing to a
#' directory where you want to run GENESPACE >=v1.0.0
#' @details Directly copies the peptide fasta files and required orthofinder
#' files. Converts the 1-indexed gff to 0-indexed bed file and re-names them.
#'
#' @return nothing
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#
#' @import data.table
#' @export
convert_input2v1 <- function(existingDir, v1Dir){
  if(dir.exists(v1Dir) & length(list.files(v1Dir)) > 0)
    stop(sprintf("v1Dir %s exists and is not empty. Cannot overwrite\n", v1Dir))
  if(!dir.exists(dirname(v1Dir)))
    stop(sprintf("parent directory of v1Dir %s does not exist.\n", v1Dir))
  if(!dir.exists(existingDir))
    stop(sprintf("existingDir %s does not exist.\n", existingDir))

  pepDir <- file.path(existingDir, "peptide")
  if(!dir.exists(pepDir))
    stop(sprintf("peptide directory %s does not exist.\n", pepDir))

  gffDir <- file.path(existingDir, "gff")
  if(!dir.exists(gffDir))
    stop(sprintf("gff directory %s does not exist.\n", gffDir))

  ofDir <- file.path(existingDir, "orthofinder")
  if(!dir.exists(ofDir))
    stop(sprintf("orthofinder directory %s does not exist.\n", ofDir))

  if(!dir.exists(v1Dir))
    dir.create(v1Dir)

  v1BedDir <- file.path(v1Dir, "bed")
  if(!dir.exists(v1BedDir))
    dir.create(v1BedDir)

  v1PepDir <- file.path(v1Dir, "peptide")
  if(!dir.exists(v1PepDir))
    dir.create(v1PepDir)

  v1ResDir <- file.path(v1Dir, "results")
  if(!dir.exists(v1ResDir))
    dir.create(v1ResDir)

  for(i in list.files(gffDir, full.names = T)){
    x <- fread(
      i, na.strings = c("", "NA"), showProgress = F,
      colClasses = c("character", "numeric", "numeric",
                     "character", "character","numeric"))
    genomeID <- gsub(".gff.gz$",  "", basename(i))
    x <- x[,c("chr", "start", "end", "id")]
    start <- NULL
    x[,start := start - 1]
    bedout <- file.path(v1BedDir, sprintf("%s.bed", genomeID))
    fwrite(x, file = bedout, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  nu <- file.copy(list.files(pepDir, full.names = T), v1PepDir)

  gids <- gsub(".fa$", "", list.files(v1PepDir))
  ofFiles <- copy_of2results(
    orthofinderDir = ofDir,
    resultsDir = v1ResDir)
}
