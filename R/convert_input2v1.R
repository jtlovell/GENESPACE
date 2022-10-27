#' @title Build genespace pangenome
#'
#' @description
#' \code{convert_input2v1} Convert orthogroup and synteny information into a
#' pangenome database. Predict locations of orthogroups that are missing a
#' node in the reference.
#'
#' @param existingDir A
#' @param v1Dir A
#' @details T...
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
    stop(sprtinf("v1Dir %s exists and is not empty. Cannot overwrite\n",
                 v1Dir))
  if(!dir.exists(dirname(v1Dir)))
    stop(sprtinf("parent directory of v1Dir %s does not exist.\n",
                 v1Dir))
  if(!dir.exists(existingDir))
    stop(sprtinf("existingDir %s does not exist.\n",
                 existingDir))

  pepDir <- file.path(existingDir, "peptide")
  if(!dir.exists(pepDir))
    stop(sprtinf("peptide directory %s does not exist.\n",
                 pepDir))

  gffDir <- file.path(existingDir, "gff")
  if(!dir.exists(gffDir))
    stop(sprtinf("gff directory %s does not exist.\n",
                 gffDir))

  ofDir <- file.path(existingDir, "orthofinder")
  if(!dir.exists(ofDir))
    stop(sprtinf("orthofinder directory %s does not exist.\n",
                 ofDir))

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
    x[,start := start - 1]
    bedout <- file.path(v1BedDir, sprintf("%s.bed", genomeID))
    fwrite(x, file = bedout, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  nu <- file.copy(list.files(pepDir, full.names = T), v1PepDir)

  gids <- gsub(".fa$", "", list.files(v1PepDir))
  ofFiles <- copy_of2results(
    orthofinderDir = ofDir,
    resultsDir = v1ResDir,
    genomeIDs = gids)
}
