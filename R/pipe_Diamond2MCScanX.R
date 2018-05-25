#' @title Run the first two steps of gene-space alignments
#'
#' @description
#' \code{pipe_Diamond2MCScanX} For all mapping files in a directory, parse with genome data
#'
#' @param inputFileMatrix A dataframe with all the necessary metadata. See details.
#' @param nmapsPerHaplotype see ploidy
#' @param dbs_radii numeric vector of length equal to `dbs_mappingsInRadius`.
#' @param dbs_mappingsInRadius numeric vector of length equal to `dbs_radius`.
#' @param mcscan.dir Path to output directory for mcscan input files
#' @param plotit Logical, should a plot be drawn?
#' @param path_to_diamond The location of the Diamond program executable
#' @param sensitive.mode Diamond sensitivity mode.
#' @param topPerc Genes this percent from the maximum hit are retained
#' @param minScore The minimum mapping score to be considered
#' @param ortherDiamondOpts Other diamond options to pass the command
#' @param nthreads Number of parallel threads to run blast with
#' @param onlyParseMap Should the mapping step be ignored (e.g. has it already
#' beed done and the mapping files are in the correct directory?)
#' @param verbose Logical, should status updates be printed?
#' @param ... Not currently in use
#' @details This is the main engine for BLAST mapping and mapping parsing. The
#' parameters in `inputFileMatrix` is used to specify the genomes and genome files
#' that are fed into `align_peptideByDiamond` and subsequently `parse_diamondBlast`.
#' The `inputFileMatrix` must contain the following column names (exactly):
#' id1,id2,ploidy1,ploidy2,abbrev1,abbrev2,pep1,pep2,gff1,gff2,blast1,blast2. These are
#' defined in the help files for `align_peptideByDiamond` and  `parse_diamondBlast`.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @import dbscan
#' @export
pipe_Diamond2MCScanX = function(inputFileMatrix,
                                path_to_diamond,
                                nmapsPerHaplotype = 1,
                                plotit = T,
                                dbs_radii = c(100,50,15),
                                dbs_mappingsInRadius = c(10,10,5),
                                mcscan.dir,
                                verbose = T,
                                sensitive.mode = "--more-sensitive",
                                topPerc = 50,
                                minScore = 200,
                                ortherDiamondOpts = "--quiet",
                                nthreads= 1,
                                onlyParseMap = F){

  if(verbose){
    cat("Running the pipeline for:\n", apply(inputFileMatrix[,1:2], 1,
                                          function(x) paste(paste(x, collapse = " vs. "),"\n")),
        "\n\n####################\n")
  }
  system(paste("rm -r", mcscan.dir))
  system(paste("mkdir",mcscan.dir))
  out = lapply(1:nrow(inputFileMatrix), function(i){
    x = inputFileMatrix[i,]

    if(!onlyParseMap){
      diamond.calls = align_peptideByDiamond(
        path_to_diamond = path_to_diamond,
        id1 = x$id1,
        id2 = x$id2,
        pep1 = x$pep1,
        pep2 = x$pep2,
        blast1 = x$blast1,
        blast2 = x$blast2,
        nthreads = nthreads,
        sensitive.mode = "--more-sensitive",
        topPerc = topPerc,
        minScore = minScore,
        ortherDiamondOpts = "--quiet",
        verbose = T)
    }

    allmap = parse_diamondBlast(
      id1 = x$id1,
      id2 = x$id2,
      abbrev1 = x$abbrev1,
      abbrev2 = x$abbrev2,
      ploidy1 = x$ploidy2,
      ploidy2 = x$ploidy1,
      blast1 = x$blast1,
      blast2 = x$blast2,
      gff1 = x$gff1,
      gff2 = x$gff2,
      nmapsPerHaplotype = nmapsPerHaplotype,
      mcscan.dir = mcscan.dir,
      dbs_radii = dbs_radii,
      topPerc = topPerc,
      dbs_mappingsInRadius = dbs_mappingsInRadius,
      verbose = verbose)
    return(allmap)
  })
  return(out)
}
