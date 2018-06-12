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
pipe_MCScanX2Blocks = function(MCScanX.path,
                               mcs_mapping.dir,

                               MCScanX.params = "-a -s 5 -m 5",
                               buffer = 1,

                               tsp.method = "Concorde",
                               Concorde.path = NULL,
                               max.jump = 5,
                               ref.id = "PhalliiHAL",
                               min.block.size = 5,

                               plotit = T,
                               chr1toplot = "Chr03",
                               chr2toplot = "Chr_03",
                               altGenome2plot = "Sviridis",

                               verbose = T){

  if(verbose)
    cat("#################\n Running MCScanX\n")

  mcs.out = run_MCScanX(ref.id = ref.id,
                        MCScanX.path = MCScanX.path,
                        mcs_mapping.dir = mcscan.dir,
                        MCScanX.params = MCScanX.params,
                        buffer = buffer)

  if(plotit){
    plot_blocksAndMapping(map = mcs.out$map,
                          blk= mcs.out$block,
                          ref.id = ref.id,
                          altGenome2plot = altGenome2plot,
                          chr1toplot = chr1toplot,
                          chr2toplot = chr2toplot,
                          main = "Raw MCScanX blocks")
  }

  if(verbose)
    cat("\n#################\n Merging adjacent blocks\n")
  merged = merge_overlappingBlocks(map = mcs.out$map, blk = mcs.out$block,
                                   verbose = verbose)

  if(plotit){
    plot_blocksAndMapping(map = map,
                          blk= blk,
                          ref.id = ref.id,
                          altGenome2plot = altGenome2plot,
                          chr1toplot = chr1toplot,
                          chr2toplot = chr2toplot,
                          main = "Merged MCScanX blocks")
  }

  if(verbose)
    cat("\n#################\n Splitting within blocks via TSP solver\n")

  tsped <- split_blocksByTSP(map = merged$map,
                                     Concorde.path = Concorde.path)

  if(plotit){
    plot_blocksAndMapping(map = tsped$map,
                          blk= tsped$block,
                          ref.id = ref.id,
                          altGenome2plot = altGenome2plot,
                          chr1toplot = chr1toplot,
                          chr2toplot = chr2toplot,
                          main = "Merged MCScanX blocks")
  }
}
