#' @title Run the exonerate program
#'
#' @description
#' \code{pipe_exonerateOrphans} A simple wrapper to run orthofinder from R.
#'
#' @param map Map object
#' @param dir.list List of directories, produced by check_environment
#' @param buffer numeric, the amount of space around the focal region to feed
#' to blastx and exonerate.
#' @param verbose Logical, should updates be printed?
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
pipe_pseudogenes <- function(map,
                             buffer = 1e3,
                             diamond.sensitive = F,
                             dir.list,
                             blast.ncores = 6,
                             genomeIDs,
                             verbose = T){
  assembly.dir = dir.list$assembly
  tmp.dir = dir.list$tmp
  peptide.dir = dir.list$peptide

  if(verbose)
    cat("Pulling fasta regions for each link ... \n\t")
  map.wregs = getfasta_map(map = map,
                           assembly.dir = assembly.dir,
                           buffer = 1e3,
                           tmp.dir = tmp.dir,
                           genomeIDs = genomeIDs)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Running blastx for each chromosome-genome combination ... \n")
  blast.regions = blastx_map(map = map.wregs,
                             tmp.dir = tmp.dir,
                             n.cores = blast.ncores,
                             diamond.sensitive = diamond.sensitive,
                             peptide.dir = dirs$peptide,
                             verbose = verbose)

  if(verbose)
    cat("\tDone!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Formatting and adding buffer to blastx hits ... ")
  exmap =  blast.regions[complete.cases(blast.regions),]
  exmap$start2 <- exmap$blast.coord.start2
  exmap$end2 <- exmap$blast.coord.end2
  exmap <- exmap[,colnames(synteny.track),with = F]

  exmap.buf <- convert_map2bufferBed(map = exmap,
                                     genomeIDs = genomeIDs,
                                     assembly.dir = dirs$assembly,buffer = 1e3)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Pulling fasta regions around each blast hit ... ")
  exmap.wfa <- get_unmapFastas(map = exmap.buf,
                               assembly.dir = dirs$assembly,
                               tmp.dir = dirs$tmp)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Writing peptide fasta for each gene of interest ... ")
  pep.files <- write_peptideByGene(geneIDs = exmap.wfa$id1,
                                   genomeIDs = genomeIDs,
                                   peptide.dir = dirs$peptide,
                                   tmp.dir = dirs$tmp)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Running exonerate for each pair ... completed:\n\t")
  setnames(pep.files,1,"id1")
  exmap.all <- merge(exmap.wfa, pep.files, by = "id1")

  exin <- with(exmap.all,
               data.table(chr = chr2,
                          start = bed.start2,
                          end = bed.end2,
                          reg.name = id2,
                          fa.file = fa.file,
                          pep.file = pep.file))
  exonerate.out <- run_exonerate(exin)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  if(verbose)
    cat("Formatting exonerate output ... ")
  exon.sum <- with(exonerate.out,
                   calc_weightedId(gff.cds =  gff.cds,
                                   gff.region = gff.region))
  setnames(exon.sum, c("id2","chr2","strand2","start2","end2","weighted.id","cds.length"))
  id1.out <- exmap.all[,c("genome1","genome2",
                          "id1","id2","block.id","og","n",
                          "chr1","start1","end1")]
  ex.out <- merge(id1.out, exon.sum, by = "id2", all = T)
  if(verbose)
    cat("Done!\n")
  ########################################################

  ########################################################
  return(list(blast.regions = blast.regions,
              ex.summary = ex.out,
              ex.gff = exonerate.out$gff.cds,
              ex.cds = exonerate.out$cds))
}
