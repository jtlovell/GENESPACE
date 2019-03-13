#' @title Cull blast to syntenic regions
#'
#' @description
#' \code{cull_syntenicBlast} An internal function, designed to be called
#' by extend_blocks and find_syntenicOrthogs.
#'
#' @param map the map data.table or data.frame. This is used to infer
#' syntenic regions.
#' @param genomeIDs the genomeIDs
#' @param gff gff data.table
#' @param dir.list Directory list
#' @param n.cores.blast Number of cores to use for blast
#' @param n.cores The number of parallel processes to run.
#' @param bed.buffer Buffer around found region in bp
#' @param best.hit.buffer Buffer to use to look for a better blast hit
#' @param max.window.bp The maximum size a bed region can be to run blast.
#' Smaller values speed things up, but lose regions.
#' @param verbose logical, should updates be printed?
#' @param ... Additional arguments passed to frNN
#'
#' @details Internal function
#'
#' @return A culled blast dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
find_pseudogenes <- function(map,
                             gff,
                             genomeIDs,
                             dir.list,
                             n.cores = 1,
                             n.cores.blast = 1,
                             bed.buffer = 1e3,
                             best.hit.buffer = 5e5,
                             max.window.bp = 1e5,
                             verbose = T){
  if(verbose)
    cat("##### Step. 1: Tracking Synteny ...\n")
  track.map <- track_synHits(map = map,
                             gff = gff,
                             max.window.bp = max.window.bp)

  if(verbose)
    cat("##### Step. 2: Writing fastas for each region ...\n")
  nomap <- convert_map2bufferBed(map = track.map$nomap,
                                 buffer = bed.buffer,
                                 genomeIDs = genomeIDs,
                                 assembly.dir = dir.list$assembly)
  map.wfa <- write_fa2file(map = nomap,
                           tmp.dir = dir.list$tmp,
                           n.cores = n.cores,
                           assembly.dir = dir.list$assembly,
                           peptide.dir = dir.list$peptide)

  if(verbose)
    cat("##### Step. 3: blastX against predicted regions ...\n")
  blast.regions <- blastx_map(map = map.wfa,
                              n.cores = n.cores.blast,
                              tmp.dir = dir.list$tmp,
                              diamond.sensitive = F,
                              clean = T,
                              verbose = T)

  if(verbose)
    cat("##### Step. 4: Choose best blast hit ...\n")
  exmap =  blast.regions[!is.na(blast.regions$score),]
  exmap$start2 <- exmap$blast.coord.start2
  exmap$end2 <- exmap$blast.coord.end2
  exmap <- exmap[,colnames(track.map$nomap),with = F]

  map.best <- choose_besthits(map = exmap,
                              peptide.dir = dir.list$peptide,
                              genomeIDs = genomeIDs,
                              buffer = best.hit.buffer)

  exmap.buf <- convert_map2bufferBed(map = map.best,
                                     genomeIDs = genomeIDs,
                                     assembly.dir = dir.list$assembly,
                                     buffer = bed.buffer)

  if(verbose)
    cat("##### Step. 5: Write peptide and fasta regions ...\nPulling fasta sequence ...")
  exmap.wfa <- get_unmapFastas(map = exmap.buf,
                               assembly.dir = dir.list$assembly,
                               tmp.dir = dir.list$tmp)
  if(verbose)
    cat("Done!\nWriting peptide sequence ...")
  pep.files <- write_peptideByGene(geneIDs = exmap.wfa$id1,
                                   genomeIDs = genomeIDs,
                                   peptide.dir = dir.list$peptide,
                                   tmp.dir = dir.list$tmp)
  if(verbose)
    cat("Done!\nFormatting for exonerate ...")
  setnames(pep.files,1,"id1")
  exmap.all <- merge(exmap.wfa, pep.files, by = "id1")

  exin <- with(exmap.all,
               data.table(chr = chr2,
                          start = bed.start2,
                          end = bed.end2,
                          reg.name = id2,
                          fa.file = fa.file,
                          pep.file = pep.file))

  if(verbose)
    cat("Done!\n##### Step. 6: Running exonerate ...\n\tCompleted:\n\t")
  exonerate.out <- run_exonerate(exin, clean = T)

  if(verbose)
    cat("##### Step. 7: Processing output ...\n")
  exon.sum <- with(exonerate.out,
                   calc_weightedId(gff.cds =  gff.cds,
                                   gff.region = gff.region))
  setnames(exon.sum, c("id2","chr2","strand2","start2","end2","weighted.id","cds.length"))
  id1.out <- exmap.all[,c("genome1","genome2",
                          "id1","id2","block.id","og",


                          "chr1","start1","end1")]
  ex.out <- merge(id1.out, exon.sum, by = "id2", all = T)

  if(verbose)
    cat("##### Pipeline complete!\n")
  return(list(exonerate.summary = ex.out,
              exonerate.cds = exonerate.out$cds,
              exonerate.gff = exonerate.out$gff.cds,
              blast.regions = blast.regions,
              best.hits = map.best,
              track.synteny = track.map))
}
