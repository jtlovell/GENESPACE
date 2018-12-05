#' @title Orthofinder utility functions
#' @description
#' \code{of_utilities} Six utilities functions meant for internal calls in compareGeneSpace
#' @name of_utilities
#' @param gff.dir character, directory containing the gff3 formatted annotation files
#' @param blast.dir character, directory containing the orthofinder output
#' @param mcscanx.input.dir character, directory where MCScanX temporary files should be stored
#' @param gff.file character, gff file name and path
#' @param genomeIDs character, genome identifiers
#' @param abbrevs character, genome abbreviations
#' @param blk data.table containing the block information
#' @param map data.table containing the map information
#' @param blast data.table containing blast hits
#' @param gff data.table containing the parsed gff annotation data
#' @param orthogroups orthogroup object from import_ofResults
#' @param gene.index gene index from import_ofResults
#' @param species.mappings species mapping object from import_ofResults
#' @param MCScanX.params character, parameters to be passed to MCScanX
#' @param mcscanx.input.dir directory for mcscan temporary files to be stored
#' @param n.mappingWithinRadius numeric, number of hits required to be in the radius
#' @param eps.radius numeric, size of the radius
#' @param pairs.only logical, should only pairs of hits in orthofinder output be retained
#' @param min.propMax numeric, minimum proportion of max score for a gene
#' @param min.score numeric, minimum score for a hit to be retained
#' @param max.hitsPerGene numeric, maximum number of hits to be retained per gene
#' @param str2drop character, string in attribute column of gff file to be dropped
#' @param str2parse character, string in attribute column of gff file to use as the separator
#' @param whichAttr numeric, which attribute should be returned in the
#' gff attribute column
#' @param verbose logical, should updates be reported?
#' @note \code{of_utilities} is a generic name for the functions documented.
#' \cr
#' If called, \code{of_utilities} returns its own arguments.
#'

#' @title Fast split of data.table
#' @description
#' \code{split.data.table} Much faster than base split.
#' @param x data.table
#' @param f factor
#' @param by factor
#' @param drop drop
#' @param flatten non-recursive unlisting
#' @rdname utilities
#' @import data.table
#' @export
concatenate_orthoGroups = function(og.inblk, n.iter = 1, verbose){
  og = og.inblk$orthogroup.datatable
  og$unique = paste0(og$block.id, "_", og$og)
  spl.gen = split(og, og$genome)
  if(verbose)
    cat("Building species-specific databases ... completed:")
  spl.out = rbindlist(lapply(names(spl.gen), function(i){
    if(verbose)
      cat(paste0(i,", "))
    print(i)
    x = spl.gen[[i]]
    spl.id = split(x, x$id)
    og.uniq = lapply(spl.id, function(y) unique(y$unique))
    og.uout = data.table(id = rep(names(spl.id),sapply(og.uniq,length)), unique = unlist(og.uniq))
    og.uout$genome = i
    return(og.uout)
  }))

  spl.gen = split(og, og$block.id)
  spl.out2 = rbindlist(mclapply(names(spl.gen), mc.cores = 6, function(i){
    if(verbose)
      cat(paste0(i,", "))
    x = spl.gen[[i]]
    spl.id = split(x, x$og)
    gene.comb = rbindlist(lapply(spl.id, function(y) expand.grid(y$id, y$id)))
    return(gene.comb)
  }))
  splt = split(spl.out, spl.out$unique)
  outl = list()
  outu = vector()
  ns = 0
  u = unique(spl.out$id)
  for(i in u){
    if(which(u == i) %% 50 == 0) cat("\tCompleted",which(u == i),"/", length(u),"\n")
    if(!i %in% outu){
      ns = ns+1
      t1 = spl.out[spl.out$id == i,]
      t2 = spl.out[spl.out$unique %in% t1$unique,]
      t3 = spl.out[spl.out$id %in% t2$id,]

      t4 = spl.out[spl.out$unique %in% t3$unique,]
      t5 = spl.out[spl.out$id %in% t4$id,]
      out = unique(t5$id)
      outl[[i]] = out
      outu = c(outu, out)
    }
  }


  graph_from_data_frame(d, directed = TRUE, vertices = NULL)

  if(verbose)
    cat("Connecting orthogroups across pairwise comparisons ... completed\n")

}
make_blockDir = function(blk,
                         block.dir,
                         verbose = T){s
  if(file.exists(block.dir)){
    system(paste("rm -rf", block.dir))
  }
  system(paste("mkdir", block.dir))

  if(verbose)
    cat("Building directories for", nrow(blk), "blocks\n")
  blk$unique = with(blk, paste(genome1, genome2, block.id, sep = "_"))
  blk$dir = sapply(1:nrow(blk), function(i){
    x = file.path(block.dir,blk$unique[i])
    system(paste("mkdir",x))
    system(paste("mkdir",file.path(x,"genome")))
    for(j in c("assembly","cds", "gff", "peptide", "transcript")){
      system(paste("mkdir",file.path(x,"genome",j)))
    }
    for(j in c("tmp","blast")){
      system(paste("mkdir",file.path(x,j)))
    }
    return(x)
  })
  if(verbose)
    cat("Done!\n")
  return(blk)
}


chop_assemblyByBlock = function(blk,
                                assembly.dir,
                                verbose = T){

  if(verbose)
    cat("Chopping up assembly fasta by block breakpoints ... \n")
  for(i in 1:nrow(blk)){
    if(verbose)
      if(i %% 100 == 0) cat("\tCompleted:",i,"/", nrow(blk),"\n")
    x = blk[i,]
    g1 = x$genome1
    g2 = x$genome2
    f = x$dir
    bed1 = x[,c("chr1","start1","end1"),with = F]
    bed1$unique = with(bed1, paste(g1,chr1,start1,end1, sep = "."))
    bed2 = x[,c("chr2","start2","end2"), with = F]
    bed2$unique = with(bed2, paste(g2,chr2,start2,end2, sep = "."))

    bed1f = file.path(f,"genome","bed1.bed")
    bed2f = file.path(f,"genome","bed2.bed")
    fa1 = file.path(assembly.dir,paste0(g1,".fa"))
    fa2 = file.path(assembly.dir,paste0(g2,".fa"))
    fafo1 = file.path(f,"genome","assembly","seq1.fa")
    fafo2 = file.path(f,"genome","assembly","seq2.fa")


    write.table(bed1, file=bed1f,
                quote=F, sep="\t", row.names=F, col.names=F)
    write.table(bed2, file=bed2f,
                quote=F, sep="\t", row.names=F, col.names=F)
    system(paste("bedtools getfasta -fi",
                 fa1, "-bed", bed1f,"-name -fo",fafo1))
    system(paste("bedtools getfasta -fi",
                 fa2, "-bed", bed2f,"-name -fo",fafo2))
  }
  if(verbose)
    cat("Done!\n")
}

chop_annotationByBlock = function(sgff,
                                  blk){


  flist = list(peptide = peptide.dir,
               cds = cds.dir,
               transcript = transcript.dir)

  for(i in names(flist)){
    if(verbose)
      cat(i,"fastas ... Importing ...")

    fastas = lapply(list.files(flist[[i]], full.names = T), function(j){
      if(i == "peptide"){
        Biostrings::readAAStringSet(j)
      }else{
        Biostrings::readDNAStringSet(j)
      }
    })
    names(fastas)<-gsub(".fa","",list.files(flist[[i]], full.names = F),fixed = F)

    if(verbose)
      cat("Parsing and writing ...")
    test = lapply(1:nrow(blk), function(j){
      x = blk[j,]
      g = sgff[[x$unique]]
      sg = split(g$id, g$genome)
      t1 = fastas[[names(sg)[1]]][sg[[1]]]
      t2 = fastas[[names(sg)[2]]][sg[[2]]]
      p1 = file.path(x$dir,"genome",i,"seq1.fa")
      p2 = file.path(x$dir,"genome",i,"seq2.fa")
      Biostrings::writeXStringSet(t1, filepath = p1)
      Biostrings::writeXStringSet(t2, filepath = p2)
    })

    if(verbose)
      cat("Done\n")
  }
}


of_inblock = function(blk, ncores = 6){
  if(verbose)
    cat("Running orthofinder in blocks\n")
  test = mclapply(1:nrow(blk), mc.cores = ncores, mc.preschedule = F, function(i){
    x = blk[i,]
    run_orthofinder(
      peptide.dir = file.path(x$dir, "genome", "peptide"),
      og.threads = 1,
      blast.threads = 1,
      og.silent = T,
      verbose = F,
      tmp.dir = file.path(x$dir, "tmp"),
      blast.dir = file.path(x$dir, "blast"))
  })
  if(verbose)
    cat("Done!\n")
}
