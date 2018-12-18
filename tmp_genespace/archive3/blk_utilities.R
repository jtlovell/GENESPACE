#' @title Orthofinder utility functions
#' @name of_utilities
#' @aliases run_exonerate
#' @aliases run_getfasta
#' @aliases run_diamondBlastx
#' @aliases make_blastRegion
#' @aliases parse_blastLoc
#' @aliases pipe_exonerate
#' @aliases add_buffer
#' @aliases cull_blastByDBS
#' @aliases cull_blast2og
#' @aliases merge_gffBlast
#' @aliases import_blast
#' @aliases import_ofResults
#' @aliases import_gff
#' @aliases parse_gff
#'
#' @description
#' \code{blk_utilities} Several utilities functions meant for internal calls in compareGeneSpace
#' @name of_utilities
#' @param gff.dir character, directory containing the gff3 formatted annotation files
#' @param blast.dir character, directory containing the orthofinder output
#' @param gff.file character, gff file name and path
#' @param genomeIDs character, genome identifiers
#' @param abbrevs character, genome abbreviations
#' @param blk data.table containing the block information
#' @param map data.table containing the map information
#' @param blast data.table containing blast hits
#' @param gff data.table containing the parsed gff annotation data
#' @param ogff list of data.tables, split gff by block.
#' @param orthogroups orthogroup object from import_ofResults
#' @param gene.index gene index from import_ofResults
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
#' @param out.dir directory where should output be written.
#' @param species.mappings the species mapping from parse_orthofinder
#' @param of.speciesIDs orthofinder species ids, from parse_orthofinder
#' @param of.geneIDs orthofinder gene ids, from parse_orthofinder
#' @param x data.table
#' @param f factor
#' @param by factor
#' @param drop drop
#' @param flatten non-recursive unlisting
#' @param verbose logical, should updates be reported?
#' @param ... additional arguments passed to data.table
#' @note \code{of_utilities} is a generic name for the functions documented.
#' \cr
#' If called, \code{blk_utilities} returns its own arguments.
#'

#' @title Fast split of data.table
#' @description
#' \code{run_getfasta} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
run_getfasta = function(bed, bed.file, ass.file, fa.file){
  write.table(bed,
              file=bed.file,
              quote=F,
              sep="\t", row.names=F, col.names=F)
  system(paste("bedtools getfasta -fi",
               ass.file, "-bed", bed.file,"-name -fo",fa.file))
}

#' @title Fast split of data.table
#' @description
#' \code{run_diamondBlastx} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
run_diamondBlastx = function(db.file,
                             pep.fa,
                             fa.file,
                             blast.file,
                             n.cores = 1,
                             max.target.seqs = 10000,
                             min.score = 20,
                             diamond.blastx.param = "--quiet"){
  system(paste("diamond makedb --quiet",
               "--in", pep.fa,
               "-d", db.file))
  system(paste("diamond blastx",
               "--threads",n.cores,
               "--max-target-seqs", max.target.seqs,
               "--min-score", min.score,
               diamond.blastx.param,
               "-d", db.file,
               "-q", fa.file,
               "-o",blast.file))
}

#' @title Fast split of data.table
#' @description
#' \code{parse_blastLoc} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
parse_blastLoc = function(bl){
  genome.info = lapply(bl$V1, function(x) strsplit(x,".", fixed = T)[[1]])
  bl$genome = sapply(genome.info, function(x) x[1])
  bl$block.id = sapply(genome.info, function(x) x[2])
  bl$chr = sapply(genome.info, function(x) x[3])
  bl$genome.start = as.numeric(sapply(genome.info, function(x) x[4]))
  bl$genome.end = bl$genome.start+bl$V8
  bl$genome.start = bl$genome.start+bl$V7
  bl$length = bl$genome.end - bl$genome.start
  bl$start = with(bl, ifelse(length < 0, genome.end, genome.start))
  bl$end = with(bl, ifelse(length > 0, genome.end, genome.start))
  bl$length = bl$end - bl$start
  setnames(bl, "V2", "id")
  bl$neg.score = bl$V12 * (-1)
  bl$neg.len = bl$length * (-1)
  setkey(bl, id, neg.score, neg.len)
  return(bl)
}

#' @title Fast split of data.table
#' @description
#' \code{make_blastRegion} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
make_blastRegion = function(bl,max.dist2besthit = 5e3){
  bl[,rank := frank(neg.score, ties.method = "random"),
     by = list(id)]
  bs = bl[bl$rank == 1,c("id","start","end")]
  setnames(bs, c("id","best.start","best.end"))
  setkey(bl,id)
  setkey(bs,id)
  bo = merge(bl, bs)
  bo$dist2best = with(bo, abs(best.start - start))
  bo = bo[bo$dist2best <= max.dist2besthit,]


  bi = bo[,list(start = min(start),
                end = max(end),
                mean.score = mean(V12),
                mean.length = mean(abs(neg.len))),
          by = list(genome, chr, id)]
  bi$length = with(bi, end-start)


  return(bi)
}

#' @title Fast split of data.table
#' @description
#' \code{add_buffer} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
add_buffer <- function(bl, fais, buffer = 1e3){
  bo = merge(fais, bl)
  bo$start <- bo$start - buffer
  bo$end <- bo$end + buffer
  bo$end[bo$end>bo$chr.length]<-bo$chr.length[bo$end>bo$chr.length]
  bo$start[bo$start<=0]<-1
  bo$reg.name = with(bo, paste0(id,"_",genome))
  return(bo)
}

#' @title Fast split of data.table
#' @description
#' \code{pipe_exonerate} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @export
pipe_exonerate <- function(hit.reg,
                           assembly.dir,
                           stas,
                           genomeIDs,
                           tmp.dir,
                           n.cores = 1,
                           min.score = 20,
                           verbose = T){

  exon.out = mclapply(1:nrow(hit.reg), mc.cores = n.cores, function(i){
    if(verbose)
      if(i %% 100 == 0)
        cat(paste0(i,", "))
    target.gene = hit.reg$id[i]
    genome = hit.reg$genome[i]
    bed = hit.reg[i,c("chr","start","end","reg.name"),with = F]
    colnames(bed)[4]<-"id"
    tmp <- run_exonerate(assembly.dir = dirs$assembly,
                         stas = stas,
                         target.gene = target.gene,
                         id = i,
                         genome = genome,
                         cds.genome = cds.genome,
                         bed = bed,
                         min.score = min.score,
                         tmp.dir = tmp.dir)
    return(tmp)
  })

  exon_cds <- do.call(c, lapply(exon.out, function(x) x$cds.seq))
  gff.gene <- rbindlist(lapply(exon.out, function(x) x$gff))
  gff.gene$ex.start = with(gff.gene, as.numeric(start)+as.numeric(exon.start))
  gff.gene$ex.end = with(gff.gene, as.numeric(start)+as.numeric(exon.end))
  gff.out <- gff.gene[,list(start = min(ex.start),
                            end = max(ex.end)),
                      by = list(id, chr, strand)]

  if(verbose)
    cat("\nDone!\n")

  return(list(gff = gff.out,
              exonerate_cds = exon_cds))

}

#' @title Fast split of data.table
#' @description
#' \code{run_exonerate} Much faster than base split.
#' @rdname blk_utilities
#' @import data.table
#' @import Biostrings
#' @export
run_exonerate <- function(ass.fasta.file,
                          pep.fasta.file,
                          bed,
                          min.score = 20,
                          tmp.dir){


  out1 = system(paste("exonerate --model protein2genome",
                      "--score",min.score,
                      "--softmasktarget --alignmentwidth 1000000",
                      "--bestn 1",
                      "--dpmemory 2000 --maxintron 10000",
                      "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                      '--ryo ">\n%tas\n"',
                      "--query",pep.fasta.file,
                      "--target", ass.fasta.file),
                intern = T)

  if(length(grep("^>", out1))==0){
    out1 = DNAStringSet("")
    names(out1) <- bed$id
    gff1 = data.frame(exon.start = NA,
                      exon.end = NA,
                      strand = NA,
                      align.info = NA)
    gffo = data.frame(bed,gff1)
  }else{
    end = which(out1 == "")[1]+1
    out1 = out1[1:end]
    gff1 = do.call(rbind,lapply(out1[grepl(bed$id, out1) & grepl("exon\t", out1)],
                                function(x) strsplit(x,"\t")[[1]][c(4:5,7,9)]))
    gff1 = data.frame(gff1,
                      stringsAsFactors = F)
    colnames(gff1)<-c("exon.start","exon.end","strand","align.info")
    gffo = data.frame(bed,gff1)

    wh = grep("^>",out1)+1
    out1 = DNAStringSet(paste(out1[wh:(length(out1)-2)], collapse = ""))
    names(out1)<-bed$id
  }
  return(list(cds.seq = out1, gff = gffo))
}

