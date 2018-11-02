#' @title Find and align pairwise best hits within blocks
#'
#' @description
#' \code{make_alignmentsInBlocks} Main function to find orthologs and align them within syntenic blocks
#'
#' @param gff.dir The path to the directory containing gff files
#' @param peptide.dir The path to the directory containing peptide fasta files
#' @param cds.dir The path to the directory containing cds fastas
#' @param tmp.dir The path to the temporary directory
#' @param results.dir The path to the directory where results should be written.
#' @param blk The block object (data.frame)
#' @param buffer Numeric, the physical distance (bp) from the end of the blocks and the last
#' best blat hit to extend the search.
#' @param genes2ignore A set of gene IDs that should not be considered. Typically
#' generate from find_genesWithoutOrthos.
#' @param ncores How many parallel processes should be run at once? Should be
#' limmited to the number of physical cores on the machine.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details This is the second step in the pipeline and generates unannotated sequences
#' that are most similar (and deemed orthologous by orthofinder).
#' @return paths do directories containing new CDS files. These contain
#' both the original annotated CDSs and new best match CDSs.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
blatExonerate_inblock = function(blk,
                                   gff.dir,
                                   peptide.dir,
                                   tmp.dir,
                                   results.dir,
                                   genes2ignore = NULL,
                                   buffer = 50000,
                                   ncores = 8,
                                   verbose = T){

  #### - Prep directories
  blk.dir = file.path(tmp.dir,"byblock")
  if(file.exists(blk.dir)){
    system(paste("rm -r", blk.dir))
  }
  system(paste("mkdir",blk.dir))

  # - Gff files
  gff.files = list.files(gff.dir, full.names = T)
  names(gff.files) = gsub(".gff3$","",basename(gff.files))

  parse_gff = function(gff){
    g = suppressWarnings(
      data.table::fread(gff,showProgress = F, verbose = F))
    g = g[g$V3 == "gene",c(9,1,4,5,7)]
    g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
    data.table::setnames(g, c("id","chr","start","end","strand"))
    return(g)
  }

  gff = rbindlist(lapply(names(gff.files), function(i){
    tmp = parse_gff(gff.files[[i]])
    tmp$genome = i
    tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
    return(tmp)
  }))

  # -- Peptide fastas
  peptide.fastas = lapply(list.files(peptide.dir, full.names = T), function(i){
    Biostrings::readAAStringSet(i)
  })
  names(peptide.fastas)<-gsub(".fa", "", list.files(peptide.dir))

  # -- CDS fastas
  cds.fastas = lapply(list.files(cds.dir, full.names = T), function(i){
    Biostrings::readAAStringSet(i)
  })
  names(cds.fastas)<-gsub(".fa", "", list.files(cds.dir))


  of_inblock = function(gff,peptide.fastas, blockline, genes2ignore, buffer){
    x = blockline
    od = file.path(blk.dir,
                   x$block.id)
    if(file.exists(od)){
      system(paste("rm -r", od))
    }
    system(paste("mkdir", od))

    p1 = file.path(od, "pep1.fa")
    p2 = file.path(od, "pep2.fa")

    bed1 = gff[gff$genome == x$genome1 &
                  gff$end >= (x$start1 - buffer) &
                  gff$start <= (x$end1 - buffer) &
                  gff$chr == x$chr1,]
    bed2 = gff[gff$genome == x$genome2 &
                  gff$end >= (x$start2 - buffer) &
                  gff$start <= (x$end2 - buffer) &
                  gff$chr == x$chr2,]
    b1 = bed1$id
    b2 = bed2$id

    bed1 = data.frame(chr = x$chr1,
                      start = min(bed1$start)-1,
                      end = max(bed1$end)-1,
                      id = paste0(x$genome1,".",x$block.id),
                      stringsAsFactors = F)
    bed2 = data.frame(chr = x$chr2,
                      start = min(bed2$start)-1,
                      end = max(bed2$end)-1,
                      id = paste0(x$genome2,".",x$block.id),
                      stringsAsFactors = F)
    bed = list(bed1,bed2)
    names(bed)<-c(x$genome1, x$genome2)
    t1 = peptide.fastas[[x$genome1]][b1]
    t2 = peptide.fastas[[x$genome2]][b2]
    t1 = t1[!names(t1) %in% genes2ignore]
    t2 = t2[!names(t2) %in% genes2ignore]
    Biostrings::writeXStringSet(t1, filepath = p1)
    Biostrings::writeXStringSet(t2, filepath = p2)

    system(paste("orthofinder -f",
                 od,
                 "-t 1 -S diamond -og 1>/dev/null 2>&1"))
    of_res = readLines(list.files(od, pattern = "^Orthogroups.txt$", recursive = T, full.names = T))
    og = lapply(of_res, function(x) strsplit(x," ")[[1]])
    og = lapply(og, function(x) x[-1])
    ol = sapply(og, length)
    sing = data.table(id = unlist(og[ol == 1]))
    out = merge(sing, gff)
    outl = split(out$id, out$genome)
    system(paste("rm -r", od))
    return(list(single.genes = outl,
                bed = bed))
  }

  blat_singleGenes = function(single.genes, cds.fastas, blockline){
    od = file.path(blk.dir,
                   blockline$block.id)
    if(file.exists(od)){
      system(paste("rm -r", od))
    }
    system(paste("mkdir", od))

    genomes = names(single.genes[[2]])
    g1 = genomes[1]
    g2 = genomes[2]
    genes = single.genes[[1]]
    c1 = file.path(od, "cds1.fa")
    c2 = file.path(od, "cds2.fa")
    t1 = cds.fastas[[g1]][genes[[g1]]]
    t1 = t1[!names(t1) %in% genes2ignore]
    t2 = cds.fastas[[g2]][genes[[g2]]]
    t2 = t2[!names(t2) %in% genes2ignore]
    Biostrings::writeXStringSet(t1, filepath = c1)
    Biostrings::writeXStringSet(t2, filepath = c2)

    bed1 = single.genes[[2]][[g1]]
    bed2 = single.genes[[2]][[g2]]


    bedf1 = file.path(od,paste0(bed1$id,".bed"))
    bedf2 = file.path(od,paste0(bed2$id,".bed"))
    bedf1a = file.path(od,paste0(bed1$id,".bed1"))
    bedf2a = file.path(od,paste0(bed2$id,".bed1"))
    write.table(bed1, file=bedf1,
                quote=F, sep="\t", row.names=F, col.names=F)
    write.table(bed2, file=bedf2,
                quote=F, sep="\t", row.names=F, col.names=F)

    faf1 = file.path(assembly.dir,paste0(g1,".fa"))
    faf2 = file.path(assembly.dir,paste0(g2,".fa"))
    fab1 = file.path(od,paste0(bed1$id,".fa"))
    fab2 = file.path(od,paste0(bed2$id,".fa"))

    system(paste("bedtools getfasta -fi",
                 faf1, "-bed", bedf1,"-name -fo",fab1))
    system(paste("bedtools getfasta -fi",
                 faf2, "-bed", bedf2,"-name -fo",fab2))

    blat1 = file.path(od,paste0(g1,"_",g2,".blast"))
    blat2 = file.path(od,paste0(g2,"_",g1,".blast"))
    system(paste("blat", fab1, c2, "-t=dna -q=dna -out=blast8 -maxIntron=20000 -minScore=100",
                 blat1,"1>/dev/null 2>&1"))
    system(paste("blat", fab2, c1, "-t=dna -q=dna -out=blast8 -maxIntron=20000 -minScore=100",
                 blat2,"1>/dev/null 2>&1"))

    blat2bed = function(blat, bed, blk.buffer = 50000, min.prop = .5, merge = T){
      d = fread(blat)
      d[, prop := V12/max(V12),
        by = list(V1)]
      d <- d[d$prop>=min.prop,]
      b.start = tapply(apply(d[,9:10],1,min),d$V1,min)-blk.buffer
      b.end = tapply(apply(d[,9:10],1,max),d$V1,min)+blk.buffer
      b.strand = tapply(ifelse(d$V10 > d$V9,"+","-"),d$V1, function(x) names(sort(table(x),decreasing=TRUE))[1])

      b.start = b.start+bed$start
      b.start[b.start<0]<-0

      wg = bed$end - bed$start
      b.end[b.end>wg]<-wg
      b.end = b.end+bed$start
      bo = data.frame(chr = bed$chr,
                       start = b.start,
                       end = b.end,
                       id = names(b.start),
                       score = "1",
                       strand = b.strand)
      rownames(bo)<-NULL
      return(bo)
    }

    blbed1 = blat2bed(blat = blat1, bed = bed1)
    blbed2 = blat2bed(blat = blat2, bed = bed2)

    write.table(blbed1, file=bedf1,
                quote=F, sep="\t", row.names=F, col.names=F)
    write.table(blbed2, file=bedf2,
                quote=F, sep="\t", row.names=F, col.names=F)

    blbeds = list(blbed1,blbed2)
    names(blbeds)<-c(g1,g2)

    system(paste("bedtools getfasta -fi",
                 faf1, "-bed", bedf1,"-name -s -fo",fab1))
    system(paste("bedtools getfasta -fi",
                 faf2, "-bed", bedf2,"-name -s -fo",fab2))

    cds1 = cds.fastas[[g1]]

    ss1 <- do.call(c,lapply(1:nrow(blbed1), function(i){
      x = blbed1[i,]
      write.table(x, file=bedf1,
                  quote=F, sep="\t", row.names=F, col.names=F)
      system(paste("bedtools getfasta -fi",
                   faf1, "-bed", bedf1,"-name -s -fo",fab1))
      Biostrings::writeXStringSet(cds.fastas[[g2]][as.character(x$id)], filepath = c1)

      # t2 = cds.fastas[[g2]][genes[[g2]]]
      # Biostrings::writeXStringSet(t1, filepath = c1)
      out1 = system(paste("exonerate --model est2genome",
                          "--forcescan target --softmasktarget",
                          "--gappedextension --refine region",
                          "--dpmemory 2000 --maxintron 10000 -n 1 ",
                          "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                          '--ryo ">\n%tas\n"',
                          "--query",c1, "--target", fab1), intern = T)
      gff1 = do.call(rbind,lapply(out1[grepl(x$id, out1) & grepl("exon\t", out1)],
                                  function(x) strsplit(x,"\t")[[1]][c(4:5,7)]))
      if(nrow(gff1)>1){
        cords = paste(apply(gff1[,1:2], 1, function(x) paste(x, collapse = "-")), collapse= ",")
      }else{
        cords = paste(gff1[,1:2], collapse = "-")
      }

      n = paste0(bed1$id,";",x$id,";",x$chr,";",cords,";",gff1[1,3])

      wh = grep("^>",out1)+1
      out1 = DNAStringSet(paste(out1[wh:(length(out1)-2)], collapse = ""))
      names(out1)<-n
      return(out1)
    }))

    ss2 = do.call(c,lapply(1:nrow(blbed2), function(i){
      x = blbed2[i,]
      write.table(x, file=bedf2,
                  quote=F, sep="\t", row.names=F, col.names=F)
      system(paste("bedtools getfasta -fi",
                   faf2, "-bed", bedf2,"-name -s -fo",fab2))
      Biostrings::writeXStringSet(cds.fastas[[g1]][as.character(x$id)], filepath = c2)

      # t2 = cds.fastas[[g2]][genes[[g2]]]
      # Biostrings::writeXStringSet(t1, filepath = c1)
      out1 = system(paste("exonerate --model est2genome",
                          "--forcescan target --softmasktarget",
                          "--gappedextension --refine region",
                          "--dpmemory 2000 --maxintron 10000 -n 1 ",
                          "--showalignment no --showtargetgff yes --showvulgar no --showcigar no",
                          '--ryo ">\n%tas\n"',
                          "--query",c2, "--target", fab2), intern = T)
      gff1 = do.call(rbind,lapply(out1[grepl(x$id, out1) & grepl("exon\t", out1)],
                                  function(x) strsplit(x,"\t")[[1]][c(4:5,7)]))
      if(nrow(gff1)>1){
        cords = paste(apply(gff1[,1:2], 1, function(x) paste(x, collapse = "-")), collapse= ",")
      }else{
        cords = paste(gff1[,1:2], collapse = "-")
      }

      n = paste0(bed2$id,";",x$id,";",x$chr,";",cords,";",gff1[1,3])

      wh = grep("^>",out1)+1
      out1 = DNAStringSet(paste(out1[wh:(length(out1)-2)], collapse = ""))
      names(out1)<-n
      return(out1)
    }))

    return(list(ss1,ss2))
  }

  if(!null(genes2ignore)){
    gff_c<-gff[!gff$id %in% genes2ignore,]
  }
  single.genes = of_inblock(blockline = merged$block[i,],
                            gff = gff_c, buffer = 10000,
                            genes2ignore = genes2ignore,
                            peptide.fastas = peptide.fastas)
  init.bed = blat_singleGenes(single.genes,
                              cds.fastas,
                              blockline = merged$block[i,])

  pipe_inblock = function(){

    final.gff = exonerate_singleGene()
    gff2bed12
    pull_unannotatedSeq()
    combine_cdsAndUnannotated()
    of_inblockFinal()
  }



}
