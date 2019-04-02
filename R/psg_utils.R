#' @title Pseudogene utility functions
#' @description
#' \code{psg_utils} Five utilities functions meant for internal calls in compareGeneSpace
#' @name syn_utils
#'
#' @param map map results data.table
#' @param assembly.dir path to assembly fastas
#' @param tmp.dir path to temp directory
#' @param peptide.dir path to peptide directory
#' @param buffer numeric, the number of basepairs outside of the range to look at
#' @param genomeIDs character, indicating genomeIDs to consider.
#' @param clean logical, should the intermediate files be removed?
#' @param diamond.sensitive logical, should diamond blastx be run in the
#' sensitive mode?
#'
#' @param verbose logical, should updates be printed?
#' @param ... not currently in use
#'
#' @note \code{psg_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{psg_utils} returns its own arguments.
#'
#' @title choose best hits
#' @description
#' \code{choose_besthits} for any unmapped region, choose the best peptide
#' @rdname psg_utils
#' @import data.table
#' @importFrom dbscan frNN dbscan
#' @importFrom Biostrings width readAAStringSet
#' @export
choose_besthits <- function(map,
                            genomeIDs,
                            peptide.dir,
                            buffer = 1e4,
                            verbose = T){
  if(verbose)
    cat("Determining physical clusters of blast hits within orthogroups ... ")
  map.all <- map
  if("n" %in% colnames(map))
    map <- map[map$n > 0,]

  map[,max.width := max(abs(end2-start2))+buffer,
      by = list(og, genome2, chr2)]
  map$mean2 <- rowMeans(map[,c("start2", "end2")])

  rdbs <- function(x, eps){
    nn <- frNN(data.frame(x,x),
               eps = eps)
    return(dbscan(nn, minPts = 0)$cluster)
  }

  map[,clus := rdbs(x = mean2, eps = max.width[1]),
      by = list(og, genome2, chr2)]
  if(verbose)
    cat("Done!\n")

  map$uniq.og <- with(map, paste0(og, "_", clus))
  spl = split(map, "uniq.og")

  pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readAAStringSet(file.path(peptide.dir,
                              paste0(i, ".fa")))))
  pep.fastas <- pep.fastas[unique(map$id1)]
  pep.len <- data.table(id1 = names(pep.fastas),
                        n.res = width(pep.fastas))
  map <- merge(map, pep.len, by = "id1")

  if(verbose)
    cat("Choosing best hits in unique orthogroup regions ... ")
  map$map.width = with(map, abs(end2 - start2))
  map$n.res1 <- map$n.res * (-1)
  map$map.width1 <- map$map.width * (-1)

  map$g1 <- as.numeric(factor(map$genome1, levels = genomeIDs))
  map$g2 <- as.numeric(factor(map$genome1, levels = genomeIDs))
  map$gdif <- with(map, abs(g1-g2))
  setkey(map, gdif, n.res1, map.width1)

  map.best <- map[, head(.SD[1]), list(uniq.og, genome2)]
  map.best <- map.best[,colnames(map.all),with = F]
  if(verbose)
    cat("Done!\n")
  return(rbind(map.all[map.all$n == 0,], map.best))
}

#' @title get fasta from id2 data in map
#' @description
#' \code{getfasta_map} get fasta from id2 data in map
#' @rdname psg_utils
#' @import data.table
#' @export
getfasta_map <- function(map,
                         genomeIDs,
                         assembly.dir,
                         tmp.dir,
                         buffer = 1e3,
                         verbose = T,
                         clean = T){

  if(verbose)
    cat("Adding", buffer, "bp buffer to unmapped region boundaries ... ")

  bed.in <- convert_map2bufferBed(map = map,
                                  assembly.dir = assembly.dir,
                                  buffer = buffer,
                                  genomeIDs = genomeIDs)
  if(verbose)
    cat("Done!\n")


  # -- Get fasta sequence of regions
  if(verbose)
    cat("Extracting fasta sequence of regions ...\n")

  spl.bed <- split(bed.in, "genome2")
  getfa.files <- rbindlist(lapply(names(spl.bed), function(i){
    x = spl.bed[[i]]
    if(verbose)
      cat(paste0("\t",i,": n. regions = ", nrow(x),"... "))
    splx = split(x, "chr2")
    out <- rbindlist(lapply(names(splx), function(j){
      y = splx[[j]]
      bed.file = file.path(tmp.dir, paste0(i,".",j,".bed"))
      ass.file = file.path(assembly.dir, paste0(i,".fa"))
      fa.file = file.path(tmp.dir, paste0(i,".",j,".fa"))
      get_fasta(bed = with(y,
                           data.frame(chr = chr2,
                                      start = bed.start2,
                                      end = bed.end2,
                                      name = id2)),
                bed.file = bed.file,
                ass.file = ass.file,
                fa.file = fa.file)
      if(clean){
        unlink(bed.file)
      }
      return(data.table(y,
                        fa.loc = fa.file))
    }))
    if(verbose)
      cat("Done!\n")
    return(out)
  }))
  return(getfa.files)
}

#' @title convert map to a buffer bed-format
#' @description
#' \code{convert_map2bufferBed} Add buffer and convert to 0-index
#' @rdname psg_utils
#' @import data.table
#' @export
convert_map2bufferBed <- function(genomeIDs,
                                  assembly.dir,
                                  buffer,
                                  map){

  assem.fas <- file.path(assembly.dir,paste0(genomeIDs,".fa"))
  if(!all(file.exists(
    file.path(assembly.dir,
              paste0(genomeIDs,".fa.fai"))))){
    if(verbose)
      cat("Assembly fasta indices not found, generating ... ")
    for (i in assem.fas)
      system(paste("samtools faidx", i))
    if(verbose)
      cat("Done!\n")
  }

  indexs <- rbindlist(sapply(genomeIDs, simplify = F, USE.NAMES = T, function(i){
    tmp <- fread(file.path(assembly.dir, paste0(i,".fa.fai")))
    out <- with(tmp,
                data.table(genome2 = i,
                           chr2 = V1,
                           chr2.start = 1,
                           chr2.end = V2,
                           stringsAsFactors = F))
    return(out)
  }))

  bed <- merge(indexs, map, by = c("genome2", "chr2"))
  bed$start.tmp <- with(bed, start2 - buffer)
  bed$end.tmp <- with(bed, end2 + buffer)
  bed$bed.start2 = with(bed,
                        ifelse(start.tmp >= chr2.start,
                               start.tmp, chr2.start) - 1)
  bed$bed.end2 = with(bed,
                      ifelse(end.tmp <= chr2.end,
                             end.tmp, chr2.end) - 1)
  bed$start.tmp <- NULL
  bed$end.tmp <- NULL
  bed$chr2.start <- NULL
  bed$chr2.end <- NULL
  return(bed)
}

#' @title blastx a map
#' @description
#' \code{blastx_map} pull all fasta sequences and peptides for a
#' chromsome-genome combination and blastx them
#' @rdname psg_utils
#' @import data.table
#' @export
blastx_map <- function(map,
                       n.cores = 6,
                       tmp.dir,
                       diamond.sensitive = F,
                       clean = T,
                       verbose = T){

  spl = split(map, "reg.id")

  if(diamond.sensitive){
    diamond.param =  paste("--threads",n.cores,
                           "--sensitive",
                           "--quiet")
  }else{
    diamond.param =  paste("--threads",n.cores,
                           "--quiet")
  }

  if(verbose)
    cat("Running", length(spl), "blastx searches ...\n")

  bl.regs <- lapply(names(spl), function(i){
    if(which(names(spl) == i) %% 10 == 0)
      cat("\t\tCompleted:",which(names(spl) == i),"/", length(spl),"\n")
    x = spl[[i]]
    if(verbose)
      cat(paste0("\t",i),"... ")
    db.file = file.path(tmp.dir, i)

    fa.file = x$reg.fa[1]
    pep.file = x$pep.fa[1]
    blast.file = file.path(tmp.dir, paste0(i,".blast.txt"))
    run_diamondBlastx(db.file = db.file,
                      pep.fa = pep.file,
                      fa.file = fa.file,
                      blast.file = blast.file,
                      top = 1000,
                      diamond.blastx.param = diamond.param)

    rl <- length(readLines(blast.file))
    if (rl < 1) {
      if(verbose)
        cat("Found no hits among", nrow(x),"links\n")
      bl <- x
      bl$blast.coord.start2 <- NA
      bl$blast.coord.end2 <- NA
      bl$blast.start <- NA
      bl$blast.end <- NA
      bl$score <- NA
      out.bl <- bl
    }else{
      bl <- fread(blast.file)
      bl <- bl[,c(1,2,7,8,12)]
      setnames(bl, c("id2","id1","st","en","score"))
      bl$blast.start = apply(bl[,c("st","en")],1,min)
      bl$blast.end = apply(bl[,c("st","en")],1,max)
      bl$en <- NULL
      bl$st <- NULL
      out.bl <- merge(x, bl, by = c("id1","id2"), all.x = T)
      out.bl$blast.coord.start2 <- with(out.bl, blast.start + bed.start2)
      out.bl$blast.coord.end2 <- with(out.bl, blast.end + bed.start2)

      out.bl$score1 <- out.bl$score * (-1)
      setkey(out.bl, score1)
      out.bl <- out.bl[, head(.SD[1]), list(id1, id2)]
      out.bl$score1 <- NULL
      if(verbose)
        cat("Found",sum(!is.na(out.bl$score)),
            "hits out of", nrow(x),
            "links. Median score = ",
            median(out.bl$score,na.rm = T),"\n")
      setkey(out.bl, start1)
    }
    if(which(names(spl) == i) == length(spl))
      cat("Done!\n")
    if(clean){
      unlink(c(paste0(db.file,".dmnd"), pep.file, fa.file, blast.file))
    }
    return(out.bl)
  })
  return(rbindlist(bl.regs))
}

#' @title get unmapped fastas
#' @description
#' \code{get_unmapFastas} for each entry in map, pull out the fasta for id2
#' @rdname psg_utils
#' @import data.table
#' @export
get_unmapFastas <- function(map,
                            assembly.dir,
                            tmp.dir,
                            clean = T){
  spl = split(map, "id2")
  out <- rbindlist(lapply(names(spl), function(j){
    y = spl[[j]]
    bed.file = file.path(tmp.dir, paste0(j,".bed"))
    ass.file = file.path(assembly.dir, paste0(y$genome2[1],".fa"))
    fa.file = file.path(tmp.dir, paste0(j,".fa"))
    get_fasta(bed = with(y,
                         data.frame(chr = chr2,
                                    start = bed.start2,
                                    end = bed.end2,
                                    name = id2)),
              bed.file = bed.file,
              ass.file = ass.file,
              fa.file = fa.file)
    if(clean){
      unlink(bed.file)
    }
    y$fa.file <- fa.file
    return(y)
  }))
  return(out)
}

#' @title write_fa2file
#' @description
#' \code{write_fa2file} swrite_fa2file
#' @rdname psg_utils
#' @import data.table
#' @export
write_fa2file <- function(map,
                          chunk.size = 100,
                          n.cores = 6,
                          tmp.dir,
                          assembly.dir,
                          peptide.dir,
                          verbose = T){
  pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readAAStringSet(file.path(peptide.dir,
                              paste0(i, ".fa")))))
  map$id2 <- with(map, paste(genome2,genome1,1:nrow(map), sep = "."))
  spl <- split(map, "genome2")
  reg.map <- rbindlist(lapply(names(spl), function(j){
    y = spl[[j]]
    if(verbose)
      cat("\tWriting", nrow(y),j,"hits")
    y$reg.id = with(y, paste0(genome2, ".", genome1, ".reg",
                              ceiling((1:nrow(y))/chunk.size)))
    sply <- split(y, "reg.id")
    if(verbose)
      cat(" split into", length(sply), paste0(chunk.size, "-hit blocks ... "))
    yo <- rbindlist(mclapply(names(sply), mc.cores = n.cores, function(k){
      z = sply[[k]]
      bed.file = file.path(tmp.dir, paste0(k,".bed"))
      ass.file = file.path(assembly.dir, paste0(j,".fa"))
      fa.file = file.path(tmp.dir, paste0(k, ".fa"))
      get_fasta(bed = with(z,
                           data.frame(chr = chr2,
                                      start = bed.start2,
                                      end = bed.end2,
                                      name = id2)),
                bed.file = bed.file,
                ass.file = ass.file,
                fa.file = fa.file)
      z$reg.fa <- fa.file
      unlink(bed.file)
      peps <- pep.fastas[unique(z$id1)]
      pep.file = file.path(tmp.dir, paste0(k, ".pep.fa"))
      writeXStringSet(peps, filepath = pep.file)
      z$pep.fa <- pep.file
      return(z)
    }))
    if(verbose)
      cat("Done!\n")
    return(yo)
  }))
  return(reg.map)
}

#' @title write peptide by gene
#' @description
#' \code{write_peptideByGene} split up multi-fasta into single fasta seqs
#' @rdname psg_utils
#' @import data.table
#' @export
write_peptideByGene <- function(geneIDs,
                                genomeIDs,
                                peptide.dir,
                                tmp.dir){
  pep.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readAAStringSet(file.path(peptide.dir,
                              paste0(i, ".fa")))))
  geneIDs <- unique(geneIDs)
  pep.fastas <- pep.fastas[geneIDs]
  out <- sapply(geneIDs, function(i){
    out.loc <- file.path(tmp.dir, paste0(i,".fa"))
    writeXStringSet(pep.fastas[i], filepath = out.loc)
    return(out.loc)
  })
  return(data.table(id = geneIDs,
                    pep.file = out))
}

#' @title calculate weighted Id
#' @description
#' \code{calc_weightedId} calculate weighted across exons from exonerate output
#' @rdname psg_utils
#' @import data.table
#' @export
calc_weightedId <- function(gff.cds, gff.region){
  gff.cds$length <- with(gff.cds, as.numeric(exon.end) - as.numeric(exon.start))
  gff.cds$identity <- sapply(gff.cds$align.info, function(x)
    as.numeric(gsub(" ","",gsub("identity","",strsplit(x,";", fixed = T)[[1]][3]))))

  gff.cds[,tot.len := sum(length),
          by = id]
  gff.cds$identity <- as.numeric(gff.cds$identity)
  gff.cds$tot.len <- as.numeric(gff.cds$tot.len)
  gff.cds$length <- as.numeric(gff.cds$length)
  gff.cds$wt.id <- with(gff.cds, identity * (length/tot.len))

  out <- gff.cds[,list(weighted.id = sum(wt.id),
                       cds.length = max(tot.len)),
                 by = id]
  setkey(out, id)
  setkey(gff.region, id)
  ret <- merge(gff.region, out)
  return(ret)
}

#' @title get fasta
#' @description
#' \code{get_fasta} call bedtools getfasta
#' @rdname psg_utils
#' @import data.table
#' @export
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
