#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{run_orthoFinderInBlock} Utility function to build metadata
#'
#' @param map The map object (data.frame or data.table)
#' @param blk The block object (data.frame)
#' @param buffer Numeric, the overlapping distance between two blocks.
#' 0 indicates that blocks that overlap by >=0 should be merged.
#' @param verbose Logical, should updates be printed.
#' @param ... Not currently in use
#' @details Primarily used in the run_MCScanX pipeline.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
run_orthoFinderInBlock = function(blk,
                                  gff.dir,
                                  peptide.dir,
                                  tmp.dir,
                                  results.dir,
                                  buffer = 50000,
                                  ncores = 4,
                                  onlyParseOF = F){

  if(!onlyParseOF){

    output.dir = file.path(tmp.dir,"peptide_ortho_inblk")

    if(file.exists(output.dir)){
      system(paste("rm -r", output.dir))
    }
    system(paste("mkdir",output.dir))

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


    peptide.fastas = lapply(list.files(peptide.dir, full.names = T), function(i){
      Biostrings::readAAStringSet(i)
    })
    names(peptide.fastas)<-gsub(".fa", "", list.files(peptide.dir))

    blk.genes = mclapply(1:nrow(blk), mc.cores = ncores, mc.preschedule = F, function(i){

      x = blk[i,]
      outdir = file.path(tmp.dir,
                         x$block.id)
      if(file.exists(outdir)){
        system(paste("rm -r", outdir))
      }
      system(paste("mkdir", outdir))
      p1 = file.path(outdir, "pep1.fa")
      p2 = file.path(outdir, "pep2.fa")
      i1 = x[,c("genome1","start1","end1")]
      b1 = gff$id[gff$genome == x$genome1 &
                    gff$end >= (x$start1 - buffer) &
                    gff$start <= (x$end1 - buffer) &
                    gff$chr == x$chr1]
      b2 = gff$id[gff$genome == x$genome2 &
                    gff$end >= (x$start2 - buffer) &
                    gff$start <= (x$end2 - buffer) &
                    gff$chr == x$chr2]

      t1 = peptide.fastas[[x$genome1]][b1]
      t2 = peptide.fastas[[x$genome2]][b2]
      Biostrings::writeXStringSet(t1, filepath = p1)
      Biostrings::writeXStringSet(t2, filepath = p2)

      system(paste("orthofinder -f",
                   outdir,
                   "-t 1 -S diamond -og 1>/dev/null 2>&1"))
      ortho.file = list.files(outdir, pattern = "^Orthogroups.txt$", recursive = T, full.names = T)
      system(paste("mv",
                   ortho.file,
                   file.path(output.dir, paste0(x$block.id,".txt"))))
      system(paste("rm -r", outdir))
    })
  }


  fs = list.files(output.dir,
                  pattern = ".txt$", full.names = T)
  names(fs) = gsub(".txt","",
                   list.files(output.dir,
                              pattern = ".txt$", full.names = F))
  og.list = lapply(fs, function(x){
    og = readLines(x)
    og = lapply(og, function(x) strsplit(x," ")[[1]])
    ons = sapply(og, function(x) x[1])
    names(og)<-ons
    og = lapply(og, function(x) x[-1])
    ol = sapply(og, length)
    oln = as.numeric(names(table(ol)))
    out = sapply(oln, USE.NAMES = T, simplify = F, function(i){
      do.call(rbind, og[ol == i])
    })
    names(out)<-names(table(ol))
    return(out)
  })

  names(og.list)<-names(fs)

  p.single =
    unlist(sapply(og.list, function(x)
      length(unlist(x[["1"]]))/length(unlist(x))))
  p.single = data.frame(block.id = names(p.single),
                        p.single = p.single,
                        stringsAsFactors = F)
  p.2 =
    unlist(sapply(og.list, function(x)
      length(unlist(x[["2"]]))/length(unlist(x))))
  p.2 = data.frame(block.id = names(p.2),
                   p.2 = p.2,
                   stringsAsFactors = F)
  p = merge(p.single, p.2, by = "block.id", all = T)
  p$p.multi = 1-rowSums(p[,-1])

  sing = rbindlist(lapply(names(og.list), function(x) data.table(id = og.list[[x]][[1]],
                                                       block = x)))
  setnames(sing, c("id","block.id"))
  sing = merge(gff, sing, by = "id")

  return(list(orthogroups = og.list, proptype = p, single.genes = sing))
}

