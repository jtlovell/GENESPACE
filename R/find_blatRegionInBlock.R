find_blatRegionInBlock = function(blk,
                                  peptide.dir,
                                  gff.dir,
                                  tmp.dir,
                                  map.buffer = 1e5,
                                  bed.buffer = 1000,
                                  prop.best = .2,
                                  min.score = 50){

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

}
