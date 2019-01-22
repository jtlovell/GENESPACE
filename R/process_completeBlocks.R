process_completeBlocks <- function(complete.regions,
                                   orthogroup.list,
                                   gff,
                                   cds.dir,
                                   verbose = T){
  if (verbose)
    cat("Merging complete block list with gff annotation ... ")
  # -- combine network list into data.table
  cblk = complete.regions
  setnames(cblk, 4:5, c("blk.start","blk.end"))

  # -- make original orthogroup data.table
  lens = sapply(og.list, length)
  og.dt = data.table(og.id = rep(names(og.list),
                                 lens),
                     id = unlist(og.list))

  # -- merge orthogroup data.table with gff annotations
  setkey(og.dt, "id")
  setkey(gff, "id")
  gff.og = merge(og.dt,gff)
  if (verbose)
    cat("Done!\nQC-ing orthogroup networks ... ")
  # -- merge orthogroup / gff within complete orthogroup data.table
  setkey(gff.og, og.id, genome, chr)
  setkey(cblk, og.id, genome, chr)
  comp.gff = merge(gff.og, cblk, all = T)

  # -- drop unnecessary columns
  comp.gff$strand <- NULL
  comp.gff$order <- NULL

  # -- drop rows where gene is not in block
  in.reg <- apply(comp.gff[,5:8],1,function(x)
    x[1] <= x[4] & x[2] >= x[3])
  in.reg[is.na(in.reg)] <- T
  comp.gff2 <- comp.gff[in.reg,]

  if(verbose)
    cat("Done!\nAdding in CDS length information ... ")
  # -- add in CDS length
  cds.fastas <- do.call(c, lapply(genomeIDs, function(i)
    readDNAStringSet(file.path(cds.dir,
                               paste0(i,".fa")))))
  cds.md <- data.table(id = names(cds.fastas),
                       length = width(cds.fastas),
                       stringsAsFactors = F)
  setkey(cds.md, id)
  setkey(comp.gff2, id)
  comp.gff3 <- merge(comp.gff2, cds.md, all.x = T)

  # -- return complete genome set
  comp.gff3<-comp.gff3[!is.na(comp.gff3$genome),]
  if(verbose)
    cat("Done!\n")
  return(comp.gff3)
}
