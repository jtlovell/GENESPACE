rerun_pairwiseOF <- function(dirs, 
                             gff,
                             genomeIDs, 
                             of.cores = 6,
                             verbose = T){
  
  gff <- subset(gff, genome %in% genomeIDs)
  # -- Step 1. Reformat peptides etc.
  if(verbose)
    cat("Preparing new orthofinder-formatted species ID database ... \n")
  make_newOFdb(tmp.dir = dirs$tmp, 
               cull.blast.dir = dirs$cull.blast, 
               peptide.dir = dirs$peptide, 
               genomeIDs = genomeIDs, 
               verbose = verbose)
  if(verbose)
    cat("\tDone!\n")
  #######################################################
  
  #######################################################

  #######################################################
  
  #######################################################
  if(verbose)
    cat("Importing new and old orthofinder gene and species IDs ... ")
  old.ids <- read_speciesIDs(of.dir = dirs$cull.score.blast, genomeIDs = genomeIDs)
  new.ids <- read_speciesIDs(of.dir = dirs$cull.blast, genomeIDs = genomeIDs)
  id.db <- merge(old.ids, new.ids, by = "genome")
  setnames(id.db,2:3,c("n.old","n.new"))
  #
  map.db <- make_mapDB(id.db = id.db,
                       blast.dir = dirs$cull.score.blast, 
                       cull.blast.dir = dirs$cull.blast)
  #
  old.genes <- read_geneIDs(of.dir = dirs$cull.score.blast, gff = gff)
  new.genes <- read_geneIDs(of.dir = dirs$cull.blast, gff = gff)
  genes <- merge(old.genes[,c("genome","id","gene.num")],
                 new.genes[,c("genome","id","gene.num")], 
                 by = c("genome","id"))
  g1 <- with(genes,  data.table(V1 = gene.num.x, new1 = gene.num.y, key = "V1"))
  g2 <- with(genes,  data.table(V2 = gene.num.x, new2 = gene.num.y, key = "V2"))
  if (verbose) 
    cat("\tDone!\n")
  #######################################################
  
  #######################################################
  if(verbose)
    cat("Reading blast file, replacing IDs and renaming ... \n")
  in.blast.files <- map.db$filename
  out.blast.files <- map.db$new.filename
  
  d <- lapply(1:nrow(map.db), function(i){
    if(verbose)
      cat(paste0("\t",map.db$genome1[i]),"-->",map.db$genome2[i],"... ")
    bl <-readRename_blastGenes(gene.dict1 = g1,
                               gene.dict2 = g2,
                               blast.file.in = in.blast.files[i],
                               blast.file.out = out.blast.files[i],
                               verbose = verbose)  
    if(verbose)
      cat("Done!\n")
    return(bl)
  })
  
  com <- paste("orthofinder", "-b", dirs$cull.blast, 
               "-a", of.cores, 
               "-og 1>/dev/null 2>&1")
  system(com)
  og.dt <- read_ogs(of.dir = dirs$cull.blast, gff = gff)
  return(og.dt)
}
  