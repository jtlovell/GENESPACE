add_homo2alignment <- function(peptide.dir,
                               of.dir,
                               results.dir){
  if(verbose)
    cat("Copying alignments to", aln.dir, "... ")
  aln.dir <- file.path(results.dir, "alignments")
  og.aln.dir <- list.files(file.path(of.dir, "Orthofinder/"), full.names =T)
  og.aln.dir <- file.path(og.aln.dir[length(og.aln.dir)],
                          "MultipleSequenceAlignments")
  dir.create(aln.dir)
  cpd <- file.copy(
    from = list.files(
      og.aln.dir,
      full.names = T),
    to = aln.dir,
    overwrite = T,
    recursive = F)

  if(verbose)
    cat("Done!\nWriting homologs to fasta for each subgraph ... ")
  genes2add <- pull_clusGenes(map = assblks)
  genelist <- lapply(split(subset(genes2add, !is.aligned), by = "og"),
                     function(x)  unique(x$gn))
  peptide.dir <- dir.locs$peptide

  pep.fs <- list.files(path = peptide.dir,
                       pattern = ".fa$",
                       full.names = T)
  pep.fastas <- do.call(c, lapply(pep.fs, function(x){
    tmp <- readAAStringSet(x)
    n <- gsub(".fa$","",basename(x))
    names(tmp) <- paste0(n,"_", names(tmp))
    return(tmp)
  }))
  pep.fastas.cull <- pep.fastas[unique(unlist(genelist))]
  pep.toadd <- lapply(genelist, function(x) pep.fastas.cull[x])
  wrote.fa <- data.table(t(sapply(names(pep.toadd), function(i){
    writeXStringSet(pep.toadd[[i]], file = file.path(aln.dir,paste0("toadd.",i,".fa")))
    c(file.path(aln.dir,paste0(i,".fa")),
      file.path(aln.dir,paste0("toadd.",i,".fa")))
  })))
  setnames(wrote.fa, c("orig.fa","add.fa"))
  wrote.fa[, out.fa := gsub(".fa$",".algn.fa",orig.fa)]
  wrote.fa[, com := paste("mafft --add",
                          add.fa,
                          "--keeplength --quiet",
                          orig.fa,">",
                          out.fa)]
  library(parallel)
  test <- lapply(1:nrow(wrote.fa), function(i){
    if(i %% 100 == 0)
      cat(i,", ")
    system(wrote.fa$com[i])
  })
  return(wrote.fa)
}
