write_pep4blastx <- function(){

  if(verbose)
    cat("Writing peptide files ... ")
  pf = pep.in[[1]]
  spl = split(1:length(pep.fastas), ceiling(1:length(pep.fastas)/100))
  spl.pep = lapply(spl, function(i) pep.fastas[i])
  spl.pep2 = lapply(1:length(spl.pep), function(j){
    print(j)
    x = spl.pep[[j]]
    n = length(x)
    return(lapply(1:n, function(i) x[i]))
  })

  sp = unlist(spl.pep2, recursive = F)
  pep.spl.files = lapply(1:length(sp), function(i){
    if(i %% 1000 == 0) cat(i,"/",length(sp),"\n")
    x = sp[[i]]
    f = file.path(tmp.dir,paste0(names(x),".spl.pep.fa"))
    writeXStringSet(x, filepath = f)
    return(f)
  })

  db.spl.files = lapply(1:length(pep.spl.files), function(i){
    if(i %% 100 == 0) cat(i,"/",length(pep.spl.files),"\n")
    x = pep.spl.files[[i]]
    db.file = gsub(".spl.pep.fa",".spl.db",x)
    system(paste("diamond makedb --quiet",
                 "--in", x,
                 "-d", db.file))
    return(db.file)
  })


  ########################################################
  if(verbose)
    cat("Writing peptide fastas to file ... \n")
  locs<-sapply(names(pep.out), USE.NAMES = T, simplify = F, function(i){
    if(verbose)
      cat("\tRunning",i,"...")
    po = pep.out[[i]]
    pep.files <- sapply(names(po), function(j){
      pep.file = file.path(tmp.dir, paste0(i,j,".pep.fa"))
      writeXStringSet(po[[j]], filepath = pep.file)
      return(pep.file)
    })
    if(verbose)
      cat("Done\n")
    return(pep.files)
  })

  pep.finfo = rbindlist(lapply(names(locs), function(i){
    x = locs[[i]]
    out<- rbindlist(lapply(names(x), function(j)
      data.table(genome = i,
                 block.id = j,
                 pep.file = x[[j]])))
    return(out)
  }))

}
