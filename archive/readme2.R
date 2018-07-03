##############################

## get into py 2.7
conda create -n py27 python=2.7 anaconda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install orthofinder
conda install diamond

# Part 1 - BLAST via diamond
## 1.1 Set up the BLAST scheme
library(data.table)
MCScanX.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX/MCScanX"

blast.scheme = rbind(c("PhalliiHAL","Phallii",2,2),
               c("PhalliiHAL","Sviridis",2,2),
               c("PhalliiHAL","Sbicolor",2,2),
               c("PhalliiHAL","Pvirgatum",2,4),
               c("Phallii","Sviridis",2,2),
               c("Phallii","Sbicolor",2,2),
               c("Phallii","Pvirgatum",2,4),
               c("Sviridis","Sbicolor",2,2),
               c("Sviridis","Pvirgatum",2,4),
               c("Sbicolor","Pvirgatum",2,4),
               c("PhalliiHAL","PhalliiHAL",2,2),
               c("Phallii","Phallii",2,2),
               c("Sviridis","Sviridis",2,2),
               c("Sbicolor","Sbicolor",2,2),
               c("Pvirgatum","Pvirgatum",4,4))
blast.scheme = cbind(blast.scheme, apply(blast.scheme[,1:2],1, function(x)
  file.path("/Users/jlovell/Documents/comparative_genomics/blast",
            paste0("Blast",paste(x, collapse = "_"),".txt"))))
blast.scheme = cbind(blast.scheme, apply(blast.scheme[,1:2],1, function(x)
  file.path("/Users/jlovell/Documents/comparative_genomics/blast",
            paste0("Blast",paste(x, collapse = "_"),".culled.txt"))))
blast.scheme = cbind(blast.scheme, apply(blast.scheme[,1:2],1, function(x)
  file.path("/Users/jlovell/Documents/comparative_genomics/mcs_mapping",
            paste0(paste(x, collapse = "_"),".gff"))))
blast.scheme = cbind(blast.scheme, apply(blast.scheme[,1:2],1, function(x)
  file.path("/Users/jlovell/Documents/comparative_genomics/mcs_mapping",
            paste0(paste(x, collapse = "_"),".blast"))))


peptide.files = list.files("/Users/jlovell/Documents/comparative_genomics/peptide",
                           pattern = ".pep.fa", full.names = T)
names(peptide.files)<-gsub(".pep.fa","",basename(peptide.files))

db.files = file.path("/Users/jlovell/Documents/comparative_genomics/diamond_db",names(peptide.files))
names(db.files)<-names(peptide.files)

gff.files = list.files("/Users/jlovell/Documents/comparative_genomics/gff",
                       pattern = "gff3$", full.names = T)
names(gff.files)<-names(peptide.files)

abbrevs = c("Ph","Pa","Sv","Sb","Pv")
names(abbrevs) = c("PhalliiHAL","Phallii","Sviridis","Sbicolor","Pvirgatum")



align_peptideByDiamond(peptide.files, db.files,
                       blast.scheme)
blast.list = parse_blast(blast.scheme, nmapsPerHaplotype = 2)
merged.blasts = merge_gffAndBlast(blast.list, gff.files)
culled.blasts = cull_blastByDens(merged.blasts)
prep_MCScanX(blast.scheme,culled.blasts, abbrevs)
mcscan.raw = run_MCScanX(blast.scheme, MCScanX.path)
parsed.map = parse_MCScanX(mcscan.raw,gff.files)
parsed.blocks = make_blocks(parsed.map)
merged = merge_overlappingBlocks(map = parsed.blocks$map,
                                 blk = parsed.blocks$block,buffer = 0,
                                 verbose = TRUE)

merged = merge_adjacentBlocks(map = parsed.blocks$map,
                              blk = parsed.blocks$block,buffer = 5,
                              verbose = TRUE)
initm = parsed.blocks
while(nrow(initm$block) > nrow(merged$block)){
  initm = merged
  merged = merge_adjacentBlocks(map = merged$map,
                                blk = merged$block,buffer = 5,
                                verbose = TRUE)
}

blast.scheme = cbind(blast.scheme, apply(blast.scheme[,1:2],1, function(x)
  file.path("/Users/jlovell/Documents/comparative_genomics/of_input",
            paste0(paste(x, collapse = "_"),".blast"))))





og = readLines("/Users/jlovell/Documents/comparative_genomics/peptide/Results_Jun18_2/Orthogroups.txt")
og = lapply(og, function(x) strsplit(x," ")[[1]])
ons = sapply(og, function(x) x[1])
names(og)<-ons
og = lapply(og, function(x) x[-1])

system.time(test <- lapply(gff$id[1:20], function(y) which(sapply(og, function(x) y %in% x))))

a = og
al = unlist(a)

# a <- list(1:3, 4:5, 6:9)
# b <- c(2, 3, 5, 8)
g <- rep(seq_along(a), sapply(a, length))


gff.of = gff
b = gff.of$id
gff.of$of.index = g[match(b, al)]


#
prep_OrthoFinder = function(blast.scheme, blast.dt,
                            peptide.files,
                            abbrev,gff.files,
                            of_wd = "/Users/jlovell/Documents/comparative_genomics/of_input"){
  num = (1:length(abbrevs)-1)
  nom = names(abbrevs)
  species_ids = paste0(num,": ", nom)
  cat(species_ids, file = file.path(of_wd,"SpeciesIDs.txt"), sep = "\n")

  index = cbind(expand.grid(nom,nom),
                expand.grid(num,num))
  index2 = index
  rownames(index)<-apply(index[,1:2],1, paste, collapse = "_")
  index$unom = paste0(index[,3],"_", index[,4])
  index$unum = paste0(index[,1],"_", index[,2])


  ##########################
  # gff
  gff = rbindlist(lapply(names(gff.files), function(i){
    tmp = parse_gff(gff.files[[i]])
    tmp$genome = i
    tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
    return(tmp)
  }))

  gff$genome = factor(gff$genome, levels = nom)
  setkey(gff, genome, order)

  sequence_ids = data.frame(genome = paste0(gff$genome, "_", gff$order,":"),
                            id = gff$id, stringsAsFactors = F)
  for(i in 1:length(nom)){
    sequence_ids$genome<-gsub(nom[i],num[i],sequence_ids$genome)
  }

  seqo = apply(sequence_ids, 1, function(x) paste(x, collapse = " "))
  cat(seqo,  file = file.path(of_wd,"SequenceIDs.txt"), sep = "\n")

  sequence_ids$genome = gsub(":","",sequence_ids$genome)
  rownames(sequence_ids)<-sequence_ids$id
  ##########################
  # fasta
  for(i in nom){
    nu = num[which(nom == i)]
    fa = Rsamtools::scanFa(peptide.files[i],as = "AAStringSet")
    names(fa)<-sequence_ids[names(fa),]$genome
    fao = file.path(of_wd,paste0("Species",nu,".fa"))
    Biostrings::writeXStringSet(fa, fao)
  }



  blast = lapply(blast.scheme[,5], fread)
  bn = gsub(".txt$","",basename(blast.scheme[,5]))
  bn =gsub("^Blast","", bn)
  names(blast)<-paste0(blast.scheme[,1],"_", blast.scheme[,2])

  ind_no = index[!rownames(index) %in% names(blast),]
  for(i in 1:nrow(ind_no)){
    bl = blast[[which(blast.scheme[,1] == ind_no[i,2] &
                 blast.scheme[,2] == ind_no[i,1])]]
    bl = bl[,c(2,1,3:6,9:10,7:8,11:12), with = F]
    blast[[rownames(ind_no)[i]]]<-bl
  }
  for(i in 1:nrow(index)){
    id = rownames(index)[i]
    tmp = blast[[i]]
    fn = file.path(of_wd,paste0("Blast",index[i,5],".txt"))
    write.table(tmp, file = fn, sep = "\t", quote = F,
                row.names = F, col.names = F)

  }

  genome1 = sapply(bn, function()

  all.genes$genome = factor(all.genes$genome, levels = nom)
  setkey(all.genes, genome)
  all.genes[, rank:=frank(chr, ties.method = "dense"),by=genome]
  dt1[, rank := frank(date),
      by = list(id)]

  sequence_ids = all.genes
  for(i in 1:length(nom)){
    sequence_ids$genome<-gsub(nom[i],num[i],sequence_ids$genome)
  }
  sequence_ids$genome = paste0()
}



plot_blocksAndMapping(map = parsed.blocks$map,
                      blk= parsed.blocks$block,
                      ref.id = "Phallii",
                      altGenome2plot = "PhalliiHAL",
                      chr1toplot = paste0("Chr01"),
                      chr2toplot =  paste0("Chr01"),
                      main = "Raw MCScanX blocks")


plot_blocksAndMapping(map = merged$map,
                      blk= merged$block,
                      ref.id = "PhalliiHAL",
                      altGenome2plot = "Pvirgatum",
                      chr1toplot = paste0("Chr0",1:9),
                      chr2toplot =  paste0("Chr0",1:9,"N"),
                      main = "Raw MCScanX blocks")

spled <-split_byDensity(map = merged$map,
                        quantile = .999)

merged = merge_adjacentBlocks(map = parsed.blocks$map,
                            blk = parsed.blocks$block,buffer = 5,
                            verbose = TRUE)
initm = parsed.blocks
while(nrow(initm$block) > nrow(merged$block)){
  initm = merged
  merged = merge_adjacentBlocks(map = merged$map,
                                blk = merged$block,buffer = 5,
                                verbose = TRUE)
}
merged = merge_adjacentBlocks(map = merged$map,
                            blk = merged$block,buffer = 1,
                            verbose = TRUE)
merged = merge_adjacentBlocks(map = merged$map,
                              blk = merged$block,buffer = 1,
                              verbose = TRUE)

plot_blocksAndMapping(map = spled$map,
                      blk= spled$block,
                      ref.id = "PhalliiHAL",
                      altGenome2plot = "Sviridis",
                      chr1toplot = paste0("Chr08"),
                      chr2toplot =  paste0("Chr_08"),
                      main = "Raw MCScanX blocks")

merge_adjacentBlocks =   function(map, blk, buffer = 1.5,
                                  verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  if(verbose)
    cat("Looking for overlapping blocks ...\n")

  whichSegmentsIntersect <- function(segment1, segment2) {
    getBoundingBox <- function(P0, P1) {
      llx <- apply(cbind(P0[,1], P1[,1]),1, min)
      lly <- apply(cbind(P0[,2], P1[,2]),1, min)
      urx <- apply(cbind(P0[,1], P1[,1]),1, max)
      ury <- apply(cbind(P0[,2], P1[,2]),1, max)
      return(cbind(llx, lly, urx, ury))
    }
    box1 <- getBoundingBox(segment1[,1:2], segment1[,3:4])
    box2 <- getBoundingBox(segment2[,1:2], segment2[,3:4])
    out = apply(box2, 1, function(k){
      chk1 <- box1[1] <= k[3]
      chk2 <- box1[3] >= k[1]
      chk3 <- box1[2] <= k[4]
      chk4 <- box1[4] >= k[2]
      return(chk1 && chk2 && chk3 && chk4)
    })
    return(out)
  }



  blk$block.id = as.character(blk$block.id)
  blk$uniq = paste0(blk$genome1,"_",blk$genome2,"_",blk$chr1, "_", blk$chr2)
  spl = split(blk, blk$uniq)
  nr = sapply(spl, nrow)
  spl = spl[nr>1]
  nm = 1
  for(i in names(spl)){
    x = spl[[i]]
    x = x[order(x$start1, x$start2),]

    has.ovl = sapply(1:nrow(x), function(y){
      k = x[,c("rankstart1","rankstart2","rankend1","rankend2")]
      z = x[y,c("rankstart1","rankstart2","rankend1","rankend2")]
      z[,1]<-max(z[,1]-buffer,0)
      z[,2]<-max(z[,2]-buffer,0)
      z[,3]<-max(z[,3]+buffer,0)
      z[,4]<-max(z[,4]+buffer,0)
      return(whichSegmentsIntersect(z, k))
    })



    colnames(has.ovl)<-x$block.id
    rownames(has.ovl)<-x$block.id
    al = any(apply(has.ovl,1, all))
    diag(has.ovl)<-FALSE
    an = any(has.ovl)
    if(!al & an){
      has.ovl = has.ovl[rowSums(has.ovl)>0, colSums(has.ovl)>0]
      t1 = apply(has.ovl, 1, function(x) rownames(has.ovl)[which(x)])
      t2 = apply(has.ovl,2, function(x) rownames(has.ovl)[which(x)])
      mlist = list()
      for(j in colnames(has.ovl)){
        tmp = unique(unlist(c(j, t1[sapply(t1, function(x) any(grepl(j,x)))],
                t2[sapply(t2, function(x) any(grepl(j,x)))],
                names(t1)[sapply(t1, function(x) any(grepl(j,x)))],
                names(t2)[sapply(t2, function(x) any(grepl(j,x)))])))
        mlist[[j]] = tmp[order(tmp)]
      }
      mlist = mlist[!duplicated(mlist)]

      for(j in 1:length(mlist)){
        map$block.id[map$block.id %in% mlist[[j]]]<-mlist[[j]][1]
      }
    }else{
      if(al){
        map$block.id[map$block.id %in% x$block.id]<-x$block.id[1]
      }
    }



  }
  out = make_blocks(map)
  map = data.frame(out[["map"]], stringsAsFactors = F)
  blk = data.frame(out[["block"]], stringsAsFactors = F)
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}


###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
split_byDensity<-function(map, max.dist = 5,
                          quantile = 0.999, verbose = T){

  split_it = function(map, radius,quantile){
    spl = split(map, map$block.id)
    out = do.call(rbind, lapply(names(spl), function(i){
      x = spl[[i]]
      x$x_a = frank(x$start1,  ties.method = "dense")
      x$x_b = frank(x$start2, ties.method = "dense")
      x$was_split = F
      if(nrow(x)>radius*2){
        nn = dbscan::frNN(x[,c("x_a","x_b")], eps = radius)

        xbar = mean(unlist(nn$dist))
        xsd = sd(unlist(nn$dist))
        qthr = qnorm(quantile, mean = xbar, sd = xsd)

        dbs = dbscan::dbscan(nn, minPts = qthr)
        x$newclust = dbs$cluster
        if(all(x$newclust==0)){
          x$newclust<-1
        }
        x = x[x$newclust!=0,]
        x$block.id = paste0(x$block.id, "_",x$newclust)
        x$was_split = length(unique(dbs$cluster))>1
        x$newclust = NULL
      }
      return(x)
    }))
  }

  ntospl = 1
  map$was_split = TRUE
  while(any(map$was_split)){
    if(verbose)
      cat("checking", sum(map$was_split[!duplicated(map$block.id)]),"blocks\n")
    map = split_it(map,
                   radius = sqrt((max.dist^2)*2),
                   quantile = quantile)
  }
  map$block.id = as.character(as.numeric(as.factor(map$block.id)))
  out = make_blocks(map)
  if(verbose)
    cat("split into",nrow(out$block),"... Done")
  return(list(map = out$map, block = out$block))
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
plot_blocksAndMapping = function(map,
                                 blk,
                                 ref.id,
                                 altGenome2plot,
                                 chr1toplot,
                                 chr2toplot,
                                 main = NULL){
  tpb = blk[blk$genome1 == ref.id & blk$genome2 == altGenome2plot,]
  tp = map[map$genome1 == ref.id & map$genome2 == altGenome2plot,]
  tpb$s1 = with(tpb, ifelse(orient == "+", rankstart1, rankend1))
  tpb$e1 = with(tpb, ifelse(orient == "+", rankend1, rankstart1))
  tpb$rankstart1 = tpb$s1
  tpb$rankend1 = tpb$e1
  tpb$s1 = NULL
  tpb$e1 = NULL
  if(!is.null(chr1toplot)){
    tpb = tpb[tpb$chr1 %in% chr1toplot,]
    tp = tp[tp$chr1 %in% chr1toplot,]
  }
  if(!is.null(chr2toplot)){
    tpb = tpb[tpb$chr2 %in% chr2toplot,]
    tp = tp[tp$chr2 %in% chr2toplot,]
  }
  sb1 = split(tpb, tpb$chr1)
  st1 = split(tp, tp$chr1)
  for(i in names(sb1)){
    sb2 = split(sb1[[i]], sb1[[i]]$chr2)
    st2 = split(st1[[i]], st1[[i]]$chr2)
    for(j in names(sb2)){
      t2 = st2[[j]]
      b2 = sb2[[j]]
      cols = rep_len(c("red3","salmon","darkorange3","gold",
                       "grey50","lightgreen","forestgreen","darkgreen",
                       "cyan","dodgerblue3","violet","purple"), length.out = nrow(b2))
      with(t2, plot(rank1, rank2,
                    col = cols[as.numeric(as.factor(block.id))],
                    pch = 16, cex = .5,
                    ylab = paste(altGenome2plot,j,"gene order"),
                    xlab = paste(ref.id,i,"gene order"),
                    main = main))
      with(b2, segments(x0 = rankstart1, x1 = rankend1,
                        y0 = rankstart2, y1 =rankend2,
                        col = "black", lwd = 1.5))
      with(b2, text(x = rowMeans(b2[,c("rankstart1","rankend1")]),
                    y = rowMeans(b2[,c("rankstart2","rankend2")]),
                    labels = block.id,
                    col = "black", cex = .5, adj = c(1,-1)))
    }
  }
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
merge_overlappingBlocks = function(map, blk, buffer = 1.5,
                                   verbose = T){
  if(verbose)
    cat("Parsing",nrow(blk), "blocks and", nrow(map),"mappings\n")
  if(verbose)
    cat("Looking for overlapping blocks ...\n")

  inrect = function(left,top,right, bottom, points, buffer = 1){
    apply(points,1, function(k){
      k[1]>(left+buffer) &
        k[1]<(right-buffer) &
        k[2]>(bottom+buffer) &
        k[2]<(top-buffer)
    })
  }

  blk$block.id = as.character(blk$block.id)
  blk$uniq = paste0(blk$genome1,"_",blk$genome2,"_",blk$chr1, "_", blk$chr2)
  spl = split(blk, blk$uniq)
  nr = sapply(spl, nrow)
  spl = spl[nr>1]
  for(i in names(spl)){
    x = spl[[i]]
    x = x[order(x$start1, x$start2),]

    has.ovl = sapply(1:nrow(x), function(y)
      rowSums(cbind(inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("start1","start2")],
                           buffer = buffer),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("end1","start2")],
                           buffer = buffer),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("start1","end2")],
                           buffer = buffer),
                    inrect(right = x$end1[y], left = x$start1[y],
                           top = x$end2[y],bottom = x$start2[y],
                           points = x[,c("end1","end2")],
                           buffer = buffer))))
    colnames(has.ovl)<-x$block.id
    rownames(has.ovl)<-x$block.id
    has.ovl = has.ovl>=1

    if(sum(has.ovl)>0){
      if(sum(colSums(has.ovl)>0)==1){
        cn = colnames(has.ovl)[colSums(has.ovl)>0]
        wh = rownames(has.ovl)[which(has.ovl[,colSums(has.ovl)>0])]
        mlist = list(wh)
        names(mlist) = cn
      }else{
        mlist = apply(has.ovl[,colSums(has.ovl)>0],2,function(y) names(which(y)))
      }
      for(j in 1:length(mlist)){
        map$block.id[map$block.id %in% mlist[[j]]]<-names(mlist)[j]
      }
    }
  }
  out = make_blocks(map)
  map = data.frame(out[["map"]], stringsAsFactors = F)
  blk = data.frame(out[["block"]], stringsAsFactors = F)
  if(verbose)
    cat("Done! Returning",nrow(blk), "blocks and", nrow(map),"mappings\n")
  return(list(block = blk, map = map))
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
make_blocks<-function(map){
  out.blk = with(map,
                 data.table(
                   block.id = as.character(tapply(block.id, block.id, function(x) x[1])),
                   genome1 = tapply(genome1, block.id, function(x) x[1]),
                   genome2 = tapply(genome2, block.id, function(x) x[1]),
                   chr1 = tapply(chr1, block.id, function(x) x[1]),
                   chr2 = tapply(chr2, block.id, function(x) x[1]),
                   start1 = tapply(start1, block.id, min),
                   start2 = tapply(start2, block.id, min),
                   end1 = tapply(end1, block.id, max),
                   end2 = tapply(end2, block.id, max),
                   rankstart1 = tapply(rank1, block.id, min),
                   rankstart2 = tapply(rank2, block.id, min),
                   rankend1 = tapply(rank1, block.id, max),
                   rankend2 = tapply(rank2, block.id, max),
                   stringsAsFactors = F))
  orient = sapply(split(map, map$block.id), function(x) cor(x$start1, x$start2))
  orient = data.table(block.id = names(orient), orient = ifelse(orient>0,"+","-"))
  out.blk = merge(out.blk, orient, by = "block.id")
  if(!any(is.na(as.numeric(out.blk$block.id)))){
    tord = order(as.numeric(out.blk$block.id))
  }else{
    tord = order(out.blk$block.id)
  }
  out.blk = out.blk[tord,]
  map = data.frame(map[order(map$genome1, map$chr1, map$start1),], stringsAsFactors = F)
  blk = data.frame(out.blk, stringsAsFactors = F)
  return(list(block = blk, map = map))
}
###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######

parse_MCScanX = function(mcscan.raw, gff.files){
  gff = rbindlist(lapply(names(gff.files), function(i){
    tmp = parse_gff(gff.files[[i]])
    tmp$genome = i
    tmp$order = frank(tmp[,c("chr","start")], ties.method = "random")
    return(tmp)
  }))

  gff1 = data.table(gff)
  setnames(gff1, paste0(colnames(gff1),1))
  gff2 = data.table(gff)
  setnames(gff2, paste0(colnames(gff2),2))

  fac = sapply(as.character(mcscan.raw$V1),function(x) strsplit(x,"-")[[1]][1])
  m = data.table(mcscan.raw[,2:3])
  setnames(m, c("id1","id2"))
  m$block.id = as.numeric(as.factor(fac))
  out = merge(gff1, merge(m, gff2, by = "id2"), by = "id1")
  setkey(out, block.id, start1)
  out$rank1 = frank(out[,c("chr1","start1")], ties.method = "dense")
  out$rank2 = frank(out[,c("chr2","start2")], ties.method = "dense")
  return(out)
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######

run_MCScanX<-function(blast.scheme, MCScanX.path,
                       MCScanX.params = "-a -s 5 -m 25", verbose = T){

  if(verbose)
    cat("concatenating gff and blast files\n")
  outgff = file.path(dirname(blast.scheme[1,7]),"all.gff")
  ingff= paste(blast.scheme[,7], collapse = " ")
  system(paste("cat", ingff,">", outgff))

  outblast = file.path(dirname(blast.scheme[1,7]),"all.blast")
  inblast= paste(blast.scheme[,8], collapse = " ")
  system(paste("cat", inblast,">", outblast))
  mcscan.input = file.path(dirname(blast.scheme[1,7]),"all")

  if(is.null(MCScanX.params)){
    com = paste(file.path(MCScanX.path,"MCScanX"),mcscan.input)
  }else{
    com = paste(file.path(MCScanX.path,"MCScanX"),MCScanX.params,mcscan.input)
  }
  if(verbose)
    cat("running MCScanX via:\n\t", com,"\n")

  system(com)
  system(paste("rm",outgff))
  system(paste("rm", outblast))
  if(verbose)
    cat("MCScanX - Done!")
  tmp = read.delim(paste0(mcscan.input,".collinearity"),
                   sep = "\t", header  =F,
                   comment.char = "#", strip.white = T,
                   stringsAsFactors = F)
  return(tmp)
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######

prep_MCScanX = function(blast.scheme,culled.blasts, abbrevs){
  spl = split(culled.blasts, paste(culled.blasts$species_id1, culled.blasts$species_id2))

  trsh = lapply(1:nrow(blast.scheme), function(i){
    x = spl[[paste(blast.scheme[i,1:2], collapse = " ")]]
    gff1 = x[,c("chr1","id1","start1","end1")]
    chr1 =  as.numeric(gsub("[^0-9\\.]", "", gff1$chr1))
    id1 = x$species_id1[1]
    abbrev1 = abbrevs[id1]
    gff1$chr1 = paste0(abbrev1,chr1)

    gff2 = x[,c("chr2","id2","start2","end2")]
    chr2 =  as.numeric(gsub("[^0-9\\.]", "", gff2$chr2))
    id2 = x$species_id2[1]
    abbrev2 = abbrevs[id2]
    gff2$chr2 = paste0(abbrev2,chr2)
    names(gff2) = names(gff1)
    gff = rbind(gff1, gff2)

    write.table(gff, file = blast.scheme[i,7],
                quote = F, sep = "\t", row.names = F, col.names = F)

    blast = x[,c("id1","id2","perc.iden","align.length",
                 "n.mismatch", "n.gapOpen", "q.start", "q.end",
                 "s.start", "s.end", "eval", "score")]

    write.table(blast, file = blast.scheme[i,8],
                quote = F, sep = "\t", row.names = F, col.names = F)

  })
}


###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
merge_gffAndBlast = function(blast.list,
                             gff.files){
  gff.list = lapply(gff.files, parse_gff)
  for(i in 1:length(gff.list)) gff.list[[i]]$species_id = names(gff.files)[i]
  blasts = rbindlist(blast.list)

  gff1 = rbindlist(gff.list)
  gff2 = rbindlist(gff.list)
  setnames(gff1, c(paste0(colnames(gff1),1)))
  setnames(gff2, c(paste0(colnames(gff2),2)))

  out = merge(gff1,merge(gff2, blasts, by = "id2"), by = "id1")
  return(out)
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
cull_blastByDens = function(merged.blasts,
                            radii = c(100, 40, 10,20),
                            mappings = c(5, 5, 5,10),
                            verbose = T){
  spl = split(merged.blasts, paste(merged.blasts$species_id1, merged.blasts$species_id2))


  outlist = lapply(spl, function(x){
    tmp = loop_dbs(map = x,
                   radii = radii,
                   mappings = mappings,
                   plotit =T)
    if(verbose)
      cat("Completed", x$species_id1[1], "vs.", x$species_id2[1],
          " ... retained", nrow(tmp), "BLAST hits\n")
    return(tmp)
  })
  return(rbindlist(outlist))
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
loop_dbs = function(map,
                    radii,
                    mappings,
                    plotit){
  run_dbs = function(map, eps_radius, mappings){
    x=data.frame(map)
    x$x_a = frank(x[,c("chr1","start1")],
                  ties.method = "dense")
    x$x_b = frank(x[,c("chr2","start2")],
                  ties.method = "dense")

    nn = dbscan::frNN(x[,c("x_a","x_b")], eps = eps_radius)
    dbs = dbscan::dbscan(nn, minPts = mappings)
    return(data.frame(rank1 = x$x_a,
                      rank2 = x$x_b,
                      cluster = dbs$cluster,
                      stringsAsFactors = F))
  }
  if(length(radii)!=length(mappings))
    stop("radii and mappings must be same length\n")
  for(i in 1:length(radii)){
    dclus = run_dbs(map, eps = radii[i], mappings = mappings[i])
    mo = cbind(map, data.table(dclus))
    if(plotit){
      with(mo, plot(rank1, rank2,
                    col = ifelse(cluster == 0,rgb(0,0,0,.2),"darkred"),
                    pch = ".",
                    xlab = paste(species_id1[1], "rank"),
                    ylab = paste(species_id2[1], "rank"),
                    main = paste("radius =", radii[i],
                                 "min mapping =", mappings[i])))
    }
    map = map[mo$cluster !=0,]
  }
  return(map)
}


###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######

merge_byOFid = function(blast, gff, of.id, genomes){
  gff1 = merge(of.id, gff, by = "id")
  gff2 = merge(of.id, gff, by = "id")

  setnames(gff1, c(paste0(colnames(gff1),1)))
  setnames(gff1,3,"q.id")

  setnames(gff2, c(paste0(colnames(gff2),2)))
  setnames(gff2,3,"s.id")
  out = merge(gff1,merge(gff2, blast, by = "s.id"), by = "q.id")
  gr = rev(genomes)[-1]
  wh1 = sapply(out$species_id1, function(x) which(genomes == x))
  wh2 = sapply(out$species_id2, function(x) which(genomes == x))
  wh = which(wh1>wh2)

  toswitch = out[wh,]
  tokeep = out[-wh,]
  setnames(toswitch, colnames(tokeep)[c(8:14,1:7,15:18,21:22,19:20,23:24)])
  toswitch = toswitch[,colnames(tokeep), with = F]
  out = rbind(toswitch, tokeep)

  spl = split(out, paste(out$species_id1, out$species_id2))

  lapply(spl, cullblast)
  return(out)
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
parse_gff = function(gff){
  g = suppressWarnings(
    data.table::fread(gff,showProgress = F, verbose = F))
  g = g[g$V3 == "gene",c(9,1,4,5,7)]
  g$V9 = sapply(g$V9, function(x) gsub("Name=","",strsplit(x,";")[[1]][2]))
  data.table::setnames(g, c("id","chr","start","end","strand"))
  return(g)
}

###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ######
parse_blast = function(blast.scheme,
                       min.score = 100,
                       min.percMax = 50,
                       nmapsPerHaplotype = 1,
                       verbose = T,
                       ignoreSelfBlast = T){
  min.propMax = min.percMax/100

  if(ignoreSelfBlast){
    blast.scheme = blast.scheme[blast.scheme[,1]!=blast.scheme[,2],]
  }
  blasts = lapply(1:nrow(blast.scheme), function(i){
    x = blast.scheme[i,5]
    subject.ploidy = as.numeric(blast.scheme[i,4])
    query.ploidy = as.numeric(blast.scheme[i,3])

    suppressWarnings(d <- data.table::fread(x,showProgress = F, verbose = F))
    n = nrow(d)
    setnames(d,c("id1","id2","perc.iden","align.length","n.mismatch",
                 "n.gapOpen","q.start","q.end","s.start","s.end","eval","score"))
    d <- d[d$score >= min.score,]
    d <- d[!duplicated(d[,c("id1","id2"), with = F]),]

    d <- d[order(d$id1, d$id2, -d$score), ]


    d[, rank := frank(score, ties.method = "dense"),
      by = list(id1)]
    d <- d[d$rank <= subject.ploidy*nmapsPerHaplotype,]

    d[, rank := frank(score, ties.method = "dense"),
      by = list(id2)]
    d <- d[d$rank <= query.ploidy*nmapsPerHaplotype,]

    d[, prop := score/max(score),
      by = list(id2)]
    d <- d[d$prop>=min.propMax,]
    d[, prop := score/max(score),
      by = list(id1)]
    d <- d[d$prop>=min.propMax,]

    d$rank = NULL
    d$prop = NULL
    if(verbose)
      cat("Done",blast.scheme[i,1],"vs.",blast.scheme[i,2],": initial hits =",n, "culled to", nrow(d),"\n")
    return(d)
  })
  return(blasts)
}

###### ###### ###### ###### ######
###### ###### ###### ###### ######
align_peptideByDiamond<- function(peptide.files,
                                  db.files,
                                  blast.scheme,
                                  topPerc = 50,
                                  minScore = 100,
                                  ortherDiamondOpts = "--quiet",
                                  nthreads= 1, verbose = T){

  if(verbose)
    cat("Making Diamond databases ...")
  dbcom = lapply(1:length(peptide.files), function(x)
    paste("diamond makedb --in", peptide.files[x],
          "-d", db.files[x], "--quiet",
          "-p", nthreads))
  trsh = lapply(dbcom, system)
  if(verbose)
    cat("Done!\n")

  if(verbose)
    cat("Running Diamond BLASTs ...")
  blcom = lapply(1:nrow(blast.scheme), function(i){
    qid = blast.scheme[i,1]
    sid = blast.scheme[i,2]
    out = blast.scheme[i,5]

    paste("diamond blastp",
          "--query",peptide.files[qid],
          "--db",db.files[sid],
          "--top",topPerc,
          "--min-score",minScore,
          ifelse(is.na(ortherDiamondOpts),"",ortherDiamondOpts),
          "-p",nthreads,
          "--out", out)
  })
  trsh = lapply(blcom, system)
  if(verbose)
    cat("Done!\n")
}



