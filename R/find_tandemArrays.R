#' @title Run MCScanX's detect_collinear_tandem_arrays
#'
#' @description
#' \code{find_pairwiseTandemArrays} Run MCScanX's detect_collinear_tandem_arrays
#'
#'
#' @param blast the blast dataset to screen for syntenic hits
#' @param gff gff data.table
#' @param mcscan.dir Location to write files / run MCScanX
#' @param orthogroup.dt Data.table, produced by upstream functions
#' @param DCTA.path path to detect_collinear_tandem_arrays
#' @param verbose logical, should updates be printed?
#' @param ... Not currently in use
#'
#' @details Internal function
#'
#' @return A culled b.last dataset
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
find_pairwiseTandemArrays <- function(gff,
                                      blast,
                                      mcscan.dir,
                                      orthogroup.dt,
                                      DCTA.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX/downstream_analyses/detect_collinear_tandem_arrays",
                                      verbose = T,
                                      ...){
  #######################################################
  #######################################################
  prep_gff4mcscan <- function(gff, genomeIDs){
    gff <- gff[gff$genome %in% genomeIDs,]
    gff[,chr.num := frank(chr,
                          ties.method = "dense"),
            by = genome]
    gff$genome <- factor(gff$genome,
                         levels = genomeIDs)
    gff$rank.start = frank(gff, genome, chr, start,
                           ties.method = "dense")
    gff$rank.end = frank(gff, genome, chr, end,
                         ties.method = "dense")
    gff[,genome.num := frank(genome,
                             ties.method = "dense")]
    lets = paste0(letters,letters)[1:length(unique(gff$genome))]
    gff$genome.abbrev = paste0(lets[gff$genome.num],gff$chr.num)
    gff = gff[,c("genome.abbrev","id","start","end")]
    return(gff)
  }
  #######################################################
  #######################################################
  prep_blast4mcscan <- function(blast, genomeIDs){
    blast <- blast[with(blast,
                        genome2 %in% genomeIDs &
                          genome1 %in% genomeIDs),]
    blast <- blast[,c("id1","id2","perc.iden","align.length",
                            "n.mismatch","n.gapOpen", "q.start", "q.end",
                            "s.start", "s.end","eval","score")]
    return(blast)
  }
  #######################################################
  #######################################################
  propagate_connectedDCTAs <- function(dcta.dt){
    aids2check = unique(dcta.dt$array.id[duplicated(dcta.dt$og) |
                                           duplicated(dcta.dt$og, fromLast = T)])
    while(length(aids2check)>0){
      x = aids2check[1]
      tid = dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% x]]
      tod = vector()
      while(length(tod)<length(tid)){
        tod = tid
        tid = dcta.dt$array.id[dcta.dt$id %in% dcta.dt$id[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)]<-tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  propagate_connectedOGs <- function(dcta.dt){
    ud = dcta.dt[!duplicated(dcta.dt[,c("og","array.id")]),]
    dup = duplicated(ud$og)

    aids2check <- ud$array.id[which(dup)]

    while(length(aids2check)>0){
      x = aids2check[1]
      tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% x]]
      tod = vector()
      while(length(tod)<length(tid)){
        tod = tid
        tid = dcta.dt$array.id[dcta.dt$og %in% dcta.dt$og[dcta.dt$array.id %in% tid]]
      }
      dcta.dt$array.id[dcta.dt$array.id %in% unique(tid)]<-tid[1]
      aids2check <- aids2check[!aids2check %in% tid]
    }
    dcta.dt <- dcta.dt[!duplicated(dcta.dt),]
    return(dcta.dt)
  }
  #######################################################
  #######################################################
  load_DCTAs <- function(dtca.output, orthogroup.dt){
    raw = fread(dtca.output, header = T)
    r2 = apply(raw, 1, function(x)
      strsplit(
        paste(x[1], x[3],
              x[2], x[4],
              sep = ","),
        ",")[[1]])
    for(i in 1:length(r2))
      r2[[i]]<-data.table(id = r2[[i]],
                          array.id = i)
    r2 = rbindlist(r2)
    r2 <- r2[!duplicated(r2),]
    r3 = propagate_connectedDCTAs(r2)
    r3 <- r3[!duplicated(r3),]
    r3 = merge(orthogroup.dt, r3, by = "id")
    r4 = propagate_connectedOGs(r3)
    r4 <- r4[!duplicated(r4),]
    r4$array.id <- frank(r4$array.id, ties.method = "dense")
    return(r4)
  }
  #######################################################
  #######################################################
  if (verbose)
    cat("Running tandem array searches for each pair of genomes ... \n")
  genome.cmbns = blast[,c("genome1","genome2")]
  genome.cmbns = genome.cmbns[!duplicated(genome.cmbns),]
  genome.cmbns = genome.cmbns[genome.cmbns$genome1 != genome.cmbns$genome2,]

  pw.cta = rbindlist(apply(genome.cmbns, 1, function(x){
    if (verbose)
      cat(x[1],"<-->",x[2], "... ")
    gff.in <- prep_gff4mcscan(gff = gff, genomeIDs = x)
    blast.in <- prep_blast4mcscan(blast = blast, genomeIDs = x)

    gff.file = file.path(mcscan.dir,"xyz.gff")
    blast.file = file.path(mcscan.dir,"xyz.blast")
    collin.file = file.path(mcscan.dir,"xyz.collinearity")
    dcta.file = file.path(mcscan.dir,"xyz.TCA.txt")

    write.table(gff.in,
                file = gff.file,
                row.names = F,
                col.names = F,
                quote = F, sep= "\t")
    write.table(blast.in,
                file = blast.file,
                row.names = F,
                col.names = F,
                quote = F, sep= "\t")
    if (verbose)
      cat("running MCScanX ... ")

    system(paste("MCScanX -a", file.path(mcscan.dir,"xyz"), "&> /dev/null"))

    if (verbose)
      cat("running detect_collinear_tandem_arrays ... \n\t")

    com = paste(dcta.path,
                "-b", blast.file,
                "-g", gff.file,
                "-c", collin.file,
                "-o", dcta.file, "&> /dev/null")
    system(com)

    dtca <- load_DCTAs(dcta.file,
                       orthogroup.dt = orthogroup.dt[,c("id","og")])
    dtca$unique.genome <- paste(x[1], x[2])
    if (verbose)
      cat("Found",
          length(unique(dtca$array.id)),"tandem arrays across",
          nrow(dtca), "genes\n")
    return(dtca)
  }))

  blast$unique.genome <- with(blast, paste(genome1, genome2))
  pw.cta$og1 <- pw.cta$og
  pw.cta$tandemarray.id <- with(pw.cta, paste(og,unique.genome,array.id))
  pw.cta <- pw.cta[,c("tandemarray.id","og1","unique.genome")]
  pw.cta <- pw.cta[!duplicated(pw.cta),]
  setkey(pw.cta, og1, unique.genome)
  setkey(blast, og1, unique.genome)
  out = merge(blast, pw.cta, all = T)
  if (verbose)
    cat("\tDone!\n")
  return(out)
}
