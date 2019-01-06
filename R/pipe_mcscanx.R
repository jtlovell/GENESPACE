#' @title Run MCScanX program
#'
#' @description
#' \code{pipe_mcscanx} Wrapper to run MCScanX program
#'
#' @param blast blast results data.table
#' @param gff concatenated gff data.table
#' @param mcscan.dir path to the directory where MCScanX will be run
#' @param mcscan.param parameters to supply MCScanX
#' @param silent.mcs logical, should MCScanX be run silently?
#' @param verbose Should updates be printed?
#' @param ... Not currently in use
#' @details ...
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @import data.table
#' @export
pipe_mcscanx <- function(blast,
                         gff,
                         mcscan.dir,
                         mcscan.param,
                         silent.mcs = T,
                         verbose = T){

  #######################################################
  #######################################################
  prep_mcs <- function(blast,
                       gff,
                       mcscan.dir,
                       mcscan.param,
                       silent.mcs){

    if(length(unique(gff$genome)) == 1){
      ga = gff
      gb = gff
      ga$genome <- "a"
      gb$genome <- "b"
      gff <- rbind(gb, ga)
      gff$genome = factor(gff$genome,
                          levels = c("a","b"))
    }else{
      gff$genome = factor(gff$genome, levels = genomeIDs)
    }

    setkey(gff, genome)

    gff[,chr.num := frank(chr, ties.method = "dense"),
        by = genome]
    gff$rank.start = frank(gff, genome, chr, start, ties.method = "random")
    gff$rank.end = frank(gff, genome, chr, end, ties.method = "random")
    gff[,genome.num := frank(genome, ties.method = "dense")]

    lets = paste0(letters,letters)[1:length(unique(gff$genome))]
    gff$genome.abbrev = paste0(lets[gff$genome.num],1)

    gff.in = gff[,c("genome.abbrev","id","rank.start","rank.end")]


    blast.in = blast[,c("id1","id2","perc.iden","align.length",
                        "n.mismatch","n.gapOpen", "q.start", "q.end",
                        "s.start", "s.end","eval","score")]

    print(gff.in)
    print(blast.in)
    write.table(gff.in,
                file = file.path(mcscan.dir,"xyz.gff"),
                row.names = F,
                col.names = F,
                quote = F, sep= "\t")
    write.table(blast.in,
                file = file.path(mcscan.dir,"xyz.blast"),
                row.names = F,
                col.names = F,
                quote = F, sep= "\t")
    if(silent.mcs){
      com <- paste("MCScanX", mcscan.param, file.path(mcscan.dir,"xyz"),"&> /dev/null")
    }else{
      com <- paste("MCScanX", mcscan.param, file.path(mcscan.dir,"xyz"))
    }

    print(com)
    system(com)
    return(file.path(mcscan.dir,"xyz"))
  }
  #######################################################
  #######################################################
  parse_mcs <- function(mcs.file){
    mcscan.raw <- read.delim(paste0(mcs.file, ".collinearity"),
                             sep = "\t", header = F,
                             comment.char = "#", strip.white = T,
                             stringsAsFactors = F)

    fac <- as.numeric(as.factor(sapply(as.character(mcscan.raw$V1), function(x)
      strsplit(x, "-")[[1]][1])))
    genes <- with(mcscan.raw, c(V2, V3))
    out <- data.table(id1 = mcscan.raw$V2,
                      id2 = mcscan.raw$V3,
                      block.id = fac)
    out<-data.table(out)
    setkey(out, id1, id2)
    return(out)
  }
  #######################################################
  #######################################################
  run_mcs <- function(blast,
                      gff,
                      mcscan.dir,
                      mcscan.param,
                      silent.mcs){
    mcs.file = prep_mcs(blast,
                        gff,
                        mcscan.dir,
                        mcscan.param,
                        silent.mcs = silent.mcs)
    mcs.parsed = parse_mcs(mcs.file)
    blast = data.table(blast)
    setkey(blast, id1, id2)
    setkey(mcs.parsed, id1, id2)
    map.out <- data.table(merge(mcs.parsed, blast))
    setkey(map.out, block.id, chr1, start1)
    return(map.out)
  }
  #######################################################
  #######################################################

  #######################################################
  gff.tmp <- gff
  gff.tmp$genome <- paste0(gff.tmp$genome,"xxxx")
  gff.tmp$id <- paste0(gff.tmp$id,"xxxx")
  gff <- data.table(rbind(gff, gff.tmp))

  genomeIDs <- c(genomeIDs, paste0(genomeIDs, "xxxx"))

  bl.dif <- blast[blast$genome1 != blast$genome2,]
  bl.same <- blast[blast$genome1 == blast$genome2,]
  bl.same$id2 <- paste0(bl.same$id2, "xxxx")
  bl.same$genome2 <- paste0(bl.same$genome2, "xxxx")
  blast <- rbind(bl.dif, bl.same)
  #######################################################
  if(verbose)
    cat("Parsing",nrow(blast),"BLAST hits by MCScanX\n\t")

  blast$unique <- with(blast, paste(genome1, genome2))
  if("block.id" %in% colnames(blast))
    blast$block.id<-NULL

  spl = split(blast, "unique")
  #######################################################
  out <- rbindlist(lapply(spl, function(x){
    genomes = c(x$genome1[1],x$genome2[1])

    if(verbose)
      cat(genomes[1],"-->", genomes[2],
          paste0("(initial hits = ",nrow(x),")"))

    gff.x <- gff[gff$genome %in% genomes,]

    tmp <- run_mcs(blast = x,
                  gff = gff.x,
                  mcscan.dir = mcscan.dir,
                  mcscan.param = mcscan.param,
                  silent.mcs = silent.mcs)

    tmp$block.id <- with(tmp, paste0(unique, block.id))

    if(verbose)
      cat(" culled hits =",nrow(tmp),"\n\t")

    return(tmp)
  }))
  #######################################################
  out$block.id <- as.numeric(as.factor(out$block.id))
  out$genome2 <- gsub("xxxx", "",  out$genome2)
  out$id2 <- gsub("xxxx", "",  out$id2)
  #######################################################
  if (verbose)
    cat("Done!\n")
  return(make_blocks(out))
}


