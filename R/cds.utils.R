#' @title selection utility functions
#' @description
#' \code{cds_utils} Five utilities functions meant for internal calls in compareGeneSpace
#' @name syn_utils
#'
#' @param cds.msaFile xxx
#' @param outDir xxx
#' @param GroupID xxx
#' @param round2digits xxx
#' @param base1 xxx
#' @param base2 xxx
#' @param proteinFasta_file xxx
#' @param orthology_file xxx
#' @param orthogroup_columnName xxx
#' @param aligner xxx
#' @param pal2nal_path xxx
#' @param n.cores xxx
#' @param process_nOrthogroups xxx
#' @param clean_tempFiles xxx
#' @param write.output xxx
#' @param outDir xxx
#' @param file_name xxx
#' @param tempDir xxx
#' @param msaFile xxx
#' @param file_prefix xxx
#' @param nSpecies xxx
#' @param clean_tempFiles xxx
#' @param problem_OrthosFile xxx
#' @param cdml_ctlDir xxx
#' @param cdmlDir xxx
#' @param outDir xxx
#' @param cdml_file xxx
#' @param nSpecies xxx
#' @param problem_OrthosFile xxx
#' @param group_columnName xxx
#' @param cds xxx
#' @param protein xxx
#' @param geneIDs xxx
#' @param speciesNames xxx
#' @param file_prefix xxx
#' @param orthoDT xxx
#' @param orthoDT_row xxx
#' @param pal2nal_tool xxx
#' @param orthogroup_columnName xxx
#' @param if.fuzzy.codon xxx
#' @param outDir_cds xxx
#' @param outDir_protein xxx
#' @param outDir_msa xxx
#' @param outDir_codon xxx
#' @param s xxx
#' @param as.character_vector xxx
#' @param remove.improperCDS xxx
#' @param if.fuzzy.codon xxx
#' @param aa_file xxx
#' @param out_msa xxx
#' @param mafft.params xxx
#' @param n.cores xxx
#' @param verbose logical, should updates be printed?
#' @param ... not currently in use
#'
#' @note \code{cds_utils} is a generic name for the functions documented.
#' \cr
#' If called, \code{cds_utils} returns its own arguments.
#'


#' @title Converts MSA file to AXT format
#' @description
#' \code{write_axt} Converts MSA file to AXT format
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readDNAMultipleAlignment DNAStringSet
#' @export
write_axt <- function(cds.msaFile,
                      outDir = NULL){

  cds.msa <- DNAStringSet(
    readDNAMultipleAlignment(cds.msaFile, format = "fasta"))

  axt <- paste(unlist(lapply(seq(cds.msa), function(x) {
    toString(cds.msa[[x]])
    })), collapse = "\n")

  axt <- paste0(paste(names(cds.msa), collapse = "-"),"\n", axt)

  if (is.null(outDir))
    outDir <- dirname(cds.msaFile)

  out_axt <- sprintf("%s/%s.axt", outDir, basename(cds.msaFile))
  writeLines(axt, con = out_axt)
  return(out_axt)
}

#' @title Calculate 4dtv from pair of CDS sequences
#' @description
#' \code{calculate_4dtv} Calculate 4dtv value for a given pair of
#' CDS sequences in MSA format
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readDNAMultipleAlignment DNAStringSet
#' @importFrom IRanges successiveViews
#' @export
calculate_4dtv <- function(cds.msaFile,
                           GroupID = NULL,
                           round2digits = 5,
                           verbose = FALSE){

  fd_AA_codonsDT <- data.table(L = c("CTT", "CTC", "CTA", "CTG"),
                               V = c("GTT", "GTC", "GTA", "GTG"),
                               S = c("TCT", "TCC", "TCA", "TCG"),
                               P = c("CCT", "CCC", "CCA", "CCG"),
                               T = c("ACT", "ACC", "ACA", "ACG"),
                               A = c("GCT", "GCC", "GCA", "GCG"),
                               R = c("CGT", "CGC", "CGA", "CGG"),
                               G = c("GGT", "GGC", "GGA", "GGG"))

  fd_AA_codonsDT <- melt(t(fd_AA_codonsDT))
  fd_AA_codonsDT$Var2 = NULL
  fd_AA_codonsDT <- data.table(fd_AA_codonsDT)
  setnames(fd_AA_codonsDT, c("Var1", "value"), c("AA", "CDS"))


  if (!file.exists(cds.msaFile))
    stop("** cds.msaFile not found **")

  tryCatch({
    cds.msa <- DNAStringSet(
      readDNAMultipleAlignment(cds.msaFile, format = "fasta"))

    cod <- data.table(
      data.frame(successiveViews(cds.msa[[1]],
                                 rep.int(3L, length(cds.msa[[1]]) %/% 3L))),
      data.frame(successiveViews(cds.msa[[2]],
                                 rep.int(3L, length(cds.msa[[2]]) %/% 3L))))
    setnames(cod, c("cds1", "cds2"))

    cod$base1 <- substring(cod$cds1, 3, 3)
    cod$base2 <- substring(cod$cds2, 3, 3)

    FD <- subset(cod,
                 cds1 %in% fd_AA_codonsDT$CDS &
                   cds2 %in% fd_AA_codonsDT$CDS)
    FD <- merge(FD, fd_AA_codonsDT,
                by.x = "cds1", by.y = "CDS")
    setnames(FD, "AA", "AA1")

    FD <- merge(FD, fd_AA_codonsDT,
                by.x = "cds2", by.y = "CDS")
    setnames(FD, "AA", "AA2")
    FD <- subset(FD, AA1 == AA2)

    FD$Tv <- is.transversion(toupper(FD$base1), toupper(FD$base2))
    FD$Ts <- is.transition(toupper(FD$base1), toupper(FD$base2))

    FDsites <- nrow(FD)

    Tv <- nrow(FD[FD$Tv == TRUE,])
    Ts <- nrow(FD[FD$Ts == TRUE,])

    FDTV <- round(Tv/FDsites, round2digits)
    FDTS <- round(Ts/FDsites, round2digits)

    ## correct raw 4dtv by ~HKY substitution model
    FDTVc <- ifelse(FDTV < 0.5, (-0.5) * log(1 - 2 * FDTV), 99.0)
    FDTVc <- round(FDTVc, round2digits)
    if (is.null(GroupID)) {
      return(data.table(FDsites = FDsites,
                        Ts = Ts,
                        Tv = Tv,
                        FDTS = FDTS,
                        FDTV = FDTV,
                        FDTVc = FDTVc))
    }else{
      return(data.table(GroupID = GroupID,
                        FDsites = FDsites,
                        Ts = Ts,
                        Tv = Tv,
                        FDTS = FDTS,
                        FDTV = FDTV,
                        FDTVc = FDTVc))
    }

  }, warning = function(w){
    if (verbose)
      print(paste("WARNING: ", w))
  }, error = function(e){
    if (verbose)
      message(sprintf("\nERROR: %s\tGroupID = %s", e, GroupID))
    if (is.null(GroupID)) {
      return(data.table(FDsites = NA,
                        Ts = NA,
                        Tv = NA,
                        FDTS = NA,
                        FDTV = NA,
                        FDTVc = NA))
    }else{
      return(data.table(GroupID = GroupID,
                        FDsites = NA,
                        Ts = NA,
                        Tv = NA,
                        FDTS = NA,
                        FDTV = NA,
                        FDTVc = NA))
    }

  })
}

#' @title Check for base transversion
#' @description
#' \code{is.transversion} Check for base transversion
#' @rdname cds_utils
#' @import data.table
#' @export
is.transversion <- function(base1,
                            base2){

  transversion <- list(A = c("T","C"),
                       C = c("A","G"),
                       G = c("T","C"),
                       T = c("A","G"))

  if (length(base1) != length(base2) || length(base1) == 0L)
    stop("base1 and base2 must be character vectors of equal length")
  out <- sapply(seq(base1), function(x){
    return(base1[x] %in% transversion[[base2[x]]])
  })
  return(out)
}

#' @title Check for base transition
#' @description
#' \code{is.transition} Check for base transition
#' @rdname cds_utils
#' @import data.table
#' @export
is.transition <- function(base1,
                          base2){

  transition <- list(A = c("G"),
                     C = c("T"),
                     G = c("A"),
                     T = c("C"))

  if (length(base1) != length(base2) || length(base1) == 0L)
    stop("base1 and base2 must be character vectors of equal length")
  out <- sapply(seq(base1), function(x){
    return(base1[x] %in% transition[[base2[x]]])
  })
  return(out)
}

#' @title Pipline to calculate 4DTV values
#' @description
#' \code{pipe_4dtv} Pipline to calculate 4DTV values for paralogs/orthologs
#'  between species in given orthology file
#' @rdname cds_utils
#' @import data.table
#' @export
pipe_4dtv <- function(cdsFasta_file,
                      proteinFasta_file = NULL,
                      orthology_file = NULL,
                      orthogroup_columnName = "V1",
                      aligner = "mafft",
                      pal2nal_path = NULL,
                      n.cores = 6,
                      process_nOrthogroups = NULL,
                      clean_tempFiles = TRUE,
                      write.output = TRUE,
                      outDir = ".",
                      file_name = basename(tempfile()),
                      tempDir = file.path(outDir, basename(tempfile())),
                      verbose = FALSE){

  aligner <- match.arg(aligner, choices = "mafft")
  ## -- check programs
  if (Sys.which(aligner) == "")
    stop(sprintf("%s - not found", aligner))
  if (!file.exists(pal2nal_tool <- file.path(pal2nal_path,"pal2nal.pl")))
    stop("** pal2nal not found **")

  if (file.exists(orthology_file)){
    orthoDT <- fread(orthology_file, header = T)
  } else
    stop("** orthology_file must be a data.frame/data.table object OR path to orthology file **")

  if (ncol(orthoDT) >3L)
    stop("** ortholog_file must have 2 columns each representing one of the homolog pairs \n\t and an optional column representing group ID**")
  if (ncol(orthoDT) == 2L && !all(orthogroup_columnName %in% colnames(orthoDT)) ){
    orthoDT[[orthogroup_columnName]] <- paste(orthoDT[[1]], orthoDT[[2]], sep = "__")
  }

  if (!is.null(process_nOrthogroups) && is.numeric(process_nOrthogroups)){
    if(process_nOrthogroups > nrow(orthoDT))
      process_nOrthogroups = nrow(orthoDT)
    orthoDT <- orthoDT[c(1:process_nOrthogroups), ]
  }

  if (nrow(orthoDT) >= 1L){
    all.genes <- unique(unlist(lapply(orthoDT[, -orthogroup_columnName, with = F], "["), use.names = F))
  } else
    stop("** ortholog_file is empty **")

  ## -- create temp directories
  suppressWarnings(dir.create(tempDir))
  suppressWarnings(dir.create(cdsDir <- file.path(tempDir, "cds")))
  suppressWarnings(dir.create(pepDir <- file.path(tempDir, "pep")))
  suppressWarnings(dir.create(msaDir <- file.path(tempDir, "msa")))
  suppressWarnings(dir.create(codonDir <- file.path(tempDir, "codon")))
  suppressWarnings(dir.create(fdtvDir <- file.path(tempDir, "fdtv")))

  cds = protein = NULL
  if(class(cdsFasta_file) != "DNAStringSet"){
    if(file.exists(cdsFasta_file))
      cds <- readDNAStringSet(cdsFasta_file) else
        stop("** cdsFasta_file must be a Biostrings::DNAStringSet object or path to CDS fasta file **")
  } else
    cds <- cdsFasta_file

  if(!all(all.genes %in% names(cds)))
    stop("** Not all GeneIDs in orthology_file were not found in cdsFasta_file **")

  cds <- cds[match(all.genes, names(cds))]
  if(verbose)
    cat("Calculating 4DTV values for ", length(cds)/2, " pairs.\n")

  if(!is.null(proteinFasta_file)){
    if(class(proteinFasta_file) != "AAStringSet"){
      if(file.exists(proteinFasta_file))
        protein <- readAAStringSet(proteinFasta_file) else
          stop("** proteinFasta must be a Biostrings::AAStringSet object or path to AA fasta file **")
    } else
      protein <- proteinFasta_file
  } else{
    if(verbose)
      cat("Translating CDS sequences .. ")
    protein <- translate_cds(cds,
                            as.character_vector = F,
                            remove.improperCDS = F,
                            if.fuzzy.codon = "solve",
                            n.cores = n.cores)
    if(verbose)
      cat("DONE !\n")
  }

  if(!all(all.genes %in% names(protein)))
    stop("** Not all GeneIDs in orthology_file were not found in proteinFasta_file **")

  protein <- protein[match(all.genes, names(protein))]

  fdtvDT <- rbindlist(mclapply(seq(nrow(orthoDT)), mc.cores = n.cores, function(i){
    if(i %% 200 == 0)
      cat(".")

    file_prefix = as.character(orthoDT[i, orthogroup_columnName, with = F])
    geneIDs <- as.character(orthoDT[i, -orthogroup_columnName, with = F])
    speciesNames <- colnames(orthoDT[i, -orthogroup_columnName, with = F])

    alignFiles <- pipe_alignCDS(cds = cds,
                                protein = protein,
                                file_prefix = file_prefix,
                                geneIDs = geneIDs,
                                speciesNames = speciesNames,
                                orthoDT = NULL,
                                orthoDT_row = NULL,
                                orthogroup_columnName = orthogroup_columnName,
                                pal2nal_tool = pal2nal_tool,
                                if.fuzzy.codon = "solve",
                                outDir_cds = cdsDir,
                                outDir_protein = pepDir,
                                outDir_msa = msaDir,
                                outDir_codon = codonDir,
                                n.cores = 1,
                                verbose = FALSE)


    if(any(is.null(alignFiles)) | any(is.na(alignFiles)) | any(alignFiles == ""))
      stop("** Making alignment files from CDS failed **")

    fdtvInfo <- calculate_4dtv(alignFiles$out_codon,
                               GroupID = file_prefix,
                               round2digits = 5)

    tempFiles <- unlist(alignFiles, use.names = F)

    ## -- clean temporary files
    if(clean_tempFiles){
      unlink(tempFiles)
    }
    if(!is.null(fdtvInfo))
      return(fdtvInfo)
  }), fill = TRUE)

  ## -- write 4dtv output
  if(write.output)
    write.csv(fdtvDT, file = sprintf("%s/%s_4dtv.csv", outDir, file_name))

  ## -- delete the temporary directory
  if(clean_tempFiles){
    unlink(tempDir, recursive = T, force = T)
  }
  return(fdtvDT)
}

#' @title calculate kaks
#' @description
#' \code{calculate_kaks} calculate kaks
#' @rdname cds_utils
#' @import data.table
#' @export
calculate_kaks <- function(msaFile,
                           file_prefix,
                           nSpecies,
                           clean_tempFiles = FALSE,
                           problem_OrthosFile = NULL,
                           cdml_ctlDir = ".",
                           cdmlDir = ".",
                           outDir = "."){
  curDir = getwd()
  cdml.ctlFile <- sprintf("%s/%s.cdml.ctl", cdml_ctlDir, file_prefix)
  out_cdml <- sprintf("%s/%s.cdml", cdmlDir, file_prefix)
  make_codemlCtl(input.msa = msaFile,
                 output.cdml = out_cdml,
                 out_cdmlCtl = cdml.ctlFile)
  run_cmd <- paste0("codeml ", basename(cdml.ctlFile))

  tryCatch({
    setwd(cdml_ctlDir)
    rCdml <- system(run_cmd, intern=TRUE)

    kaks <- parse_Codeml(cdml_file = out_cdml,
                         nSpecies = nSpecies,
                         problem_OrthosFile = problem_OrthosFile,
                         outDir = outDir,
                         verbose = F)
  }, warning = function(w){
    print(paste("WARNING:  ",w))
  }, error = function(e){
    setwd(curDir)
    cat(paste("\nERROR:  ",e))
    stop("\n** codeml run failed **\n")
  } )
  if(clean_tempFiles){
    unlink(c(out_cdml, cdml.ctlFile))
    cdml.ctlFile = out_cdml = NULL
  }
  return(list(kaks = kaks, out_cdml = out_cdml, cdml.ctlFile = cdml.ctlFile))
}

#' @title Parse codeml output
#' @description
#' \code{parse_Codeml} Parse codeml output to obtain dN/dS and other values
#' @rdname cds_utils
#' @import data.table
#' @export
parse_Codeml <- function(cdml_file,
                         nSpecies = NULL,
                         problem_OrthosFile = NULL,
                         group_columnName = "og",
                         outDir = ".",
                         verbose = F){

  ## -- helper function
  rm_spaces<-function(s){
    s<-gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", s, perl=TRUE)
    return(s)
  }
  if(is.null(nSpecies))
    stop("*** Must provide number of pairs analyzed through paml/codeml ***")

  if(is.null(problem_OrthosFile))
    problem_OrthosFile <- sprintf("%s/%s", outDir, "problematic_Orthologs.kaks.txt")

  codeml_df = NULL

  cdmlInfo <- readLines(cdml_file)
  if(nSpecies==2){
    s <- cdmlInfo[length(cdmlInfo)]
    ## could use this - but errors out when encounters missing/integer values
    ## Ka <- as.numeric(regmatches(s, regexec("\\sdN\\s*=\\s*(\\d+\\.\\d+)", s))[[1]][[2]])
    s <- rm_spaces(s)
    s <- gsub("= ","=", s)
    s <- gsub(" =","=", s)
    s <- unlist(strsplit(s,"[= ]"))
    codeml_df <- data.table(t(as.numeric(s[c(1:length(s))%%2==0])))
    setnames(codeml_df, s[c(1:length(s))%%2!=0])
    setnames(codeml_df, old = c("dN/dS", "dN", "dS"), new = c("KaKs", "Ka", "Ks"))
    codeml_df[[group_columnName]] <- sub("\\.cdml", "", basename(cdml_file))
    out_columns <- c(group_columnName, "Ka", "Ks", "KaKs")
    codeml_df <- codeml_df[, ..out_columns]
  }
  if(nSpecies > 2){
    combs <- grep("\\.\\.\\.", cdmlInfo)
    if(length(combs)!=length(combn(nSpecies, 2, simplify=F))){
      write(basename(cdml_file),file = problem_OrthosFile, append = T, ncolumns = 1)
      if(verbose)
        message("\nProblem with this file:\n", cdml_file)
    } else {
      for(i in combs){
        locs <- gregexpr("[()]", cdmlInfo[i], perl = T)
        id1 <- substr(cdmlInfo[i],locs[[1]][1]+1, locs[[1]][2]-1)
        id2 <- substr(cdmlInfo[i],locs[[1]][3]+1, locs[[1]][4]-1)
        gComparison <- ifelse(id1>id2, paste0(id1, ".vs.", id2), paste0(id2, ".vs.", id1))
        s <- cdmlInfo[i+4]
        s <- rm_spaces(s)
        s <- gsub("= ","=", s)
        s <- gsub(" =","=", s)
        s <- unlist(strsplit(s,"[= ]"))
        tmp <- data.table(t(as.numeric(s[c(1:length(s))%%2==0])))
        setnames(tmp, s[c(1:length(s))%%2!=0])
        setnames(tmp, old = c("dN/dS", "dN", "dS"), new = c("KaKs", "Ka", "Ks"))
        tmp <- tmp[, c("Ka", "Ks", "KaKs"), with = T]
        setnames(tmp, paste(gComparison, colnames(tmp), sep="@"))
        if(is.null(codeml_df))
          codeml_df <- tmp else
            codeml_df <- cbind(codeml_df, tmp)
      }
      codeml_df[[group_columnName]] <- sub("\\.cdml","", basename(cdml_file))
      out_columns <- c(group_columnName, colnames(codeml_df)[!colnames(codeml_df) %in% group_columnName])
      codeml_df <- codeml_df[, ..out_columns]
    }
  }
  return(codeml_df)
}

#' @title Generate codeml input control file
#' @description
#' \code{make_codemlCtl} Generate codeml input control file
#' @rdname cds_utils
#' @import data.table
#' @export
make_codemlCtl <- function(input.msa,
                           output.cdml,
                           out_cdmlCtl){
  codeml_template.ctl <- "seqfile = input.msa\n    outfile = output.cdml\n       noisy = 0\n     verbose = 0\n     runmode = -2\n   cleandata = 1\n     seqtype = 1\n   CodonFreq = 2\n       model = 2\n     NSsites = 0\n       icode = 0\n       Mgene = 0\n   fix_kappa = 0\n       kappa = 2\n   fix_omega = 0\n       omega = 1\n   fix_alpha = 1\n       alpha = .0\n      Malpha = 0\n       ncatG = 4\n       clock = 0\n       getSE = 0\nRateAncestor = 0\n      method = 0"
  codeml_template.ctl <- sub("input.msa", input.msa, codeml_template.ctl)
  codeml_template.ctl <- sub("output.cdml", output.cdml, codeml_template.ctl)
  writeLines(codeml_template.ctl, con = out_cdmlCtl)
}

#' @title Pipline to calculate KaKS values
#' @description
#' \code{pipe_kaks} Pipline to calculate KaKS values for
#' orthologs between species in given orthology file
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet translate
#' @export
pipe_kaks <- function(cdsFasta_file,
                      proteinFasta_file = NULL,
                      orthology_file = NULL,
                      orthogroup_columnName = "V1",
                      aligner = "mafft",
                      pal2nal_path = NULL,
                      calculate_4dtv = FALSE,
                      n.cores = 6,
                      process_nOrthogroups = NULL,
                      clean_tempFiles = TRUE,
                      write.output = TRUE,
                      outDir = ".",
                      file_name = basename(tempfile()),
                      tempDir = file.path(outDir, basename(tempfile()))){

  aligner <- match.arg(aligner, choices = "mafft")
  ## -- check programs
  if(Sys.which(aligner) == "")
    stop(sprintf("%s - not found", aligner))
  if(Sys.which("codeml") == "")
    stop("codeml - not found")
  if(!file.exists(pal2nal_tool <- file.path(pal2nal_path,"pal2nal.pl")))
    stop("** pal2nal not found **")

  curDir <- getwd()

  if(file.exists(orthology_file)){
    orthoDT <- fread(orthology_file, header = T)
  } else
    stop("** orthology_file must be a data.frame/data.table object OR path to orthology file **")

  if(!all(orthogroup_columnName %in% colnames(orthoDT)))
    stop("** orthogroup_columnName must be a column in orthoDT **")

  if(!is.null(process_nOrthogroups) && is.numeric(process_nOrthogroups) && process_nOrthogroups >= 1L){
    orthoDT <- orthoDT[c(1:process_nOrthogroups), ]
  }

  if(nrow(orthoDT) >= 1L){
    all.genes <- unique(unlist(lapply(orthoDT[, -orthogroup_columnName, with = F], "["), use.names = F))
  } else
    stop("** ortholog_file is empty **")
  ## -- create temp directories
  suppressWarnings(dir.create(tempDir))
  suppressWarnings(dir.create(cdsDir <- file.path(tempDir, "cds")))
  suppressWarnings(dir.create(pepDir <- file.path(tempDir, "pep")))
  suppressWarnings(dir.create(msaDir <- file.path(tempDir, "msa")))
  suppressWarnings(dir.create(codonDir <- file.path(tempDir, "codon")))
  suppressWarnings(dir.create(ctlDir <- file.path(tempDir, "ctl")))
  suppressWarnings(dir.create(cdmlDir <- file.path(tempDir, "cdml")))
  problem_OrthosFile <- sprintf("%s/problematic_Orthologs.kaks.txt", outDir)

  cds = protein = NULL
  if(class(cdsFasta_file) != "DNAStringSet"){
    if(file.exists(cdsFasta_file))
      cds <- readDNAStringSet(cdsFasta_file) else
        stop("** cdsFasta_file must be a Biostrings::DNAStringSet object or path to CDS fasta file **")
  } else{
    cds <- cdsFasta_file
  }
  if(!all(all.genes %in% names(cds)))
    stop("** Not all GeneIDs in orthology_file were not found in cdsFasta_file **")

  cds <- cds[match(all.genes, names(cds))]

  if(!is.null(proteinFasta_file)){
    if(class(proteinFasta_file) != "AAStringSet"){
      if(file.exists(proteinFasta_file))
        protein <- readAAStringSet(proteinFasta_file) else
          stop("** proteinFasta must be a Biostrings::AAStringSet object or path to AA fasta file **")
    } else{
      protein <- proteinFasta_file
    }
  } else{
    protein <-translate_cds(cds,
                            as.character_vector = F,
                            remove.improperCDS = F,
                            if.fuzzy.codon = "solve",
                            n.cores = n.cores)
  }
  if(!all(all.genes %in% names(protein)))
    stop("** Not all GeneIDs in orthology_file were not found in proteinFasta_file **")

  protein <- protein[match(all.genes, names(protein))]


  kaksDT <- rbindlist(mclapply(seq(nrow(orthoDT)), mc.cores = n.cores, function(i){
    if(i %% 200 == 0)
      cat(".")

    file_prefix = as.character(orthoDT[i, orthogroup_columnName, with = F])
    geneIDs <- as.character(orthoDT[i, -orthogroup_columnName, with = F])
    speciesNames <- colnames(orthoDT[i, -orthogroup_columnName, with = F])

    alignFiles <- pipe_alignCDS(cds = cds,
                                protein = protein,
                                file_prefix = file_prefix,
                                geneIDs = geneIDs,
                                speciesNames = speciesNames,
                                orthoDT = NULL,
                                orthoDT_row = NULL,
                                orthogroup_columnName = orthogroup_columnName,
                                pal2nal_tool = pal2nal_tool,
                                if.fuzzy.codon = "solve",
                                outDir_cds = cdsDir,
                                outDir_protein = pepDir,
                                outDir_msa = msaDir,
                                outDir_codon = codonDir,
                                n.cores = 1,
                                verbose = FALSE)

    tempFiles <- unlist(alignFiles, use.names = F)
    if(any(is.null(alignFiles)) | any(is.na(alignFiles)) | any(alignFiles == ""))
      stop("** Making alignment files from CDS failed **")

    fdtvInfo = NULL
    if(calculate_4dtv && length(speciesNames) == 2){
      suppressWarnings(dir.create(fdtvDir <- file.path(tempDir, "fdtv")))
      fdtvInfo <- calculate_4dtv(alignFiles$out_codon,
                                 GroupID = file_prefix,
                                 round2digits = 5)
    }

    ## -- kaks calculation
    ks <- calculate_kaks(msaFile = alignFiles$out_codon,
                         file_prefix = file_prefix,
                         nSpecies = length(speciesNames),
                         clean_tempFiles = clean_tempFiles,
                         problem_OrthosFile = problem_OrthosFile,
                         cdml_ctlDir = ctlDir,
                         cdmlDir = cdmlDir,
                         outDir = outDir)
    tempFiles <- c(tempFiles, ks$cdml.ctlFile, ks$out_cdml)
    kaks <- ks$kaks
    ## -- clean temporary files
    if(clean_tempFiles){
      unlink(tempFiles)
    }
    if(nrow(kaks) == 1L){
      if(!is.null(fdtvInfo)){
        kaks <- cbind(kaks, round(fdtvInfo[, c("FDTVc"), with = F], 5))
      }
      return(kaks)
    }
  }))
  setwd(curDir)

  ## -- delete the temporary directory
  if(clean_tempFiles)
    unlink(tempDir, recursive = T, force = T)

  ## -- write kaks output
  if(write.output)
    write.csv(kaksDT, file = sprintf("%s/%s_kaks.csv", outDir, file_name))
  return(kaksDT)
}

#' @title Aligns amino acid sequences
#' @description
#' \code{alignAA} Aligns amino acid sequences using multiple sequence aligner, mafft
#' @rdname cds_utils
#' @import data.table
#' @export
align_pep <- function(pep.fa,
                      out_msa,
                      mafft.params = "--retree 1 --quiet"){
  cmd <- sprintf("mafft %s --clustalout %s > %s",
                 mafft.params,
                 aa_file,
                 out_msa)
  msatmp <- system(cmd)
}

#' @title Translates CDS sequences
#' @description
#' \code{translate_cds} Translates CDS sequences using Biostrings::translate
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings translate AAStringSet
#' @export
translate_cds <- function(s,
                          as.character_vector = TRUE,
                          remove.improperCDS = TRUE,
                          if.fuzzy.codon = "solve"){

  aa <- lapply(seq(s), function(x) tryCatch({
    if(as.character_vector)
      as.character(translate(s[[x]], genetic.code = GENETIC_CODE, no.init.codon=FALSE,
                                         if.fuzzy.codon = if.fuzzy.codon)) else
                                           translate(s[[x]], genetic.code = GENETIC_CODE, no.init.codon=FALSE,
                                                                 if.fuzzy.codon = if.fuzzy.codon)
  }, error = function(e) {
    if(remove.improperCDS)
      return(NA) else {
        stop("Could not translate given coding sequences.")
      }
  }, warning = function(w) {
    if(remove.improperCDS)
      return(NA) else {
        stop("Could not translate given coding sequences.")
      }
  }) )
  if(as.character_vector)
    return(unlist(aa)) else {
      aa <- AAStringSet(aa)
      names(aa) <- names(s)
      return(aa)
    }
}

#' @title Prepares CDS and AA fasta for alignment
#' @description
#' \code{pipe_alignCDS} Prepares CDS and AA fasta for alignment
#' @rdname cds_utils
#' @import data.table
#' @importFrom Biostrings readAAStringSet readDNAStringSet writeXStringSet translate
#' @export
pipe_alignCDS <- function(cds,
                          protein = NULL,
                          geneIDs = NULL,
                          speciesNames = NULL,
                          file_prefix = NULL,
                          orthoDT = NULL,
                          orthoDT_row = 1,
                          pal2nal_tool = NULL,
                          orthogroup_columnName = "og",
                          if.fuzzy.codon = "solve",
                          outDir_cds = ".",
                          outDir_protein = ".",
                          outDir_msa = ".",
                          outDir_codon = ".",
                          n.cores = 6,
                          verbose = FALSE){

  if(is.null(pal2nal_tool) || !file.exists(pal2nal_tool))
    stop("** path to pal2nal perl script is missing **")

  if(class(cds) != "DNAStringSet")
    if(file.exists(cds))
      cds <- readDNAStringSet(cds) else
        stop("** cds must be a Biostrings::DNAStringSet object or path to CDS fasta file **")

    if(!is.null(protein)){
      if(class(protein) != "AAStringSet")
        if(file.exists(protein))
          protein <- readAAStringSet(protein) else
            stop("** proteinFasta must be a Biostrings::AAStringSet object or path to AA fasta file **")
    }
    if(any(is.null(geneIDs) || is.null(speciesNames) || is.null(file_prefix))){
      if(is.null(orthoDT) || !is.data.frame(orthoDT))
        stop("** orthoDT must be a data.frame or data.table object or path to orthology file **")

      if(!all(orthogroup_columnName %in% colnames(orthoDT)))
        stop("** orthogroup_columnName must be a column in orthoDT **")
      if(nrow(orthoDT) < orthoDT_row)
        stop(sprintf("** orthoDT must be of at least %d rows **", orthoDT_row))

      file_prefix = as.character(orthoDT[orthoDT_row, orthogroup_columnName, with = F])
      geneIDs <- as.character(orthoDT[orthoDT_row, -orthogroup_columnName, with = F])
      speciesNames <- colnames(orthoDT[orthoDT_row, -orthogroup_columnName, with = F])
    }

    cds_file <- sprintf("%s/%s.cds.fa", outDir_cds, file_prefix)
    sub_cds <- cds[match(geneIDs,names(cds))]
    names(sub_cds) <- speciesNames
    writeXStringSet(sub_cds, cds_file, format="fasta", append = F)

    pep_file <- sprintf("%s/%s.pep.fa", outDir_protein, file_prefix)
    if(is.null(protein)){
      sub_protein <- translate(sub_cds,
                                           if.fuzzy.codon = if.fuzzy.codon)
    } else {
      sub_protein <- protein[match(geneIDs,names(protein))]
      names(sub_protein) <- speciesNames
    }
    writeXStringSet(sub_protein, pep_file, format="fasta", append = F)

    out_msa <- sprintf("%s/%s.%s", outDir_msa, file_prefix, "clu")
    alignAA(pep_file,
            out_msa,
            mafft.params = "--retree 1 --quiet")

    out_codon <- sprintf("%s/%s.aln", outDir_codon, file_prefix)
    paltmp <- system(sprintf("perl %s %s %s -output fasta -nogap > %s", pal2nal_tool, out_msa, cds_file, out_codon))

    return(list(cds_file = cds_file, pep_file = pep_file, out_msa = out_msa, out_codon = out_codon))
}
