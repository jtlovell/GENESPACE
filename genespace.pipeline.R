#
#
# pipe_genespace(parameter.list = NULL){
#   if(!is.null(parameter.list)){
#     attach(parameter.list)
#   }
#   cat("Running the GENESPACE pipeline on pre-computed blast results\n")
#   if(!exists("wd"))
#     stop("Must specify working.directory (wd)\n")
#
#
#   cat("Setting working directory as", wd,"\n")
#   setwd(wd)
#
#   cat("Using the following parameters:\n")
#   #############################################################
#   if(!exists("MCScanX.path"))
#     stop("Need to specify MCScanX.path\n")
#   if(!file.exists(file.path(MCScanX.path,"MCScanX"))
#     stop("cannot find MCScanX program in", MCScanX.path)
#   cat("\t1. Path to MCScanX exectuable:",MCScanX.path,"\n")
#   #############################################################
#   if(!exists("genomeIDs"))
#     stop("Must specify genomeIDs as a character vector")
#   if(length(genomeIDs) <= 1)
#     stop("Must specify > 1 genome\n")
#   if(!all(genomeIDs %in% dir(file.path(wd, "raw_annotations"))))
#      stop("Cannot find all genomeIDs in",file.path(wd, "raw_annotations"))
#   cat("\t2. Genomes to test:",paste(genomeIDs, collapse = ", "),"\n")
#   #############################################################
#   if(!exists("ploidy"))
#     stop("must specify ploidy of genomes")
#
#   # -- 1. Set up the required variables
#
#   ## Path to directory that contains the MCScanX executable
#   MCScanX.path = "/Users/jlovell/Documents/comparative_genomics/programs/MCScanX"
#
#   ## Working directory, that contains at least the raw_annotation and raw_assembly subdirectories
#   wd = "/Users/jlovell/Desktop/genespace_runs/grasses2"
#   setwd(wd)
#
#   ## Genome IDs, which must match subdirectories in the raw_annotation folder
#   genomeIDs = c("Sbicolor", "Zmays", "Bstacei")
#
#   ## Ploidy of the genomes
#   ploidy = c(2,4,2)
#   names(ploidy) <- genomeIDs
#
#   ## Method to define synteny, either 'global' or 'pairwise'
#   synteny.method = "pairwise"
#
#   ## Minimum number of hits in a block
#   min.blockSize = 10
#
#   ## MCScanX parameters
#   MCScanX.m.param = 50
#   MCScanX.s.param = 5
#
#   ## Maximum number of blast duplicates to retain for each haploid genome
#   max.hits.per.haplotype = 4
#
#   ## N. genes away from known syntenic hits to call syntenic
#   ngenes.synteny.search = 200
#
#   ## N. genes away from known syntenic hits to call syntenic
#   ngenes.synteny.buffer = 50
#
#   ## Set number of parallel processes to run
#   n.parallel.procs = 6
#
#   ## Plotting diagnostics options
#   dotplot.palette = colorRampPalette(c("red3","darkorange","blue3"))
#   dotplot.minHitsByRefChr = 100
# }
