#' @title Accessory function to run OrthoFinder
#' @description
#' \code{run_orthofinder} Call OrthoFinder from R
#'
#' @param gsParam a genespace param list, made by init_genespace.
#' @param verbose Logical, should updates be printed to the console?
#'
#' @details Simple directory parser to find and check the paths to all
#' annotation and assembly files.
#'
#' @return A list containing paths to the raw files. If a file is not found,
#' path is returned as null and a warning is printed.
#'
#' @examples
#' \dontrun{
#' # coming soon
#' }
#'
#' @export
run_orthofinder <- function(gsParam, verbose = TRUE){
  ##############################################################################
  # 1. Get things set up
  # -- 1.1 combine genomeIDs and outgroups
  genomeIDs <- gsParam$genomeIDs
  if(!is.na(gsParam$outgroup)[1]){
    genomeIDs <- c(genomeIDs, gsParam$outgroup)
    genomeIDs <- genomeIDs[!duplicated(genomeIDs)]
  }


  # -- 1.2 get paths
  ofDir <- gsParam$paths$orthofinder
  tmpDir <- gsParam$paths$tmp
  pepDir <- gsParam$paths$peptide

  # -- 1.3 check to make sure that the peptide fastas are in there
  pepf <- file.path(pepDir, sprintf("%s.fa", genomeIDs))
  if(!all(file.exists(pepf)))
    stop(sprintf(
      "something is wrong with the peptide files. could not find: %s",
      paste(pepf[!file.exists(pepf)], collapse = "\n")))

  # -- 1.4 check that the orthofinder directory does not exist or is empty
  if(dir.exists(ofDir)){
    fs <- list.files(path = ofDir)
    if(length(fs) == 0){
      unlink(ofDir, recursive = T)
    }else{
      stop(sprintf(
        "orthofinder directory %s exists, remove this before proceeding\n",
        ofDir))
    }
  }

  # -- 1.5 rename the required parameters
  onewayBlast <- gsParam$params$onewayBlast
  diamondUltraSens <- gsParam$params$diamondUltraSens
  path2orthofinder <- gsParam$shellCalls$orthofinder
  nCores <- gsParam$params$nCores
  runOfInR <- !is.na(path2orthofinder)

  ############################################################################
  # 2. set up the directory structure and make sure things look good
  # -- copy peptides over to tmp directory
  # this allows for more genomeIDs in /peptide than just those in gsParam
  if(verbose)
    cat(strwrap(sprintf(
      "Copying files over to the temporary directory: %s",
      tmpDir), indent = 8, exdent = 16), sep = "\n")

  # -- if the tmpDir exists, remove and re-create it
  if(dir.exists(tmpDir))
    unlink(tmpDir, recursive = T)
  dir.create(tmpDir)

  # -- copy over the peptide files
  nu <- file.copy(pepf, tmpDir)

  ############################################################################
  # 3. Get the orthofinder command
  ofComm <- sprintf(
    "-f %s -t %s -a 1 %s %s -X -o %s",
    tmpDir, nCores,
    ifelse(onewayBlast, "-1", ""),
    ifelse(diamondUltraSens, "-S diamond_ultra_sens", ""),
    ofDir)

  # -- strip out extra spaces that may exist
  ofComm <- gsub("  ", " ", gsub("  ", " ", ofComm))

  ############################################################################
  # 4. If orthofinder is available, run it from R
  if(runOfInR){
    if(verbose)
      cat(strwrap(sprintf(
        "Running the following command in the shell: `%s %s`.This can take a
        while. To check the progress, look in the `WorkingDirectory` in the
        output (-o) directory", path2orthofinder, ofComm),
        indent = 8, exdent = 16), sep = "\n")

    outp <- system2(
      path2orthofinder,
      ofComm,
      stdout = TRUE, stderr = TRUE)
    if(verbose)
      cat(paste(c("\t", outp), collapse = "\n\t"))
  }else{
    ############################################################################
    # 5. If not, print how to do it and stop.
    stop(cat(strwrap(
      "Could not find a valid path to the orthofinder program from R. To run
        orthofinder, ensure that the orthofinder program is in the $PATH, then
        call the following from the shell: \n", indent = 0, exdent = 8),
      sprintf("orthofinder %s", ofComm)),
      "Once OrthoFinder has been run, re-call run_genespace",sep = "\n")
  }
}

