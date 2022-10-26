# This is **beta** version 1.0.0 of the GENESPACE R package. **USE WITH CAUTION** 

The manuscript describing GENESPACE is now pubished in eLife. [Find the article here](https://elifesciences.org/articles/78526). Please cite this if you use GENESPACE. 

Updates since v0.9.4 ...

- **major changes!! backwards compatibility to < v1 is broken in v1.0.1**: 
    - The user must now specify pre-formatted bed files instead of raw gff3 formatted annotations. This allows GENESPACE to catch issues with data formatting before running orthofinder etc. We are truly sorry to alter the specification and input so much, but it is unavoidable -- we want GENESPACE to integrate with existing bioinformatic pipelines, and too much internal customization was required, which produced hard-to-resolve errors.
    - **Resolved issues with plot visualization in interactive R sessions**: xy dotplots and riparian plots are no longer written to the plot window. These are now saved as vector pdfs in the new /dotplot and /riparian directories respectively.
    Better organization of results output. Previously all output went into /results. /results now only contains the parsed orthofinder results and annotated bed files. All other results are placed in /syntenicHits, /dotplots, /riparian, and /pangenome. 
    - All GENESPACE functions are now integrated into `run_genespace()`, which completes all steps internally (the user can still run each step independently too). 
    - The only change to the actual algorithm is a modification to how large tandem arrays are handled and masked. This will only affect a small number of genomes. 

- **minor changes**: 
    - better integration if running orthofinder separately
    - code cleanup and better documentation, both in line and in help files
    - much better data QC and handing of errors in parameter specification
    - ggplot2-integrated graphics and R object output for downstream interactive viz through ggplotly. 
    - internal position calculation for riparian plots and inferred syntenic positions transferred to a new function `integrate_synteny()`. This speeds up riparian plotting by 10x or more. 
    - `plot_riparian()` function generalization that allows the user to specify a block coordinate matrix and nothing else (allowing plotting of data from other packages). 
    - parallel processing switched from parallel to biocparallel
    - replace `orthofinderMethod = "fast"` with `onewayOF = TRUE` to match the new `orthofinder -1` specification. 

- **upcoming changes for v1.1.0**:
    - An executable Rscript for `run_genespace` and parameterization via .json to allow GENESPACE to be called directed from the shell 
    - `synteny()` will be able to be called from the shell on a single blast8-formatted text file
    - a new function `query_genespace()`, which allows positional or orthology based queries of pan-genome, pairwise hits, raw orthogroups/orthologues, or riparian source data. 
    
## 1. Overview

GENESPACE conducts synteny- and orthology-constrained comparative genomics, which can be used to (1) make multi-genome graphical comparisons, (2) define syntenic block breakpoints, and (3) build a pan-genome annotation. 

For the purposes of this readme, we assume that both OrthoFinder and GENESPACE can be installed on the same machine and run together. However, this may not always be the case. If you need to run OrthoFinder separately, see section 7.1 Additional considerations: running OrthoFinder separately. 

#### 1.1 Dependencies

(see 2.2-2.5 below for details) There are three dependencies outside of R (OrthoFinder, DIAMOND2 (included with OrthoFinder), MCScanX) required to run GENESPACE. See below for installation instructions. In R, four additional dependency packages are required: Biostrings, dbscan, igraph and data.table. 

#### 1.2 Input data

(see 3.1-3.4 below for details) GENESPACE operates on a directory containing gene .bed and peptide .fasta files where each entry is the primary transcript of a gene and the headers of the peptide exactly match the name (4th) bed column. You can make these yourself, or use GENESPACE's built in conversion functions 

#### 1.3 Running GENESPACE

For most purposes, GENESPACE v1 requires just two function calls in R. 

```
gsParams <- set_params(...)
run_genespace(gsParam)
```

This will run all steps (or return information to the user on what needs to be done to have a successful run). However, the various sub-functions are also available to the user and can be run independently. Here is how each step works:

1. `gsParams <- run_orthofinder(gsParams)` call orthofinder from R (or if orthofinder cannot be called, return the code to run it). 
2. `gsParams <- annotate_bed(gsParams)` calculate tandem array members and add orthofinder-derived data. 
3. `gsParams <- synteny(gsParams)` annotate blast files with synteny and orthogroup information. If `orthofinderInBlk = TRUE`, this re-calls annotate_bed on re-calculated orthogroups. 
4. `gsParams <- pangenome(gsParams)` combine synteny and orthology information into a single matrix
5. `gsParams <- integrate_synteny(gsParams)`
6. `gsParams <- plot_riparian(gsParams)` make multi-genome synteny plots

#### 1.5 Notes and considerations when parameterizing your run

- If any of your genomes are polyploid or have a history of whole-genome duplications, the phylogenetic context of the genomes is VERY important. This topic is covered in [the paper](https://elifesciences.org/articles/78526), but also visited in [this issue posted on github](https://github.com/jtlovell/GENESPACE/issues/28#issuecomment-1241280788).  
- When making riparian plots, it is far better to have blocks colored by a HAPLOID genome than a polyploid. GENESPACE can handle polyploid references in plot_riparian, but, since blocks get split by reference chromosome of origin, this can produce a much larger file with less clear chromosomes-of-origin. 
- The hardest part of GENESPACE is getting the dependencies installed and your annotations formatted correctly. Don't be dissuaded if you have having trouble with these steps! Open an issue and we can help!

## 2. Installation

#### 2.1 Install R

GENESPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

#### 2.2 Install third party programs

GENESPACE also requires third-party software that can be installed as follows.

`OrthoFinder` (which includes `DIAMOND2`) is most simply installed via conda (in the shell, not R). 

```
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources. See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

#### 2.3 Open an interactive R session

You can open R using the the GUI, Rstudio, or as an interactive session in the terminal. If you are using a conda environment or specifying the path to orthofinder in the `$PATH`, you need to enter an R interactive session from that shell environment, which can be done by either calling `R` (for command line interface) `open -na rstudio` if using Rstudio. 

#### 2.4 Install GENESPACE R package

Once in R, GENESPACE can be installed directly from github via:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")
```

#### 2.5 Install R dependencies

The above command will install the CRAN-sourced dependencies (`data.table`, `dbscan` and `R.utils`). The bioconductor dependencies (`rtracklayer` and `Biostrings`) need to be installed separately via:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```

## 3. Getting the raw annotation files formatted

This can be the trickiest part of running GENESPACE. If you are having trouble, open an issue. 

#### 3.1 GENESPACE-readable annotations

For each genome, GENESPACE needs .bed formatted coordinates of each gene and peptides that exactly match the "name" (4th) bed column. 

The bed file needs to have the first four tab-separated fields (chromosome ID, start, end, name) in a bed file. Other fields are allowed, but will be ignored by GENESPACE. A header is not allowed. For example:

```
chr1  1500  2000  gene1
chr1  2500  2200  gene2
```

The fasta file is just the standard format with the header exactly matching the names. The order of the sequenced does not matter, but each header must uniquely match the names in the bed file.

```
>gene1
MGASGRGAGEQQSPSSTG
>gene2
MLVMSECKGRDRSPSSSM
```

These files must be stored in a directory, lets call it `genomeRepo` with a single uniquely-named subdirectory for each genome. For example:

```
/genomeRepo
└─ genome1
    └─ peptides.fa
    └─ genes.bed
└─ genome2
    └─ peptides.fa
    └─ genes.bed
  ...
```

Here, the directory names "genome1" and "genome2" will be used as the genomeIDs in GENESPACE unless alternative IDs are given (see below). The files themselves can be named whatever, as long as there is a single file ending in ".fa" that contains peptide sequences and a bed-formatted annotation ending in .bed in each genome subdirectory. If there are other files in the subdirectories, that is fine as long as none end in .fa or .bed. The parameter checking step in GENESPACE will ensure these conditions and will return an error if they are not met. 

#### 3.2 Parsing a single pair of annotation files for one genome

Here, we assume that the starting data is a gff3-formatted annotation file and a fasta file contain peptide sequences. Both the gff3 and fasta can be gzip compressed or not. Other compression formats are not supported. 

The first step is to convert the gff3 into the simpler .bed format and ensure that the peptide fasta headers exactly match the 'id' field (4th bed column). GENESPACE offers a flexible method `match_annotations`. For a single genome, the user provides path names to the peptide and gff3 annotations for a single genome. This approach will require some scripting in R to ensure that all annotations of interest are parsed correctly, but also give the user flexibility to store the annotations in whatever directory structure works for them. 

```
match_annotations(
  path2fasta = "path/to/rawAnnotations/[fileName].fa",
  path2gff3 = "path/to/rawAnnotations/[fileName].gff3",
  genomeID = "nameForGenome",
  path2genomeRepo = "/place/to/store/genomes",
  ...)
```

#### 3.3 Parse many annotations at once. 

If you are analyzing many genomes that are all formatted identically, you can automate the annotation process in `match_annotations` by specifying the path to your genome repo and the names of the subdirectories (genomeIDs) that contain the annotation files. 

```
match_annotations(
  genomeIDs = c("genome1", "genome2"), # names of the subdirectories
  path2rawAnnotationRepo = "/place/to/find/raw/annotationFiles"
  path2genomeRepo = "/place/to/store/genomes",
  ...)
```

#### 3.4 Parameterizing annotation matching 

There are two presets that parse NCBI- and phytozome-formatted annotation (`preset = "NCBI"`, `preset = "phytozome"`). For all other formats, the user will need to specify parameters to parse the attributes gff3 field and the fasta header:

- `whichAttributeField`: the name (or index) of the ';' separated field in the attributes column contains the geneID. 
- `whichHeaderField`: the name (or index) of the `headerDelim`-separated field that contains the geneID in the fasta header
- `stripAttribute`: character string to strip off the geneID in the gff3 attribute column
- `stripHeader`: character string to strip off the geneID in the fasta header

Parsing can take some troubleshooting, which can be aided by setting `troubleshoot = TRUE`. This prints the first 10 lines of the raw and parsed gff and fasta headers. 


## 4. Setting and checking GENESPACE parameters

#### 4.1 Required GENESPACE parameters

For simplicity, the user should set the GENESPACE parameters once before starting the run. This produces an R object and a yaml file. This is accomplished by the function `set_params`, which requires the following user-specified input:

- `path2orthofinder` - the file path to the OrthoFinder installation
- `path2mcscanx` - the file path to the MCScanX installation
- `path2genomeRepo` - same as above, file path to where the parsed annotations are stored. 
- `genomeIDs` - the names of the genomes and subdirectories in `path2genomeRepo` to use in the GENESPACE run
- `wd` - file path to the working directory for the GENESPACE run. 

#### 4.2 Common optional parameters

There are many ways to customize the GENESPACE run. Check the helpfile for `set_params` for details about them all. Here are the other most common:

- `useHOGs` - should hierarchical orthogroups be used instead of global orthogroups. HOGs tend to split arrays and better reflect orthologs. Global orthogroups better capture paralogous regions. Default is TRUE unless any of the genomes have ploidy > 1.  
- `genomeNames` - replaces genomeIDs as the names for genomes in all GENESPACE output. 
- `ploidy` - integer vector specifying the ploidy of each genome. If not specified, assumes all genomes are haploid
- `outgroup` - genomeIDs that should only be used in the orthofinder run but not considered for synteny etc. If not specified, no outgroup is used. 
- `blkSize` - default is 5. This is the minimum number of anchored blast hits needed to call a syntenic block. 
- `nGaps` - default is 5. Passed to `MCScanX -m` as the maximum number of gaps in a syntenic block. 
- `synBuff` - default is 100. This is the maximum euclidean distance (in gene-rank order units) from a syntenic anchor for a gene to be considered in the syntenic buffer. Also the minimum distance between two genes to split a tandem array into two clusters.   
- `nCores` - the number of parallel processes to run
- `fastOrthofinder` - should the 'fast' OrthoFinder method be used in lieu of the default?
- `orthofinderInBlk` - should OrthoFinder be re-run within each large syntenic region? Default is TRUE if any genome ploidy > 1. 

#### 4.3 `set_params` methods

`set_params` first checks that all user-specified input is OK:

1. Make sure that OrthoFinder > v2.5.4, DIAMOND > v2.0 and MCScanX are installed and the paths are specified correctly. 
2. Checks that the genomeRepo contains each of the genomeID subdirectories and that each subdirectory contains exactly one peptide fasta and one .bed annotation with matching headers and names. 
3. Checks that the working directory exists or can be created
4. Checks that all optional arguments specify valid parameters

Any errors in parameter specifications or annotation parsing are caught up front, and if problems were found, informative errors are given. Once the parameters are checked, `set_params` converts the parameters into a list and writes the list to a .yaml formatted text file stored in the working directory. The parameter list needs to be stored as an R object e.g. `gsParams <- set_params(...)`, then is provided to downstream functions. If you want to re-visit a run after closing R or removing the gsParam object, you just need to run `gsParams <- set_params(wd = "path/to/the/run/directory")`. This will automatically find the yaml file and re-check and populate the parameter list, returning it as an R object. If the user wants to modify the parameters later on, you can simply re call and specify the change desired. For example to decrease the number of parallel processes: `gsParams <- set_params(gsParams, nCores = 1)`. 


## 5. Running orthofinder

#### 5.1 Running orthofinder outside of R

Running OrthoFinder is by far the most time and CPU consuming step of GENESPACE. Large runs can take days on a small machine. In these cases, it may be useful to run OrthoFinder outside of R on a cluster, but only if you want to run default OrthoFinder and not `fastOrthofinder`. 

To do this, you have two options. If you are not planning on running OrthoFinder within blocks, just set `set_params(path2orthofinder = NA)` and call `run_orthofinder()` normally. If you want to run OrthoFinder within blocks, you will still need a valid `path2orthofinder` specified in `set_params`. This this case, you should call `run_orthofinder(onlyReturnCode = TRUE)`. Either method will just return the code needed to run OrthoFinder and won't actually run it. You are then free to paste that code into a shell script and run on another machine (assuming you have access to the working directory).

#### 5.2 Running OrthoFinder from R

GENESPACE can call OrthoFinder from R assuming a valid `path2orthofinder` via `run_orthofinder()`. If the default method is used, this is a simple function where peptides are copied from the genomeRepo and fed to OrthoFinder. The `fastOrthofinder` method is more involved and runs four steps: 

1. OrthoFinder is run through the file-formatting step (-op), which returns DIAMOND2 blastp calls.
2. blastp calls where the query genome has >= the number of genes in the target are run, but replacing the default --more-sensitive with --fast. 
3. Each inter-genomic blast8 output file is "mirrored" and re-named so that the target --> query and query --> target
4. OrthoFinder -b is run on the resulting blast files. 

#### 5.3 Parse OrthoFinder results

GENESPACE uses four sets of data from OrthoFinder:

1. Species number --> species name and geneID --> gene name dictionaries
2. Gene ID --> orthogroup dictionary
3. N0 phylogenetically hierarchical orthogroup table. 
4. Pairwise ortholog table. 

Each of these datasets are copied into /results as SpeciesIDs.txt, SequenceIDs.txt, N0.tsv,  Orthogroups.tsv, and the Orthologues directory (Orthologues/Orthologues_genomeID/genomeID1__v__genomeID2.tsv). **NOTE** if you run OrthoFinder outside of R in a separate directory, you need to copy these files over manually.  

If you call `run_orthofinder` and a previous run is detected, it will NOT overwrite the results ever (in previous GENESPACE versions you could set `overwrite = TRUE`). If a previous run exists, GENESPACE will just look for the above set of five files and copy them over to /results. If the files exist in the /results directory `run_orthofinder` will return a warning telling you that, if you want to re-run orthofinder, you will need to manually delete the contents of /orthofinder and /results directories. 

## 6. Prepare data for synteny

#### 6.1 Combine and annotate the bed files

The genome ID ('genome' column) is added to each bed file, then the files are concatenated into a long-format text file stored in results/annotBed.txt. The orthofinder geneIDs (column = "ofID"), phylogenetically hierarchical orthogroups ("globHOG") and global orthogroups ("globOG") and then added to this file. If useHOGs = TRUE, the globHOG column is copied into the 'og' column. Otherwise the globOG column is used. 

#### 6.2 Infer positions and representatives of tandem arrays

We define tandem arrays as orthogroup (or HOG if useHOGs = TRUE) members that are <= synBuff genes apart along a chromosome. Euclidean distances between genes are calculated within all chromosomes and orthogroups and tandem arrays are extracted. The most physically central gene (gene rank order, then physical position, then peptide length as tiebreakers) are flagged as the array representative gene (arrayRep = TRUE) and the order of array representative genes are re-calculated. This process of clustering and re-ranking genes is repeated until now new genes are added to an array. 

#### 6.3 QC the genomes

With the orthogroups integrated, we count the number of array representatives per genome-by-chromosome combination. Chromosomes (or contigs, scaffolds etc) with n. array representatives < blkSize are not considered for synteny and flagged ignore = TRUE in the annotBed file. If any genome has > 5% of all genes on such small chromosomes, an error will be returned. This error can be turned off by specifying `annotate_bed(allowBadGenomes = TRUE)` or resolved by ignoring that genome by running `set_params(genomeIDs = c("genomes", "other", "than", "bad", "ones"))`. 

## 7. Run synteny

Synteny runs in a pairwise manner - BLAST hits from each unique pair of genomes where the query has >= the number of genes as the target are concatenated and read into R. The following steps are accomplished ...

#### 7.1 If self-blast (intragenomic)

Self hits for array representative genes are the anchors. Array representative hits within synBuff of self hits are in the synteny buffer. 

If ploidy is > 1, the hits in the synteny buffer are masked and the synteny search is accompished via steps 7.2-7.6

#### 7.2 Flag 'potential anchor hits'

Defined hierarchically, based on the following

1. both the query and the target are array representatives
2. if onlyOgAnchors, both the query and target genes are in the same orthogroup
3. the top n scoring hits for each query gene where n is the ploidy of the target genome
4. re-calculate gene rank order positions for only the potential anchor hits. 

#### 7.3 Pull 'initial anchors'

To get high-quality initial anchors that extend to the very ends of collinear regions, we run MCScanX twice, first globally on potential anchor hits, then only on hits near global collinear hits. 

1. Feed potential anchor hits to MCScanX_h with -m = nGaps and -s = blkSize
2. Pull all hits near collinear hits with dbscan, where minPts = blkSize and radius = sqrt(2) * synBuff and group them into 'initial syntenic regions'.
3. Re-calculate gene rank order positions within each syntenic region and re-run MCScanX_h on potential anchor hits separately in each region. Collinear hits are the 'initial anchors'. 

#### 7.4 Finalize syntenic anchor hits

1. Pull all arrayRep [and optionally orthogroup-constrained] hits within sqrt(2) * synBuff gene rank order of initial anchors, retaining dbscan cluster identity as 'synRegion'. 
2. Re-calculate gene-rank order position within each synRegion and re-run MCScanX_h. 
3. Re-calculate gene-rank order on collinear hits and cluster with dbscan where minPts = blkSize and radius = sqrt(2) * blkSize. Collinear hits in blocks with > blkSize hits are the 'final anchors' and flagged in the blast files as isAnchor = TRUE. 
4. Re-calculate gene-rank order for all final anchor hits. Split blocks that have no duplicate members but have overlapping coordinates through run-length equivalent decoding. These are the final blocks and flagged in the blast files as 'blkID'. 

#### 7.5 Make syntenic regions and hits in buffer

1. Re-calculate gene-rank order for all hits. Pull hits within synBuff * sqrt(2) of the final anchors. These are 'inBuffer' syntenic hits, flagged in the blast files as inBuffer = TRUE.
2. Re-calculate gene-rank order on all inBuffer hits and cluster iwith dbscan using a radius of synBuff * sqrt(2). These are the syntenic regions, flagged in the blast files as 'regID'.

#### 7.6 Make syntenic orthogroups

The syntenic orthogroup blast hits across all pairs of genomes are agregated and clustered in igraph into cyclic undirected graph. This is decoded into an integer vector and added to the annotated bed file in the column labeled 'synOG'. 

#### 7.7 Modifications if nSecondaryHits > 0. 

If the user wants to search for secondary (likely paralagous) hits, and specifies nSecondaryHits > 0, GENESPACE repeats 7.2-7.7, except that all inBuffer hits are masked, and the following parameters are replaced: onlyOgAnchors <-- onlyOgAnchorsSecond, blkSize <-- blkSizeSecond, nGaps <-- nGapsSecond. 

#### 7.8 Modification if orthofinderInBlk = TRUE

If the user wants to use the orthofinderInBlk method (recommended if any genomes are polyploid or nSecondaryHits > 0), the following steps are run. 

1. For each syntenic region, all blast hits (even non-array representatives) are fed to orthofinder
2. The orthogroup (or HOG) graph is extracted and compiled across regions
3. This vector is added to the annotated bed file in the 'inBlkOG' column. 
4. The 'og' column is re-made by merging the synOG and inBlkOG subgraphs. 
5. Steps 6.1-7.7 are re-run. 

## 8 Aggregate syntenic information across genomes

#### 8.1 Calculate coordinates of syntenic blocks

For each syntenic block, we calculate the positions of the bounding genes and return the coordinates as a text file in the results directory. 

#### 8.2 Determine the expected position of each gene against each genome. 

We seek to know the syntenic position of each array representative gene. This lets us quickly phase hits and blocks by ancestral chromosome and build the coordinate system for a pan-genome annotation. For each pair of genomes, we do the following:

1. Pull all arrayRep genes and split them by syntenic block. 
2. Within each block, re-calculate the gene-rank order positions of the anchor hits. 
3. Iteratively prune/cluster anchor hits so that each microblock is has a 1:1 ratio of positions between the query and target genes (perfect diagonal).
4. Calculate the syntenic position of each arrayRep gene using linear interpolation. 
5. Use the nearest neighbor position of genes within a block, but not a microblock. Flag these as 'medium-confidence' position estimates. 
6. Use the nearest neighbor positions of genes within a region, but not a block. Flag these as 'low-confidence' position estimates. 
7. Aggregate into the annotated bed format like above, but with repeated entries for each syntenic position across chromsomes and genomes. 

#### 8.3 Split syntenic blocks by chromosome of origin for each genome. 

1. Hits are read in and the expected chromosomal positions of both query and target genes are inferred. 
2. Runs of > synBuffer hits that map to a single chromosome are assigned to that chromosome. 
3. If the same gene hits multiple chromosomes, runs of those chromosomes are both included. 








#### 1.6 If you want to run OrthoFinder on a server, but can't get GENESPACE installed there ...

GENESPACE (and R in general) have a number of installation requirements that can be problematic when running on a server. If you don't want to run orthofinderInBlk (i.e. all of your genomes are haploid), but want to do a BIG orthofinder run that can't be accomplished on a machine where GENESPACE can be installed, you can now more easily integrate the results of orthofinder in GENESPACE with the following:

1. Ensure that your bed 'name' (4th) column matches exactly your peptide headers
2. Run orthofinder on a server using the re-named peptides 
3. Make a working directory on the machine where you will run genespace
4. Move the .bed and peptide fasta files into the GENESPACE working directory, labeled bed/genomeID.bed and peptide/genomeID.fa, respectively. 
5. Move the following orthofinder results files into the /results subdirectory of the working directory:
    - /Orthogroups/orthogroups.tsv
    - /WorkingDirectory/SpeciesIDs.txt
    - /WorkingDirectory/SequenceIDs.txt
    - /WorkingDirectory/Blast*.txt.gz
    - Phylogenetic_Hierarchical_Orthogroups/N0.tsv
    - The full /Orthologues directory
6. Run GENESPACE as above

