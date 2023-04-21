## GENESPACE: an R package for synteny- and orthology-constrained comparative genomics. 

This is version 1 of GENESPACE. GENESPACE is an analytical pipeline to understand patterns of synteny and orthology across multiple genomes. The manuscript describing GENESPACE is now pubished in eLife. [Find the article here](https://elifesciences.org/articles/78526). Please cite this if you use GENESPACE. 

There are currently two tutorials that (1) illustrate [what GENESPACE does and how to use it](https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genespaceGuide.html) and (2) demostrate how to [customize your riparian plots](https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/riparianGuide.html). This README is primarily to show how to get GENESPACE up and running.

Beginning in v1.2.0, GENESPACE also supports some basic genome-level visualization and quality control methods. Find details on these new functions in the [tutorial on sliding windows and contig mapping](https://htmlpreview.github.io/?https://github.com/jtlovell/tutorials/blob/main/genomeVizGuide.html). These are new functions and may have some bugs. Please report any issues that arise. 

**NOTE** v0.9.x is no longer supported. Please upgrade to v1.1+. There are some significant changes to the structure of GENESPACE in V1 (but few changes to the underlying algorithms). There have been a few minor changes since the v1.1.4 release on 1-March 2023. Either install from the master branch or most recent release. The next planned major update will be v1.2.x during the spring or summer 2023, which will include Rscript integration so that GENESPACE can be called directly from the command line without an interactive R session. 

##################
## 1. Quick start

#### 1.1 Input files and dependencies

To run GENESPACE, you need:

1. A valid installation of OrthoFinder, MCScanX and R (see part 2: software installation)
2. A peptide (protein) fasta and a bed-formatted gene model annotation file for 2 or more genomes (see part 3: formatting your annotations)

#### 1.2 Running GENESPACE

With these data in hand, you initialize a GENESPACE run, then run the pipeline from R:

```
library(GENESPACE)
gpar <- init_genespace(
  wd = "/path/to/your/workingDirectory", 
  path2mcscanx = "/path/to/MCScanX/")
gpar <- run_genespace(gsParam = gpar) 
```

`run_genespace` completes the full pipeline, which includes:

1. Tandem array discovery
2. Syntenic block coordinate calculation
3. Synteny-constrained orthogroups
4. Pairwise dotplots
5. Syntenic position interpolation of all genes
6. Pan-genome annotation construction
7. Multi-genome riparian plotting

See the tutorial for more details

#### 1.3 Running with example data

We recommend doing a test run with real data before applying GENESPACE to your own datasets. To help facilitate this, run these commands in R (assuming valid installation):

```
###############################################
# -- change paths to those valid on your system
genomeRepo <- "~/path/to/store/rawGenomes"
wd <- "~/path/to/genespace/workingDirectory"
path2mcscanx <- "~/path/to/MCScanX/"
###############################################

# -- download raw data from NCBI for human and chicken genomes
dir.create(genomeRepo)
rawFiles <- download_exampleData(filepath = genomeRepo)

# -- parse the annotations to fastas with headers that match a gene bed file
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi",
  genespaceWd = wd)

# -- initalize the run and QC the inputs
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = path2mcscanx)

# -- accomplish the run
out <- run_genespace(gpar)
```

##################

## 2. Software installation

Installation for v1 is identical to v0.9.4. Use the same conda (or other) installation environment as before. To install from scratch, you need R, a few 3rd party programs and a few R packages. Here is how to install these:

#### 2.1 Install R

GENESPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

#### 2.2 Install orthofinder

`OrthoFinder` (which includes `DIAMOND2`) is most simply installed via conda (in the shell, not R). 

```
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources. See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  
Regardless of how `OrthoFinder` is installed, ensure that you have `OrthoFinder` version >= 2.5.4 and `DIAMOND` version >= 2.0.14.152.

#### 2.3 Install MCScanX

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

#### 2.4 Install GENESPACE

Once the above 3rd party dependencies are installed, get into R. If you made a conda environment, its useful to open R directly from that environment so that OrthoFinder stays in the path. 

```
conda activate orthofinder
open -na rstudio # if using rstudio, otherwise, simply `R`
```

Once in R, the easiest way to install GENESPACE uses the package `devtools` (which may need to be installed separately):

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")
```

#### 2.5 Install R dependencies

If they are not yet installed, `install_github` will install a few dependencies directly (ggplot2, igraph, dbscan, R.utils, parallel). However, you will need to install the bioconductor packages separately:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)
```

##################

## 3. Formatting your annotations

This can be the trickiest part of running GENESPACE. If you are having trouble, open an issue. 

#### 3.1 GENESPACE-readable annotation format

For each genome, GENESPACE needs:

1. bed formatted coordinates of each gene (chr, start, end, name), other fields are allowed, but will be ignored by GENESPACE
2. peptide sequences in fasta format, where the header exactly matches the "name" (4th) bed column

Example bed file:

```
cat /bed/genome1.bed | head -2
chr1  1500  2000  gene1
chr1  2500  2200  gene2
```

Example peptide fasta file: 

```
cat /peptide/genome1.fa | head -4
>gene1
MGASGRGAGEQQSPSSTG
>gene2
MLVMSECKGRDRSPSSSM
```

#### 3.2 Manual parsing of raw annotation files

It's possible to format all files using custom scripts. If you opt for this route, re-name and place the genespace-formatted annotations in your working directory. For example for two genomes ("genome1", "genome2"):

```
/workingDirectory
└─ peptide
    └─ genome1.fa
    └─ genome2.fa
└─ bed
    └─ genome1.bed
    └─ genome2.bed
```

#### 3.3 Naming and storing raw annotation files

Its usually preferable to have a static repository containing the raw versions of all of your genome annotations, lets call this directory `/genomeRepo`. Each unique genome annotation should be stored in its own subdirectory with an informative name. Here, each genome is given a species id, genotype id and genome version. The annotation files of interest are the gff3-formated protein coding gene features and the primary transcript translated cds (peptides). These are the two most common types of annotation files available. Download them from your favorite repositories and place them in a directory structure like below. The file names within each genome directory can be anything as long as they have file names that can be found through regular expression matching (see `parse_annotations(gffString, faString)`).

```
/genomeRepo
└─ species1_genoX_v1.0_NCBI
    └─ peptidesGenoX.fa
    └─ genesGenoX.gff3
└─ species2_genoY_v1.0_NCBI
    └─ peptidesSpecies2.fa
    └─ genes.gff3
└─ species3_genoW_v1.0_phytozome
    └─ peptides.fa
    └─ genesSpecies2v1.gff3
└─ species4_genoZ_v1.0_otherRepo
    └─ peptides.fa
    └─ genes.gff3
  ...
```

#### 3.4 Partially-automated annotation parsing

Assuming the directory structure from 2.4, we want to parse the raw annotations into the directory structure outlined in 2.3. The convenience function `parse_annotations` does this. Lets take the simple case where we have translatedCDS and gene.gff3 files from NCBI. In this case, we can specify the "ncbi" present which will automatically parse these files correctly:

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = c("species1_genoX_v1.0_NCBI", "species2_genoY_v1.0_NCBI"),
  genomeIDs = c("species1", "species2"),
  gffString = "gff3",
  faString = "fa",
  presets = "ncbi", 
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

This will parse the two NCBI-formatted annotations into matched bed and fasta files and place those in the /bed and /peptide subdirectories of the GENESPACE working directory. Importantly, it will re-name the annotations with the genomeIDs. This allows for informative naming conventions in the raw file repo, but simpler names for plotting and analysis. 

#### 3.5 Non-standard annotation formats

Most repositories have their own annotation formats. You may want to combine annotations from NCBI with those from phytozome (plant genome repo) or non-standard repositories. To do this, you simply re-run parse_annotations with different genomes. For example add a phytozome-formatted annotation:

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = "species3_genoW_v1.0_phytozome",
  genomeIDs = "species3",
  gffString = "gff3",
  faString = "fa",
  presets = "phytozome",
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

Then add a non-standard annotation where there is only a simple geneID header on the peptide and the gene names are in the "GeneID" ';'-delimited 9th gff column. 

```
parsedPaths <- parse_annotations(
  rawGenomeRepo = "/genomeRepo", 
  genomeDirs = "species4_genoZ_v1.0_otherRepo",
  genomeIDs = "species4",
  gffString = "gff3",
  faString = "fa",
  headerEntryIndex = 1, 
  gffIdColumn = "GeneID",
  genespaceWd = "/path/to/GENESPACE/workingDir")
```

#### 3.6 Additional parameterization

There are lots of parameters that can be combined to parse every annotation format. See the help file for details. The main parameters to get your annotation parsed correctly are: 

  - `gffIdColumn`: the name (or index) of the ';' separated field in the attributes column contains the geneID.
  - `headerEntryIndex`: the name (or index) of the `headerSep`-separated field that contains the geneID in the fasta header
  - `gffStripText`: character string to strip off the geneID in the gff3 attribute column
  - `headerStripText`: character string to strip off the geneID in the fasta header

Parsing can take some troubleshooting, which can be aided by setting `troubleshoot = TRUE`. This prints the first 10 lines of the raw and parsed gff and fasta headers. For some annotation types, you may need to give a fairly complex regular expression to gffStripText or headerStripText. For example, if the fasta name has some additional unwanted info after the last "." in the name: `headerStripText = '[.][^.]+$'`. 


#### 3.7 Convert a previous run and re-run GENESPACE

Earlier versions of GENESPACE used 'simplified gff3' annotation files. To make GENESPACE compatible with existing bioinformatic software, we have opted to use the bed format for v1 onward. If you have a previous GENESPACE run, you can easily convert the input to v1 as follows (this also copies over orthofinder results so those can be re-used too): 

```
wd <- "/place/to/store/new/run"
convert_input2v1(
  v1Dir = wd, 
  existingDir = "/path/to/existing/run")
```
