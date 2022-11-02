# This is **beta** version 1.0.0 of the GENESPACE R package. **USE WITH CAUTION** 

## 1. Quick start 

There are some big changes here compared to v0.9.4 (previous release), see section 2 below. v1.0.0 really is a beta release. There are some known issues that will be resolved before the first full v1.0.1 release, and probably some bugs will come up. If you use this beta version (thanks!), please report any problems you run into and we will address them ASAP. 

The manuscript describing GENESPACE is now pubished in eLife. [Find the article here](https://elifesciences.org/articles/78526). Please cite this if you use GENESPACE. 

Run v1 GENESPACE as follows:

### 1.1 Install v1

Get into R like you would for previous versions (see below for 3rd party installation instructions). Then install v1 from the dev branch as follows:

```
detach("package:GENESPACE", unload = TRUE)
devtools::install_github("jtlovell/GENESPACE@dev")
library(GENESPACE)
```

### 1.2 Convert a previous run and re-run GENESPACE

If you have a previous GENESPACE run, you can convert the input directly to v1, then run the entire pipe:

First, do the conversion

```
wd <- "/place/to/store/new/run"
convert_input2v1(
  v1Dir = wd, 
  existingDir = "/path/to/existing/run/workingDirectory")
```

Then intialize the run ... lots of parameters to customize here if you'd like

```
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = "~/Desktop/programs/MCScanX/")
```

The the whole pipe (including riparian plotting and pan-genome annotation generation) can be run as:

```
gpar <- run_genespace(gsParam = gpar)
```

### 1.3 If you are starting from scratch ... 

Previously GENESPACE operated on raw gff and peptide annotations. The internal conversion to a simpler format caused most of the reported errors. To simplify this and reduce errors, GENESPACE now only takes processed .bed and .fasta annotations. See below, but in short, these are standard 4-column bed files where the 4th column in the gene ID, which exactly matches the header in the fasta file. 

Before running GENESPACE you need to either build these files yourself (and store them in /bed and /peptide) subdirectories. Alternatively, with the right parameterization, `parse_annotations()` can do the formatting for you.

Below is an example of how to run the full GENESPACE pipeline between human and chicken genomes from NCBI:

First off, make a directory to store the raw genome annotation files

```
genomeRepo <- file.path("~/Downloads/testGenespace/rawRepo")
if(!dir.exists(genomeRepo))
  dir.create(genomeRepo)
dir.create(file.path(genomeRepo, "human"))
dir.create(file.path(genomeRepo, "chicken"))
```

Then download human annotation from NCBI

```
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
  destfile = file.path(genomeRepo, "human", "GCF_000001405.40_GRCh38.p14_genomic.gff.gz"))
```

And chicken from NCBI

```
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz"))
download.file(
  url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz",
  destfile = file.path(genomeRepo, "chicken", "GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.gff.gz"))
```

Then parse the annotations into the simplified GENESPACE bed and fasta. See 4.3 (automatic parsing of annotations) for more details. 

```
wd <- "~/Downloads/testGenespace/verts"
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo, 
  genomeDirs = c("human", "chicken"),
  genomeIDs = c("human", "chicken"),
  presets = "ncbi", 
  genespaceWd = wd)
```

Now, we can initialize the GENESPACE run

```
gpar <- init_genespace(
  wd = wd, 
  path2mcscanx = "/path/to/MCScanX/")
```

Finally, with all the parameters set and the annotations formatted, we can run the entire GENESPACE pipe with one function call:

```
gpar <- run_genespace(gsParam = gpar)
```

## 2. Changes since v0.9.4 ...

### 2.1 major changes!!

GENESPACE is now structured differently. This breaks backwards compatibility to < v1, but was required to reduce the issues with raw data input. Some big changes:

  1. The user must now specify pre-formatted bed files instead of raw gff3 formatted annotations. This allows GENESPACE to catch issues with data formatting before running orthofinder etc. We are truly sorry to alter the specification and input so much, but it is unavoidable -- we want GENESPACE to integrate with existing bioinformatic pipelines, and too much internal customization was required, which produced hard-to-resolve errors.
  2. Resolved issues with plot visualization in interactive R sessions: xy dotplots and riparian plots are no longer written to the plot window. These are now saved as vector pdfs in the new /dotplot and /riparian directories respectively.
  3. Better organization of results output. Previously all output went into /results. /results now only contains the parsed orthofinder results and annotated bed files. All other results are placed in /syntenicHits, /dotplots, /riparian, and /pangenome. 
  4. All GENESPACE functions are now integrated into `run_genespace()`, which completes all steps internally (the user can still run each step independently too). 
  5. The only change to the actual algorithm is a modification to how large tandem arrays are handled and masked. This will only affect a small number of genomes. 

### 2.2 minor changes

  1. code cleanup and better documentation, both in line and in help files
  2. much better data QC and handing of errors in parameter specification
  3. ggplot2-integrated graphics and R object output for potential downstream interactive viz.
  4. internal position calculation for riparian plots and inferred syntenic positions transferred to a new function `integrate_synteny()`. This speeds up riparian plotting by 10x or more. 
  5. `plot_riparian()` function generalization that allows the user to specify a block coordinate matrix and nothing else (allowing plotting of data from other packages). 
  6. replace `orthofinderMethod = "fast"` with `onewayOF = TRUE` to match the new `orthofinder -1` specification. 
  7. a new function `query_genespace()`, allows positional or orthology based queries of pan-genome, pairwise hits, or raw orthogroups/orthologues. 

### 2.3 upcoming changes for v1.1.0 (planned for spring or summer 2023)

  1. An executable Rscript for `run_genespace` and parameterization via .json to allow GENESPACE to be called directed from the shell 
  2. `synteny()` will be able to be called from the shell on a single blast8-formatted text file
  3. Zoom-in on specific regions of the riparian plot

## 3. Installation

Installation for v1 is identical to v0.9.4. Use the same conda (or other) installation environment as before. To install from scratch, you need R, a few 3rd party programs and a few R packages. Here is how to install these:

### 3.1 Install R

GENESPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

### 3.2 Install third party programs

GENESPACE also requires third-party software that can be installed as follows.

`OrthoFinder` (which includes `DIAMOND2`) is most simply installed via conda (in the shell, not R). 

```
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources. See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

### 3.3 Open an interactive R session

You can open R using the the GUI, Rstudio, or as an interactive session in the terminal. If you are using a conda environment or specifying the path to orthofinder in the `$PATH`, you need to enter an R interactive session from that shell environment, which can be done by either calling `R` (for command line interface) `open -na rstudio` if using Rstudio. 

### 3.4 Install GENESPACE R package

Once in R, GENESPACE can be installed directly from github via:

```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")
```

### 3.5 Install R dependencies

The above command will install the CRAN-sourced dependencies (`data.table`, `dbscan` and `R.utils`). The bioconductor dependencies (`rtracklayer` and `Biostrings`) need to be installed separately via:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```

## 4. Getting your annotation files formatted

This can be the trickiest part of running GENESPACE. If you are having trouble, open an issue. 

### 4.1 GENESPACE-readable annotation format

For each genome, GENESPACE needs .bed formatted coordinates of each gene and peptides that exactly match the "name" (4th) bed column. 

The bed file needs to have the first four tab-separated fields (chromosome ID, start, end, name) in a bed file. Other fields are allowed, but will be ignored by GENESPACE. A header is not allowed. For example:

```
chr1  1500  2000  gene1
chr1  2500  2200  gene2
```

The fasta file is just the standard format with the header exactly matching the names in the bed file. The order of the sequenced does not matter, but each header must uniquely match the names in the bed file.

```
>gene1
MGASGRGAGEQQSPSSTG
>gene2
MLVMSECKGRDRSPSSSM
```

### 4.2 Manual parsing of raw annotation files

It's possible to format all files using custom scripts. If you opt for this route, place the genespace-formatted annotations in your working directory. For example for two genomes ("genome1", "genome2")

```
/workingDirectory
└─ peptide
    └─ genome1.fa
    └─ genome2.fa
└─ bed
    └─ genome1.bed
    └─ genome2.bed
```

### 4.3 Getting and storing raw annotation files

Its usually preferable to have a static repository containing the raw versions of all of your genome annotations, lets call this directory `/genomeRepo`. Each unique genome annotation should be stored in its own subdirectory with an informative name. Here, each genome is given a species id, genotype id and genome version. The annotation files of interest are the gff3-formated protein coding gene features and the primary transcript translated peptides. These are the two most common types of annotation files available. Download them from your favorite repositories and place them in a directory structure like below:

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

The names within each genome directory can be anything as long as they have file names that can be found through regular expression matching (see `parse_annotations(gffString, faString)`).

### 4.4 Partially-automated annotation parsing

Assuming the directory structure from 4.3, we want to parse the raw annotations into the directory structure outlined in 4.2. The convenience function `parse_annotations` does this. Lets take the simple case from 1.3 where we have translatedCDS and gene.gff3 files from NCBI. In this case, we can specify the "ncbi" present which will automatically parse these files correctly:

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

### 4.5 Non-standard annotation formats

Most repositories have their own annotation formats. You may want to combine annotations from NCBI with those from phytozome (plant genome repo) or other non-standard repositories. To do this, you simply re-run parse_annotations with different genomes. For example add a phytozome-formatted annotation:

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

The add a non-standard annotation where there is only a simple geneID header on the peptide and the gene names are in the "GeneID" ';'-delimited 9th gff column. 

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


There are lots of parameters that can be combined to parse every annotation format. See the help file for details. The main parameters to get your annotation parsed correctly are: 

  - `gffIdColumn`: the name (or index) of the ';' separated field in the attributes column contains the geneID. 
  - `headerEntryIndex`: the name (or index) of the `headerSep`-separated field that contains the geneID in the fasta header
  - `gffStripText`: character string to strip off the geneID in the gff3 attribute column
  - `headerStripText`: character string to strip off the geneID in the fasta header

Parsing can take some troubleshooting, which can be aided by setting `troubleshoot = TRUE`. This prints the first 10 lines of the raw and parsed gff and fasta headers. 


## 5. Setting and checking GENESPACE parameters

### 5.1 Initalizing the GENESPACE run with defaults

The only inputs required for genespace are a directory containing matched annotations in /peptide and /bed subdirectories and the path to the MCScanX installation.

```
gpar <- init_genespace(
    wd = "/path/to/GENESPACE/workingDir", 
    path2mcscanx = "/path/to/MCScanX")
```

This initializes a run where each genome is named following the file names in the peptide and /bed subdirectories and assumes that all genomes are haploid. 

### 5.2 Specifying custom paths to diamond and orthofinder

Installation of orthofinder (and diamond) represents one of the most common issues getting GENESPACE up and running. If you cannot install into a conda environment and open R from there (or put orthofinder in the path), you will to specify a custom path to the orthofinder and diamond program installations through `init_genespace(..., path2orthofinder = "path/to/orthofinder", path2diamond = "path/to/diamond")`. 

### 5.3 Specifying ploidy

Genome ploidy represents the most common required customization. If all of your genomes have the same ploidy, they can be set by specifying a single value. For example for all diploid genomes (likely representing inbred tetraploid genomes) `init_genespace(..., ploidy = 2)`. 

However, you may only have a couple higher ploidy genomes. In this case, its preferable to specify both the genomeIDs and ploidy to ensure you have the right ploidy matching the right genome. For example: `init_genespace(..., ploidy = c(1, 2, 1), genomeIDs = c("species1", "species2", "species3"))`, where two genomes are haploid and one is an outbred diploid. 

### 5.4 Changing other parameters.

There are many ways to customize the GENESPACE run. Check the helpfile for `init_genespace` for details about them all. Here are the other most common:

- `useHOGs` - should hierarchical orthogroups be used instead of global orthogroups. HOGs tend to split arrays and better reflect orthologs. Global orthogroups better capture paralogous regions. Default (NA) chooses whichever matches the ploidy of the genomes best. 
- `outgroup` or `ignoreTheseGenomes` - genomeIDs that should only be used in the orthofinder run but not considered for synteny etc. If not specified, all genomes are used. `outgroup` is a bit of a missnomer and the argument is kept here for now for better backwards compatibility to <v1. 
- `blkSize` - default is 5. This is the minimum number of anchored blast hits needed to call a syntenic block. 
- `nGaps` - default is 5. Passed to `MCScanX -m` as the maximum number of gaps in a syntenic block. 
- `synBuff` - default is 100. This is the maximum euclidean distance (in gene-rank order units) from a syntenic anchor for a gene to be considered in the syntenic buffer. Also the minimum distance between two genes to split a tandem array into two clusters.   
- `nCores` - the number of parallel processes to run
- `orthofinderInBlk` - should OrthoFinder be re-run within each large syntenic region? Default is TRUE if any genome ploidy > 1. 



