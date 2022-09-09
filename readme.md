###################

**Version 1 of GENESPACE is coming soon** Check back for updates in the fall 2022. 
V1 will include the following changes:

- Use of bed files instead of gff3 (will still have accessory functions for conversion)
- ggplot2 graphics and options to output to plotly for interactive viz
- better error catching and checking for genome quality
- optional use of phylogenetically hierarchical orthogroups 
- optional use of within-block orthologs (instead of orthogroups)
- slightly modified tandem-array collapsing, to better integrate with synteny

Thanks to the users who offered suggestions (and found bugs)!

###################

This is v0.9.3 of the GENESPACE R package. This is new software, which we are actively working to make more user friendly. Please open an issue above or email John directly (jlovell[at]hudsonalpha[dot]org) if you run into problems or the help files are not sufficient. Thanks for using GENESPACE! 

[**Click here for a more detailed vignette**](https://htmlpreview.github.io/?https://github.com/jtlovell/GENESPACE/blob/master/doc/genespaceOverview.html) including pipeline details, output descriptions, etc. Below is a quick overview of installation and running the package.

**NOTE!** The hardest part of GENESPACE is getting the dependencies installed and your annotations formatted correctly. Don't be dissuaded if you have having trouble with these steps! Open an issue and we can help!

## 1. Overview

GENESPACE conducts synteny- and orthology-constrained comparative genomics, which can be used to (1) make multi-genome graphical comparisons, (2) define syntenic block breakpoints and (3) build a pan-genome annotation. GENESPACE outputs a synteny-constrained and -anchored orthogroup pan-genome annotation among multiple genomes. This simple text file allows for extraction and exploration of regional gene-level variation, a necessary step to integrate comparative and quantitative genomic goals.

## 2. Installation

#### 2.1 Install R

GENESPACE is meant to be run interactively in the R environment for statistical computing. So, you need to have R installed. See [CRAN](https://www.r-project.org/) for the most recent release. 

#### 2.2 Install third party programs

GENESPACE also requires third-party software that can be installed as follows.

`orthofinder` (which includes `diamond`) is most simply installed via conda (in the shell, not R). 

```{bash, eval = FALSE}
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources.  See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 

#### 2.3 Open an interactive R session

You can open R using the the GUI, Rstudio, or as an interactive session in the terminal. If you are using a conda environment or specifying the path to orthofinder in the `$PATH`, you need to enter an R interactive session from that shell environment, which can be done by either calling `R` (for command line interface) `open -na rstudio` if using Rstudio. 

#### 2.4 Install GENESPACE R package

Once in R, GENESPACE can be installed directly from github via:

```{r, eval = FALSE}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE", upgrade = F)
```

#### 2.5 Install R dependencies

The above command will install the CRAN-sourced dependencies (`data.table`, `dbscan` and `R.utils`). The bioconductor dependencies (`rtracklayer` and `Biostrings`) need to be installed separately via:

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))
```

## 3. Run GENESPACE on the tutorial data

#### 3.1 Get started

First, in R, require the GENESPACE package and make a file path to a directory where you want to put your run. The directory should either be empty or non-existent (GENESPACE will make it for you, assuming that the parent directory exists).  

```{r}
library(GENESPACE)
runwd <- file.path("~/Desktop/testGenespace")
```

#### 3.2 Get some example data

To illustrate all steps of the pipeline, lightly subset NCBI-formatted annotations of human/chimpanzee chromosomes 3-4 and rhesus chromosomes 2 & 5 are provided in the extData of the GENESPACE R package. These data can be added to the above directory with the correct subdirectory structure for GENESPACE via:

```{r}
make_exampleDataDir(writeDir = runwd)
```

This creates a subdirectory called `/rawGenomes`. For downstream flexibility (e.g. multiple genome versions for one species, metadata or assembly data, etc), the raw genome directory structure follows: `/rawGenomes/[speciesID]/[versionID]/annotion`. When working with your own data, place the raw annotation files in this same directory structure with separate directories for each species, separate subdirectories for each genome version, and the annotation files in a subdirectory called "annotation".

#### 3.3 Initialize the GENESPACE run

All elements of GENESPACE require a list of parameters, specified to functions as `gsParam`. This contains paths to files, working directories, program executables the basic parameterization of the run. See the [**vignette**](https://htmlpreview.github.io/?https://github.com/jtlovell/GENESPACE/blob/master/doc/genespaceOverview.html) for argument descriptions. 

```{r}
gids <- c("human","chimp","rhesus")
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1,3),
  wd = runwd, 
  gffString = "gff", 
  pepString = "pep",
  path2orthofinder = "orthofinder", 
  path2mcscanx = "~/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
```

#### 3.4 Format the raw annotations for GENESPACE

To ensure accurate merges between peptide fasta headers and gff3 annotations, GENESPACE uses a simplified gff3-like text file with a column `id` that exactly matches the peptide fasta header. GENESPACE has built in functions to parse NCBI (`parse_ncbi`) and phytozome (`parse_phytozome`); however, if using non-standard formatted annotations, this can be a tricky step. While the example data was originally downloaded from NCBI, we stripped much of the NCBI formatting to make the files smaller. They now can be parsed with the generic `parse_annotations`:

```{r parse annotations}
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=")
```

#### 3.5 Run orthofinder

GENESPACE requires orthofinder to be run. Here, since orthofinder is in the path, we can run it straight from R, using the 'fast' search method

```{r orthofinder}
gpar <- run_orthofinder(gsParam = gpar)
```

#### 3.6 Run the GENESPACE synteny search

The main engine of GENESPACE is `synteny`. This is a complex function that parses pairwise blast hits into syntenic regions and blocks. See the [**vignette**](https://htmlpreview.github.io/?https://github.com/jtlovell/GENESPACE/blob/master/doc/genespaceOverview.html) for details. Here, we will just use the defaults:

```{r synteny}
gpar <- synteny(gsParam = gpar)
```

This function populates the results directory directory with dotplots and annotated blast hits. 

#### 3.7 Make multi-species synteny plots

GENESPACE visualizes multi-species synteny with a 'riparian' plot. The default specification orders chromosomes by maximum synteny to a reference genome and colors the synteny by their reference chromosome of origins. There are lots of ways to customize this plot. 

```{r riparian, fig.width=5}
ripdat <- plot_riparianHits(gpar)
```

#### 3.8 Build a pangenome annotation

The main output of GENESPACE is a synteny-anchored pan-genome annotation, where every unique synteny-constrained orthogroup is placed in a position on the reference genome gene order. This is constructed by `pangenome`. See the [**vignette**](https://htmlpreview.github.io/?https://github.com/jtlovell/GENESPACE/blob/master/doc/genespaceOverview.html) for more details. 

```{r build pangenome}
pg <- pangenome(gpar)
```


## Legal

GENESPACE R Package (GENESPACE) Copyright (c) 2021, HudsonAlpha Institute for Biotechnology. All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.

NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative  works, and perform publicly and display publicly, and to permit others to do so.


*** License Agreement ***

GENESPACE R Package (GENESPACE) Copyright (c) 2021, HudsonAlpha Institute for Biotechnology. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the HudsonAlpha Institute for Biotechnology nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.
