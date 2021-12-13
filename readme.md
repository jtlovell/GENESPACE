NOTE: **Significant updates to beta v0.9.1 (at /dev as of 12:09MT 13-Dec 2021).**
Preparing source code at /dev for release as v0.9.1. Updates here will be documented upon final release by 16-Dec 2021

NOTE: **Significant updates to beta v0.8.6b (pushed 8:00pm USA MT 23-November 2021).** *If you are running within-block orthofinder, please feel free to test out v0.8.6b. If you find problems, please open an issue on github. Either way, we highly suggest re-running after upgrading to v0.8.6 once the stable version is released in early December 2021.* Updates include:

1. Iterative array parsing - eliminates nearly-syntenic duplicates in the pangenome and improves synteny-constrained and within-block orthogroup assignment.
2. Batch parallelization of synteny. Modest improvement of speed in single threaded mode. Significant and linear improvement in speed with multi-threaded mode. 
3. Smaller file size and clearer dotplot graphics
4. More streamlined functions: easier to follow help files and documentation. 



For a full tutorial, see: https://github.com/jtlovell/GENESPACE/blob/dev/doc/genespaceOverview.pdf

## Overview
GENESPACE is a comparative genomics framework implemented in the R environment for statistical computing. The premise is that, when analyzing high-quality genome assemblies and annotations, we can improve the confidence of evolutionary inference by combining two sources of evidence for common ancestry: synteny (i.e. collinearity of gene order) and coding sequence similarity (homology). In addition to providing a second line of evidence beyond sequence similarity, combining synteny and homology have several benefits:

- exclude paralogous regions
- control for variable ploidy among genomes
- search for orthologs within single-copy regions
- infer expected gene position across genomes: a priori hypothesis about presence-absence variation
- track genes in regions of interest (e.g. QTL intervals)

GENESPACE outputs a synteny-constrained and -anchored orthogroup pan-genome annotation among multiple genomes. This simple text file allows for extraction and exploration of regional gene-level variation, a necessary step to integrate comparative and quantitative genomic goals. 

For details and a complete example run, see the [vignette](https://github.com/jtlovell/GENESPACE/blob/dev/doc/genespaceOverview.pdf) and source code in `/doc` directory. 

## Installation

#### Install required programs

GENESPACE requires third-party software can be installed as follows.

`orthofinder` (which includes `diamond`) is most simply installed via conda (in the terminal, not R). 

```{bash, eval = FALSE}
conda create -n orthofinder
conda activate orthofinder
conda install -c bioconda orthofinder 
```

If conda is not available on your machine, you can install orthofinder from a number of other sources.  See [orthofinder documentation](https://github.com/davidemms/OrthoFinder) for details.  

`MCScanX` should be installed from [github](https://github.com/wyp1125/MCScanX). 
**NOTE** Due to nested dependencies in the `orthofinder` install, `orthofinder` must be in the path when called, while the file path to `MCScanX` does not need to be in the path and can be specified in R. 

#### Install GENESPACE

If you are planning to run orthofinder from within R, which is needed if using orthofinderInBlk (see below for details) and recommended when using any genomes with ploidy > 1, enter R from the terminal with either `R` (for command line interace), or if Rstudio is installed on your machine `open -na rstudio`. 

Once in R, GENESPACE can be installed directly from github via:

```{r, eval = FALSE}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE", upgrade = F)
```

#### Install R dependencies

The above command will install the CRAN-sourced dependencies (`data.table`, `dbscan` and `R.utils`). The bioconductor dependencies (`Rtracklayer` and `Biostrings`) need to be installed separately via:

```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "Rtracklayer"))
```


## Legal

GENESPACE R Package (GENESPACE) Copyright (c) 2021,HudsonAlpha Institute for Biotechnology. All rights reserved.

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
