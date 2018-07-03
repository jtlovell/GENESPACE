# Overview

A common goal in comparative genomics is to compare gene sequences of many species with reference to the ancestral state. While simple in highly conserved, single copy genes and dipliod/haploid organisms, such alignments are tricky when extending to the biological complexity of most genomes. 

Here, we develop a pipeline to accomplish alignemnts in "gene space" (within, or near gene models) to gain a more complete picture of genomic polymorphism among species. The approach can handle any ploidy, ancient or recent whole-genome duplicates, inversions (and other small-scale rearangements) and most other known complexity. Such flexibility is affordered by accomplishing basic comparative genomes tasks within collinear blocks. Once the blocks are constructed, any basic task can be accomplished, knowing that one is comparing evolutionarily related, diploid sequences. 

## The workflow employed here follows:

1. Conduct all pairwise and intra-specific peptide-peptide BLAST alignments and pipe these into the program orthofinder. Orthofinder extracts networks of genes that appear to share an evolutionary origin. Cull BLAST results to genes that are found in orthofinder groups and exist in dense patches in physical genomic space. 
2. Feed highly-culled BLAST results into the multiple-colinear block generation program MCScanX. Parse the block coordinates by position, aiming to retain larger, denser collinear blocks. 
3. Re-run #1, within each block, finding the best network of genes. For those with >0 genes in the network, output results as before. For those without an orthologous sequence in the block, blat against the genomic reference and find the best hit. If above specified score thresholds, refine the sequence via blast, then exonerate. Write new CDS to a new fasta file. 
4. After compiling all CDS fastas (both real and found by exonerate), re-run #1-2. Finally, parsing orthofinder alignments and returning these to file. 

# Getting started

## Installation of dependencies

The following programs must be installed and added to your path. All except MCScanX are available on conda

- MCScanX
- Diamond
- orthofinder
- exonerate
- mafft
- bedtools

Also a few R packages are also needed:

- data.table
- dbscan
- Biostrings (from bioconductor)

Make sure that these are installed by opening R and running:

`R check_environment()`

## Input files

The program requires four directories to be populated with properly formatted files:

`peptide.dir` the directory containing fasta files with the translated peptides for each gene
`cds.dir` the directory containing fasta files with the coding sequences for each gene
`gff.dir` the directory containing gene annotations in gff3 format
`assembly.dir` the genomic DNA sequence assembly for each genome

These all must be named according to the species, so A. thaliana will have:
Athaliana.fa, Athaliana.fa, Athaliana.gff3, Athaliana.fa

## Running it

If running the pipeline (as a single step), one can either use the R interface, or call the program from the command line. 

`compareGeneSpace --genome.dir "~/Downloads/genome" `



Since the pipeline works conditional on previous steps, parallelization must be done either:

1. Stepwise, so each step is parallelized, or
2. Run on a single node. 

For 5 <1.5Gb genomes, the whole program takes ~4h on an 8 core node, so it is not overly slow. However, the time increases exponentially as the number of genomes increases. If running on a single machine, consider analyzing <=8 genomes at a time. Say, if 32 genomes are to be compared, one could use a single genome as the tester and run 8 sets of 5-genome comparisons. Then, concatenate the results at the end. 








Here we extend basic annotation-annotation alignments by searching for un-annotated sequences
that appear to be orthologous to CDS regions annotated in related species. To accomplish this goal, 
we constrain our search within collinear blocks that share an evolutionary origin between species. We search
for, then add the un-annotated sequences into the fasta files, then fina

We start by employing the standard blast-orthofinder workflow to find pairwise orthology networks between species.
