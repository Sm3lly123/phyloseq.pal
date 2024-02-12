# phyloseq.pal

This README is a WIP!

## Overview

`phyloseq.pal` provides a collection of functions to simplify plotting taxa abundance, conducting alpha diversity analysis, and performing statistical tests on richness metrics with phyloseq objects.

## Installation

```         
devtools::install_github("ellenelizabethsmith/phyloseq.pal")
```

## Functions

1.  `getgloms`. Useful if you will be using agglomerated objects many times in your analysis! Returns a named list of phyloseq objects.

```         
# Example Usage:
library(phyloseq.pal)
library(phyloseq)
data(GlobalPatterns)
gloms <- getgloms(GlobalPatterns, ranks = c("Phylum", "Class", "Order", "Family", "Genus","Species")) 
```

2.  `plot_taxa_abundance` To simplify generating stacked bar plots to visualize the abundance of taxa in a phyloseq object. It includes options for ordering taxa by abundance and plotting relative or absolute abundance. You can add one more more faceting variables seperated by a `+`

```         
#Pre-agglomerated object:
plot_taxa_abundance(gloms$Genus, rank = "Genus", x = "X.SampleID", wrap = "SampleType", n = 20, byabundance = TRUE, abs = FALSE, size = 10)

#Non-agglomerated
plot_taxa_abundance(GlobalPatterns, rank = "Family", x = "X.SampleID", wrap = "SampleType", n = 10, byabundance = TRUE, abs = TRUE, size = 10)
```

3.  `get_percentage_assigned` Returns a table with the percentage of ASVs which have a taxonomic assigned at each rank.

```         
get_percentage_assigned(GlobalPatterns)
```

4.  `ps_from_ampliseq`For creating a phyloseq object from the dada2 output directory from [nf-core/ampliseq](https://nf-co.re/ampliseq/2.8.0). Note, newer versions of ampliseq may change file names and break this. Metadata addition is optional.

```         
metadata <- read.csv("my_metadata.csv")
ps <- ps_from_ampliseq("path/to/top/level/outfolder",metadata)
```

5. `do_microshades_plot` Convenience wrapper for using [KarstensLab/microshdes](https://github.com/KarstensLab/microshades) with a phyloseq object. Uses top 5 most abundant in "rankhigh". Facet is optional.

```
do_microshades_plot(ps = GlobalPatterns, rankhigh = "Phylum", ranklow = "Genus", x = "SampleID", facet = SampleType)
```
