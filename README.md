# cscplots

### Overview

**cscplots** is an R package providing a library of functions for use as 
pipeline tools for identifying and analysing cancer stem cells in single-cell 
RNA-Seq datasets. 

Tools for plots to identify putative cancer stem cells from single-cell scRNA-Seq data in expression matrix format e.g. containing read counts or TPMs, and to analyse the cancer stem cells using preferential expression analysis. This will assist in finding transcription driving the tumours that the cell populations are sampled from, e.g. drug resistance, immune evasion, and other molecular features underlying tumour evolution and hardiness.

The idea behind **cscplots** is to provide a set of functions that can be used 
together within R to create a pipeline that takes as input an expression matrix, 
and through successive steps identify and then analyses cancer stem cells. 
Three demonstration pipelines are provided in tutorials 1, 2 and 3, to 
familiarise the use with the steps necessary, and the additional code that is 
used to join the analysis and plotting together in the necessary way.

### Getting started 

If you are just getting started with **cscplots** it is recommended to start 
with the tutorial [vignettes](https://p-carter.github.io/cscplots/articles/index.html).

### Resources

* [p-carter.github.io/cscplots](https://p-carter.github.io/cscplots) (online documentation, vignettes)
* [Open an issue](https://github.com/p-carter/cscplots/issues) (GitHub issues for bug reports, feature requests)

### Installation

* Install latest development version from GitHub (requires [devtools](https://github.com/hadley/devtools) package):

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("p-carter/cscplots", dependencies = TRUE, build_vignettes = FALSE)
```

This installation won't include the vignettes (they take some time to build), but all of the vignettes are 
available online at [p-carter.github.io/cscplots/articles](https://p-carter.github.io/cscplots/articles/).

## Pipeline schema

![alt text](https://github.com/p-carter/cscplots/blob/master/cscplots_pipeline.jpg?raw=true)
