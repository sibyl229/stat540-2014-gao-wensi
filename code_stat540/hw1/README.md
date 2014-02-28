Homework 1 README
========================================================

### Intro

This analysis is based on publicly-available expression study of mouse brain tissue with single gene mutation (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7191). 

S1P2, when mutated, results in seizures. While mutated S1P3, a related gene, does not. In this study, gene expression of two brain regions (hippocampus and neocortex) from three mouse strains (wild type, S1P2 mutant, and S1P3 mutant) are measured. Additional information of gender of the mice and processing date is available in the data being analyzed.

Expression is measured on the Affymetrix MG_U74Av2 platform.

#### Report
[Report quick preview](assign1.md)
[Report Source](assign1.rmd)
[Report preview (won't work on private repo)](http://htmlpreview.github.io/?https://raw.github.com/sibyl229/stat540-2014-gao-wensi-hw/master/code_stat540/hw1/assign1.html)

#### Data
[Data](../../data/mouseBrain/)

### Dependencies

* plyr
* ggplot2
* lattice
* xtable
* RColorBrewer
* gplots
* preprocessCore
* **reshape**
* limma
* **MIfuns**