Assignment 1
========================================================

### Intro ###

This analysis is based on publicly-available expression study of mouse brain tissue with single gene mutation (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7191). 

S1P2, when mutated, results in seizures. While mutated S1P3, a related gene, does not. In this study, gene expression of two brain regions (hippocampus and neocortex) from three mouse strains (wild type, S1P2 mutant, and S1P3 mutant) are measured. Additional information of gender of the mice and processing date is available in the data being analyzed.

Expression is measured on the Affymetrix MGU74Av2 platform.

```r
library(plyr)
library(ggplot2)
library(xtable)
```







```r
design <- read.table("../data/mouseBrain/GSE7191-design.txt", row.names = 1, 
    header = TRUE)
expDat <- read.table("../data/mouseBrain/GSE7191-data.txt", row.name = "probe", 
    header = TRUE)
```


#### Q1) Basic Characteristics ####

a) Dimensionality of the data

```r
dim(expDat)  # [probes, samples]
```

```
## [1] 12422    50
```

```r
dim(design)  # [samples, design factors]
```

```
## [1] 50  4
```

```r
fa <- colnames(design)  # design factors
```


b) Factors in experimental design

```r
with(design, table(Genotype))
```

```
## Genotype
##   S1P2_KO   S1P3_KO Wild_type 
##        20        10        20
```

```r
with(design, table(BrainRegion))
```

```
## BrainRegion
## hippocampus   neocortex 
##          25          25
```

```r
with(design, table(Sex))
```

```
## Sex
## female   male 
##     26     24
```

```r
with(design, table(DateRun))
```

```
## DateRun
## 01/16/04 03/11/04 07/23/04 08/14/03 08/21/03 09/11/03 10/23/03 12/18/03 
##        7        4        4        8        8        7        7        5
```



```r
# maxLevel <- max(laply(design, nlevels))
faCounts <- ldply(design, function(x) {
    faFeq <- table(x)
    iFaCounts <- data.frame(faFeq)
    colnames(iFaCounts) <- c("factor.level", "counts")
    return(iFaCounts)
})
html_print(faCounts)
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Fri Feb 21 17:24:07 2014 -->
<TABLE border=1>
<TR> <TH> .id </TH> <TH> factor.level </TH> <TH> counts </TH>  </TR>
  <TR> <TD> DateRun </TD> <TD> 01/16/04 </TD> <TD align="right"> 7 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 03/11/04 </TD> <TD align="right"> 4 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 07/23/04 </TD> <TD align="right"> 4 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 08/14/03 </TD> <TD align="right"> 8 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 08/21/03 </TD> <TD align="right"> 8 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 09/11/03 </TD> <TD align="right"> 7 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 10/23/03 </TD> <TD align="right"> 7 </TD> </TR>
  <TR> <TD> DateRun </TD> <TD> 12/18/03 </TD> <TD align="right"> 5 </TD> </TR>
  <TR> <TD> Genotype </TD> <TD> S1P2_KO </TD> <TD align="right"> 20 </TD> </TR>
  <TR> <TD> Genotype </TD> <TD> S1P3_KO </TD> <TD align="right"> 10 </TD> </TR>
  <TR> <TD> Genotype </TD> <TD> Wild_type </TD> <TD align="right"> 20 </TD> </TR>
  <TR> <TD> BrainRegion </TD> <TD> hippocampus </TD> <TD align="right"> 25 </TD> </TR>
  <TR> <TD> BrainRegion </TD> <TD> neocortex </TD> <TD align="right"> 25 </TD> </TR>
  <TR> <TD> Sex </TD> <TD> female </TD> <TD align="right"> 26 </TD> </TR>
  <TR> <TD> Sex </TD> <TD> male </TD> <TD align="right"> 24 </TD> </TR>
   </TABLE>

You can also embed plots, for example:


```r
plot(cars)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


