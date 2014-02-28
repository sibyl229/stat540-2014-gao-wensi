Assignment 1
========================================================

Wen Si (Sibyl) Gao

### Intro ###

This analysis is based on publicly-available expression study of mouse brain tissue with single gene mutation (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7191). 

S1P2, when mutated, results in seizures. While mutated S1P3, a related gene, does not. In this study, gene expression of two brain regions (hippocampus and neocortex) from three mouse strains (wild type, S1P2 mutant, and S1P3 mutant) are measured. Additional information of gender of the mice and processing date is available in the data being analyzed.

Expression is measured on the Affymetrix MG_U74Av2 platform.

```r
library(plyr)
library(ggplot2)
library(lattice)
library(xtable)
library(RColorBrewer)
library(gplots, warn.conflicts = FALSE)
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
```

```r
library(preprocessCore)  # to help normalization
library(reshape, warn.conflicts = FALSE)
library(limma)
library(MIfuns, warn.conflicts = FALSE, quietly = TRUE)  # to help print ftable
```

```
## MIfuns 5.1
```









```r
design <- read.table("../../data/mouseBrain/GSE7191-design.txt",
                     row.names=1, header=TRUE)
expDat <- read.table("../../data/mouseBrain/GSE7191-data.txt",
                     row.name="probe", header=TRUE)
design$Sid <- factor(rownames(design))

# order expression data columns (samples) by the rows in the design table
expDat <- expDat[, row.names(design)]
```



```r
# Reorder Genotype factor levels so that wildtype has level 1
design$Genotype <- 
  factor(design$Genotype, levels=c("Wild_type","S1P2_KO","S1P3_KO"))
```



```r
# create new factor for DateRun
# levels are order by increasing date, and labeled by "Day-1", etc
newDateRun <- with(design, {
  ## Reorder DateRun levels
  dates <- as.Date(levels(DateRun),"%m/%d/%y")
  newDateLevels <- levels(DateRun)[order(dates)]
  nfDateRun <- factor(DateRun, newDateLevels)
  
  ## Re-label DateRun level labels
  levels(nfDateRun) <- paste("Day", 1:nlevels(DateRun), sep="-") 
  return(nfDateRun)
})
design$DateRun <- newDateRun

# more informative column names for expDat
colnames(expDat) <- 
  with(design, paste(DateRun, 
                     substr(Genotype,1,4), 
                     substr(BrainRegion,1,4), 
                     substr(Sex,1,1), Sid,
                     sep="_"))

# define sample order
sampleOrd <- with(design, order(DateRun, Genotype, BrainRegion, Sex))
expDat <- expDat[,sampleOrd]
design <- design[sampleOrd,]
```



### Q1) Basic Characteristics ###

#### a) Dimensionality of the data

```r
dim(expDat)  # [probes, samples]
```

```
[1] 12422    50
```

```r
dim(design)  # [samples, design factors]
```

```
[1] 50  5
```

```r
fa <- colnames(design) # design factors
```



#### b) Factors in experimental design

Looking at each factor individually, samples are mostly evenly distributed across different levels of a factor. The exception is that there is smaller sample size for S1P3 knockout group.


```r
with(design, table(Genotype))
```

```
Genotype
Wild_type   S1P2_KO   S1P3_KO 
       20        20        10 
```

```r
with(design, table(BrainRegion))
```

```
BrainRegion
hippocampus   neocortex 
         25          25 
```

```r
with(design, table(Sex))
```

```
Sex
female   male 
    26     24 
```

```r
with(design, table(DateRun))
```

```
DateRun
Day-1 Day-2 Day-3 Day-4 Day-5 Day-6 Day-7 Day-8 
    8     8     7     7     5     7     4     4 
```






Cross-tabulating genotype and brain region also finds no issue with the factorial experimental design. As sample sizes are relatively even across all genotype-brain region cobminations (besides the issue of fewer S1P3 knockouts):


```r
x <- with(design, table(BrainRegion, Genotype))
html_print(addmargins(x))
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Fri Feb 28 02:18:40 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH> Wild_type </TH> <TH> S1P2_KO </TH> <TH> S1P3_KO </TH> <TH> Sum </TH>  </TR>
  <TR> <TD align="right"> hippocampus </TD> <TD align="right"> 10 </TD> <TD align="right"> 10 </TD> <TD align="right"> 5 </TD> <TD align="right"> 25 </TD> </TR>
  <TR> <TD align="right"> neocortex </TD> <TD align="right"> 10 </TD> <TD align="right"> 10 </TD> <TD align="right"> 5 </TD> <TD align="right"> 25 </TD> </TR>
  <TR> <TD align="right"> Sum </TD> <TD align="right"> 20 </TD> <TD align="right"> 20 </TD> <TD align="right"> 10 </TD> <TD align="right"> 50 </TD> </TR>
   </TABLE>


Even as the gender of the mouse is taken into account. The design seem okay:


```r
x <- with(design, table(Sex, BrainRegion, Genotype))
html_print(ftable(x))
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Fri Feb 28 02:18:40 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH>        </TH> <TH>             </TH> <TH> Genotype </TH> <TH> Wild_type </TH> <TH> S1P2_KO </TH> <TH> S1P3_KO </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> Sex    </TD> <TD> BrainRegion </TD> <TD>          </TD> <TD>           </TD> <TD>         </TD> <TD>         </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> female </TD> <TD> hippocampus </TD> <TD>          </TD> <TD>         5 </TD> <TD>       5 </TD> <TD>       3 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD>        </TD> <TD> neocortex   </TD> <TD>          </TD> <TD>         5 </TD> <TD>       5 </TD> <TD>       3 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> male   </TD> <TD> hippocampus </TD> <TD>          </TD> <TD>         5 </TD> <TD>       5 </TD> <TD>       2 </TD> </TR>
  <TR> <TD align="right"> 5 </TD> <TD>        </TD> <TD> neocortex   </TD> <TD>          </TD> <TD>         5 </TD> <TD>       5 </TD> <TD>       2 </TD> </TR>
   </TABLE>



Problem seem to appear when examining DateRun as a factor. DateRun could potentially confound both Genotype and gender of the mouse. It is often the case that a single day contains only samples of the same genotype and samples of the same gender. As a result, one may not know if different expression readings is due to different genotypes, or different dates the sample is run. It's the same case with gender of the mouse. If DateRun (the batch effect) is an important factor of the expression readings, one has to correct for this effect. 


```r
x <- with(design, table(Genotype, DateRun))
html_print(ftable(x))
```

<!-- html table generated in R 3.0.2 by xtable 1.7-1 package -->
<!-- Fri Feb 28 02:18:40 2014 -->
<TABLE border=1>
<TR> <TH>  </TH> <TH>           </TH> <TH> DateRun </TH> <TH> Day-1 </TH> <TH> Day-2 </TH> <TH> Day-3 </TH> <TH> Day-4 </TH> <TH> Day-5 </TH> <TH> Day-6 </TH> <TH> Day-7 </TH> <TH> Day-8 </TH>  </TR>
  <TR> <TD align="right"> 1 </TD> <TD> Genotype  </TD> <TD>         </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> <TD>       </TD> </TR>
  <TR> <TD align="right"> 2 </TD> <TD> Wild_type </TD> <TD>         </TD> <TD>     4 </TD> <TD>     8 </TD> <TD>     7 </TD> <TD>     0 </TD> <TD>     1 </TD> <TD>     0 </TD> <TD>     0 </TD> <TD>     0 </TD> </TR>
  <TR> <TD align="right"> 3 </TD> <TD> S1P2_KO   </TD> <TD>         </TD> <TD>     4 </TD> <TD>     0 </TD> <TD>     0 </TD> <TD>     7 </TD> <TD>     1 </TD> <TD>     0 </TD> <TD>     4 </TD> <TD>     4 </TD> </TR>
  <TR> <TD align="right"> 4 </TD> <TD> S1P3_KO   </TD> <TD>         </TD> <TD>     0 </TD> <TD>     0 </TD> <TD>     0 </TD> <TD>     0 </TD> <TD>     3 </TD> <TD>     7 </TD> <TD>     0 </TD> <TD>     0 </TD> </TR>
   </TABLE>



#### c) Explore differenitial readings on one probe

I decided to examining the readings of the probe for S1P2 (S1pr2). The probe for S1P2 is 99372_at on this platform.


```r
extractByProbe <- function(eDat, probeId, eDesign=design){
  (theProbe <- which(row.names(eDat) == probeId))
  pDat <- data.frame(eDesign, gExp = unlist(eDat[theProbe, ]))
  return(pDat)
}

plotDiffExp <- function(pDat){
  p <- ggplot(pDat, 
              aes(x=Genotype, y=gExp, 
                  group=BrainRegion, color=BrainRegion, shape=Sex)) + 
      geom_point(alpha=0.6, size=3) + 
      geom_smooth(method=loess, se=F)
  return(p)
}

pDat <- extractByProbe(expDat, '99372_at')
plotDiffExp(pDat) + facet_grid(.~Sex)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


Clearly something isn't right here. The probe reading of S1P2 in the S1P2 knockout group is supposed to decrease, yet the reading is the same or even higher there than that in the wildtype. 

#### d) Differential readings on one probe (numeric)


```r
# compute mean expression at each combination of Genetype, BrainRegion and Sex
pMean <- aggregate(gExp ~ Genotype + BrainRegion + Sex, pDat, FUN = mean)

# display the result in wide form
x <- xtabs(gExp~BrainRegion+Sex+Genotype,pMean)
ftable(x)
```

```
                   Genotype Wild_type S1P2_KO S1P3_KO
BrainRegion Sex                                      
hippocampus female              6.623   6.724   6.445
            male                6.539   6.491   6.676
neocortex   female              6.645   6.608   6.398
            male                6.587   6.518   6.851
```


### Q2 Examine the sample correlation matrix

#### a)


```r
# creating some color palette
jPurplesFun <- colorRampPalette(brewer.pal(n = 9, "Purples"))
jGreysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
```



```r
#fig.width=12, fig.height=12}
sampleCor <- cor(expDat)

myheatmap <- function(sampCor, ...){
  heatmap.2(sampCor, 
            Rowv = FALSE, dendrogram="none",
            symm=TRUE, margins=c(10,10),
            trace="none", scale="none", col = jGreysFun(256))
}
myheatmap(sampleCor)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 



b) Outlier sample

To show how the outlier stands out visually, I computed the mean correlation for each row. And display the distribution of the means using a density plot.

```r
meanCor <- rowMeans(sampleCor)
ggplot(data.frame(meanCor=meanCor),
       aes(x=meanCor)) +
  geom_bar(alpha=0.5)
```

```
stat_bin: binwidth defaulted to range/30. Use 'binwidth = x' to adjust this.
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 

```r

outlierIdx <- which(meanCor==min(meanCor))
(outlierName <- names(outlierIdx))
```

```
[1] "Day-6_S1P3_hipp_f_GSM172976"
```


Numerically, I computed the rank sum test to show that the correlation involving outlier sample has a different mean from those not involving the outlier sample. I suppose I have to assume independence between all these correlation coefficeints which obviously isn't true, but maybe I can be sloppy here... 




#### Q2c: Examine the outlier in the context of its experimental group.

Which group is the outlier in, in terms of Genotype, BrainRegion, and Sex? 


```r

group.def <- c("Genotype", "BrainRegion", "Sex")
thisGroup <- design[outlierIdx,]
```



```r
isSameGroup <- function(x, y, group.def=TRUE){
  return(all(x[group.def] == y[group.def]))
}
members <- alply(design, 1, isSameGroup,
                   y=design[outlierIdx,], group.def=group.def)
members <- unlist(members)
colnames(expDat)[members]
```

```
[1] "Day-5_S1P3_hipp_f_GSM172974" "Day-6_S1P3_hipp_f_GSM172975"
[3] "Day-6_S1P3_hipp_f_GSM172976"
```


Scatter plot these samples against each other and comment on the plot. Try to address overplotting! Hint: `splom()` from `lattice` is helpful.


```r
sampLabels <- design$Sid
sampLabels[outlierIdx] <- "Outlier"
grpSampCor <- expDat[, members]
colnames(grpSampCor) <- sampLabels[members]
splom(grpSampCor, panel = panel.smoothScatter)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21.png) 


### 3  Normalization 

#### a) Prior to normalization


```r
# elongate <- function(hDat){
#   sampleName <- factor(colnames(eDat))
#   ldply(hDat, summarize){
# #     print(dim(x))
# #     print(names(x))
# #     return(data.frame(gExp=x, sampleName=x))
#   })
# }
```



```r
# set.seed(540)
# hSize <- 2000
# theseProbes <- sample(1:nrow(expDat), hSize)
# expDat2 <- expDat[theseProbes,]
# expDat2 <- expDat

sampleBoxplot <- function(eDat){
  
  # create long format of the expression data, 
  # better suited for plotting
  longExpDat <- 
    melt(cbind(eDat, probe=rownames(eDat)), 'probe')
  
  longExpDat <- 
    rename(longExpDat, 
           c("variable"="SampleName", "value"="gExp"))

  p <- ggplot(longExpDat, 
            aes(SampleName, gExp, color=SampleName==outlierName))
  p <- p + geom_boxplot()
  return(p)
}

sampleBoxplot(expDat)
```

![plot of chunk unnamed-chunk-23](figure/unnamed-chunk-23.png) 




```r

myNormalize <- function(eDat){
  neDat <- normalize.quantiles(as.matrix(eDat))
  neDat <- data.frame(neDat)
  rownames(neDat) <- rownames(eDat)
  colnames(neDat) <- colnames(eDat) 
  return(neDat)
}

nExpDat <- myNormalize(expDat)
sampleBoxplot(nExpDat)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 



#### b) With normalization alone

```r
myheatmap(cor(nExpDat))
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25.png) 



#### c-e) With outlier removed and quantile normalization

```r
nrExpDat <- expDat[, colnames(expDat) != outlierName]
nrExpDat <- myNormalize(nrExpDat)
nrDes <- design[-c(outlierIdx),]
myheatmap(cor(nrExpDat))
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26.png) 



```r
sampleBoxplot(nrExpDat)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27.png) 




```r
plotDiffExp(extractByProbe(nrExpDat, '99372_at', nrDes))
```

![plot of chunk unnamed-chunk-28](figure/unnamed-chunk-28.png) 

```r
#plotDiffExp(miniDat <- extractByProbe(nrExpDat, '92352_at', nrDes))
```

### 4 Differential expression across genotypes (within neocortex brain region)

a)


```r
isNc <- nrDes$BrainRegion == "neocortex"
ncDes <- subset(nrDes, subset=isNc)
ncDat <- subset(nrExpDat, select=isNc)
```



```r
ncDesMat <- model.matrix(~Genotype, ncDes)
ncFit <- lmFit(ncDat, ncDesMat)
ncEbFit <- eBayes(ncFit)
koHits <- topTable(ncEbFit, number=50, adjust.method="BH",
                   coef = grep("KO", colnames(coef(ncEbFit))))
```



#### b)

```r
hitsExp <- ncDat[rownames(koHits), order(ncDes$Genotype)]
heatmap.2(as.matrix(hitsExp), 
          Colv = NA, Rowv = NA, scale="none", trace="none",
          margins = c(12, 10), col = jGreysFun(256))
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31.png) 


#### c)

```r
pMax <- 1e-3
# unadjTop <- topTable(ncEbFit, number=Inf, 
#                      adjust.method="none", p.value=pMax,
#                      coef = grep("KO", colnames(coef(ncEbFit))))
koFullTable <- topTable(ncEbFit, number=Inf,
                        coef = grep("KO", colnames(coef(ncEbFit))))
hitsByPval <- subset(koFullTable, P.Value < pMax)
```

The number of hits with p-value smaller than 0.001 is 


```

Error in base::parse(text = code, srcfile = NULL) : 
  2:0: unexpected end of input
1: nrow(hitsByPval
   ^

```

.

Computing FDR using q-value (adj.P.Val)

```r
q <- tail(koHits, n=1)$adj.P.Val # FDR of hits
eFP <- nrow(koHits) * q # expected # FD
# P <- nrow(koHits)
# FDR <- eFP / P
```


d)

```r
# selecting interesting and boring probes
intrProbes <- rownames(head(koFullTable,3))
borProbes <- rownames(tail(koFullTable,3))

multiProbeExtract <- function(probesL, xDat, xDes, status){
  multiProbeDat <- ldply(probesL, function(jProbe){
    data <- extractByProbe(xDat, jProbe, xDes)
    data$status <- factor(status)
    data$probe <- factor(jProbe)
    return(data)
  })
}

intDat <- multiProbeExtract(intrProbes, ncDat, ncDes, status="interesting")

borDat <- multiProbeExtract(borProbes, ncDat, ncDes, status="boring")

selectNcDat <- rbind(intDat, borDat)
```



```r
ggplot(selectNcDat, 
         aes(x=Genotype, y=gExp, group=probe, color=status)) + 
         geom_point() + geom_smooth(method="loess") + 
  facet_wrap(~ probe, ncol=3)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35.png) 




e) 

```r
s1p3KOHits <- topTable(ncEbFit, number=Inf,
                       adjust.method="BH", p.value=0.1,
                       coef = grep("S1P3_KO", colnames(coef(ncEbFit))))
nrow(s1p3KOHits)
```

```
[1] 61
```



Q5 (5 points) Differential expression analysis for Genotype * BrainRegion.

You should be using data with the worst outlier removed, after quantile normalization and for both brain regions.

Q5a: Fit a 3x2 full factorial model.

Test for any effect of Genotype and/or BrainRegion, i.e. test your model against a model with just an intercept. How many probes have a BH-adjusted p-value, a.k.a. q-value, less than 1e-3?



```r
desMat <- model.matrix(~Genotype*BrainRegion, nrDes)
fit <- lmFit(nrExpDat, desMat)
ebFit <- eBayes(fit)
hits <- topTable(ebFit, number=Inf, adjust.method="BH", p.value=1e-3,
                   coef = colnames(coef(ebFit))[-1])
nrow(hits)
```

```
[1] 1524
```



Q5b: Test the null hypothesis that BrainRegion doesn't matter, i.e. that all terms involving BrainRegion are zero.
How many probes have a BH-adjusted p-value less than 0.1?


```r
fullTableBrn <- topTable(ebFit, number=Inf, adjust.method="BH",
  coef = grep("BrainRegion", colnames(coef(ebFit))))
nrow(subset(fullTableBrn, P.Value < 0.1))
```

```
[1] 4377
```


Q5c: Highlight some probes where BrainRegion does and does not matter.

Using the results from Q5b, plot and describe some results for a couple of hits and non-hits.


```r
# selecting interesting and  boring probes
interestProbesBrn <- rownames(head(fullTableBrn,3))
borProbesBrn <- rownames(tail(fullTableBrn,3))

# intDatBrn <- ldply(interestProbesBrn, function(probe){
#   data <- extractByProbe(nrExpDat, probe, nrDes)
#   data$status <- factor("Interesting")
#   data$probe <- factor(probe)
#   return(data)
# })
# 

# borDatBrn <- ldply(borProbesBrn, function(probe){
#   data <- extractByProbe(nrExpDat, probe, nrDes)
#   data$status <- factor("Boring")
#   data$probe <- factor(probe)
#   return(data)
# })

intDatBrn <- multiProbeExtract(interestProbesBrn, nrExpDat, nrDes, status="interesting")

borDatBrn <- multiProbeExtract(borProbesBrn, nrExpDat, nrDes, status="boring")

selectDatBrn <- rbind(intDatBrn, borDatBrn)
# for (pb in interestBrn){
#   print(plotDiffExp(extractByProbe(nrExpDat, pb, nrDes)))
# }
```



```r
ggplot(selectDatBrn, 
       aes(x=Genotype, y=gExp, group=BrainRegion, color=BrainRegion)) +
  geom_point() + geom_smooth(method="loess") + 
  facet_wrap(~ probe, ncol=3)
```

![plot of chunk unnamed-chunk-40](figure/unnamed-chunk-40.png) 





Q5d: Test the null hypothesis that Genotype doesn't matter, i.e. that all terms involving Genotype are zero.

How many probes have a BH-adjusted p-value less than 0.1? Compare these results to those obtained for BrainRegion in Q5b, either based on this count or based on all the p-values. What do you conclude about the relative magnitude of the influence of brain region vs. genotype on gene expression?


```r
hitsGn <- topTable(ebFit, number=Inf, adjust.method="BH", p.value=0.1,
  coef = grep("Genotype", colnames(coef(ebFit))))
nrow(hitsGn)
```

```
[1] 141
```


