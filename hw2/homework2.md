Homework 02
======================================================================









## Q1) Microarray Analysis

The six samples in this study, 3 replicates for each of the two conditions, were analysed on the Affymetrix Yeast Genome Array 2.0 platform. 

### a) (1pt) Load Microarray Data

Load the normalized data.   

```r
mcDat <- read.table("../data/yeast/GSE37599-data.tsv", 
                    header=TRUE, sep="\t", row.names=1)
```

  
What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent and how you're interpreting column names.

The dataset has 10928 rows, each of which represents a probe. It has 6 columns, each of which represents a sample. 

First letter in the column name refers to the condition (batch vs chemostat) and the number refers to which replicate the sample belongs to.



### b) (2pt) Identify Sample Swap

The labels on two of the samples have been swapped, that is one of the batch samples has been labelled as chemostat and vice-versa. Produce the plots described below and explain how they allow you to identify the swapped samples.
  
i. (High volume) scatter plot matrix. 

```r
# png("figure/sampleSwap.scatter.png")
# splom(mcDat, panel = panel.smoothScatter)
# dev.off()
```

![Scatter plot across all pair of samples](figure/sampleSwap.scatter.png)

b1 and c2 are swapped, because the scatter plots show that b1 and c2 are each better correlated (showing less scatter) with samples of the Opposite condition than samples of its labeled condition.

ii. A heatmap of the first 100 genes (you can try more but it gets slow).


```r
# creating some color palette
jPurplesFun <- colorRampPalette(brewer.pal(n = 9, "Purples"))
jGreysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))



myheatmap <- function(sampCor, ...){
  heatmap.2(sampCor, 
            Rowv = FALSE, dendrogram="none",
            symm=TRUE, margins=c(10,10),
            trace="none", scale="none", col = jPurplesFun(256))
}

# cooment out because it's slow
# png("figure/sampleSwapSample.heatmap.png")
# myheatmap(as.matrix(mcDat[1:100,]))
# dev.off()
```

![100 gene expression heatmap](figure/sampleSwapSample.heatmap.png)

The heatmap of expression level of 100 gene across samples shows that the expression pattern is similar for (b2, b3, c2) and (b1, c1, c3). This indicates that b1 and c2 are swapped.

iii. Compute the Pearson correlation of the samples and plot the results using a heatmap.


```r
sampleCor <- cor(mcDat)
myheatmap(sampleCor)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7.png) 


Darker color here means stronger correlation between the corresponding pair of samples. b1 are more correlated with c1 and c3, opposite of b1's labeled group.
Similarly for c2, more correlated to samples of the oppositely labeled samples.

iv. Scatterplot the six data samples with respect to the first two principal components and label the samples.


```r
createDesign <- function(mcDat, chemoLabels){
  sampleLabels <- colnames(mcDat)
  conditions <- factor(rep("batch", length(sampleLabels)),
                       levels=c("batch", "chemostat"))
  conditions[sampleLabels %in% chemoLabels] <- "chemostat"
  mcDes <- data.frame(condition=conditions)
  rownames(mcDes) <- sampleLabels  
  return(mcDes)
} 

sampleLabels <- colnames(mcDat)
chemoSmpls <- grep("c", sampleLabels, value=TRUE)
design <- createDesign(mcDat, 
                       chemoLabels=chemoSmpls)

mcPC <- prcomp(mcDat, scale=TRUE)$rotation[,1:2]
mcPC <- cbind(design, mcPC)
p <- ggplot(mcPC, aes(x=PC1, y=PC2, color=condition, 
                       label=rownames(mcPC)))
p <- p + geom_point() 
p <- p + geom_text()
p
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8.png) 


The first two principal components of b1 and c2 are quite dissimilar to those of other samples labeled with the same condtion, while more similar to the principal components of the samples of the labeled the opposite codition.

### c) (2pt) Microarray Differential Expression

Fix the label swap identified in question 1b. We want to swap b1 <--> c2. Revisit one or more elements of question 1b to sanity check before proceeding. 


```r
# correcting the labels (relying on copy-on-modify property of R)
oldSampLabels <- sampleLabels
sampleLabels[grep("b1", oldSampLabels)] <- "c2"
sampleLabels[grep("c2", oldSampLabels)] <- "b1"

print(oldSampLabels)
```

```
[1] "b1" "b2" "b3" "c1" "c2" "c3"
```

```r
print(sampleLabels)
```

```
[1] "c2" "b2" "b3" "c1" "b1" "c3"
```

```r

# relabel samples by the correct labels
colnames(mcDat) <- sampleLabels  
# sort the data samples to match with the design
mcDat <- mcDat[,rownames(design)]

sampleCor <- cor(mcDat)
myheatmap(sampleCor)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9.png) 


After fixing the sample swap, samples of the same condition are more similar to each other (shown by the high correlation - in dark color), than samples of differnet conditions.

Now use this data to do a differential expression analysis with `limma`.


```r
desMat <- model.matrix(~condition, design)
mcEBFit <- eBayes(lmFit(mcDat, desMat))
mcTT <- topTable(mcEBFit, number=Inf, adjust.method="BH",
                 coef = grep("condition", colnames(desMat)))
```



Package these results in a data frame with six columns:

* probe.id - The array probe id.

* gene.id - The id of the gene which the probe overlaps (see below).

* p.value - The raw p-value for the probe.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the limma results table.

* test.stat - The test statistics which for limma is the moderated t statistic. This is the column called "t" in the limma results table.

>  The gene id can be retrieved using the `yeast2.db` package from Bioconductor. In particular, the `yeast2ORF` object available after loading `yeast2.db` contains the mapping between probe IDs and yeast gene ids. Assuming you have a character vector of probes called `probe.ids`, the gene IDs can be retrieved using `gene.ids <- unlist(mget(probe.ids, yeast2ORF))`.


```r
# gene names sorted by toptable result above
genes <- unlist(mget(rownames(mcTT), yeast2ORF))
expect_equal(names(genes), rownames(mcTT))

mcTT <- mcTT[, c("P.Value", "adj.P.Val", "logFC", "t")]
colnames(mcTT) <- c("p.value", "q.value", "log.fc", "test.stat")
mcTT <- cbind(probe.id=rownames(mcTT), gene.id=genes, mcTT)
```


Remove any rows with probes which don't map to genes. You'll be able to find these because they will have `NA` as their gene id. Work with this data.frame to answer the questions below.


```r
mcTT <- mcTT[!is.na(mcTT$gene.id),]
```



i. How many probes did we start with and how many remain after removing probes without gene ids?

We started with 10928 probes. And 5705 remain after removing probes without gene ids.



ii. Illustrate the differential expression between the batch and the chemostat samples for the top hit (i.e., probe with the lowest p- or q-value).


```r
expect_equal(rownames(design), colnames(mcDat))

prepareDataMC <- function(probeId){
  probeDat <- mcDat[rownames(mcDat)==probeId,]
  
  newDat <- cbind(design, t(probeDat))
  
  # in case of multiple probes selected
  newDat <- melt(newDat, "condition", variable_name="probe.id")
  newDat <- rename(newDat, c("value"="expression"))  # rename column
  return(newDat)
}

topHitPid <- mcTT[1,"probe.id"]
topHitDat <- prepareDataMC(topHitPid)

p <- ggplot(topHitDat, aes(x=condition, y=expression, color=condition)) 
p <- p + geom_point()
p
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13.png) 


Top hit shows little within group (condition) differnces, but very large between group (conditions) differences. This makes sense that limma picks it as significant, because the small within group difference is unlikely to produces by chance the large difference across groups.


iii. How many probes are identified as differentially expressed at a false discovery rate (FDR) of 1e-5 (note: this is a FDR cutoff used in the original paper)?


```r
fdr <- 1e-5
nHits <- sum(mcTT$q.value <= fdr)
```


725 probes are identified as differentially expressed at FDR of 10<sup>-5</sup>.

iv. Save your results for later with `write.table()`.

> When using write.table to save the data, you need to pass the arguments `row.names = TRUE, col.names = NA` to make sure the headers are written correctly.


```r
write.table(mcTT, "results/condition.topTable.GSE37599.tsv",
            row.names = TRUE, col.names = NA, sep="\t")
```


## Q2) RNA-Seq Analysis

> We have aligned the RNA-Seq library using the [Stampy](http://www.well.ox.ac.uk/project-stampy) aligner and generated count data. The data file is available as [stampy.counts.tsv](../../examples/yeastPlatforms/data/stampy.counts.tsv). In this question you will use this data to do a differential expression analysis using different packages from Bioconductor.


### a) (1pt) Load RNA Count Data and Sanity Check


```r
rcDat <- read.table("../data/yeast/stampy.counts.tsv", 
                    header=TRUE, row.names=1)
rcDat <- rcDat[, rownames(design)]  # make sure that sample (column) order matches the design
```


i) What are dimensions of the dataset? In addition to reporting number of rows and columns, make it clear what rows and columns represent. What is the difference between the rows of this dataset versus rows of the array data in question 1a?

The dataset has 6542 rows, each representing a gene (I guess Stampy aligns to the genome?), and 6 columns, each representing a sample. Here we have genes for the rows instead of probes with the microarray data, because, first, RNA-seq doesn't use probes, and, second, with that aligner RNA transcripts are mapped to genes.


ii) Do a sanity check to make sure there is no sample swap by plotting a heatmap of the sample correlations.


```r
sampleCor <- cor(rcDat)
myheatmap(sampleCor)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17.png) 


Since samples are strongly correlated within same condition, and less correlated between different conditions, no smaple seems to have occured.

### b) (2pt) `edgeR` Differential Expression Analysis

Now you will use `edgeR` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  Recall that `edgeR` needs to estimate the dispersion parameter in the negative binomial model using an empirical Bayes method. Estimate the dispersion parameters using `estimateGLMCommonDisp`, `estimateGLMTrendedDisp` and `estimateGLMTagwiseDisp`. Plot the tagwise dispersion against log2-CPM (counts per million).  


```r
dge.glm <- DGEList(counts=rcDat, group=design$condition)
desMat <- model.matrix(~condition, design)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, desMat, verbose=TRUE)
```

```
Disp = 0.00551 , BCV = 0.0742 
```

```r
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp, desMat)
```

```
Loading required package: splines
```

```r
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, desMat)
#plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18.png) 






ii)  Use the glm functionality of `edgeR`, i.e. use the `glmFit` function, to identify differentially expressed genes between conditions. 


```r
edger.fit <- glmFit(dge.glm.tag.disp, desMat)
edger.lrt <- glmLRT(edger.fit, coef=grep("condition", colnames(desMat)))
edger.results <- topTags(edger.lrt, n=Inf, adjust.method="BH", sort.by="PValue")
```


Package these results in a data.frame called 'edger.results' with five columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.

* test.stat - The test statistic, which for `edgeR` is a likelihood ratio. This is the column called "LR" in the `edgeR` results table.


```r
edger.results <- data.frame(edger.results[,c("PValue", "FDR", "logFC", "LR")])
colnames(edger.results) <- c("p.value", "q.value", "log.fc", "test.stat")
edger.results <- cbind(gene.id=rownames(edger.results), edger.results)
```


Save your results for later with `write.table()` in file called `stampy.edger.results.tsv`.


```r
write.table(edger.results, "results/stampy.edger.results.tsv",
            row.names=TRUE, col.names = NA, sep="\t")
```



iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?


```r
fdr <- 1e-5
edger.nHits <- sum(edger.results$q.value <= 1e-5)
```


2669 genes are differentially expressed between conditions at FDR of 10<sup>-5</sup>.

iv) How many genes are differentially over-expressed in chemostat compared to batch medium samples at a false discovery rate (FDR) of 1e-5?


```r
edger.OX.nHits <- sum(
  edger.results$q.value <= 1e-5 & 
  edger.results$log.fc > 0)
```


1515 are differentially over-expressed in chemostat compared to batch with 10<sup>-5</sup>.

### c) (2pt) `DESeq` Differential Expression Analysis

Now you will use `DESeq` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  `DESeq` also needs to estimate the dispersion. Use `estimateSizeFactors` and `estimateDispersions` to normalize the data. Plot the estimated dispersions against the mean normalized counts.


```r

#reading in the same count table data and grouping information
deSeqDat <- newCountDataSet(rcDat, conditions=design$condition)
#head(counts(deDat))

# account for differences in library coverage and
deSeqDat <- estimateSizeFactors(deSeqDat)  
#sizeFactors(deSeqDat)
deSeqDat <- estimateDispersions(deSeqDat)   # estimate variance

#plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24.png) 


ii)  Use the negative binomial test of `DESeq`, i.e. use the `nbinomTest` function, to identify differentially expressed genes between conditions. Note that the output of this function does not return results ordered by p-values or logged fold-changes. You can manually reorder the results if you want (not required for this homework).


```r
## this takes a minute or so for JB
deseq.results <- nbinomTest(deSeqDat, 
                      levels(design$condition)[1], 
                      levels(design$condition)[2])
plotMA(deseq.results)
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25.png) 


Package these results in a data.frame called 'deseq.results' with four columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "log2FoldChange" in the `deseq` results table.

Save your results for later with `write.table()` in file called `stampy.deseq.results.tsv`.


```r
deseq.results <- deseq.results[,c("id", "pval", "padj", "log2FoldChange")]
colnames(deseq.results) <- c("gene.id", "p.value", "q.value", "log.fc")
rownames(deseq.results) <- deseq.results$gene.id
write.table(deseq.results, "results/stampy.deseq.results.tsv",
            row.names=TRUE, col.names = NA, sep="\t")
```


iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?


```r
deseq.nHits <- sum(deseq.results$q.value <= fdr)
```


2198 genes are differentially expressed between conditions at a false discovery rate (FDR) of 10<sup>-5</sup>.

iv) How many differentially expressed genes are identified by both 'edgeR' and 'DESeq'?


```r
getHitsGeneId <- function(results, fdr=1e-5){
  isDiff <- (results$q.value <= fdr)
  return(results$gene.id[isDiff])
}

edgerHitGenes <- getHitsGeneId(edger.results)
deseqHitGenes <- getHitsGeneId(deseq.results)
bothEDHitGenes <- intersect(edgerHitGenes, deseqHitGenes)
```


2176 differentially expressed genes are identified by both 'edgeR' and 'DESeq'.


### d) (2pt) `voom` Differential Expression Analysis

Now you will use `voom+limma` to identify differentially expressed genes between the batch medium vs. chemostat conditions.

i)  `voom` normalizes the counts before it converts counts to log2-cpm. Use `calcNormFactors` to normalize counts.


```r
norm.factor <- calcNormFactors(rcDat)
lib.size <- colSums(rcDat)*norm.factor
```



ii)  Use `voom' to convert count data into logged CPM data and then use 'limma' to identify differentially expressed genes between conditions. 


```r
dat.voomed <- voom(rcDat, desMat, plot=TRUE,lib.size=lib.size)
```

![plot of chunk unnamed-chunk-30](figure/unnamed-chunk-30.png) 

```r

voom.fit <- eBayes(lmFit(dat.voomed, desMat))
voom.results <- topTable(voom.fit, number=Inf, adjust.method="BH",
                   coef = grep("condition", colnames(desMat)))
```


Package these results in a data.frame called 'voom.limma.results' with five columns:

* gene.id - The id of the gene which reads were aligned to.

* p.value - The raw p-value for the gene.

* q.value - The BH corrected p-value, aka the q-value.

* log.fc - The log fold change which is the column called "logFC" in the `edgeR` results table.

* test.stat - The test statistic, which is the column called "t".

Save your results for later with `write.table()` in file called `stampy.limma.results.tsv`.


```r
voom.results <- voom.results[, c("P.Value", "adj.P.Val", "logFC", "t")]
colnames(voom.results) <- c("p.value", "q.value", "log.fc", "test.stat")
voom.results <- cbind(gene.id=rownames(voom.results), voom.results)
```


iii) How many genes are differentially expressed between conditions at a false discovery rate (FDR) of 1e-5?


```r
voomHitGenes <- getHitsGeneId(voom.results)
```




```

Error in base::parse(text = code, srcfile = NULL) : 
  2:0: unexpected end of input
1: length(voomHitGenes
   ^

```

 genes are differentially expressed between conditions at a false discovery rate (FDR) of 10<sup>-5</sup>.

iv)  What fraction of the genes identified using `voom+limma` are also found by `edger` and `DESeq` methods? For example if the DE analysis using `voom+limma` found 1000 genes and both `edgeR` and `DESeq`  found 500 of these, the fraction of genes found would be $\frac{500}{1000}=0.5$.


```r
allHitGenes <- intersect(voomHitGenes, bothEDHitGenes)
r <- length(allHitGenes) / length(voomHitGenes)
```


99.8885% of the genes identified using `voom+limma` are also found by `edger` and `DESeq` methods.

### e) (3pt) Comparison of Differential Expression Analyses
 
Now that we have the results of the differential expression analysis performed by three popular methods, we are going to compare and illustrate the results.

i) In previous questions, we noticed that different methods identified different differentially expressed genes. Create a Venn diagram showing all genes identified as differentially expressed by `edgeR`, `DESeq`, and `voom+limma`. Check your if your answers to questions 2c-iv, and 2d-iv are correct.


```r

# Put the things you want to plot in a list. The names in the list will be
# put on the plot.
de.genes <- list("edgeR"=edgerHitGenes, 
                 "DEseq"=deseqHitGenes, 
                 "voom+limma"=voomHitGenes)

# Start a new plot
plot.new()

# Draw the Venn diagram. Note the argument `filename=NULL` tells it to
# create a plot object instead of outputting to file.
venn.plot <- venn.diagram(de.genes, filename = NULL, fill = c("blue", 
    "green", "yellow"))

# Draw the plot on the screen.
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-34](figure/unnamed-chunk-34.png) 



> The Venn diagram can be drawn using the `VennDiagram` package. It's a little obtuse if you want to plot to screen (or embed image using knitr), but the following code should get you started. Also be aware there is an argument called `force.unique`, which defaults to TRUE, that determines how elements that appear more than once in a set are handled when forming the Venn counts (in particular useful in question 3). 


ii) Using the function `plotSmear` function from `edgeR`, you can look at a scatterplot of observed differential expression (y-axis) against overall abundance (x-axis), both axes logarithmically transformed -- to check that putative DE genes seem plausible. Create a smear plot. Within this plot, identify the set of genes which are differentially expressed at an FDR of 1e-5 using all three methods (i.e., the q-values estimated by `edgeR`, `DESeq`, and `voom` are below 1e-5). Explain how you interpret this plot. Do you see reasonable results?


```r
plotSmear(dge.glm, de.tags=allHitGenes)
```

![plot of chunk unnamed-chunk-35](figure/unnamed-chunk-35.png) 


The red dots marking the genes that 3 methods all identify as significantly differentially expressed. These genes have in common that they are expressed in a relatively large amount (CPM) and their expression fold change is noticeable between the conditions.

This makes sense because if a gene's overall abundance is small, a small fluctuation by chance would result in a big fold change. On the other hand, gene of great abundance changing by a big number may only mean a small fold change. These genes don't appear in the hits agreed upon by all 3 methods.

> Use the output of `DGEList` as the object of `plotSmear`. Use de.tags to highlight genes selected by all methods.




iii) There are two genes identified by `edgeR` and `voom+limma` but not by `DESeq`. Illustrate the logged counts of them. Compare the (log) counts of these two genes with those of two genes identified by the three methods (see example below)


```r
bothEVHitGenes <- intersect(edgerHitGenes, voomHitGenes)
xHits <- setdiff(bothEVHitGenes, deseqHitGenes)

prepareData <- function(data, rid){
  jDat <- data[rownames(data) %in% rid,]
  
  newDat <- cbind(design, t(jDat))
  
  # in case of multiple probes selected
  newDat <- melt(newDat, "condition", variable_name="rid")
  newDat <- rename(newDat, c("value"="expression"))  # rename column
  return(newDat)
}

xHitsDat <- prepareData(rcDat, xHits)
```



```r
p <- ggplot(xHitsDat, aes(y=log2(expression), x=condition, group=rid))
p + geom_point() + geom_smooth(method="loess") + facet_grid(. ~ rid)
```

![plot of chunk unnamed-chunk-37](figure/unnamed-chunk-37.png) 


This two hits are missed probably because the variance of the expression within the same condition is too large, or that the expression fold change between conditions is too small.

I tried use expression scaled by library size and then log2 transform.


```r
p <- ggplot(prepareData(dat.voomed$E, xHits), 
            aes(y=expression, x=condition, group=rid))
p + geom_point() + geom_smooth(method="loess") + facet_grid(. ~ rid)
```

![plot of chunk unnamed-chunk-38](figure/unnamed-chunk-38.png) 



## Q3) Compare DEA results between RNA-Seq and array

In question 1, you performed a DEA of array data using `limma`. In question 2, you used different methods to perform DEA of RNA-Seq data. In particular, in this question, you will focus on the DEA using `edgeR` (question 2b) to compare the results of RNA-Seq DEA with those of the array DEA . 

i) Use a Venn diagram to display the overlap and non-overlap of the __genes__ identified as differentially expressed at an FDR of 1e-5 by these analyses (i.e., array vs `edgeR` differentially expressed genes).



```r
arrayHitGenes <- getHitsGeneId(mcTT)

plot.new()  # Start a new plot

# note multiple probes mapped to the same gene show up as one in the Venn diagram
venn.plot <- venn.diagram(list(array=arrayHitGenes,
                               edgeR=edgerHitGenes),
                          filename = NULL, fill = c("blue", "yellow"),
                          force.unique = TRUE) 
grid.draw(venn.plot)
```

![plot of chunk unnamed-chunk-39](figure/unnamed-chunk-39.png) 


RNA-seq (with EdgeR) finds vast majority to hits found by microarray (with limma). This is probably due to the larger dynamic range of count data than signal intensity, making RNA-seq better at detecting expression changes. (more sensitive and powerful?)


ii) As expected, more genes were identified as differentially expressed using RNA-Seq data. In this question, you will examine the difference between the q-values from both analyses (i.e., array and `edgeR`) by overlaying density plots of the q-values from each analysis.


```r
# prepare results for comparision between array and RNAseq (edgeR)
arrVsEdgeR <- rbind(cbind(mcTT[, c("gene.id","q.value")],
                          platform="Microarray+limma"),
                    cbind(edger.results[, c("gene.id","q.value")],
                          platform="RNA-seq+edgeR"))

arrEdgerBothGenes <- intersect(mcTT$gene.id, edger.results$gene.id)
```


  * To respond to this question, make two plots. One plot that includes the densities of q-values of the genes analyzed by both platforms (i.e., genes shared by both data frames).
  

```r
p <- ggplot(subset(arrVsEdgeR, gene.id %in% arrEdgerBothGenes),
            aes(x=q.value, color=platform)) + geom_density()
p 
```

![plot of chunk unnamed-chunk-41](figure/unnamed-chunk-411.png) 

```r

p <- ggplot(subset(arrVsEdgeR, gene.id %in% arrEdgerBothGenes),
            aes(x=sqrt(q.value), color=platform)) + geom_density()
p 
```

![plot of chunk unnamed-chunk-41](figure/unnamed-chunk-412.png) 


Since I'm mostly interested in the location of the peak in the density plot of q-value, it seems to make sense to square root transform (which is monotonic) the q.value, so I can see separation of the peaks for the two platforms better. 

It seems that RNA-seq has peak of q.value density plot shifted slightly to the left, regardless of using genes appearing in both (plots above) or either (plots below) analysis. This probably indicates that RNA-seq is more powerful test for differential analysis.
  
  * Another plot that includes the densities of q-values of ALL genes analyzed by at least one of the platforms.
  

```r
p <- ggplot(arrVsEdgeR, aes(x=q.value, color=platform)) + 
  geom_density() + xlim(0,0.1)
p 
```

![plot of chunk unnamed-chunk-42](figure/unnamed-chunk-421.png) 

```r

p <- ggplot(arrVsEdgeR, aes(x=sqrt(q.value), color=platform)) + 
  geom_density()
p 
```

![plot of chunk unnamed-chunk-42](figure/unnamed-chunk-422.png) 




Make some observations about the strengths of these two platforms.

iii) We provide a data set with array expression and count data for 5 interesting genes; below is also code to load it and a figure depicting it. Consult the array and the `edgeR` DEA results from your previous analyses for these genes. For each gene, state its status with respect to these analyses, i.e. where it falls in those Venn diagrams. Comment on the results and plots.



```r
jDat <- dget("../data/yeast/featGenesData-q3-DPUT.txt")
jDat <- melt(jDat, measure.vars=c("log.count", "arrayExp"),
             variable_name="platform")
jDat <- rename(jDat, c("value"="expression"))

p <- ggplot(jDat, aes(x=expression, y=gene.id, color=cond))
p + geom_point(alpha=0.5) + facet_grid(.~platform)
```

![plot of chunk unnamed-chunk-43](figure/unnamed-chunk-43.png) 



 * "YGL209W"
 I expect count to find no significance between the two conditions,
 while array to find it significantly under-expressed in chemostat condition. It likely belongs to the tiny region for microarray (in the Venn diagram) that's not part of the overlapping region.
 
 * "YCL042W" 
  I expect count to find it significantly over-expressed in chemostat condition,  
 while array to find no significance between the two conditions. It likely belongs to the big region for edgeR (in the Venn diagram) that's not part of the overlapping region.
 
 * "YBL025W" 
 I expect both count and array to find no significance between the two conditions, meaning that it does NOT belong to the Venn diagram at all.
 
 * "YDR384C" 
 I expect both count and array to find it significantly over-expressed in chemostat condition, meaning that it belongs to the overlapping region of the Venn diagram.
 
 * "YDR345C"
 I expect both count and array to find it significantly under-expressed in chemostat condition, meaning that it belongs to the overlapping region of the Venn diagram.
 
Because microarray and RNA-seq are quite different techniques, they don't necessarily detect the same change in expression. For example a microarray might not be able to detect small expression change or change beyond saturation of the signal. Also some ways of using RNA-seq ignores the effect of expression changes across different isoform of the same gene, which is an effect that might be picked up my some microarray. It's possible for the two platform (microarray and RNAseq) to give very different results (over whether a gene is diffentially expressed). To determine which technique to use may depend on the question.
