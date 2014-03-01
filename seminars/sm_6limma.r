library(limma)
library(ggplot2)

prDat <- read.table("../data/GSE4051_data.tsv")
preDes <- readRDS("../data/GSE4051_design.rds")

m <- 1000
n <- 3
x <- matrix(rnorm(m * n), nrow = m)
obsVars <- apply(x, 1, var)


wtDes <- subset(prDes, gType == "wt")