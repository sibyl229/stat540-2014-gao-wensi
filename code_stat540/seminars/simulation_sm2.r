

sampSize <- 10
numSamps <- 4

# generate random samples from normal distribution
set.seed(540)
allRands <- rnorm(sampSize*numSamps)
x <- matrix(allRands, nrow=sampSize)

# modify rpw and column names
rownames(x) <- sprintf("obs%02d", 1:sampSize)
colnames(x) <- sprintf("samp%d", 1:numSamps)

# compute mean for each colunm
## MARGIN argument (2nd) 1 means row dimension, 2 means column dimension
apply(x, 2, mean)

