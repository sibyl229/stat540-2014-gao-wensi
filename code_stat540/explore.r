eDat <- read.table("data/GSE4051_MINI.txt", header=TRUE, row.names = 1)

weeDat <- eDat[eDat$poisonFang > 7.5,]

sample(weeDat)
sampleRowNum <- sample(1:nrow(eDat), size=3, replace=TRUE)
eDat[sampleRowNum,]