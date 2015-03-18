# Basic DESeq2 commands for a read counts table
# Using the selex dataset from ALDEx2

library(DESeq2)
library(ALDEx2)

#load selex datafame
data("selex")

#format into matrix of read counts
countsTable <- data.matrix(selex,rownames.force=NA)

#make table of conditions in same order as columns of counts table
#i.e. the first 7 samples are "N" condition and the next 7 samples are "S" condition
coldata <- data.frame(condition=factor(c(rep("N",7), rep("S", 7))))

#DESeq2 formats data from a matrix of read counts
dds <- DESeqDataSetFromMatrix(countData=countsTable, colData=coldata, design = ~condition)

#See documentation for fit types
ddsoutput <- DESeq(dds, fitType="local")

#get results
res <- results(ddsoutput)

#summarize results
summary(res)

plotMA(ddsoutput)
