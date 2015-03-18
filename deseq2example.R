# Basic DESeq2 commands for a read counts table
# Using the selex dataset from ALDEx2

#DESEq2 is a Bioconductor package available here:
#http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html

# Please see the Bioconductor page for the manual and install instructions


#load libraries for DESeq2 and ALDEx2
#ALDEx2 contains the selex dataset
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

#Look 
plotMA(ddsoutput)


#----------------------------------------------------------------------------------------
# Another way to load data (from a tab-delimited, plaintext table)
# This would replace data("selex")

#The pasilla dataset is available from:
#http://bioconductor.org/packages/2.11/data/experiment/html/pasilla.html

d<-read.table("pasilla.txt", sep="\t", quote="", check.names=F, header=T, row.names=1)
