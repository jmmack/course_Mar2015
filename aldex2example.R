#example to run ALDEx2
#ALDEx2 is a Bioconductor package available here:
#http://www.bioconductor.org/packages/release/bioc/html/ALDEx2.html

# Please see the Bioconductor page for the manual and install instructions

library(ALDEx2)

#read in a table of counts data
d<-read.table("table.txt", sep="\t", quote="", check.names=F, header=T, row.names=1)

#alternatively, you can load the example dataset "selex"
#data(selex)

#Make a vector of conditions. This must be in the same order and the same number as the columns (samples) of the input table
conds<-c(rep("N", 7), rep("S", 7)


#get the clr values
#this is the main ALDEx function for all downstream analyses
#mc.samples=128 is often sufficient
x <- aldex.clr(d, mc.samples=1024, verbose=TRUE)

#perform t-test (both Welches and Wilcoxon, plus a Benjamini-Hochberg multiple test correction)
x.tt <- aldex.ttest(x, conds, paired.test=FALSE)

#estimate effect size and the within and between condition values
#include indiv. samples or not
x.effect <- aldex.effect(x, conds, include.sample.summary=TRUE, verbose=TRUE)


#merge the data
x.all <- data.frame(x.tt, x.effect)

#write a .txt with your results
write.table(x.all, file="aldex_ttest.txt", sep="\t", quote=F, col.names=NA)
