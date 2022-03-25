# set the working directory
setwd("/YOUR/FILE/PATH/")

# import gene count data
tribolium_counts <- read.csv("TriboliumCounts.csv", row.names="X")

# select a subset of the gene count data
#tribolium_counts[ , 1:6]

#Add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))

#BiocManager::install("edgeR")

library("edgeR")

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)

#Plot the library sizes before normalization and write to a jpg file
jpeg("exactTest_librarySizes_beforeNorm.jpg")
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
dev.off() 

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects before normalization
jpeg("exactTest_MDS_beforeNorm.jpg")
plotMDS(list, col=rep(1:12, each=3))
dev.off()

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Calculate normalized factors
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Write the normalized counts to a file
write.table(normList, file="tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)

#View normalization factors
list$samples
dim(list)

#Plot the library sizes after normalization
jpeg("exactTest_librarySizes_afterNorm.jpg")
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
dev.off()

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
jpeg("exactTest_MDS_afterNorm.jpg")
plotMDS(list, col=rep(1:12, each=3))
dev.off()

#Calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
jpeg("exactTest_logCPM.jpg")
heatmap(logcpm)
dev.off()

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion

#View dispersion estimates and biological coefficient of variation
jpeg("exactTest_BCV.jpg")
plotBCV(list)
dev.off()

##############
#Perform an exact test for treat_4h vs ctrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

#Create results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table))$table

#Create a table of DE genes filtered by FDR
resultsTbl_4h.keep <- resultsTbl_4h$FDR <= 0.05
resultsTbl_4h_filtered <- resultsTbl_4h[resultsTbl_4h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_4h_filtered, file="exactTest_4h_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_4h$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("exactTest_4h_DE.jpg")
plotMD(tested_4h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("exactTest_4h_smear.jpg")
plotSmear(tested_4h)
dev.off()

##############
#Perform an exact test for treat_24h vs ctrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

#Create a table of DE genes filtered by FDR
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table))$table

#Create filtered results table of DE genes
resultsTbl_24h.keep <- resultsTbl_24h$FDR <= 0.05
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_24h_filtered, file="exactTest_24h_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_24h$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("exactTest_24h_DE.jpg")
plotMD(tested_24h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("exactTest_24h_smear.jpg")
plotSmear(tested_24h)
dev.off()

##############
#Perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

#Create a table of DE genes filtered by FDR
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table))$table

#Create filtered results table of DE genes
resultsTbl_treat.keep <- resultsTbl_treat$FDR <= 0.05
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_treat_filtered, file="exactTest_treat_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_treat$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("exactTest_treat_DE.jpg")
plotMD(tested_treat)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("exactTest_treat_smear.jpg")
plotSmear(tested_treat)
dev.off()

##############
#Perform an exact test for cntrl_4h vs nctrl_24h
tested_cntrl <- exactTest(list, pair=c("treat_24h", "treat_4h"))

#Create a table of DE genes filtered by FDR
resultsTbl_nctrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table))$table

#Create filtered results table of DE genes
resultsTbl_ctrl.keep <- resultsTbl_nctrl$FDR <= 0.05
resultsTbl_cntrl_filtered <- resultsTbl_nctrl[resultsTbl_ctrl.keep,]

#Write the results of the exact tests to a csv file
write.table(resultsTbl_cntrl_filtered, file="exactTest_cntrl_filtered.csv", sep=",", row.names=TRUE)

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested_cntrl$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
jpeg("exactTest_cntrl_DE.jpg")
plotMD(tested_cntrl)
abline(h=c(-1, 1), col="blue")
dev.off()

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
jpeg("exactTest_cntrl_smear.jpg")
plotSmear(tested_cntrl)
dev.off()
