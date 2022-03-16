setwd("/Users/bamflappy/NFCDSWorkshop_Day2")

##prefetch SRR8288561 SRR8288562 SRR8288563 SRR8288564 SRR8288557 SRR8288560 SRR8288558 SRR8288559 SRR8288565 SRR8288566 SRR8288567 SRR8288568

##fastq-dump --gzip SRR8288561; fastq-dump --gzip SRR8288562; fastq-dump --gzip SRR8288563; fastq-dump --gzip SRR8288564; fastq-dump --gzip SRR8288557; fastq-dump --gzip SRR8288560
##fastq-dump --gzip SRR8288558; fastq-dump --gzip SRR8288559; fastq-dump --gzip SRR8288565; fastq-dump --gzip SRR8288566; fastq-dump --gzip SRR8288567; fastq-dump --gzip SRR8288568

#brew install fastqc

#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh

#bash Miniconda3-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda

#source /Users/bamflappy/miniconda/bin/activate

#conda init zsh

#conda list

#conda install -c bioconda gffread

#gffread -E -F -T Tribolium_castaneum.gff3 -o Tribolium.gtf

#export PATH=$PATH:$PWD/hisat2-2.2.1

#hisat2-build Tribolium_castaneum.genome.fa TriboliumBuild

##cntrl samples 4h
#hisat2 -q -x TriboliumBuild -U SRR8288561.fastq.gz -S SRR8288561_accepted_hits.sam --summary-file SRR8288561_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288562.fastq.gz -S SRR8288562_accepted_hits.sam --summary-file SRR8288562_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288563.fastq.gz -S SRR8288563_accepted_hits.sam --summary-file SRR8288563_alignedSummary.txt

##cntrl samples 24h
#hisat2 -q -x TriboliumBuild -U SRR8288558.fastq.gz -S SRR8288558_accepted_hits.sam --summary-file SRR8288558_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288567.fastq.gz -S SRR8288567_accepted_hits.sam --summary-file SRR8288567_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288568.fastq.gz -S SRR8288568_accepted_hits.sam --summary-file SRR8288568_alignedSummary.txt

##treat samples 4h
#hisat2 -q -x TriboliumBuild -U SRR8288564.fastq.gz -S SRR8288564_accepted_hits.sam --summary-file SRR8288564_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288557.fastq.gz -S SRR8288557_accepted_hits.sam --summary-file SRR8288557_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288560.fastq.gz -S SRR8288560_accepted_hits.sam --summary-file SRR8288560_alignedSummary.txt

##treat samples 24h
#hisat2 -q -x TriboliumBuild -U SRR8288559.fastq.gz -S SRR8288559_accepted_hits.sam --summary-file SRR8288559_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288565.fastq.gz -S SRR8288565_accepted_hits.sam --summary-file SRR8288565_alignedSummary.txt
#hisat2 -q -x TriboliumBuild -U SRR8288566.fastq.gz -S SRR8288566_accepted_hits.sam --summary-file SRR8288566_alignedSummary.txt

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Rsubread")

library("Rsubread")

#?featureCounts

#cntrl samples 4h
cntrl1_fc_4h <- featureCounts(files="SRR8288561_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl2_fc_4h <- featureCounts(files="SRR8288562_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl3_fc_4h <- featureCounts(files="SRR8288563_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)

#cntrl samples 24h
cntrl1_fc_24h <- featureCounts(files="SRR8288558_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl2_fc_24h <- featureCounts(files="SRR8288567_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
cntrl3_fc_24h <- featureCounts(files="SRR8288568_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)

#treat samples 4h
treat1_fc_4h <- featureCounts(files="SRR8288564_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat2_fc_4h <- featureCounts(files="SRR8288557_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat3_fc_4h <- featureCounts(files="SRR8288560_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)

#treat samples 24h
treat1_fc_24h <- featureCounts(files="SRR8288559_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat2_fc_24h <- featureCounts(files="SRR8288565_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)
treat3_fc_24h <- featureCounts(files="SRR8288566_accepted_hits.sam", annot.ext="Tribolium.gtf", isGTFAnnotationFile=TRUE)

names(cntrl1_fc_4h)
names(treat1_fc_24h)

head(cntrl1_fc_4h$counts)
head(treat1_fc_24h$counts)

tribolium_counts <- data.frame(
  SRR8288561 = unname(cntrl1_fc_4h$counts),
  SRR8288562 = unname(cntrl2_fc_4h$counts),
  SRR8288563 = unname(cntrl2_fc_4h$counts),
  SRR8288558 = unname(cntrl1_fc_24h$counts),
  SRR8288567 = unname(cntrl2_fc_24h$counts),
  SRR8288568 = unname(cntrl2_fc_24h$counts),
  SRR8288564 = unname(treat1_fc_4h$counts),
  SRR8288557 = unname(treat2_fc_4h$counts),
  SRR8288560 = unname(treat2_fc_4h$counts),
  SRR8288559 = unname(treat1_fc_24h$counts),
  SRR8288565 = unname(treat2_fc_24h$counts),
  SRR8288566 = unname(treat2_fc_24h$counts)
)

head(tribolium_counts)

rownames(tribolium_counts) <- rownames(cntrl1_fc$counts)

head(tribolium_counts)

#Add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("cntrl_24h",3), rep("treat_4h",3), rep("treat_24h",3)))

#BiocManager::install("edgeR")

library("edgeR")

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)

#Plot the library sizes before normalization
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects before normalization
plotMDS(list, col=rep(1:12, each=3))

#There is no purpose in analysing genes that are not expressed in either 
# experimental condition, so genes are first filtered on expression levels
keep <- filterByExpr(list)
table(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Calculate normalized factors
list <- calcNormFactors(list)
normList <- cpm(list, normalized.lib.sizes=TRUE)

#View normalization factors
list$samples
dim(list)

#Plot the library sizes after normalization
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

#Draw a MDS plot to show the relative similarities of the samples
# and to view batch and treatment effects after normalization
plotMDS(list, col=rep(1:12, each=3))

#Draw a heatmap of individual RNA-seq samples using moderated
# log-counts-per-million after normalization
logcpm <- cpm(list, log=TRUE)
heatmap(logcpm)

#Produce a matrix of pseudo-counts
#Estimate common dispersion and tagwise dispersions
list <- estimateDisp(list)
list$common.dispersion

#View dispersion estimates and biological coefficient of variation
plotBCV(list)


##############
#Perform an exact test for treat_4h vs ctrl_4h
tested <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table

#Create filtered results table of DE genes
resultsTbl.keep <- resultsTbl$FDR <= 0.05
resultsTblFiltered <- resultsTbl[resultsTbl.keep,]

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
plotMD(tested)
abline(h=c(-1, 1), col="blue")

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
plotSmear(tested)


##############
#Perform an exact test for treat_24h vs ctrl_24h
tested <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

#Create results table of DE genes
resultsTbl <- topTags(tested, n=nrow(tested$table))$table

#Create filtered results table of DE genes
resultsTbl.keep <- resultsTbl$FDR <= 0.05
resultsTblFiltered <- resultsTbl[resultsTbl.keep,]

#Look at the counts-per-million in individual samples for the top genes
o <- order(tested$table$PValue)
cpm(list)[o[1:10],]

#View the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested))

#Plot log-fold change against log-counts per million, with DE genes highlighted
#The blue lines indicate 2-fold changes
plotMD(tested)
abline(h=c(-1, 1), col="blue")

#Make a mean-difference plot of two libraries of count data with smearing of points
#  with very low counts, especially those that are zero for one of the columns
plotSmear(tested)
