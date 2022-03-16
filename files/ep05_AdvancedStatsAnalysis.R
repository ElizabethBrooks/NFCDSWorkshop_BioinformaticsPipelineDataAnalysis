setwd("/Users/bamflappy/NFCDSWorkshop_Day2")

#export PATH=$PATH:$PWD/sratoolkit.2.11.2-mac64/bin

##prefetch SRR8288561 SRR8288562 SRR8288563 SRR8288564 SRR8288557 SRR8288560 SRR8288558 SRR8288559 SRR8288565 SRR8288566 SRR8288567 SRR8288568

##fastq-dump --gzip SRR8288561; fastq-dump --gzip SRR8288562; fastq-dump --gzip SRR8288563; fastq-dump --gzip SRR8288564; fastq-dump --gzip SRR8288557; fastq-dump --gzip SRR8288560
##fastq-dump --gzip SRR8288558; fastq-dump --gzip SRR8288559; fastq-dump --gzip SRR8288565; fastq-dump --gzip SRR8288566; fastq-dump --gzip SRR8288567; fastq-dump --gzip SRR8288568

#brew install fastqc

#fastqc SRR8288561.fastq.gz --extract

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
  cntrl1_4h = unname(cntrl1_fc_4h$counts),
  cntrl2_4h = unname(cntrl2_fc_4h$counts),
  cntrl3_4h = unname(cntrl3_fc_4h$counts),
  treat1_4h = unname(cntrl1_fc_24h$counts),
  treat2_4h = unname(cntrl2_fc_24h$counts),
  treat3_4h = unname(cntrl3_fc_24h$counts),
  cntrl1_24h = unname(treat1_fc_4h$counts),
  cntrl2_24h = unname(treat2_fc_4h$counts),
  cntrl3_24h = unname(treat3_fc_4h$counts),
  treat1_24h = unname(treat1_fc_24h$counts),
  treat2_24h = unname(treat2_fc_24h$counts),
  treat3_24h = unname(treat3_fc_24h$counts)
)

# export gene count data
write.csv(tribolium_counts, "TriboliumCounts_fullSet.csv", row.names=FALSE)
tribolium_counts <- read.csv("TriboliumCounts_fullSet.csv")

head(tribolium_counts)

rownames(tribolium_counts) <- rownames(cntrl1_fc_4h$counts)

head(tribolium_counts)

#Import grouping factor
targets <- read.csv(file="groupingFactors.csv")
  
#BiocManager::install("edgeR")

library("edgeR")


#############
#Setup a design matrix
group <- factor(paste(targets$treatment,targets$hours,sep="."))
#cbind(targets,Group=group)

#Create DGE list object
list <- DGEList(counts=tribolium_counts,group=group)
colnames(list) <- rownames(targets)
head(list)

#Plot the library sizes before normalization
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)
#list$samples
#Write normalized counts to file
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Verify TMM normalization using a MD plot
#Write plot to file
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1)
colors <- rep(c("blue", "darkgreen"), 2)
#Write plot without legend to file
plotMDS(list, col=colors[group], pch=points[group])
#Write plot with legend to file
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)

#The experimental design is parametrized with a one-way layout, 
# where one coefficient is assigned to each group
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#design

#Next, the NB dispersion is estimated
list <- estimateDisp(list, design, robust=TRUE)
#list$common.dispersion
#Visualize the dispersion estimates with a BCV plot
#Write plot to file
plotBCV(list)

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)
#Write plot to file
plotQLDisp(fit)

#Test whether the average across all cntrl groups is equal to the average across
#all treat groups, to examine the overall effect of treatment
con.treat_cntrl <- makeContrasts(set.treat_cntrl = (treat.4h + treat.24h)/2
                                                    - (cntrl.4h + cntrl.24h)/2,
                                                    levels=design)

#Look at genes with significant expression across all UV groups
anov.treat_cntrl <- glmTreat(fit, contrast=con.treat_cntrl, lfc=log2(1.2))
summary(decideTests(anov.treat_cntrl))
#Write plot to file
plotMD(anov.treat_cntrl)
abline(h=c(-1, 1), col="blue")


#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.24h_4h <- makeContrasts(set.24h_4h = (cntrl.24h + treat.24h)/2
                                          - (cntrl.4h + treat.4h)/2,
                                          levels=design)

#Look at genes with significant expression across all UV groups
anov.24h_4h <- glmTreat(fit, contrast=con.24h_4h, lfc=log2(1.2))
summary(decideTests(anov.24h_4h))
#Write plot to file
plotMD(anov.24h_4h)
abline(h=c(-1, 1), col="blue")


#Test whether there is an interaction effect
con.interaction <- makeContrasts(set.interaction = ((treat.4h + treat.24h)/2
                                                    - (cntrl.4h + cntrl.24h)/2)
                                                    - ((cntrl.24h + treat.24h)/2
                                                    - (cntrl.4h + treat.4h)/2),
                                                    levels=design)

#Look at genes with significant expression
anov.interaction <- glmTreat(fit, contrast=con.interaction, lfc=log2(1.2))
summary(decideTests(anov.interaction))
#Write plot to file
plotMD(anov.interaction)
abline(h=c(-1, 1), col="blue")
