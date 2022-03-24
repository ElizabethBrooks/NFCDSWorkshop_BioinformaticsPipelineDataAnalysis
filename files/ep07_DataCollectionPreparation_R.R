# set the working directory
setwd("/Users/bamflappy/NFCDSWorkshop_Day2/omics_data")

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

head(tribolium_counts)

rownames(tribolium_counts) <- rownames(cntrl1_fc_4h$counts)

head(tribolium_counts)

# export gene count data
write.csv(tribolium_counts, "/Users/bamflappy/NFCDSWorkshop_Day2/script_data/TriboliumCounts.csv", row.names=TRUE)
