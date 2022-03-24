# set the working directory
setwd("/Users/bamflappy/NFCDSWorkshop_Day2/script_data/")

# import gene count data
tribolium_counts <- read.csv("TriboliumCounts.csv", row.names="X")

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
jpeg("glm_librarySizes_beforeNorm.jpg")
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")
dev.off()

#Retain genes only if it is expressed at a minimum level
keep <- filterByExpr(list)
summary(keep)
list <- list[keep, , keep.lib.sizes=FALSE]

#Use TMM normalization to eliminate composition biases
# between libraries
list <- calcNormFactors(list)

#Retrieve normalized counts
normList <- cpm(list, normalized.lib.sizes=TRUE)

#Write the normalized counts to a file
write.table(normList, file="tribolium_normalizedCounts.csv", sep=",", row.names=TRUE)

#View normalization factors
list$samples
dim(list)

#Verify TMM normalization using a MD plot
# and write the plot to jpg file
jpeg("glm_MD_afterNorm.jpg")
plotMD(cpm(list, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)
dev.off()

#Use a MDS plot to visualizes the differences
# between the expression profiles of different samples
points <- c(0,1,2,3,15,16,17,18)
colors <- rep(c("blue", "darkgreen", "red", "black"), 2)

#Write plot without legend to file
jpeg("glm_MDS_withoutLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
dev.off()

#Write plot with legend to file
jpeg("glm_MDS_withLegend.jpg")
plotMDS(list, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
dev.off()

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
jpeg("glm_BCV.jpg")
plotBCV(list)
dev.off()

#Now, estimate and plot the QL dispersions
fit <- glmQLFit(list, design, robust=TRUE)
#head(fit$coefficients)

#Write plot to file
jpeg("glm_QLDisp.jpg")
plotQLDisp(fit)
dev.off()


########
#Test whether the average across all cntrl groups is equal to the average across
#all treat groups, to examine the overall effect of treatment
con.treat_cntrl <- makeContrasts(set.treat_cntrl = (treat.4h + treat.24h)/2
                                                    - (cntrl.4h + cntrl.24h)/2,
                                                    levels=design)

#Look at genes with significant expression across all UV groups
anov.treat_cntrl <- glmTreat(fit, contrast=con.treat_cntrl, lfc=log2(1.2))
summary(decideTests(anov.treat_cntrl))

#Write plot to file
jpeg("glm_treat_cntrl_MD.jpg")
plotMD(anov.treat_cntrl)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_treat_cntrl.filtered <- topTags(anov.treat_cntrl, n=nrow(anov.treat_cntrl$table), adjust.method="fdr")$table
write.table(tagsTbl_treat_cntrl.filtered, file="glm_treat_cntrl.csv", sep=",", row.names=TRUE)


########

#Test whether the average across all tolerant groups is equal to the average across
#all not tolerant groups, to examine the overall effect of tolerance
con.24h_4h <- makeContrasts(set.24h_4h = (cntrl.24h + treat.24h)/2
                                          - (cntrl.4h + treat.4h)/2,
                                          levels=design)

#Look at genes with significant expression across all UV groups
anov.24h_4h <- glmTreat(fit, contrast=con.24h_4h, lfc=log2(1.2))
summary(decideTests(anov.24h_4h))

#Write plot to file
jpeg("glm_24h_4h_MD.jpg")
plotMD(anov.24h_4h)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_24h_4h.filtered <- topTags(anov.24h_4h, n=nrow(anov.24h_4h$table), adjust.method="fdr")$table
write.table(tagsTbl_24h_4h.filtered, file="glm_24h_4h.csv", sep=",", row.names=TRUE)


########
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
jpeg("glm_interaction_MD.jpg")
plotMD(anov.interaction)
abline(h=c(-1, 1), col="blue")
dev.off()

#Generate table of DE genes
tagsTbl_inter.filtered <- topTags(anov.interaction, n=nrow(anov.interaction$table), adjust.method="fdr")$table
write.table(tagsTbl_inter.filtered, file="glm_interaction.csv", sep=",", row.names=TRUE)
