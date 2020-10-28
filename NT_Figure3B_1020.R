#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations on Figure 3B
#  Revision 10/20 d.bulgarelli@dundee.ac.uk  
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#load required packages (not all of these are needed for this so far but might be handy downstream)
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library("grid")
#ternary plot
source("tern_e.R")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#backdoor entry for JHI restrictions - THIS IS SRA SPECIFIC
#.libPaths("C:/Users/sr42605/R/library")

#set working directory: Senga
#setwd("C:/Users/sr42605/OneDrive - University of Dundee/RODRIGO FIGURE/data")

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#############################################################
#generate a phyloseq object using the normalised GO terms
#############################################################

#import the normalised GO terms
dat_info <- read.delim("cog_count_table.txt" , sep = "\t", header=T, row.names= 2, blank.lines.skip = FALSE)
dim(dat_info)
#remove the first column
dat_info_2 <- dat_info[, 2:13]
colnames(dat_info_2)
sum(colSums(dat_info_2))

#design file
design <- read.delim("NT_GO_annotations_map_2.txt", sep = "\t", header=TRUE, row.names=1)
design

#prune unplanted soil samples
design_2 <- design[colnames(dat_info_2), ]
design_2

#construct a DESeq dataset combining count data and sample information
#A DESeqDataSet object must have an associated design formula  The formula should be a tilde (???) followed by the variables of interest. In this case the column "Description" in the desing file depicts the variable of interest
NT_COG_cds <- DESeqDataSetFromMatrix(countData =dat_info_2, colData=design_2, design= ~ Genotype)

#execute the differential count analysis with the function DESeq 
NT_cds_test <- DESeq(NT_COG_cds, fitType="local", betaPrior=FALSE) 

#Genotype effect
ED <- results(NT_cds_test , contrast = c("Genotype",  "Elite", "Desert")) 
EN <- results(NT_cds_test , contrast = c("Genotype",  "Elite", "North"))
ND <- results(NT_cds_test , contrast = c("Genotype",  "North", "Desert"))

#FDR filtering
ED_FDR_005 <- ED[(rownames(ED)[which(ED$padj <0.05)]), ]
EN_FDR_005 <- EN[(rownames(EN)[which(EN$padj <0.05)]), ]
ND_FDR_005 <- ED[(rownames(ND)[which(ND$padj <0.05)]), ]

#Elite enriched
ED_FDR_005_elite <- ED_FDR_005[(rownames(ED_FDR_005)[which(ED_FDR_005$log2FoldChange > 0)]), ]
EN_FDR_005_elite <- EN_FDR_005[(rownames(EN_FDR_005)[which(EN_FDR_005$log2FoldChange > 0)]), ]
#Desert enriched
ED_FDR_005_desert <- ED_FDR_005[(rownames(ED_FDR_005)[which(ED_FDR_005$log2FoldChange < 0)]), ]
#North enriched
EN_FDR_005_north <- EN_FDR_005[(rownames(EN_FDR_005)[which(EN_FDR_005$log2FoldChange < 0)]), ]

#Genotype effect corrected
#Elite
Elite_effect <- intersect(row.names(EN_FDR_005_elite), row.names(ED_FDR_005_elite))
#wild genotypes
wild_effect <- intersect(row.names(ED_FDR_005_desert), row.names(EN_FDR_005_north))

#Ternary plot
#extract the means of the genotypes for plotting
NT_base_mean <- sapply(levels(NT_cds_test$Genotype), function(lvl) rowMeans( counts(NT_cds_test,normalized=TRUE)[,NT_cds_test$Genotype == lvl] ) )
#Build a ternary matrix for ternary plots
mean_North <- as.data.frame(NT_base_mean[, 4])
mean_Elite <- as.data.frame(NT_base_mean[, 3])
mean_Desert <- as.data.frame(NT_base_mean[, 2])
temat_1A <- cbind(mean_Desert, mean_North, mean_Elite)
colnames(temat_1A) <- c("Desert","North","Elite")
dim(temat_1A)
#colors and plotting
dev.off()
fig_colors <- ifelse(rownames(temat_1A) %in% Elite_effect, "magenta","darkgrey")
names(fig_colors) <- rownames(temat_1A)
fig_colors[wild_effect] <- "goldenrod4"
tern_e(temat_1A, scale = 1, prop=F, col=fig_colors, grid_color="black", labels_color="black", pch=23, main="COG enrichments")
#save this plot for figure generation


