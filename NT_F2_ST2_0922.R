#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations presented in:
#  Figure 2, Supplementary Figure 2 and Table S3
#  Revision 09/22 
#  rodrigo.alegria@um6p.ma
#  d.bulgarelli@dundee.ac.uk 
#
#############################################################

#####################################################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################


#required packages 
library(vegan)
library(phyloseq)
library(DESeq2)
library(UpSetR)
library(ggplot2)

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("/cluster/db/R_shared/NT2020_R_data/")

#import the file
JH02_JH16_data_phyloseq <- readRDS("NT_rare_25K_phyloseq_ASVs_Silva138_0820.rds")

#check total number of reads
sum(colSums(otu_table(JH02_JH16_data_phyloseq)))

#subset for Library JH02_NT (Figure 1)
JH02_NT_ASVs <- subset_samples(JH02_JH16_data_phyloseq, LibraryID =="JH02")

#check number of reads of JH02 library 
sum(colSums(otu_table(JH02_NT_ASVs)))


#remove the control samples from the dataset
JH02_NT_NoWA_ASVs <- subset_samples(JH02_NT_ASVs, Description !="Bulk.Agar")
colSums(otu_table(JH02_NT_NoWA_ASVs))

#examine metadata information
sample_data(JH02_NT_NoWA_ASVs)

########################################
#Figure 2A Betadiversity
########################################

#identify the levels
sample_data(JH02_NT_NoWA_ASVs)$Treatment
#Order the levels according to a desired sequence (e.g., first the soil then the ecotypes)
sample_data(JH02_NT_NoWA_ASVs)$Treatment <- ordered(sample_data(JH02_NT_NoWA_ASVs)$Treatment, levels=c("N0", "N25", "N100"))
sample_data(JH02_NT_NoWA_ASVs)$Genotype <- ordered(sample_data(JH02_NT_NoWA_ASVs)$Genotype, levels=c("Quarryfield", "Desert", "North", "Morex"))

#CAP Bray distance
JH02_CAP_BC <- ordinate(JH02_NT_NoWA_ASVs, "CAP", "bray", ~ Treatment * Genotype)
plot_ordination(JH02_NT_NoWA_ASVs, JH02_CAP_BC, color = "Genotype")

#assign shapes to Treatment and color to Genotype
p=plot_ordination(JH02_NT_NoWA_ASVs, JH02_CAP_BC, shape ="Treatment", color = "Genotype")
p = p + geom_point(size = 7, alpha = 0.75)
p = p + scale_colour_manual(values = c("#000000","#E69F00","#009E73","#CC79A7"))
p = p + scale_shape_manual(values = c(15,16,17))
p + ggtitle("CAP 16S data, Bray distance")

#permanova calculation
#All samples (included in the figure)
BC <- phyloseq::distance(JH02_NT_NoWA_ASVs, "bray")
adonis(BC ~ Genotype* Treatment , data= as.data.frame(as.matrix(sample_data(JH02_NT_NoWA_ASVs))), permutations = 5000)

########################################
#Figure 1B Bacterial enrichment
########################################

countData = as.data.frame(otu_table(JH02_NT_NoWA_ASVs))
class(countData)
colnames(countData)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH02_NT_NoWA_ASVs)))
class(colData)

#construct a DESeq dataset combining count data and sample information
#A DESeqDataSet object must have an associated design formula  The formula should be a tilde (???) followed by the variables of interest. In this case the column "Description" in the desing file depicts the variable of interest
NT_N0_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Description)

#execute the differential count analysis with the function DESeq 
NT_N0_cds_test <- DESeq(NT_N0_cds, fitType="local", betaPrior=FALSE) 

#N100
#define the taxa differentially enriched in the rhizosphere samples
Elite_Desert <- results(NT_N0_cds_test, contrast = c("Description", "Modern.100", "Desert.100")) 
Elite_North <- results(NT_N0_cds_test, contrast = c("Description",  "Modern.100", "North.100")) 
#North_Desert <- results(NT_N0_cds_test, contrast = c("Description",  "North.100", "Desert.100")) 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Elite_Desert_FDR_N100 <- Elite_Desert[(rownames(Elite_Desert)[which(Elite_Desert$padj <0.05)]), ]
Elite_North_FDR_N100 <- Elite_North[(rownames(Elite_North)[which(Elite_North$padj <0.05)]), ]
#North_Desert_FDR_N100 <- North_Desert[(rownames(North_Desert)[which(North_Desert$padj <0.05)]), ]

#extract base mean
Elite_Desert_FDR_N100_2 <- as.data.frame(Elite_Desert_FDR_N100[ ,1])
Elite_North_FDR_N100_2 <- as.data.frame(Elite_North_FDR_N100[ ,1])
#North_Desert_FDR_N100_2 <- as.data.frame(North_Desert_FDR_N100[ ,1])

#rename rows
rownames(Elite_Desert_FDR_N100_2) <- row.names(Elite_Desert_FDR_N100)
rownames(Elite_North_FDR_N100_2) <- row.names(Elite_North_FDR_N100)
#rownames(North_Desert_FDR_N100_2) <- row.names(North_Desert_FDR_N100)

#rename columns
colnames(Elite_Desert_FDR_N100_2) <- c("counts_ED_N100")
colnames(Elite_North_FDR_N100_2) <- c("counts_EN_N100")
#colnames(North_Desert_FDR_N100_2) <- c("counts_ND_N100")


Elite_Desert_FDR_N100_2[Elite_Desert_FDR_N100_2 > 1] <- 1
Elite_North_FDR_N100_2[Elite_North_FDR_N100_2 > 1] <- 1
#North_Desert_FDR_N100_2[North_Desert_FDR_N100_2 > 1] <- 1

#N25
#define the taxa differentially enriched in the rhizosphere samples
Elite_Desert <- results(NT_N0_cds_test, contrast = c("Description", "Modern.25", "Desert.25")) 
Elite_North <- results(NT_N0_cds_test, contrast = c("Description",  "Modern.25", "North.25")) 
#North_Desert <- results(NT_N0_cds_test, contrast = c("Description",  "North.25", "Desert.25")) 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Elite_Desert_FDR_N25 <- Elite_Desert[(rownames(Elite_Desert)[which(Elite_Desert$padj <0.05)]), ]
Elite_North_FDR_N25 <- Elite_North[(rownames(Elite_North)[which(Elite_North$padj <0.05)]), ]
#North_Desert_FDR_N25 <- North_Desert[(rownames(North_Desert)[which(North_Desert$padj <0.05)]), ]

#extract base mean
Elite_Desert_FDR_N25_2 <- as.data.frame(Elite_Desert_FDR_N25[ ,1])
Elite_North_FDR_N25_2 <- as.data.frame(Elite_North_FDR_N25[ ,1])
#North_Desert_FDR_N25_2 <- as.data.frame(North_Desert_FDR_N25[ ,1])

#rename rows
rownames(Elite_Desert_FDR_N25_2) <- row.names(Elite_Desert_FDR_N25)
rownames(Elite_North_FDR_N25_2) <- row.names(Elite_North_FDR_N25)
#rownames(North_Desert_FDR_N25_2) <- row.names(North_Desert_FDR_N25)

#rename columns
colnames(Elite_Desert_FDR_N25_2) <- c("counts_ED_N25")
colnames(Elite_North_FDR_N25_2) <- c("counts_EN_N25")
#colnames(North_Desert_FDR_N25_2) <- c("counts_ND_N25")


Elite_Desert_FDR_N25_2[Elite_Desert_FDR_N25_2 > 1] <- 1
Elite_North_FDR_N25_2[Elite_North_FDR_N25_2 > 1] <- 1
#North_Desert_FDR_N25_2[North_Desert_FDR_N25_2 > 1] <- 1

#N0
#define the taxa differentially enriched in the rhizosphere samples
Elite_Desert <- results(NT_N0_cds_test, contrast = c("Description", "Modern.0", "Desert.0")) 
Elite_North <- results(NT_N0_cds_test, contrast = c("Description",  "Modern.0", "North.0")) 
#North_Desert <- results(NT_N0_cds_test, contrast = c("Description",  "North.0", "Desert.0")) 

# extract  OTUs whose adjusted p.value in a given comparison is below 0.05. 
Elite_Desert_FDR_N0 <- Elite_Desert[(rownames(Elite_Desert)[which(Elite_Desert$padj <0.05)]), ]
Elite_North_FDR_N0 <- Elite_North[(rownames(Elite_North)[which(Elite_North$padj <0.05)]), ]
#North_Desert_FDR_N0 <- North_Desert[(rownames(North_Desert)[which(North_Desert$padj <0.05)]), ]

#extract base mean
Elite_Desert_FDR_N0_2 <- as.data.frame(Elite_Desert_FDR_N0[ ,1])
Elite_North_FDR_N0_2 <- as.data.frame(Elite_North_FDR_N0[ ,1])
#North_Desert_FDR_N0_2 <- as.data.frame(North_Desert_FDR_N0[ ,1])

#rename rows
rownames(Elite_Desert_FDR_N0_2) <- row.names(Elite_Desert_FDR_N0)
rownames(Elite_North_FDR_N0_2) <- row.names(Elite_North_FDR_N0)
#rownames(North_Desert_FDR_N0_2) <- row.names(North_Desert_FDR_N0)

#rename columns
colnames(Elite_Desert_FDR_N0_2) <- c("counts_ED_N0")
colnames(Elite_North_FDR_N0_2) <- c("counts_EN_N0")
#colnames(North_Desert_FDR_N0_2) <- c("counts_ND_N0")
dim(North_Desert_FDR_N0_2)

Elite_Desert_FDR_N0_2[Elite_Desert_FDR_N0_2 > 1] <- 1
Elite_North_FDR_N0_2[Elite_North_FDR_N0_2 > 1] <- 1
#North_Desert_FDR_N0_2[North_Desert_FDR_N0_2 > 1] <- 1

#define a list of unique Taxa
taxa_list_N100 <- unique(c(rownames(Elite_Desert_FDR_N100_2), rownames(Elite_North_FDR_N100_2)))
taxa_list_N25 <- unique(c(rownames(Elite_Desert_FDR_N25_2), rownames(Elite_North_FDR_N25_2)))
taxa_list_N0 <- unique(c(rownames(Elite_Desert_FDR_N0_2), rownames(Elite_North_FDR_N0_2))) 
taxa_list_all <- unique(c(c(taxa_list_N100, taxa_list_N25), taxa_list_N0))
length(taxa_list_all)

#create 9 new dataset for merging
#N100
Elite_Desert_FDR_N100_3 <- as.data.frame(Elite_Desert_FDR_N100_2[taxa_list_all, ])
Elite_North_FDR_N100_3 <- as.data.frame(Elite_North_FDR_N100_2[taxa_list_all, ])
#North_Desert_FDR_N100_3 <- as.data.frame(North_Desert_FDR_N100_2[taxa_list_all, ])

colnames(Elite_Desert_FDR_N100_3) <- c("ED_N100")
colnames(Elite_North_FDR_N100_3) <- c("EN_N100")
#colnames(North_Desert_FDR_N100_3) <- c("ND_N100")

row.names(Elite_Desert_FDR_N100_3) <- as.vector(taxa_list_all)
row.names(Elite_North_FDR_N100_3) <- as.vector(taxa_list_all)
#row.names(North_Desert_FDR_N100_3) <- as.vector(taxa_list_all)

#N25
Elite_Desert_FDR_N25_3 <- as.data.frame(Elite_Desert_FDR_N25_2[taxa_list_all, ])
Elite_North_FDR_N25_3 <- as.data.frame(Elite_North_FDR_N25_2[taxa_list_all, ])
#North_Desert_FDR_N25_3 <- as.data.frame(North_Desert_FDR_N25_2[taxa_list_all, ])

colnames(Elite_Desert_FDR_N25_3) <- c("ED_N25")
colnames(Elite_North_FDR_N25_3) <- c("EN_N25")
#colnames(North_Desert_FDR_N25_3) <- c("ND_N25")

row.names(Elite_Desert_FDR_N25_3) <- as.vector(taxa_list_all)
row.names(Elite_North_FDR_N25_3) <- as.vector(taxa_list_all)
#row.names(North_Desert_FDR_N25_3) <- as.vector(taxa_list_all)

#N0
Elite_Desert_FDR_N0_3 <- as.data.frame(Elite_Desert_FDR_N0_2[taxa_list_all, ])
Elite_North_FDR_N0_3 <- as.data.frame(Elite_North_FDR_N0_2[taxa_list_all, ])
#North_Desert_FDR_N0_3 <- as.data.frame(North_Desert_FDR_N0_2[taxa_list_all, ])

colnames(Elite_Desert_FDR_N0_3) <- c("ED_N0")
colnames(Elite_North_FDR_N0_3) <- c("EN_N0")
#colnames(North_Desert_FDR_N0_3) <- c("ND_N0")

row.names(Elite_Desert_FDR_N0_3) <- as.vector(taxa_list_all)
row.names(Elite_North_FDR_N0_3) <- as.vector(taxa_list_all)
#row.names(North_Desert_FDR_N0_3) <- as.vector(taxa_list_all)

#Merge the dataset
enriched_Taxa_N100 <- cbind(Elite_Desert_FDR_N100_3, Elite_North_FDR_N100_3)
enriched_Taxa_N25 <- cbind(Elite_Desert_FDR_N25_3, Elite_North_FDR_N25_3)
enriched_Taxa_N0 <- cbind(Elite_Desert_FDR_N0_3, Elite_North_FDR_N0_3)

enriched_Taxa_allN <- cbind(cbind(enriched_Taxa_N100, enriched_Taxa_N25), enriched_Taxa_N0)

#set NA to 0
enriched_Taxa_allN[is.na(enriched_Taxa_allN)] <- 0
enriched_Taxa_allN[, 1:5]
colnames(enriched_Taxa_allN)

#visualisation
upset(enriched_Taxa_allN, sets = c("ED_N100", "EN_N100", "ED_N25", "EN_N25", "ED_N0", "EN_N0"), sets.bar.color = "#0072B2",
      order.by = "freq")

#Elite effect
Elite_enriched_D <- Elite_Desert_FDR_N0[(rownames(Elite_Desert_FDR_N0)[which(Elite_Desert_FDR_N0$log2FoldChange > 1)]), ]
dim(Elite_enriched_D)
Elite_enriched_N <- Elite_North_FDR_N0[(rownames(Elite_North_FDR_N0)[which(Elite_North_FDR_N0$log2FoldChange > 1)]), ]
dim(Elite_enriched_N)

#B1K effect
Desert_enriched_E <- Elite_Desert_FDR_N0[(rownames(Elite_Desert_FDR_N0)[which(Elite_Desert_FDR_N0$log2FoldChange < -1)]), ]
dim(Desert_enriched_E)
North_enriched_E <- Elite_North_FDR_N0[(rownames(Elite_North_FDR_N0)[which(Elite_North_FDR_N0$log2FoldChange < -1)]), ] 
dim(North_enriched_E)

#########################################################
#Supplementary Table S2
#Define impact of 'agar plugs' on taxonomic composition
#########################################################

#generate phylum classification 
JH02_phylum <- tax_glom(JH02_NT_ASVs, taxrank= "Phylum", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

#subset for bulk soil samples
JH02_phylum_bulk <- subset_samples(JH02_phylum, Microhabitat=="Bulk") 

#filter taxa with no counts
JH02_phylum_bulk <- prune_taxa(taxa_sums(JH02_phylum_bulk) > 0, JH02_phylum_bulk) 

#extract the count matrix
JH02_phylum_bulk_count <- as.data.frame(otu_table(JH02_phylum_bulk))
JH02_phylum_bulk_count[, 1:3]

#generate a new design file for bulk soil samples
bulk_design <- as.data.frame(sample_data(JH02_phylum_bulk))
#make sure the order of samples is the same between count data and design file
bulk_design_2 <- bulk_design[colnames(JH02_phylum_bulk_count), ]
#replace the individual sample IDs with their attribute
colnames(JH02_phylum_bulk_count) <- bulk_design_2$Description
colnames(JH02_phylum_bulk_count)

#generate average values for the phylum counts
JH02_phylum_bulk_count_average <- as.data.frame(do.call(cbind, by(t(JH02_phylum_bulk_count),INDICES=names(JH02_phylum_bulk_count),FUN=colMeans)))
JH02_phylum_bulk_count_average[1:5, ]

#asses the normality of the distributions
shapiro.test(JH02_phylum_bulk_count_average$Bulk.0)
shapiro.test(JH02_phylum_bulk_count_average$Bulk.25)
shapiro.test(JH02_phylum_bulk_count_average$Bulk.100)
shapiro.test(JH02_phylum_bulk_count_average$Bulk.Agar)
#all below alpha 0.05 => Spearman's rank correlation then

#Stats
cor.test(JH02_phylum_bulk_count_average$Bulk.0, JH02_phylum_bulk_count_average$Bulk.Agar,
         method = "spearman")
cor.test(JH02_phylum_bulk_count_average$Bulk.25, JH02_phylum_bulk_count_average$Bulk.Agar,
         method = "spearman")
cor.test(JH02_phylum_bulk_count_average$Bulk.100, JH02_phylum_bulk_count_average$Bulk.Agar,
         method = "spearman")
#END
