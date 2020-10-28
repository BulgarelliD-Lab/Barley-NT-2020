#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations presented on Table S3
#  Revision 09/20 d.bulgarelli@dundee.ac.uk 
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
library ("forcats")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#######################
#Pre-processing
######################

#Import the Phyloseq Object
JH02_JH16_phyloseq <- readRDS("NT_rare_25K_phyloseq_genus_0820_Silva138.rds")
JH02_JH16_phyloseq 

#subset for JH02 samples 
JH02_sample <- subset_samples(JH02_JH16_phyloseq, LibraryID == "JH02")


#######################
#Table S3
######################

#generate phylum classification 
JH02_phylum <- tax_glom(JH02_sample, taxrank= "Phylum", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))

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

#End
