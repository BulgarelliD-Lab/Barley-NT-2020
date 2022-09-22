#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations presented in Figure 3
#  Revision 09/22 d.bulgarelli@dundee.ac.uk 
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

#load required packages
library(ggplot2)
library(vegan)
library(phyloseq)

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#######################
#Pre-processing: Bacteria
######################

#Import the Kraken count matrices (metagenomics) 
NT_Bacteria_kraken <- readRDS("kraken_bacteria_unfiltered_barley_db.rds")
NT_Fungi_kraken <- readRDS("kraken_fungi_unfiltered_barley_db.rds")

#Figure S4A: overall reads assigned to phylum level

#Phylum level: Bacteria
NT_Bacteria_kraken_phylum <- as.data.frame(NT_Bacteria_kraken$phylum)
#inspect the file
dim(NT_Bacteria_kraken_phylum)
NT_Bacteria_kraken_phylum[, 1:5]
#prune the sample column
NT_Bacteria_kraken_phylum <- NT_Bacteria_kraken_phylum[, 2:73]
#stats
sum(rowSums(NT_Bacteria_kraken_phylum))
mean(rowSums(NT_Bacteria_kraken_phylum))
min(rowSums(NT_Bacteria_kraken_phylum))
max(rowSums(NT_Bacteria_kraken_phylum))

#create a vector with the number of reads classified 
Bacterial_reads <- rowSums(NT_Bacteria_kraken_phylum)

#create a dataset to visualise the proportion of reads per sample 
Bacterial_reads_d <- as.data.frame(Bacterial_reads)

#convert in log scale
Bacterial_reads_d <- log10(Bacterial_reads_d)

#rename the columns in the generated datasets
colnames(Bacterial_reads_d) <- c("Bacterial_reads")

#rename the rows of the dataset
names_new <- c("Bulk.1", "Bulk.2", "Bulk.3", "Desert.1", "Desert.2", "Desert.3", "North.1", "North.2", "North.3", "Elite.1", "Elite.2", "Elite.3")
row.names(Bacterial_reads_d) <- names_new

#create the mapping file
genotype <- c("Bulk", "Bulk", "Bulk", "Desert", "Desert", "Desert", "North", "North", "North", "Elite", "Elite", "Elite")
genotype <-as.data.frame(genotype)
rownames(genotype) <- names_new
genotype 

#combine these datasets with the design file
design_bacteria <- cbind(genotype, Bacterial_reads_d)

#order levels
design_bacteria$genotype <- ordered(design_bacteria$genotype , levels=c("Bulk","Desert","North","Elite"))

#Phylum level: Fungi
NT_Fungi_kraken_phylum <- as.data.frame(NT_Fungi_kraken$phylum)
#inspect the file
dim(NT_Fungi_kraken_phylum)
NT_Fungi_kraken_phylum[, 1:5]
#prune the sample column
NT_Fungi_kraken_phylum <- NT_Fungi_kraken_phylum[, 2:9]
#stats
sum(rowSums(NT_Fungi_kraken_phylum))
mean(rowSums(NT_Fungi_kraken_phylum))
min(rowSums(NT_Fungi_kraken_phylum))
max(rowSums(NT_Fungi_kraken_phylum))

#create a vector with the number of reads classified 
Fungal_reads <- rowSums(NT_Fungi_kraken_phylum)

#create a dataset to visualise the proportion of reads per sample 
Fungal_reads_d <- as.data.frame(Fungal_reads)

#convert in log scale
Fungal_reads_d <- log10(Fungal_reads_d)

#rename the columns in the generated datasets
colnames(Fungal_reads_d) <- c("Fungal_reads")

#rename the rows of the dataset
names_new <- c("Bulk.1", "Bulk.2", "Bulk.3", "Desert.1", "Desert.2", "Desert.3", "North.1", "North.2", "North.3", "Elite.1", "Elite.2", "Elite.3")
row.names(Fungal_reads_d) <- names_new

#combine these datasets with the design file
design_fungi <- cbind(genotype, Fungal_reads_d)

#order levels
design_fungi$genotype <- ordered(design_fungi$genotype , levels=c("Bulk","Desert","North","Elite"))

###############
#Figure 3A
###############

#combine these datasets with the design file
#revise the rownames of the files
rownames(design_fungi) <- c("Bulk.1F", "Bulk.2F", "Bulk.3F", "Desert.1F", "Desert.2F", "Desert.3F", "North.1F", "North.2F", "North.3F", "Elite.1F", "Elite.2F","Elite.3F")
rownames(design_bacteria) <- c("Bulk.1B", "Bulk.2B", "Bulk.3B", "Desert.1B", "Desert.2B", "Desert.3B", "North.1B", "North.2B", "North.3B",  "Elite.1B", "Elite.2B","Elite.3B")
#add an extra column specifying the Kingdom
#bacteria
KCB <- c("Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria","Bacteria")
KCB <- as.data.frame(KCB)
row.names(KCB) <-rownames(design_bacteria)
design_bacteria_2 <- cbind(KCB, design_bacteria)
colnames(design_bacteria_2) <- c("Kingdom", "genotype", "reads")
design_bacteria_2

#fungi
KCF <- c("Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi","Fungi")
KCF <- as.data.frame(KCF)
row.names(KCF) <-rownames(design_fungi)
design_fungi_2 <- cbind(KCF, design_fungi)
colnames(design_fungi_2) <- c("Kingdom", "genotype", "reads")

#merge the dataset
design_BF <- rbind(design_bacteria_2, design_fungi_2)
design_BF

#order levels
design_BF$genotype <- ordered(design_BF$genotype , levels=c("Bulk","Desert","North","Elite"))
design_BF$Kingdom <- ordered(design_BF$Kingdom , levels=c("Bacteria", "Fungi"))

#Plot
p <- ggplot(design_BF, aes(x=genotype, y=reads, fill=Kingdom))
p + geom_dotplot(binaxis='y', stackdir='center', binpositions="all", dotsize = 2) +  scale_fill_manual(values = c("#0072B2", "#56B4E9"))

#stats (computed on the independent dataset)
#assess sample distribution
shapiro.test(design_bacteria_2$reads)
shapiro.test(design_fungi_2$reads)

#parametric test
B_anova <- aov(reads~ genotype, design_bacteria_2)
summary(B_anova)
F_anova <- aov(reads~ genotype, design_fungi_2)
summary(F_anova)

###############
#Figure 3B
###############

#Class calculation
NT_Fungi_kraken_class <- as.data.frame(NT_Fungi_kraken$class)

#inspect the file
dim(NT_Fungi_kraken_class)
NT_Fungi_kraken_class[, 1:5]

#prune the sample column
NT_Fungi_kraken_class <- NT_Fungi_kraken_class[, 2:34]

#transpose the data set
NT_Fungi_kraken_class_t <- t(NT_Fungi_kraken_class) 
NT_Fungi_kraken_class_t[, 1:5]

#create an AMF dataset using Glomeromycetes as a proxy for AMF
NT_AMF_kraken <- as.data.frame(NT_Fungi_kraken_class_t["Glomeromycetes", ])
colnames(NT_AMF_kraken) <-c("AMF")
NT_AMF_kraken_t <- t(NT_AMF_kraken)

#convert in relative abundance
NT_AMF_kraken_t_2 <- (NT_AMF_kraken_t/colSums(NT_Fungi_kraken_class_t))*100
colSums(NT_AMF_kraken_t_2)

#combine these datasets with the design file
design_AMF <- cbind(genotype, t(NT_AMF_kraken_t_2))

#order levels
design_AMF$genotype <- ordered(design_AMF$genotype , levels=c("Bulk","Desert","North","Elite"))

#plotting
p <- ggplot(design_AMF, aes(x=genotype, y=AMF, fill=genotype)) + ylim (0,2)
p + geom_dotplot(binaxis='y', stackdir='center', binpositions="all", dotsize = 2) +  scale_fill_manual(values = c("#F0E442", "#F0E442", "#F0E442", "#F0E442"))

#######################################################
#Figure 3C
#######################################################

#Family calculation
NT_Bacteria_kraken_family <- as.data.frame(NT_Bacteria_kraken$family)

#inspect the file
dim(NT_Bacteria_kraken_family)
NT_Bacteria_kraken_family[, 1:5]

#prune the sample column
NT_Bacteria_kraken_family <- NT_Bacteria_kraken_family[, 2:367]
NT_Bacteria_kraken_family[1:5, ]

#create a vector with the number of reads familyified 
row.names(NT_Bacteria_kraken_family) <- names_new
NT_Bacteria_kraken_family[, 1:5]
NT_Bacteria_kraken_family_t <- t(NT_Bacteria_kraken_family)

#The counts table 
NT_kraken_Bacteria <- otu_table(NT_Bacteria_kraken_family_t, taxa_are_rows=TRUE)

#The mapping file
#add a microhabitat column
microhabitat <- c("soil", "soil", "soil", "rhizo", "rhizo", "rhizo", "rhizo", "rhizo", "rhizo", "rhizo", "rhizo", "rhizo")
design_gm <- cbind(genotype, microhabitat)
colnames(design_gm) <- c("Sample", "Microhabitat")
NT_kraken_map <- sample_data(design_gm)

#The phyloseq object
NT_kraken_Bacteria_phyloseq <- merge_phyloseq(NT_kraken_Bacteria, NT_kraken_map)
NT_kraken_Bacteria_phyloseq

#Transform the count in relative abundance
NT_kraken_Bacteria_phyloseq_prop <- transform_sample_counts(NT_kraken_Bacteria_phyloseq,  function(x) 1e+06 * x/sum(x))

#Compute a dissimilarity matrix using BC distance
BC <- phyloseq::distance(NT_kraken_Bacteria_phyloseq_prop, "bray")
#cluster dendrogram
NT_Bacteria_kraken_family_hc <- hclust(BC , method = "ward.D2")
#plotting
plot(NT_Bacteria_kraken_family_hc, hang = -1, cex = 1)

#stats sample
adonis(BC ~ Sample, data = design_gm, permutations = 5000)

#######################################################
#Figure 3D
#######################################################

#Family calculation
NT_Fungi_kraken_family <- as.data.frame(NT_Fungi_kraken$family)

#inspect the file
dim(NT_Fungi_kraken_family)
NT_Fungi_kraken_family[, 1:5]

#prune the sample column
NT_Fungi_kraken_family <- NT_Fungi_kraken_family[, 2:177]
NT_Fungi_kraken_family[1:5, ]

#create a vector with the number of reads familyified 
row.names(NT_Fungi_kraken_family) <- names_new
NT_Fungi_kraken_family[, 1:5]
NT_Fungi_kraken_family_t <- t(NT_Fungi_kraken_family)

#The counts table 
NT_kraken_Fungi <- otu_table(NT_Fungi_kraken_family_t, taxa_are_rows=TRUE)

#The mapping file: same as before, NT_kraken_map

#The phyloseq object
NT_kraken_Fungi_phyloseq <- merge_phyloseq(NT_kraken_Fungi, NT_kraken_map)
NT_kraken_Fungi_phyloseq

#Transform the count in relative abundance
NT_kraken_Fungi_phyloseq_prop <- transform_sample_counts(NT_kraken_Fungi_phyloseq,  function(x) 1e+06 * x/sum(x))

#Compute a dissimilarity matrix using BC distance
BC <- phyloseq::distance(NT_kraken_Fungi_phyloseq_prop, "bray")
#cluster dendrogram
NT_Fungi_kraken_family_hc <- hclust(BC , method = "ward.D2")
#plotting
plot(NT_Fungi_kraken_family_hc, hang = -1, cex = 1)

#stats microhabitat
adonis(BC ~ Sample, data = design_gm, permutation = 5000)

#END