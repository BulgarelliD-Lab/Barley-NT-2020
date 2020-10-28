#####################################################################################
#JH02_DADA2 CALCULATION  for NT manuscript by Figure 1 and S3 1020
#####################################################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################


#required packages 
library("phyloseq")
library("DESeq2")
library("limma")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("plyr")
library("VennDiagram")#not working in this r version, run from JH02_phyloseq_DADA2_07_20_RA.R, r version 3.5.2
library("treemap")


#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set the working directory
setwd("/cluster/db/ralegriaterrazas/R_data_JH02_DADA2")


#import the file
JH02_JH16_data_phyloseq <- readRDS("NT_rare_25K_phyloseq_genus_0820_Silva138.rds")

#check total number of reads
sum(colSums(otu_table(JH02_JH16_data_phyloseq)))

#subset for Library JH02_NT (Figure 1)
JH02_NT_ASVs <- subset_samples(JH02_JH16_data_phyloseq, LibraryID =="JH02")

#check number of reads of JH02 library 
sum(colSums(otu_table(JH02_NT_ASVs)))


#remove the control samples from the dataset
JH02_NT_NoWA_ASVs <- subset_samples(JH02_NT_ASVs, Description !="Bulk.Agar")
sum(colSums(otu_table(JH02_NT_NoWA_ASVs)))

#examine metadata information
sample_data(JH02_NT_NoWA_ASVs)

#examine the ASVs 
otu_table(JH02_NT_NoWA_ASVs)[1:10, 1:6]

#examine taxonomy
tax_table(JH02_NT_NoWA_ASVs)[1:1, 1:6]
JH02_NT_taxa <- tax_table(JH02_NT_NoWA_ASVs)
dim(tax_table(JH02_NT_NoWA_ASVs))

#save taxonomy information for format editing in Excel
#write.table(JH02_NT_taxa, file="JH02_NT_taxa_Genus.txt", sep="\t")

JH02_NT_taxa_ordered <- read.delim ("JH02_NT_taxa_Genus.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)


#extract ASV table
JH02_NT_NoWA_ASVs_table <- otu_table(JH02_NT_NoWA_ASVs)
colnames(JH02_NT_NoWA_ASVs_table)
rownames(JH02_NT_NoWA_ASVs_table)
sample_sums(JH02_NT_NoWA_ASVs_table)

#extract mapping file
design <- read.delim("JH02_NT_JH16_simplified_mapping.txt", sep = "\t", header=TRUE, row.names=1)

design_NT <- design[colnames(otu_table(JH02_NT_NoWA_ASVs)), ]

##########################################################################################

#PCoA bray distance
JH02_NT_NoWA_ASVs_bray <- ordinate(JH02_NT_NoWA_ASVs, "PCoA", "bray")
plot_ordination(JH02_NT_NoWA_ASVs,JH02_NT_NoWA_ASVs_bray  , color = "Treatment")

#assign shapes to Treatment type and color to Genotype
p=plot_ordination(JH02_NT_NoWA_ASVs, JH02_NT_NoWA_ASVs_bray , shape ="Treatment", color = "Genotype")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("orange","magenta","green","brown"))
p + ggtitle("PCoA 16S data, Bray distance")


###########
# Figure 1A: constrained ordinatiton for Description
###########

JH02_CAP <- ordinate(JH02_NT_NoWA_ASVs, "CAP", "bray", ~ Description)
plot_ordination(JH02_NT_NoWA_ASVs, JH02_CAP, color = "Description")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_NT_NoWA_ASVs,JH02_CAP , shape ="Treatment", color = "Genotype")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("orange","magenta","green","brown"))
p + ggtitle("CAP Description 16S data, Bray distance")

#anova on the axis
anova(JH02_CAP, permutations=5000)

##########
#Figute 1C:Permanova analysis
##########

BC <- phyloseq::distance(JH02_NT_NoWA_ASVs, "bray")
adonis(BC ~ Microhabitat*Treatment, data= design_NT, permutations = 5000)

#Host effect (rhizosphere only)
#Subsetting
JH02_NT_NoWA_ASVs_rhizo <- subset_samples(JH02_NT_NoWA_ASVs, Microhabitat == "Rhizosphere")
design_rhizosphere <- design_NT[colnames(otu_table(JH02_NT_NoWA_ASVs_rhizo)), ]

#BC distance
BC <- phyloseq::distance(JH02_NT_NoWA_ASVs_rhizo, "bray")
adonis(BC ~ Genotype * Treatment , data= design_rhizosphere, permutations = 5000)


########
#Deseq2: Differential abundance analysis for Figures 1(C,D)and Figure S3
########

#ASV table
JH02_OTU_counts_integer <- otu_table(JH02_NT_NoWA_ASVs)

countData = as.data.frame(JH02_OTU_counts_integer)
class(countData)
colnames(JH02_OTU_counts_integer)

#the design file containing sample information
colData = design_NT[colnames(JH02_OTU_counts_integer), ]
class(colData)
rownames(colData)

#construct a DESeq dataset combining count data and sample information
#A DESeqDataSet object must have an associated design formula  The formula should be a tilde followed by the variables of interest. In this case the column "Description" in the desing file depicts the variable of interest
JH02_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData, design= ~ Description)

#execute the differential count analysis with the function DESeq 
JH02_cds_test <- DESeq(JH02_cds, fitType="local", betaPrior=FALSE) 


#####
#Differential abundance analysis for Figure S3
#####

#define the Genus significantly enriched in the rhizosphere samples
DESeq_OUTPUT_rhizosphere_D0 <- results(JH02_cds_test , contrast = c("Description", "Bulk.0","Desert.0")) 
DESeq_OUTPUT_rhizosphere_N0 <- results(JH02_cds_test , contrast = c("Description", "Bulk.0","North.0"))
DESeq_OUTPUT_rhizosphere_M0 <- results(JH02_cds_test , contrast = c("Description", "Bulk.0","Modern.0"))
DESeq_OUTPUT_rhizosphere_D25 <- results(JH02_cds_test , contrast = c("Description", "Bulk.25","Desert.25")) 
DESeq_OUTPUT_rhizosphere_N25 <- results(JH02_cds_test , contrast = c("Description", "Bulk.25","North.25"))
DESeq_OUTPUT_rhizosphere_M25 <- results(JH02_cds_test, contrast = c("Description", "Bulk.25","Modern.25"))
DESeq_OUTPUT_rhizosphere_D100 <- results(JH02_cds_test, contrast = c("Description", "Bulk.100","Desert.100")) 
DESeq_OUTPUT_rhizosphere_N100 <- results(JH02_cds_test , contrast = c("Description", "Bulk.100","North.100"))
DESeq_OUTPUT_rhizosphere_M100 <- results(JH02_cds_test , contrast = c("Description", "Bulk.100","Modern.100"))

#inspect the result file
DESeq_OUTPUT_rhizosphere_D0  
mcols(DESeq_OUTPUT_rhizosphere_D0   , use.names=TRUE)


#Extract Genera whose adjusted p.value in a given comparison is below 0.05.
DESeq_OUTPUT_rhizosphere_FDR005_D0 <-DESeq_OUTPUT_rhizosphere_D0[(rownames(DESeq_OUTPUT_rhizosphere_D0)[which(DESeq_OUTPUT_rhizosphere_D0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_N0 <-DESeq_OUTPUT_rhizosphere_N0[(rownames(DESeq_OUTPUT_rhizosphere_N0)[which(DESeq_OUTPUT_rhizosphere_N0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_M0 <-DESeq_OUTPUT_rhizosphere_M0[(rownames(DESeq_OUTPUT_rhizosphere_M0)[which(DESeq_OUTPUT_rhizosphere_M0$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_D25 <-DESeq_OUTPUT_rhizosphere_D25[(rownames(DESeq_OUTPUT_rhizosphere_D25)[which(DESeq_OUTPUT_rhizosphere_D25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_N25 <-DESeq_OUTPUT_rhizosphere_N25[(rownames(DESeq_OUTPUT_rhizosphere_N25)[which(DESeq_OUTPUT_rhizosphere_N25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_M25 <-DESeq_OUTPUT_rhizosphere_M25[(rownames(DESeq_OUTPUT_rhizosphere_M25)[which(DESeq_OUTPUT_rhizosphere_M25$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_D100 <-DESeq_OUTPUT_rhizosphere_D100[(rownames(DESeq_OUTPUT_rhizosphere_D100)[which(DESeq_OUTPUT_rhizosphere_D100$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_N100 <-DESeq_OUTPUT_rhizosphere_N100[(rownames(DESeq_OUTPUT_rhizosphere_N100)[which(DESeq_OUTPUT_rhizosphere_N100$padj <0.05)]), ]
DESeq_OUTPUT_rhizosphere_FDR005_M100 <-DESeq_OUTPUT_rhizosphere_M100[(rownames(DESeq_OUTPUT_rhizosphere_M100)[which(DESeq_OUTPUT_rhizosphere_M100$padj <0.05)]), ]

#Identify Genera enriched in the rhizosphere (second term of the comparison, negative fold change)
DESeq_OUTPUT_rhizosphere_enriched_D0 <- DESeq_OUTPUT_rhizosphere_D0[(rownames(DESeq_OUTPUT_rhizosphere_D0)[which(DESeq_OUTPUT_rhizosphere_D0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N0 <- DESeq_OUTPUT_rhizosphere_N0[(rownames(DESeq_OUTPUT_rhizosphere_N0)[which(DESeq_OUTPUT_rhizosphere_N0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M0 <- DESeq_OUTPUT_rhizosphere_M0[(rownames(DESeq_OUTPUT_rhizosphere_M0)[which(DESeq_OUTPUT_rhizosphere_M0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M0 <- DESeq_OUTPUT_rhizosphere_M0[(rownames(DESeq_OUTPUT_rhizosphere_M0)[which(DESeq_OUTPUT_rhizosphere_M0$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_D25 <- DESeq_OUTPUT_rhizosphere_D25[(rownames(DESeq_OUTPUT_rhizosphere_D25)[which(DESeq_OUTPUT_rhizosphere_D25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N25 <- DESeq_OUTPUT_rhizosphere_N25[(rownames(DESeq_OUTPUT_rhizosphere_N25)[which(DESeq_OUTPUT_rhizosphere_N25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M25 <- DESeq_OUTPUT_rhizosphere_M25[(rownames(DESeq_OUTPUT_rhizosphere_M25)[which(DESeq_OUTPUT_rhizosphere_M25$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_D100 <- DESeq_OUTPUT_rhizosphere_D100[(rownames(DESeq_OUTPUT_rhizosphere_D100)[which(DESeq_OUTPUT_rhizosphere_D100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_N100 <- DESeq_OUTPUT_rhizosphere_N100[(rownames(DESeq_OUTPUT_rhizosphere_N100)[which(DESeq_OUTPUT_rhizosphere_N100$log2FoldChange < 0)]), ]
DESeq_OUTPUT_rhizosphere_enriched_M100 <- DESeq_OUTPUT_rhizosphere_M100[(rownames(DESeq_OUTPUT_rhizosphere_M100)[which(DESeq_OUTPUT_rhizosphere_M100$log2FoldChange < 0)]), ]

#Identify OTUs enriched in the rhizosphere (second term of the comparison, negative fold change)
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D0), rownames(DESeq_OUTPUT_rhizosphere_enriched_D0))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N0), rownames(DESeq_OUTPUT_rhizosphere_enriched_N0))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M0), rownames(DESeq_OUTPUT_rhizosphere_enriched_M0))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D25), rownames(DESeq_OUTPUT_rhizosphere_enriched_D25))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N25), rownames(DESeq_OUTPUT_rhizosphere_enriched_N0))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M25), rownames(DESeq_OUTPUT_rhizosphere_enriched_M25))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_D100), rownames(DESeq_OUTPUT_rhizosphere_enriched_D100))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_N100), rownames(DESeq_OUTPUT_rhizosphere_enriched_N100))
DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 <- intersect(rownames(DESeq_OUTPUT_rhizosphere_FDR005_M100), rownames(DESeq_OUTPUT_rhizosphere_enriched_M100))

########
# Figure S3-A. Venn_diagram Rhizosphere effect
#######

#Identify Genera enriched at each N treatment
Enriched_0_rhizo <-unique(c(DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0, DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0))
Enriched_25_rhizo <-unique(c(DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25, DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25))
Enriched_100_rhizo <-unique(c(DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100, DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100))


#Save files
#saveRDS(Enriched_0_rhizo,"Enriched_0_rhizo.rds")
#saveRDS(Enriched_25_rhizo,"Enriched_25_rhizo.rds")
#saveRDS(Enriched_100_rhizo,"Enriched_100_rhizo.rds")


#to belaunched from r version 3.5.2

#Enriched_0R<-readRDS("Enriched_0_rhizo.rds")
#Enriched_25R<-readRDS("Enriched_25_rhizo.rds")
#Enriched_100R<-readRDS("Enriched_100_rhizo.rds")

#Construct Venn diagram
#venn.diagram(
# x = list( Enriched_0R, Enriched_25R,Enriched_100R),
#category.names = c( "N0" , "N25","N100"), filename= "Venn_Rhizosphere.png",
#output=TRUE)

##########
##Figure S3 Treemap taxonomy
#########

#Retrive taxonomy exclusive to N0% rhizo
#select enriched at N25_N100
Enriched_25_100_rhizo <-unique(c(Enriched_25_rhizo, Enriched_100_rhizo))

#select uniquely enriched at N0
Enriched_0_rhizo_unique <- Enriched_0_rhizo[!Enriched_0_rhizo %in% Enriched_25_100_rhizo]

#Prune enriched at N25_N100 from N0
Enriched_0_rhizo_unique2 <- Enriched_0_rhizo_unique[!Enriched_0_rhizo_unique %in%Enriched_25_100_rhizo]


#Retrieve tax info
Enriched_0_rhizo_unique_taxonomy <- JH02_NT_taxa_ordered[Enriched_0_rhizo_unique2,  ]

#Save taxonomy table to edited in Excel
#write.table(Enriched_0_rhizo_unique_taxonomy, "Enriched_0_rhizo_unique_taxonomy.txt")


# Build Dataset for treemap constructed in excel from Enriched_0_rhizo_unique_taxonomy according to treemap parameters

tm_subgroup_rhizo<- read.table("tm_subgroup_rhizo.txt",header= T)

tm_subgroup_rhizo_df <- data.frame(tm_subgroup_rhizo)

# treemap
treemap(tm_subgroup_rhizo_df ,
        index=c("Phylum","Class"),
        vSize="Value",
        type="index",
        
        fontsize.labels=c(10,10),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,2),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("right", "bottom"), 
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
        
) 

####
##Figure 1C Genotypes comparison , differential abundance amalysis
####

#define the Genus significantly enriched in the genotypes cpmoarison
DESeq_OUTPUT_D0_M0 <- results(JH02_cds_test , contrast = c("Description", "Desert.0","Modern.0")) 
DESeq_OUTPUT_N0_D0 <- results(JH02_cds_test , contrast = c("Description", "North.0","Desert.0"))
DESeq_OUTPUT_N0_M0 <- results(JH02_cds_test , contrast = c("Description", "North.0","Modern.0"))
DESeq_OUTPUT_D25_M25 <- results(JH02_cds_test , contrast = c("Description", "Desert.25","Modern.25")) 
DESeq_OUTPUT_N25_D25 <- results(JH02_cds_test , contrast = c("Description", "North.25","Desert.25"))
DESeq_OUTPUT_N25_M25 <- results(JH02_cds_test , contrast = c("Description", "North.25","Modern.25")) 
DESeq_OUTPUT_D100_M100 <- results(JH02_cds_test , contrast = c("Description", "Desert.100","Modern.100")) 
DESeq_OUTPUT_N100_D100 <- results(JH02_cds_test , contrast = c("Description", "North.100","Desert.100"))
DESeq_OUTPUT_N100_M100 <- results(JH02_cds_test , contrast = c("Description", "North.100","Modern.100")) 

#Genotype enriched
#D0 vs M0
# extract the Genus whose adjusted p.value in a given comparison is below 0.05.
DESeq_OUTPUT_FDR005_D0_M0 <- DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_D0_M0 <- DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M0_D0 <-  DESeq_OUTPUT_D0_M0[(rownames(DESeq_OUTPUT_D0_M0)[which(DESeq_OUTPUT_D0_M0$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_D0_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_M0), rownames(DESeq_OUTPUT_enriched_D0_M0))
DESeq_OUTPUT_enriched_FDR005_M0_D0 <- intersect(rownames(DESeq_OUTPUT_FDR005_D0_M0), rownames(DESeq_OUTPUT_enriched_M0_D0))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )

#N0 vs D0 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05.
DESeq_OUTPUT_FDR005_N0_D0 <- DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N0_D0 <- DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_D0_N0 <-  DESeq_OUTPUT_N0_D0[(rownames(DESeq_OUTPUT_N0_D0)[which(DESeq_OUTPUT_N0_D0$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N0_D0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_D0), rownames(DESeq_OUTPUT_enriched_N0_D0))
DESeq_OUTPUT_enriched_FDR005_D0_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_D0), rownames(DESeq_OUTPUT_enriched_D0_N0))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_D0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_D0_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D0 )

#N0 vs M0 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_N0_M0<- DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N0_M0 <- DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M0_N0 <-  DESeq_OUTPUT_N0_M0[(rownames(DESeq_OUTPUT_N0_M0)[which(DESeq_OUTPUT_N0_M0$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N0_M0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_M0), rownames(DESeq_OUTPUT_enriched_N0_M0))
DESeq_OUTPUT_enriched_FDR005_M0_N0 <- intersect(rownames(DESeq_OUTPUT_FDR005_N0_M0), rownames(DESeq_OUTPUT_enriched_M0_N0))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0 <- intersect(DESeq_OUTPUT_enriched_FDR005_N0_M0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N0 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0 <- intersect(DESeq_OUTPUT_enriched_FDR005_M0_N0,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M0 )


#D25 vs M25 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_D25_M25 <- DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_D25_M25 <- DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M25_D25 <-  DESeq_OUTPUT_D25_M25[(rownames(DESeq_OUTPUT_D25_M25)[which(DESeq_OUTPUT_D25_M25$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_D25_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_M25), rownames(DESeq_OUTPUT_enriched_D25_M25))
DESeq_OUTPUT_enriched_FDR005_M25_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_D25_M25), rownames(DESeq_OUTPUT_enriched_M25_D25))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )

#N25 vs D25 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_N25_D25 <- DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N25_D25 <- DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_D25_N25 <-  DESeq_OUTPUT_N25_D25[(rownames(DESeq_OUTPUT_N25_D25)[which(DESeq_OUTPUT_N25_D25$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N25_D25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_D25), rownames(DESeq_OUTPUT_enriched_N25_D25))
DESeq_OUTPUT_enriched_FDR005_D25_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_D25), rownames(DESeq_OUTPUT_enriched_D25_N25))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_D25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_D25_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D25 )

#N25 vs M25 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_N25_M25<- DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N25_M25 <- DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M25_N25 <-  DESeq_OUTPUT_N25_M25[(rownames(DESeq_OUTPUT_N25_M25)[which(DESeq_OUTPUT_N25_M25$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N25_M25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_M25), rownames(DESeq_OUTPUT_enriched_N25_M25))
DESeq_OUTPUT_enriched_FDR005_M25_N25 <- intersect(rownames(DESeq_OUTPUT_FDR005_N25_M25), rownames(DESeq_OUTPUT_enriched_M25_N25))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25 <- intersect(DESeq_OUTPUT_enriched_FDR005_N25_M25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N25 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25 <- intersect(DESeq_OUTPUT_enriched_FDR005_M25_N25,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M25 )

#D100 vs M100 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_D100_M100 <- DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_D100_M100 <- DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M100_D100 <-  DESeq_OUTPUT_D100_M100[(rownames(DESeq_OUTPUT_D100_M100)[which(DESeq_OUTPUT_D100_M100$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_D100_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D100_M100), rownames(DESeq_OUTPUT_enriched_D100_M100))
DESeq_OUTPUT_enriched_FDR005_M100_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_D100_M100), rownames(DESeq_OUTPUT_enriched_M100_D100))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )

#N100 vs D100 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_N100_D100 <- DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N100_D100 <- DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_D100_N100 <-  DESeq_OUTPUT_N100_D100[(rownames(DESeq_OUTPUT_N100_D100)[which(DESeq_OUTPUT_N100_D100$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N100_D100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_D100), rownames(DESeq_OUTPUT_enriched_N100_D100))
DESeq_OUTPUT_enriched_FDR005_D100_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_D100), rownames(DESeq_OUTPUT_enriched_D100_N100))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_D100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_D100_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_D100 )

#N100 vs M100 
# extract the Genus whose adjusted p.value in a given comparison is below 0.05
DESeq_OUTPUT_FDR005_N100_M100<- DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$padj <0.05)]), ]

#Identify Genus enriched in the first and second terms of the comparison
DESeq_OUTPUT_enriched_N100_M100 <- DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$log2FoldChange > 0)]), ]
DESeq_OUTPUT_enriched_M100_N100 <-  DESeq_OUTPUT_N100_M100[(rownames(DESeq_OUTPUT_N100_M100)[which(DESeq_OUTPUT_N100_M100$log2FoldChange < 0)]), ]

#intersect the datasets
DESeq_OUTPUT_enriched_FDR005_N100_M100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_M100), rownames(DESeq_OUTPUT_enriched_N100_M100))
DESeq_OUTPUT_enriched_FDR005_M100_N100 <- intersect(rownames(DESeq_OUTPUT_FDR005_N100_M100), rownames(DESeq_OUTPUT_enriched_M100_N100))

#further filtering for Genus significantly enriched vs. soil in the respective comparision
DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100 <- intersect(DESeq_OUTPUT_enriched_FDR005_N100_M100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_N100 )
DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100 <- intersect(DESeq_OUTPUT_enriched_FDR005_M100_N100,DESeq_OUTPUT_rhizosphere_enriched_FDR005_M100 )


#################################################################################################################################
## Venn_diagram Genotype comaparisons Figure 1C
#################################################################################################################################

#Retrive Genus enriched for each Genotype at N0%
Enriched_M_0 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_M0_D0,DESeq_OUTPUT_enriched_rhizo_FDR005_M0_N0))
Enriched_D_0 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_D0_M0, DESeq_OUTPUT_enriched_rhizo_FDR005_D0_N0))
Enriched_N_0 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_N0_D0,DESeq_OUTPUT_enriched_rhizo_FDR005_N0_M0))

Enriched_0_gen <- unique(c(Enriched_M_0,Enriched_D_0,Enriched_N_0))

#Retrive Genus enriched for each Genotype at N25%
Enriched_M_25 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_M25_D25,DESeq_OUTPUT_enriched_rhizo_FDR005_M25_N25))
Enriched_D_25 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_D25_M25, DESeq_OUTPUT_enriched_rhizo_FDR005_D25_N25))
Enriched_N_25 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_N25_D25,DESeq_OUTPUT_enriched_rhizo_FDR005_N25_M25))

Enriched_25_gen <- unique(c(Enriched_M_25,Enriched_D_25,Enriched_N_25))

#Retrive Genus enriched for each Genotype at N100%
Enriched_M_100 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_M100_D100,DESeq_OUTPUT_enriched_rhizo_FDR005_M100_N100))
Enriched_D_100 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_D100_M100, DESeq_OUTPUT_enriched_rhizo_FDR005_D100_N100))
Enriched_N_100 <-unique(c(DESeq_OUTPUT_enriched_rhizo_FDR005_N100_D100,DESeq_OUTPUT_enriched_rhizo_FDR005_N100_M100))

Enriched_100_gen <- unique(c(Enriched_M_100,Enriched_D_100,Enriched_N_100))

#Save files for its use in r version 3.5.2
#saveRDS(Enriched_0_gen,"Enriched_0_gen.rds")
#saveRDS(Enriched_25_gen,"Enriched_25_gen.rds")
#saveRDS(Enriched_100_gen,"Enriched_100_gen.rds")


Enriched_0_gen<-readRDS("Enriched_0_gen.rds")
Enriched_25_gen<-readRDS("Enriched_25_gen.rds")
Enriched_100_gen<-readRDS("Enriched_100_gen.rds")

#Construct Venn Diagram
#venn.diagram(
#  x = list( Enriched_0_gen, Enriched_25_gen,Enriched_100_gen),
#  category.names = c( "N0 ","N25","N100"), filename= "Venn_Genotype.png",
#  output=TRUE)


#######
##Treemap for taxonomy Figure 1D
######

##Retrieve Genera enriched at genotype comparison at all N treatments
Enriched_0_25_gen <-intersect(Enriched_0_gen, Enriched_25_gen)
Enriched_0_100_gen <- intersect(Enriched_0_gen, Enriched_100_gen)
Enriched_0_25_100_gen <-c(Enriched_0_25_gen,Enriched_0_100_gen) 

#Select uniquely enriched at N0%
Enriched_0_gen_unique<- Enriched_0_gen[!Enriched_0_gen %in% Enriched_0_25_100_gen]

#Retrieve taxonomy
Enriched_0_taxonomy_gen <- JH02_NT_taxa_ordered[Enriched_0_gen_unique,  ]

#write.table(Enriched_0_taxonomy_gen, "Enriched_0_taxonomy_gen.txt")


# Build Dataset in excel from Enriched_0_taxonomy_gen according to treemap parameters

tm_subgroup <- read.table("tm_subgroup.txt", header=T)

tm_subgroup_df <- data.frame(tm_subgroup)

# treemap
treemap(tm_subgroup_df ,
        index=c("Phylum","Class"),
        vSize="Value",
        type="index",
        
        fontsize.labels=c(12,12),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
        fontcolor.labels=c("white","black"),    # Color of labels
        fontface.labels=c(2,2),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
        bg.labels=c("transparent"),              # Background color of labels
        align.labels=list(
          c("right", "bottom"), 
          c("center", "center")
        ),                                   # Where to place labels in the rectangle?
        overlap.labels=0.5,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
        inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
        
) 


##################################################################################################################







