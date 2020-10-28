#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations on Figure 4
#  Revision 10/20 SRobertsonalberty002@dundee.ac.uk  
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

#backdoor entry for JHI restrictions - THIS IS SRA SPECIFIC
#.libPaths("C:/Users/sr42605/R/library")

#set working directory: Senga
#setwd("C:/Users/sr42605/OneDrive - University of Dundee/RODRIGO FIGURE/data")

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#######################
#Pre-processing
######################

#Import the Phyloseq Object
JH02_JH16_phyloseq <- readRDS("NT_rare_25K_phyloseq_genus_0820_Silva138.rds")
JH02_JH16_phyloseq 

#subset for JH16 samples 
JH16_sample <- subset_samples(JH02_JH16_phyloseq, LibraryID == "JH16")
JH16_sample
sample_sums(JH16_sample)

#######################
#Figure 4A
######################

#subset for rhizosphere samples
JH16_rhizo <- subset_samples(JH16_sample, Microhabitat == "Rhizosphere")
#extract dryweight information
design_rhizo <- as.data.frame(as.matrix(sample_data(JH16_rhizo)))
#data distribution
hist(as.numeric(design_rhizo$Dryweight))
shapiro.test(as.numeric(design_rhizo$Dryweight))
#stat
stat <-aov(as.numeric(Dryweight) ~ Description, data = design_rhizo)
summary(stat)
TukeyHSD(stat)
#observed
design_rhizo$Description <- ordered(design_rhizo$Description, levels=c("Modern.N","Modern.A", "B1K.N","B1K.A"))
p <- ggplot(design_rhizo, aes(x=Description, y=as.numeric(Dryweight), fill=Description)) + geom_boxplot() 
p + geom_jitter( size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("magenta","blue","goldenrod", "goldenrod4"))

#export as EPS size 550 X 300

#######################
#Figure 4B
######################

#CAP analysis for plotting
JH16_CAP <- ordinate(JH16_sample, "CAP", "bray", ~ Soil * Microhabitat)
plot_ordination(JH16_sample, JH16_CAP, color = "Description", shape ="Treatment")

#assign shapes to Soil and color to Ecotype
p=plot_ordination(JH16_sample, JH16_CAP, shape ="Treatment", color = "Description")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_colour_manual(values = c("goldenrod4","black","goldenrod","blue", "brown","magenta"))
p + ggtitle("CAP 16S data, Bray distance")

#Permanova calculation rhizosphere samples
BC <- phyloseq::distance(JH16_rhizo, "bray")
adonis(BC ~ Soil * Treatment, data= design_rhizo, permutations = 5000)

#export as EPS size 550 X 300

#######################
#Figure 4C & D
######################

#create a DESeq object
#extract count data 
JH16_counts_integer <- otu_table(JH16_sample)
countData = as.data.frame(JH16_counts_integer)
colnames(JH16_counts_integer)

#the design file containing sample information
colData = as.data.frame(as.matrix(sample_data(JH16_sample)[colnames(JH16_counts_integer), ]))
rownames(colData)
class(colData)

#construct a DESeq dataset combining count data and sample information t
JH16_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH16_cds_test <- DESeq(JH16_cds, fitType="local", betaPrior=FALSE) 

#Morex conditioned calculation
##################

#define the Genera significantly enriched in the rhizosphere samples
Native_Rhizo <- results(JH16_cds_test, contrast = c("Description",  "Modern.B", "Modern.N")) 

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Native_Rhizo_FDR_005 <- Native_Rhizo[(rownames(Native_Rhizo)[which(Native_Rhizo$padj <0.05)]), ]

#Identify Genera enriched in the rhizosphere (second term of the comparison, negative fold change)
Native_Rhizo_enriched <-  Native_Rhizo[(rownames(Native_Rhizo)[which(Native_Rhizo$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the native rhizosphere; Morex conditioned
Native_Rhizo_enriched_FDR005 <- intersect(rownames(Native_Rhizo_FDR_005), rownames(Native_Rhizo_enriched))

#define the Genera significantly enriched in the rhizosphere samples vs autoclaved ones
Autoclaved_Rhizo <- results(JH16_cds_test, contrast = c("Description",  "Modern.A", "Modern.N")) 

# extract  Genra whose adjusted p.value in a given comparison is below 0.05. 
Autoclaved_Rhizo_FDR_005 <- Autoclaved_Rhizo[(rownames(Autoclaved_Rhizo)[which(Autoclaved_Rhizo$padj <0.05)]), ]

#Identify Genera enriched in the rhizosphere (second term of the comparison, negative fold change)
Autoclaved_Rhizo_enriched <-  Autoclaved_Rhizo[(rownames(Autoclaved_Rhizo)[which(Autoclaved_Rhizo$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the autoclaved rhizosphere; Morex conditioned
Autoclaved_Rhizo_enriched_FDR005 <- intersect(rownames(Autoclaved_Rhizo_FDR_005), rownames(Autoclaved_Rhizo_enriched))

#diagnostic genera for Morex conditioned; those are the genera simoultaneously enriched in Native rhizosphere versus bulk soil and autoclaved samples
Modern_Native_005 <- intersect(Native_Rhizo_enriched_FDR005, Autoclaved_Rhizo_enriched_FDR005)
length(Modern_Native_005)
#15

#B1K conditioned calculation
##################

#define the Genera significantly enriched in the rhizosphere samples
Native_Rhizo_b <- results(JH16_cds_test, contrast = c("Description",  "B1K.B", "B1K.N")) 

# extract  Genera whose adjusted p.value in a given comparison is below 0.05. 
Native_Rhizo_b_FDR_005 <- Native_Rhizo_b[(rownames(Native_Rhizo_b)[which(Native_Rhizo_b$padj <0.05)]), ]

#Identify Genera enriched in the rhizosphere (second term of the comparison, negative fold change)
Native_Rhizo_b_enriched <-  Native_Rhizo_b[(rownames(Native_Rhizo_b)[which(Native_Rhizo_b$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the native rhizosphere; Morex conditioned
Native_Rhizo_b_enriched_FDR005 <- intersect(rownames(Native_Rhizo_b_FDR_005), rownames(Native_Rhizo_b_enriched))

#define the Genera significantly enriched in the rhizosphere samples vs autoclaved ones
Autoclaved_Rhizo_b <- results(JH16_cds_test, contrast = c("Description",  "B1K.A", "B1K.N")) 

# extract  Genra whose adjusted p.value in a given comparison is below 0.05. 
Autoclaved_Rhizo_b_FDR_005 <- Autoclaved_Rhizo_b[(rownames(Autoclaved_Rhizo_b)[which(Autoclaved_Rhizo_b$padj <0.05)]), ]

#Identify Genera enriched in the rhizosphere (second term of the comparison, negative fold change)
Autoclaved_Rhizo_b_enriched <-  Autoclaved_Rhizo_b[(rownames(Autoclaved_Rhizo_b)[which(Autoclaved_Rhizo_b$log2FoldChange < 0)]), ]

#Intersection: this is the list of genera signficanlty enriched in the autoclaved rhizosphere; B1K conditioned
Autoclaved_Rhizo_b_enriched_FDR005 <- intersect(rownames(Autoclaved_Rhizo_b_FDR_005), rownames(Autoclaved_Rhizo_b_enriched))

#diagnostic genera for B1K conditioned; those are the genera simoultaneously enriched in Native rhizosphere versus bulk soil and autoclaved samples
B1K_Native_005 <- intersect(Native_Rhizo_b_enriched_FDR005, Autoclaved_Rhizo_b_enriched_FDR005)
length(B1K_Native_005)
#24

#######################
#Plotting
######################

#identify the genera enriched in both conditioned soils (Morex and B1K): we aim at testing the hypothesis of convergent enrichment
enriched_genera <- intersect(Modern_Native_005, B1K_Native_005)
#venn diagram constructed in Illustrator

#overlap with JH02 enrichment
Elite_enriched_N0 <- read.delim("NT_Figure_4_input_data_2.txt")
colnames(Elite_enriched_N0)
#overlap
overlap <- intersect(enriched_genera, Elite_enriched_N0$X.ASV_ID)
overlap

#Plotting taxa
#Define the enriched family
JH16_enriched <- prune_taxa(enriched_genera, JH16_sample)

#melting phyloseq object into a dataframe
df_family_JH16_all <- psmelt(JH16_enriched)
#plotting
plot_order_JH16_all <- ggplot(df_family_JH16_all, aes(Sample, Abundance, fill = fct_reorder(Class, Abundance))) + geom_col()+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#Order the samples
df_family_JH16_all$Sample <- factor(df_family_JH16_all$Sample,levels = c("AC13","AC14","AC62","AC64","AC65","AC67","AC68", "AC15", "AC16", "AC17", "AC35", "AC36", "AC37", "AC38", "AC39", "AC60", "AC63", "AC66", "AC01", "AC02", "AC03", "AC04", "AC05", "AC25", "AC26", "AC27", "AC28", "AC29", "AC50", "AC51", "AC52", "AC53", "AC54", "AC21", "AC22", "AC42", "AC43","AC46", "AC47", "AC49", "AC71", "AC72", "AC18", "AC19", "AC20", "AC23", "AC24", "AC40", "AC41", "AC44", "AC45", "AC48", "AC70", "AC73", "AC75", "AC76", "AC06", "AC07", "AC08", "AC09", "AC10", "AC30", "AC31", "AC32", "AC33", "AC34", "AC55", "AC56", "AC57", "AC58", "AC59"))
df_family_JH16_all_plot <- ggplot(df_family_JH16_all, aes(Sample, Abundance, fill = fct_reorder(Class, Abundance))) + geom_col() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
df_family_JH16_all_plot

#add cb friendly "grey" pallette
cbp1 <- c("black", "olivedrab3", "red", "olivedrab4","deepskyblue")
df_family_JH16_all_plot_cb <- df_family_JH16_all_plot + scale_fill_manual(breaks = c ("Polyangia", "Gammaproteobacteria", "Actinobacteria", "Alphaproteobacteria", "Bacteroidia"), values=cbp1)
df_family_JH16_all_plot_cb


#Save as EPS size 1000 X 500

#end





