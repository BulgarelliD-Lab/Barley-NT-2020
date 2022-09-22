#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations on Figure 7 and Supplementary Figure 8
#  Revision 09/22 SRobertsonalberty002@dundee.ac.uk
#  d.bulgarelli@dundee.ac.uk      
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
library ("PMCMR")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()


#set working directory: Senga
#setwd("C:/Users/sr42605/OneDrive - University of Dundee/RODRIGO FIGURE/data")

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#######################
#Pre-processing
######################

#Import the Phyloseq Object
JH02_JH16_phyloseq <- readRDS("NT_rare_25K_phyloseq_ASVs_Silva138_0820.rds")
JH02_JH16_phyloseq 

#subset for JH16 samples 
JH16_sample <- subset_samples(JH02_JH16_phyloseq, LibraryID == "JH16")
JH16_sample
sample_sums(JH16_sample)

#######################
#Figure 7B
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
p + geom_jitter( size=4,shape=22, position=position_jitter(0.2))+ scale_fill_manual(values = c("#CC79A7","#0072B2","#E69F00", "#56B4E9"))

########################
#Figure 7C
########################

#CAP analysis for plotting
JH16_CAP <- ordinate(JH16_sample, "CAP", "bray", ~ Soil * Microhabitat)
plot_ordination(JH16_sample, JH16_CAP, color = "Description", shape ="Treatment")

#assign shapes to Soil and color to Ecotype
p=plot_ordination(JH16_sample, JH16_CAP, color = "Description", shape = "Description")
p = p + geom_point(size = 5, alpha = 0.75)
p = p + scale_shape_manual(values = c(15, 15, 15, 15, 15, 15))
p = p + scale_colour_manual(values = c("#56B4E9","black","#E69F00","#0072B2", "white","#CC79A7"))
p + ggtitle("CAP 16S data, Bray distance")

#Permanova calculation rhizosphere samples
BC <- phyloseq::distance(JH16_rhizo, "bray")
adonis(BC ~ Soil * Treatment, data= design_rhizo, permutations = 5000)

#Permutation: free
#Number of permutations: 5000

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Soil            1    0.1845  0.1845  1.0751 0.01446 0.3223    
#Treatment       1    3.6569  3.6569 21.3142 0.28664 0.0002 ***
#  Soil:Treatment  1    0.1664  0.1664  0.9701 0.01305 0.4033    
#Residuals      51    8.7501  0.1716         0.68586           
#Total          54   12.7579                 1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#######################
#Figure 7D
#######################

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

#define the ASVs significantly enriched in the rhizosphere samples
Native_Rhizo <- results(JH16_cds_test, contrast = c("Description",  "Modern.B", "Modern.N")) 

# extract  ASVs whose adjusted p.value in a given comparison is below 0.05. 
Native_Rhizo_FDR_005 <- Native_Rhizo[(rownames(Native_Rhizo)[which(Native_Rhizo$padj <0.05)]), ]

#Identify ASVs enriched in the rhizosphere (second term of the comparison, negative fold change)
Native_Rhizo_enriched <-  Native_Rhizo[(rownames(Native_Rhizo)[which(Native_Rhizo$log2FoldChange < 0)]), ]

#Intersection: this is the list of ASVs signficanlty enriched in the native rhizosphere; Morex conditioned
Native_Rhizo_enriched_FDR005 <- intersect(rownames(Native_Rhizo_FDR_005), rownames(Native_Rhizo_enriched))

#define the ASVs significantly enriched in the rhizosphere samples vs autoclaved ones
Autoclaved_Rhizo <- results(JH16_cds_test, contrast = c("Description",  "Modern.A", "Modern.N")) 

# extract  ASVs whose adjusted p.value in a given comparison is below 0.05. 
Autoclaved_Rhizo_FDR_005 <- Autoclaved_Rhizo[(rownames(Autoclaved_Rhizo)[which(Autoclaved_Rhizo$padj <0.05)]), ]

#Identify ASVs enriched in the rhizosphere (second term of the comparison, negative fold change)
Autoclaved_Rhizo_enriched <-  Autoclaved_Rhizo[(rownames(Autoclaved_Rhizo)[which(Autoclaved_Rhizo$log2FoldChange < 0)]), ]

#Intersection: this is the list of ASVs signficanlty enriched in the autoclaved rhizosphere; Morex conditioned
Autoclaved_Rhizo_enriched_FDR005 <- intersect(rownames(Autoclaved_Rhizo_FDR_005), rownames(Autoclaved_Rhizo_enriched))

#diagnostic ASVs for Morex conditioned; those are the ASVs simoultaneously enriched in Native rhizosphere versus bulk soil and autoclaved samples
Modern_Native_005 <- intersect(Native_Rhizo_enriched_FDR005, Autoclaved_Rhizo_enriched_FDR005)
length(Modern_Native_005)
#15

#B1K conditioned calculation
##################

#define the ASVs significantly enriched in the rhizosphere samples
Native_Rhizo_b <- results(JH16_cds_test, contrast = c("Description",  "B1K.B", "B1K.N")) 

# extract  ASVs whose adjusted p.value in a given comparison is below 0.05. 
Native_Rhizo_b_FDR_005 <- Native_Rhizo_b[(rownames(Native_Rhizo_b)[which(Native_Rhizo_b$padj <0.05)]), ]

#Identify ASVs enriched in the rhizosphere (second term of the comparison, negative fold change)
Native_Rhizo_b_enriched <-  Native_Rhizo_b[(rownames(Native_Rhizo_b)[which(Native_Rhizo_b$log2FoldChange < 0)]), ]

#Intersection: this is the list of ASVs signficanlty enriched in the native rhizosphere; Morex conditioned
Native_Rhizo_b_enriched_FDR005 <- intersect(rownames(Native_Rhizo_b_FDR_005), rownames(Native_Rhizo_b_enriched))

#define the ASVs significantly enriched in the rhizosphere samples vs autoclaved ones
Autoclaved_Rhizo_b <- results(JH16_cds_test, contrast = c("Description",  "B1K.A", "B1K.N")) 

# extract  ASVs whose adjusted p.value in a given comparison is below 0.05. 
Autoclaved_Rhizo_b_FDR_005 <- Autoclaved_Rhizo_b[(rownames(Autoclaved_Rhizo_b)[which(Autoclaved_Rhizo_b$padj <0.05)]), ]

#Identify ASVs enriched in the rhizosphere (second term of the comparison, negative fold change)
Autoclaved_Rhizo_b_enriched <-  Autoclaved_Rhizo_b[(rownames(Autoclaved_Rhizo_b)[which(Autoclaved_Rhizo_b$log2FoldChange < 0)]), ]

#Intersection: this is the list of ASVs signficanlty enriched in the autoclaved rhizosphere; B1K conditioned
Autoclaved_Rhizo_b_enriched_FDR005 <- intersect(rownames(Autoclaved_Rhizo_b_FDR_005), rownames(Autoclaved_Rhizo_b_enriched))

#diagnostic ASVs for B1K conditioned; those are the ASVs simoultaneously enriched in Native rhizosphere versus bulk soil and autoclaved samples
B1K_Native_005 <- intersect(Native_Rhizo_b_enriched_FDR005, Autoclaved_Rhizo_b_enriched_FDR005)
length(B1K_Native_005)
#35

#######################
#Plotting
######################

#identify the ASVs enriched in both conditioned soils (Morex and B1K): we aim at testing the hypothesis of convergent enrichment
enriched_ASVs <- intersect(Modern_Native_005, B1K_Native_005)
length(enriched_ASVs)
#10
write.table(enriched_ASVs, file="enriched_ASVs.txt", sep="\t")

#overlap with JH02 enrichment
Elite_enriched_N0 <- read.delim("NT_Figure_4_input_data_2.txt", row.names =1)
colnames(Elite_enriched_N0)
#overlap
overlap <- intersect(enriched_ASVs, rownames(Elite_enriched_N0))
overlap
#7
#Plotting taxa
#Define the enriched family
JH16_enriched <- prune_taxa(enriched_ASVs, JH16_sample)

#melting phyloseq object into a dataframe
df_family_JH16_all <- psmelt(JH16_enriched)
#plotting
plot_order_JH16_all <- ggplot(df_family_JH16_all, aes(Sample, Abundance, fill = fct_reorder(Class, Abundance))) + geom_col()+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#Order the samples
df_family_JH16_all$Sample <- factor(df_family_JH16_all$Sample,levels = c("AC13","AC14","AC62","AC64","AC65","AC67","AC68", "AC15", "AC16", "AC17", "AC35", "AC36", "AC37", "AC38", "AC39", "AC60", "AC63", "AC66", "AC01", "AC02", "AC03", "AC04", "AC05", "AC25", "AC26", "AC27", "AC28", "AC29", "AC50", "AC51", "AC52", "AC53", "AC54", "AC21", "AC22", "AC42", "AC43","AC46", "AC47", "AC49", "AC71", "AC72", "AC18", "AC19", "AC20", "AC23", "AC24", "AC40", "AC41", "AC44", "AC45", "AC48", "AC70", "AC73", "AC75", "AC76", "AC06", "AC07", "AC08", "AC09", "AC10", "AC30", "AC31", "AC32", "AC33", "AC34", "AC55", "AC56", "AC57", "AC58", "AC59"))
df_family_JH16_all_plot <- ggplot(df_family_JH16_all, aes(Sample, Abundance, fill = fct_reorder(Class, Abundance))) + geom_col() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
df_family_JH16_all_plot

#add cb friendly pallette
cbp1 <- c("#e8f5ff", "#8fff7d", "#829796", "#cec2be","#afbff6")
df_family_JH16_all_plot_cb <- df_family_JH16_all_plot + scale_fill_manual(breaks = c ("Polyangia", "Gammaproteobacteria", "Actinobacteria", "Alphaproteobacteria", "Bacteroidia"), values=cbp1)
df_family_JH16_all_plot_cb

#Save as EPS size 1000 X 500

#############################
#Supplementary Figure 8A and B
##############################

#Index calculations
JH16_no_plants_alpha_rare <-  estimate_richness(JH16_sample, measures = c("Observed", "Shannon"))

#generate a new dataframe for data visualisation

#Sample information
JH16_map <- as.data.frame(as.matrix(sample_data(JH16_sample)))
dim(JH16_map)
#description
design_Description <- as.data.frame(JH16_map[, 7])
rownames(design_Description) <- rownames(JH16_map)
colnames(design_Description) <- c("Description")


#treatment
design_treatment <- as.data.frame(JH16_map[, 6])
rownames(design_treatment) <- rownames(JH16_map)
colnames(design_treatment) <- c("Treatment")

#data frame Genotype_Description
design_GD <- cbind(design_Description, design_treatment)

#Observed ASVs
JH16_Observed <- as.data.frame(JH16_no_plants_alpha_rare[ ,1])
rownames(JH16_Observed) <- rownames(JH16_no_plants_alpha_rare)
colnames(JH16_Observed) <- c("Observed")

#Combine the dataset sample description and Observed ASVs
JH16_Observed_GD <- cbind(design_GD, JH16_Observed)

#Order the levels according to a defined order
JH16_Observed_GD$Description <- ordered(JH16_Observed_GD$Description, levels=c("Modern.B", "Modern.N", "Modern.A", "B1K.B", "B1K.N", "B1K.A"))
nt16 <- c("white", "#CC79A7","#0072B2","black", "#E69F00","#56B4E9")
#plotting
p5 <- ggplot(JH16_Observed_GD, aes(x=Description, y=Observed, fill=Description))
p5 + geom_jitter(size=4,shape=22, width = 0.2) + scale_fill_manual(values = nt16)


#stats
#check the distribution of the data
shapiro.test(JH16_Observed_GD$Observed)

#Observed
kruskal.test(Observed~ Description, data = JH16_Observed_GD)

#Kruskal-Wallis rank sum test

#data:  Observed by Description
#Kruskal-Wallis chi-squared = 55.617, df = 5, p-value = 9.745e-11

posthoc.kruskal.dunn.test (x=JH16_Observed_GD$Observed, JH16_Observed_GD$Description, p.adjust.method="BH")

#Pairwise comparisons using Dunn's-test for multiple	
                         #comparisons of independent samples 

#data:  JH16_Observed_GD$Observed and JH16_Observed_GD$Description 

         #Modern.B Modern.N Modern.A B1K.B   B1K.N  
#Modern.N 0.11708  -        -        -       -      
#Modern.A 5.8e-06  0.00067  -        -       -      
#B1K.B    0.75551  0.20206  5.8e-06  -       -      
#B1K.N    0.12228  0.86422  0.00016  0.21731 -      
#B1K.A    5.8e-06  0.00126  0.86422  7.9e-06 0.00034

#P value adjustment method: BH 

#Shannon
JH16_Shannon <- as.data.frame(JH16_no_plants_alpha_rare[ ,2])
rownames(JH16_Shannon) <- rownames(JH16_no_plants_alpha_rare)
colnames(JH16_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH16_Shannon_GD <- cbind(design_GD, JH16_Shannon)

#Order the levels according to a defined order
JH16_Shannon_GD$Description <- ordered(JH16_Shannon_GD$Description, levels=c("Modern.B", "Modern.N", "Modern.A", "B1K.B", "B1K.N", "B1K.A"))
#plotting
p5 <- ggplot(JH16_Shannon_GD, aes(x=Description, y=Shannon, fill=Description))
p5 + geom_jitter(size=4,shape=22, width = 0.2) + scale_fill_manual(values = nt16)

#stats
#check the distribution of the data
shapiro.test(JH16_Shannon_GD$Shannon)

#Shannon
kruskal.test(Shannon~ Description, data = JH16_Shannon_GD)

#Kruskal-Wallis rank sum test

#data:  Shannon by Description
#Kruskal-Wallis chi-squared = 60.725, df = 5, p-value = 8.606e-12

posthoc.kruskal.dunn.test (x=JH16_Shannon_GD$Shannon, JH16_Shannon_GD$Description, p.adjust.method="BH")

#Pairwise comparisons using Dunn's-test for multiple	
                         #comparisons of independent samples 

#data:  JH16_Shannon_GD$Shannon and JH16_Shannon_GD$Description 

        # Modern.B Modern.N Modern.A B1K.B   B1K.N  
#Modern.N 0.03830  -        -        -       -      
#Modern.A 3.0e-07  0.00044  -        -       -      
#B1K.B    0.60765  0.09787  4.8e-07  -       -      
#B1K.N    0.02367  0.89069  0.00034  0.06477 -      
#B1K.A    1.4e-06  0.00252  0.60765  3.7e-06 0.00226

#P value adjustment method: BH


#end

