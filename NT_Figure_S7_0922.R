####################################################
## R code for NT SFing7
######################################################
###Quantification of bacterial and fungal DNA in soil feed-back experiment
#######################################################
#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()
#############################################################

library ("ggplot2")
library("dplyr")
library("MASS")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMRplus")



###To begin with set the working directory/ folder where the data sets are stored
setwd()
getwd()


#####################################################
#Data pre-processing such a concentration and copy number calculations and log scale were performed in Excel
###################################################

##The CT values of the qPCR were calculated as follow: 

#1) Linear extrapolation of the samples CT values in the standard curve (CT standards and Ln of concentrations (fg)), with the formula: y=ax+b, where y=CT value of the sample
#Y(x) = Y(1)+ (x- x(1)/x(2)-x(1)) * (Y(2) - Y(1)) in Excel

#2) Transformed in copy numbers:
## Bacteria 16S= Genome copy # = DNA (g) / (g_to_bp const. x genome size) With DNA (g) = 20 ng g to bp const = 1.096 x 10^-21 g genome size = 4.6 x 10^6 bp 20 x 10^-9 g / (1.096 x 10^-21 g/bp x 4.6 x 10^6 bp) = 4.0 x 10^6 copies in 20 ng
##Fungi ITS=  Genome copy # = DNA (g) / (g_to_bp const. x genome size) With DNA (g) = 20 ng g to bp const = ITS copy#=20ng* 1.079 x 10^-12 genome size#=12.1*10^6  20 ng / ( 1.079 x 10^-12 ng/bp x 12.1 x 10^6 bp) = 1531874 copies in 20 ng

#3) Log transformed


###########################################################
##Data 16S copy numbers
###########################################################

#Import data file 
qJH16c <-(read.delim("JH16_Input_R_16s_ITS_copies.txt", sep = "\t"))

########################################
#Plotting Bacterial log 16S copy numbers

#order the factor
qJH16c$Sample <-  ordered(qJH16c$Sample, levels=c( "Rhizo_modern", "Heat_modern","Rhizo_desert", "Heat_desert"))

Scale <- scale_y_continuous(limit = c(0, 25))
p <-ggplot(qJH16c, aes(x=Sample, y=Log_16S_copy, fill=Sample)) + geom_boxplot() + ggtitle("")+ylab("Log Bacterial 16S rRNA copies/ul")
p + geom_jitter( size=5,shape=21, position=position_jitter(0.02))+ scale_fill_manual(values = c("#CC79A7","#0072B2","#E69F00", "#56B4E9"))+Scale+ theme(legend.position = "none",axis.text=element_text(size=16),axis.title.x = element_blank(),axis.title.y = element_text(size=16))



########################################
#Stats Bacterial 16S copy numbers

qJH16$Log_bacterial <- as.factor(qJH16c$Log_bacterial)
qJH16$Sample <- as.factor(qJH16c$Sample)
kruskal.test(Log_bacterial~Sample , data=qJH16c)
kwAllPairsDunnTest (x= qJH16$Log_bacterial, g= qJH16$Sample, p.adjust.method="BH")


###################################################################
###################################################################
###########################################################
##Data ITS copy numbers
###########################################################

##Plotting Fungal ITS copy numbers

#order the factor
qJH16c$Sample <-  ordered(qJH16c$Sample, levels=c( "Rhizo_modern", "Heat_modern","Rhizo_desert", "Heat_desert"))

Scale <- scale_y_continuous(limit = c(0, 23), )
p <-ggplot(qJH16c, aes(x=Sample, y=Log_ITS_copy, fill=Sample)) + geom_boxplot() + ggtitle("")+ylab("Log ITS fungal copies/ul")
p + geom_jitter( size=5,shape=21, position=position_jitter(0.02))+ scale_fill_manual(values = c("#CC79A7","#0072B2","#E69F00", "#56B4E9"))+Scale+ theme(legend.position = "none",axis.text=element_text(size=16),axis.title.x = element_blank(),axis.title.y = element_text(size=16))


##################################################################
#Stats Fungal ITS copy numbers

kruskal.test(Log_Fungal~Sample , data=qJH16c)

###################################################################
##End



