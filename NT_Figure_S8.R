####################################################
## R code for Figure S8
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
library("ggfortify")
library("ggthemes")
library("colorspace")
library("grDevices")
library("vegan")
library ("ape")
library("PMCMRplus")

#Data pre-procesing such a concentration calculations and change to log scale were performed in Excel

#Import data file 
qJH16 <-(read.delim("JH16_SFigure.txt", sep = "\t"))

###############plotting Bacterial DNA log concentrations

#order the factor
qJH16$Sample <-  ordered(qJH16$Sample, levels=c( "Rhizo_modern", "Heat_modern","Rhizo_desert", "Heat_desert"))

Scale <- scale_y_continuous(limit = c(0, 23))
p <-ggplot(qJH16, aes(x=Sample, y=Log_bacterial, fill=Sample)) + geom_boxplot() + ggtitle("")+ylab("Log Bacterial DNA concentration(fg/ul)")
p + geom_jitter( size=5,shape=21, position=position_jitter(0.02))+ scale_fill_manual(values = c("magenta", "blue","goldenrod","goldenrod4"))+Scale+ theme(legend.position = "none",axis.text=element_text(size=16),axis.title.x = element_blank(),axis.title.y = element_text(size=16))



##################################################################
#Stats Bacterial DNA concentrations

kruskal.test(Log_bacterial~Sample , data=qJH16)
posthoc.kruskal.dunn.test (x= qJH16$Log_bacterial, g= qJH16$Sample, p.adjust.method="BH")

###################################################################
###################################################################
##Ploting Fungal DNA log concentrations
#order the factor
qJH16$Sample <-  ordered(qJH16$Sample, levels=c( "Rhizo_modern", "Heat_modern","Rhizo_desert", "Heat_desert"))

Scale <- scale_y_continuous(limit = c(0, 23), )
p <-ggplot(qJH16, aes(x=Sample, y=Log_Fungal, fill=Sample)) + geom_boxplot() + ggtitle("")+ylab("Log Fungal DNA concentration(fg/ul)")
p + geom_jitter( size=5,shape=21, position=position_jitter(0.02))+ scale_fill_manual(values = c("magenta", "blue","goldenrod","goldenrod4"))+Scale+ theme(legend.position = "none",axis.text=element_text(size=16),axis.title.x = element_blank(),axis.title.y = element_text(size=16))


##################################################################
#Stats Fungal DNA concentrations

kruskal.test(Log_Fungal~Sample , data=qJH16)

###################################################################




