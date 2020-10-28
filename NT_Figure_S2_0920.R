#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations presented on Figure S2
#  Revision 09/20 d.bulgarelli@dundee.ac.uk 
#
#  
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

#load required packages (not all of these are needed for this so far but might be handy downstream)
library("ggplot2")
library ("PMCMR")

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#DAVIDE_DUNDEE
setwd("/cluster/db/R_shared/NT2020_R_data/")

#############################################################
#import the dry weight and Nitrogen data
#############################################################

dat_info <- read.delim("NT_Figure_S2_input_data.txt", row.names = 1)

#plant data: remove bulk soil samples
#####################################
dat_info_rhizo <- dat_info[(rownames(dat_info)[which(dat_info$Microhabitat != "Bulk")]), ]

#check the distribution of the data
shapiro.test(dat_info_rhizo$Dryweight)
shapiro.test(dat_info_rhizo$N_content)

#check the effect of the treatmenton both parameters
#Dry weight
kruskal.test(Dryweight ~ Treatment, data = dat_info_rhizo)
posthoc.kruskal.dunn.test (x=dat_info_rhizo$Dryweight, g=dat_info_rhizo$Treatment, p.adjust.method="BH")
#Nitrogen content
kruskal.test(N_content ~ Treatment, data = dat_info_rhizo)
posthoc.kruskal.dunn.test (x=dat_info_rhizo$N_content, g=dat_info_rhizo$Treatment, p.adjust.method="BH")

#order levels
dat_info_rhizo$Treatment <- ordered(dat_info_rhizo$Treatment, levels=c("N0","N25","N100"))

#Plot
p <- ggplot(dat_info_rhizo, aes(x=Treatment, y=Dryweight, fill=Treatment)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("gray","gray45","gray25"))

#Plot
p <- ggplot(dat_info_rhizo, aes(x=Treatment, y=N_content, fill=Treatment)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("gray","gray45","gray25"))


#residual nitrogen in rhizosphere
#####################################

#check the distribution of the data
shapiro.test(dat_info_rhizo$NH4)
shapiro.test(dat_info_rhizo$NO3)

#check the effect of the treatmenton both parameters
#NH4
kruskal.test(NH4 ~ Treatment, data = dat_info_rhizo)
posthoc.kruskal.dunn.test (x=dat_info_rhizo$NH4, g=dat_info_rhizo$Treatment, p.adjust.method="BH")

#NO3
kruskal.test(NO3 ~ Treatment, data = dat_info_rhizo)
posthoc.kruskal.dunn.test (x=dat_info_rhizo$NO3, g=dat_info_rhizo$Treatment, p.adjust.method="BH")

#order levels
dat_info_rhizo$Treatment <- ordered(dat_info_rhizo$Treatment, levels=c("N0","N25","N100"))

#Plot
p <- ggplot(dat_info_rhizo, aes(x=Treatment, y=NH4, fill=Treatment)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("gray","gray45","gray25"))

#Plot
p <- ggplot(dat_info_rhizo, aes(x=Treatment, y=NO3, fill=Treatment)) + geom_boxplot()
p + geom_jitter(size=5,shape=21, position=position_jitter(0.2))+ scale_fill_manual(values = c("gray","gray45","gray25"))

#End

