#############################################################
#
# Ref to the ARTICLE
# 
#  Code to compute calculations presented in Figure S6
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

#loand the package(s) needed for the analysis
library(vegan)

#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()

#set working directory: Davide
setwd("/cluster/db/R_shared/NT2020_R_data/")

#######################################################################
#analyses on soil environmental variables
#######################################################################

# upload the data file:mineral nitrogen
env_data <- read.delim("NT_Figure_S6_input_data.txt", sep = "\t", header=TRUE, row.names=1)
dim(env_data)

#generate a non-metric multidimensional scaling
env_data.mds <- metaMDS(env_data, autotransform = FALSE)

#plotting
dev.off()
#save this image with sites name
plot(env_data.mds, type= "t", display=c ("sites"))

#fitting environmental parameters
env_data.em <- envfit(env_data.mds, env_data, permutations = 5000)
env_data.em

#save the resulting data as Additional_file5_ws3.txt 

#Figure
plot(env_data.em, p.max = 0.002)

#End 

