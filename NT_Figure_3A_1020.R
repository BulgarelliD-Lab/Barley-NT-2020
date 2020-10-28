#############################################################
#r version 3.5.3
#############################################################
#Script required to generate figure 3A 
#############################################################
# Libraries and functions required
#############################################################
#clean up
#These initial commands are required to clean-up the memory and start a new session
rm(list=ls())
dev.off()

# Load required packages -----------------------------------------------------
library(tidyr)
library(ggfortify)
library(dplyr)
library(UpSetR)
library(factoextra)
library(RColorBrewer)



#retrieve R and package versions and compare to the uploaded file in gtHub for the reproducibility of the code
sessionInfo()


# Load and prepare data ------------------------------------------------------

seq_sum_up <- read.table("go_biol_proc_up_SEQ.txt", header = TRUE)
seq_sum_slim_up <- read.table("slim_biol_proc_up_SEQ.txt", header = TRUE)
seq_sum_down <- read.table("go_biol_proc_down_SEQ.txt", header = TRUE)
seq_sum_slim_down <- read.table("slim_biol_proc_down_SEQ.txt", header = TRUE)




# get a count of each KO in each taxa
seq_ko_table_up <- seq_sum_up %>% count(SEQ,KO)
seq_ko_slim_table_up <- seq_sum_slim_up %>% count(SEQ,KO)
seq_ko_table_down <- seq_sum_down %>% count(SEQ,KO)
seq_ko_slim_table_down <- seq_sum_slim_down %>% count(SEQ,KO)


# Convert from long to wide format
seq_wide_up <- spread(seq_ko_table_up, KO, n)
write.csv(seq_wide_up, "seq_wide_nt_upr.csv")

seq_wide_slim_up <- spread(seq_ko_slim_table_up, KO, n)
write.csv(seq_wide_slim_up, "seq_wide_slim_upr.csv")

seq_wide_down <- spread(seq_ko_table_down, KO, n)
write.csv(seq_wide_down, "seq_wide_nt_downr.csv")

seq_wide_slim_down <- spread(seq_ko_slim_table_down, KO, n)
write.csv(seq_wide_slim_down, "seq_wide_slim_downr.csv")

#amended files to include sample names instead of Numerical ID
#seq_wide <- read.csv("seq_wide_nt2.csv", header = TRUE)
seq_wide_up <- read.csv("seq_wide_nt_upr2.csv", header = TRUE)
seq_wide_down <- read.csv("seq_wide_nt_downr2.csv", header = TRUE)
seq_wide_slim_up <- read.csv("seq_wide_slim_upr2.csv", header = TRUE)
seq_wide_slim_down <- read.csv("seq_wide_slim_downr2.csv", header = TRUE)

row.names(seq_wide_up) <- seq_wide_up$SEQ
row.names(seq_wide_slim_up) <- seq_wide_slim_up$SEQ
row.names(seq_wide_down) <- seq_wide_down$SEQ
row.names(seq_wide_slim_down) <- seq_wide_slim_down$SEQ
# Set missing values to zero and any count greater than zero to one to give 
#  a boolean presence/abscence matrix
seq_wide_up[is.na(seq_wide_up)] <- 0
for (x in 2:ncol(seq_wide_up)) {
  seq_wide_up[,x] <- ifelse(seq_wide_up[,x]>0,1,0)
}
write.csv(seq_wide_up, "binary.seq.vs.KO.nt_upr.csv")


seq_wide_slim_up[is.na(seq_wide_slim_up)] <- 0
for (x in 2:ncol(seq_wide_slim_up)) {
  seq_wide_slim_up[,x] <- ifelse(seq_wide_slim_up[,x]>0,1,0)
}
write.csv(seq_wide_slim_up, "binary.seq.vs.KO.slim_upr.csv")

seq_wide_down[is.na(seq_wide_down)] <- 0
for (x in 2:ncol(seq_wide_down)) {
  seq_wide_down[,x] <- ifelse(seq_wide_down[,x]>0,1,0)
}
write.csv(seq_wide_down, "binary.seq.vs.KO.nt_downr.csv")

seq_wide_slim_down[is.na(seq_wide_slim_down)] <- 0
for (x in 2:ncol(seq_wide_slim_down)) {
  seq_wide_slim_down[,x] <- ifelse(seq_wide_slim_down[,x]>0,1,0)
}
write.csv(seq_wide_slim_down, "binary.seq.vs.KO.slim_downr.csv")





# Genrate an intersect plot --------------------------------------------------
#all up
row.names(seq_wide_up)
colnames(seq_wide_up)
bin.seq_wide_up <- seq_wide_up[,-1]

row.names(bin.seq_wide_up) <- row.names(seq_wide_up)
row.names(bin.seq_wide_up)

t.seq.wide_up<- t(bin.seq_wide_up)
t.seq.wide_up <-as.data.frame(t.seq.wide_up)
row.names(t.seq.wide_up)
colnames(t.seq.wide_up)
setnames_up <- colnames(t.seq.wide_up)
setnames_up

upset(t.seq.wide_up, sets=setnames_up, order.by = "freq", nintersects=50)

#slim up

row.names(seq_wide_slim_up)
colnames(seq_wide_slim_up)
bin.seq_wide_slim_up <- seq_wide_slim_up[,-1]

row.names(bin.seq_wide_slim_up) <- row.names(seq_wide_slim_up)
row.names(bin.seq_wide_slim_up)

t.seq.wide_slim_up<- t(bin.seq_wide_slim_up)
t.seq.wide_slim_up <-as.data.frame(t.seq.wide_slim_up)
row.names(t.seq.wide_slim_up)
colnames(t.seq.wide_slim_up)
setnames_slim_up <- colnames(t.seq.wide_slim_up)
setnames_slim_up

upset(t.seq.wide_slim_up, sets=setnames_slim_up, order.by = "freq", nintersects=50)

#all down
row.names(seq_wide_down)
colnames(seq_wide_down)
bin.seq_wide_down <- seq_wide_down[,-1]

row.names(bin.seq_wide_down) <- row.names(seq_wide_down)
row.names(bin.seq_wide_down)

t.seq.wide_down<- t(bin.seq_wide_down)
t.seq.wide_down <-as.data.frame(t.seq.wide_down)
row.names(t.seq.wide_down)
colnames(t.seq.wide_down)
setnames_down <- colnames(t.seq.wide_down)
setnames_down

upset(t.seq.wide_down, sets=setnames_down, order.by = "freq", nintersects=50)

#slim down

row.names(seq_wide_slim_down)
colnames(seq_wide_slim_down)
bin.seq_wide_slim_down <- seq_wide_slim_down[,-1]

row.names(bin.seq_wide_slim_down) <- row.names(seq_wide_slim_down)
row.names(bin.seq_wide_slim_down)

t.seq.wide_slim_down<- t(bin.seq_wide_slim_down)
t.seq.wide_slim_down <-as.data.frame(t.seq.wide_slim_down)
row.names(t.seq.wide_slim_down)
colnames(t.seq.wide_slim_down)
setnames_slim_down <- colnames(t.seq.wide_slim_down)
setnames_slim_down

upset(t.seq.wide_slim_down, sets=setnames_slim_down, set_size.show = TRUE, order.by = "freq", nintersects=50)


