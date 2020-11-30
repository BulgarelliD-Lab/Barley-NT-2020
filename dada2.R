#!/bin/env Rscript

#$ -j y
#$ -o job_logs/$JOB_NAME.$JOB_ID
#$ -cwd 
#$ -pe smp 24

# Script to generate ASVs from 16S data and store as phyloseq objects using dada2

# Fastq files are expected to have the sample name preceding the first '_', and have
# pairing indicated by _1 and _2 in the filenames

library("getopt")

optspec=matrix(c(
	'reads',    'r', 1, 'character', 'Path to directory of fastq files',
	'metadata', 'm', 1, 'character', 'Path to metadata file',
	'taxonomy', 't', 1, 'character', 'Path to taxonomy database',
	'name',     'n', 1, 'character', 'Name of job',
	'help',     'h', 0, 'logical',   'Display help'
),byrow=TRUE,ncol=5)

opt=getopt(optspec)

if (!is.null(opt$help)) {
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$reads)) {
	cat("Error: no reads argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$metadata)) {
	cat("Error: no metadata argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$name)) {
	cat("Error: no name argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

if (is.null(opt$taxonomy)) {
	cat("Error: no taxonomy argument provided\n")
	cat(getopt(optspec,usage=TRUE))
	q(status=1)
}

library("ggplot2")
library("gridExtra")
library("dada2")
library("phyloseq")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

set.seed(100)

read_path=opt$reads
ref_fasta=opt$taxonomy
metadata_file=opt$metadata
filt_path="filtered"

cat("dada2 analysis\n==============\n\n")
cat(paste("Reads: ",read_path,"\n"))
cat(paste("Metadata:",metadata_file,"\n"))
cat(paste("Taxonomy:",ref_fasta,"\n\n"))

fns<-sort(list.files(read_path, full.names=TRUE))
fnFs<-fns[grepl("_1",fns)]
fnRs<-fns[grepl("_2",fns)]

samdf<-read.csv(metadata_file,header = TRUE,sep="\t")
rownames(samdf)<-samdf$SampleID

filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

if (!file_test("-d",filt_path)) {
	dir.create(filt_path)
	for (i in seq_along(fnFs)) {
		cat(paste("forward: ", fnFs[[i]], "; reverse: ", fnRs[[i]]))
		# First 10bp trimmed, but 3' end of read consistently >Q20, so leave as is
  		fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]), c(filtFs[[i]], filtRs[[i]]),
                      trimLeft=10, maxN=0, maxEE=2, truncQ=2, compress=TRUE)
	}
}

message('Derepping...')
sam.names<-sapply(strsplit(basename(filtFs),"_"),'[',1)
derepFs<-derepFastq(filtFs)
derepRs<-derepFastq(filtRs)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

message('Building error model')
ddF <- dada(derepFs, err=NULL, selfConsist=TRUE, MAX_CONSIST=20)
ddR <- dada(derepRs, err=NULL, selfConsist=TRUE, MAX_CONSIST=20)

message("Check convergence on error model...")
dada2:::checkConvergence(ddF[[1]])
dada2:::checkConvergence(ddR[[1]])

plotErrors(ddF)
plotErrors(ddR)

message('applying error model...')
dadaFs<-dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE, multithread = TRUE)
dadaRs<-dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE, multithread = TRUE)

message('merging reads...')
mergers<-mergePairs(dadaFs,derepFs,dadaRs,derepRs)
seqtab.all<-makeSequenceTable(mergers)
seqtab<-removeBimeraDenovo(seqtab.all)

message('Assigning taxonomy')
# minboot defined bootstrap threshold - recommended for 50 if reads <250 bp, otherwise 80
taxtab<-assignTaxonomy(seqtab,refFasta = ref_fasta,minBoot=50)
colnames(taxtab)<-c("Kingdom","Phylum","Class","Order","Family","Genus")

# to aid debugging phyloseq object creation problems...
saveRDS(taxtab,file='taxtab.rds')
saveRDS(samdf,file='samdf.rds')
saveRDS(seqtab,file='seqtab.rds')

message('Creating phyloseq object')
ps<-phyloseq(
	tax_table(taxtab),
	sample_data(samdf),
	otu_table(t(seqtab),taxa_are_rows=TRUE) 
)

message('saving Rds file...')
saveRDS(ps,file=paste0(opt$name,'_dada2.rds'))

