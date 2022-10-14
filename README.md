# Update October 2022

This repository contains scripts and material to reproduce calculations and figures of the manuscript Alegria Terrazas, Robertson-Albertyn et al., https://www.biorxiv.org/content/10.1101/605204v3

For each composite/sub-figure and/or table we produce a separate script: e.g., NT_Figure_1_S3_1020.R refers to the script needed to reproduce figure and calculatios of both Figure 3 and Figure S3. The folder inputfiles contain files required for the reproduction of the codes.

# NT Metagenomic assembly, Kraken and dada2 code

Scripts for metagenomic analysis included in the NT_scripts directory.

The code used to carry out the metagenomic assembly and functional classification, Kraken2 classification and dada2 16S classification are contained in separate subdirectories of this directory. The contents of each of subdirectory are as follows:

*  `bin` - directory containing scripts used for analysis
*  `*.ipynb` - Juypter notebook workflow for each analysis
*  `*.html` - HTML export of Jupyter notebook for viewing without a Jupyter environment
*  `*.yaml` - anaconda yaml configuration file for replicating analysis runtime environment 

Computationally intensive workloads are carried out with the scripts in the appropriate 'bin' directory on an HPC cluster. Data reprocessing and analysis is carried out directly within the Jupyter notebooks.

## Scripts for HPC Jobs

Note that scripts are designed to work with the University of Dundee HPC cluster, running Univa Grid Engine, and modifications may be required to run in different environments. The GridEnginge directives (lines beginning `#$`) indicate the resource requirements of the script. The following are commonly used in the scripts which may need altering:

*  `$TMPDIR`:  path to local temporary directory on compute node
*  `$TMPDIR1`: path to local ramdisk on compute node (used as high-performance temporary space)
*  Job classes: 
   +  Scripts with the '-jc long' directive require >12 hours and <72 hours to run
   +  Scripts with the '-jc inf' directive require >72 hours to run
*  Memory:  Default memory allocation is 8Gb. Scripts requiring greater amounts of memory will contain the directive `#$ -mods l_hard mfree 64G` for example, requesting 64Gb RAM
*  CPU cores: The number of cores required is requested with the `#$ -pe smp` directive
*  Ramdisk: Jobs using a ramdisk will contain a `#$ -adds l_hard ramdisk 4G` directive, in the case requiring 4Gb for the ramdisk

All scripts provide usage information when run without parameters

## Conda environments

Conda .yaml environments allow the same software environment used to run the analysis to be recreated to improve reproducibility. In order to make use of these, follow the instructions in the bioconda [User guide](https://bioconda.github.io/user/install.html#install-conda).

An environment can then be created (for example, for dada2) as follows:

`conda env create -f nt_dada2.yaml`

Note for complex environments this can take some time. Once complete the environment can be activated with 

`conda activate nt_dada2`

The same versions of the software used in this analysis will then be available on your path.
