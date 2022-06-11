## Script written by Jonathan Renk
## Adapted by Travis Hattery
## 31 Jan 2022

## This script runs FarmCPU to calculate suggested p-value thresholds required for GWAS

############################################################################################################
#### Prepare Environment
############################################################################################################

#library(bigmemory)
library(biganalytics)
require(compiler)
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

############################################################################################################
#### import data
############################################################################################################

myY_full <- read.table("gwas_blup_full.txt", head = TRUE) #PHENOTYPE DATA
myGD <- read.big.matrix("GAPIT.Metab.Genotype.Numerical.txt", type="char", sep="\t", head = T) #GENOTYPE DATA
myGM <- read.table("GAPIT.Metab.Genotype.map.txt", head = TRUE) #MAP DATA

############################################################################################################
#### Generating the pvalue threshold used in FarmCPU
############################################################################################################

pvals <- FarmCPU.P.Threshold(
  Y=myY_full[,c(1,73)], #only two columns allowed, the first column is taxa name and the second is phenotype
  GD=myGD,
  GM=myGM,
  trait="hc_psl_full", #name of the trait, only used for the output file name
  theRep=100 #number of permutation times
)

#outputs p-values to working directory
