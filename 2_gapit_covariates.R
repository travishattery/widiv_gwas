## Script written by Jonathan Renk
## Adapted by Travis Hattery
## 31 Jan 2022

## This script runs GAPIT to calculate 5 PCA values required for GWAS

############################################################################################################
#### Prepare Environment
############################################################################################################

#library(bigmemory)
library(biganalytics)
require(compiler)
source("http://zzlab.net/GAPIT/gapit_functions.txt")

############################################################################################################
#### import data
############################################################################################################

myY <- read.table("gwas_blup_2016.txt", head = TRUE) #PHENOTYPE DATA
myGD <- read.big.matrix("GAPIT.Metab.Genotype.Numerical.txt", type="char", sep="\t", head = T) #GENOTYPE DATA
myGM <- read.table("GAPIT.Metab.Genotype.map.txt", head = TRUE) #MAP DATA

############################################################################################################
#### Generating the PCA file used in FarmCPU #####
############################################################################################################

myGAPIT <- GAPIT( 
  Y=myY[,c(1,49)],  #column 49 is Total SL "sl"
  GD=myGD,
  GM=myGM,
  PCA.total=5, #first 5 PC's will be sufficient
  method.bin="optimum",
  model="FarmCPU"
)

#PCA values output to working directory
