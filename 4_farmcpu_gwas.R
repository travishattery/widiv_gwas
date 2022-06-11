## Script written by Jonathan Renk
## Adapted by Travis Hattery
## 2 Feb 2022

## This script compiles all components to conduct GWAS

############################################################################################################
#### Prepare Environment
############################################################################################################

#Installation of required R packages (once per computer)
#install.packages("bigmemory")
#install.packages("biganalytics")

#Import library (each time to start R)
library(bigmemory)
library(biganalytics)
require(compiler)

## Clearing the global environment
rm(list=ls(all=TRUE))
gc()

#GAPIT and FarmCPU Web source code
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")


############################################################################################################
#### Generating the numerical format genotype files used in FarmCPU
############################################################################################################

#myY <- read.table("gwas_blup_traits_christine_v2.txt", head = TRUE)
#myG <- read.table("widiv_446g_christine_SNPs.hmp.txt", head=FALSE)

#myY <- read.table("gwas_blup_metabtraits.txt", head = TRUE)
#myG <- read.table("widiv_448metab_SNP.hmp.txt", head=FALSE)

# Converting HapMap format to numerical for FarmCPU
#myGAPIT <- GAPIT(G=myG, output.numerical=TRUE)


############################################################################################################
#### import data
############################################################################################################

setwd("<WORKING DIRECTORY>")
myY <- read.table("gwas_blup_full.txt", head = TRUE) #PHENOTYPE DATA NEEDED FOR GWAS (1/5)
myGD <- read.big.matrix("GAPIT.Metab.Genotype.Numerical.txt", type="char", sep="\t", head = T) #GENOTYPE DATA NEEDED FOR GWAS (2/5)
myGM <- read.table("GAPIT.Metab.Genotype.map.txt", head = TRUE) #MAP DATA NEEDED FOR GWAS (3/5)

############################################################################################################
#### Generating the PCA file used in FarmCPU #####
############################################################################################################

#gc()

#myGAPIT <- GAPIT( 
#  Y=myY[,c(1,49)],  #column 49 is Total SL "sl"
#  GD=myGD,
#  GM=myGM,
#  PCA.total=5, #first 5 PC's will be sufficient
#  method.bin="optimum",
#  model="FarmCPU"
#)

# Use output from previous block of script as input here
# load in the PC's to use as covariates in the model
# covariates are based on map data, so will be identical for all metabolite datasets
setwd("C:/<WORKING DIRECTORY>")
PCA <- read.csv("covariates.csv", stringsAsFactors=FALSE)
PCA <- PCA[,c(-1)] #COVARIATE DATA NEEDED FOR GWAS (4/5)

############################################################################################################
#### Generating the pvalue threshold used in FarmCPU
############################################################################################################

#pvals <- FarmCPU.P.Threshold(
#  Y=myY[,c(1,49)], #only two columns allowed, the first column is taxa name and the second is phenotype
#  GD=myGD,
#  GM=myGM,
#  trait="sl", #name of the trait, only used for the output file name
#  theRep=30 #number of permutation times
#  #theRep was set to 100, I changed back to default 30
#)

#hc_psl_
#hcsat_psl_
#hcunsat_psl_
#fa_psl_
#oh_psl_

# Use output from previous block of script as input here
setwd("C:/<WORKING DIRECTORY>")
pvals <- read.table("FarmCPU.p.threshold.optimize.oh_psl_full.txt", head = FALSE) # change this for each trait
#pvals <- read.table("FarmCPU.p.threshold.optimize.hc_isu16.txt", head = FALSE) # change this for each trait
#pvals <- read.table("FarmCPU.p.threshold.optimize.hc_umn16.txt", head = FALSE) # change this for each trait
#pvals <- read.table("FarmCPU.p.threshold.optimize.hc_isu17.txt", head = FALSE) # change this for each trait
pval <- quantile(pvals$V1, 0.01) #P-VALUE CUTOFFS NEEDED FOR GWAS (5/5)

############################################################################################################
#### Run FarmCPU for GWAS results
############################################################################################################


#sl=49, hc=50, hcsat=51, hcunsat=52, fa=59, oh=64
#hc_psl=73, hcsat_psl_74, hcunsat_psl=75, fa_psl=82, oh_psl=87

setwd("C:/<WORKING DIRECTORY>")
myFarmCPU <- FarmCPU(
  Y=myY[,c(1,87)], #phenotype CHANGE FOR EACH TRAIT
  GD=myGD, #genotype matrix
  GM=myGM, #genotypic map
  CV=PCA, #covariates
  threshold.output=0.01, #P value smaller than threshold.output will be output in GWAS table
  p.threshold=pval, 
  #p.threshold=0.05/nrow(myGM)
  MAF.calculate=TRUE, #Calculate minor allele frequency (MAF) or not, if set to TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  method.bin="optimum",
  maf.threshold=0.05, #When MAF.calculate=TRUE, the SNPs with a lower MAF (<maf.threshold) will be deleted
  maxLoop=50 #Maximum number of iterations allowed
)

# FarmCPU will output plots and excel sheets of results to the Working Directory



