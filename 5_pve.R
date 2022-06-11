## Script written by Travis Hattery
## May 2022

## RUN IN COMMAND LINE AFTER DOWNLOADING PLINK AND GCTA64

############################################################################################################
#### PLINK
############################################################################################################

# start with hapmap SNP file in TASSEL
# Filter to remove minor SNP states (binarize)
# save as PLINK to get .plk.map and .plk.ped files
# run in CMD

plink --noweb --file trait.plk --make-bed --out trait

# outputs .bed .bim and .fam

############################################################################################################
#### GCTA
############################################################################################################

# take outputs of PLINK (.bed .bim and .fam)
# run in CMD

gcta64 --bfile trait --make-grm --out trait

# outputs .grm.bin and others

# add phenotypic data to working directory
# column 1 is al "-9"
# column 2 is genotype name
# column 3 is BLUP phenotypic data
# make sure file is saved as tab-delim text, then rename to .phen
# run in CMD

gcta64 --grm trait --pheno widiv448_trait.phen --reml --out trait

# outputs trait.hsq which can be opened in Excel
# extract V(G)/Vp variance and SE

