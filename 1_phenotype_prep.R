## Script written by Jonathan Renk
## 7 June 2020
## Script adapted by Travis Hattery
## 3 FEB 2022

# Clear global environment
rm(list=ls(all=TRUE))
gc()

# Load required packages
#install.packages("car")
#install.packages("dfoptim")
#install.packages("lme4")
#install.packages("qqplotr")
library("car")
library("dfoptim")
library("lme4")
library("qqplotr")





############################################################
####### SPEARMAN CORRELATION ########
############################################################
setwd("C:/<WORKING DIRECTORY>")

input <- read.csv("widiv_3env_metab_spearman_comb.csv", header=T, stringsAsFactors=F)
traits <- colnames(input[3:110])

metab_est <- data.frame(matrix(ncol = 108, nrow = 108))
colnames(metab_est) <- traits
rownames(metab_est) <- traits
metab_p <- data.frame(matrix(ncol = 108, nrow = 108))
colnames(metab_p) <- traits
rownames(metab_p) <- traits

# Loop through all traits (i) by all traits (j) and store to grids
for(i in 1:108){
  
  for(j in 1:108){
    xname <- traits[i]
    yname <- traits[j]
    
    x <- input[,xname]
    y <- input[,yname]
    
    spearman <- cor.test(x,y,  method = "spearman")
    
    metab_est[i,j] = spearman$estimate
    metab_p[i,j] = spearman$p.value
    
    spearman
    
  }
}

write.csv(metab_est,"C:/<WORKING DIRECTORY>/spearman_metab_est.csv", row.names = TRUE)
write.csv(metab_p,"C:/<WORKING DIRECTORY>/spearman_metab_p.csv", row.names = TRUE)




#used for env-comparisons
for(i in 1:1){

#Uncomment for ISU16 vs UMN16  
  xname <- colnames(input)[i]
  yname <- colnames(data)[i+108]
#Uncomment for UMN16 vs ISU17  
#  xname <- colnames(data)[i+108]
#  yname <- colnames(data)[i+216]
#Uncomment for ISU16 vs ISU17  
#  xname <- colnames(data)[i]
#  yname <- colnames(data)[i+216]

#  ggscatter(data, x = xname, y = yname, 
#            add = "reg.line", conf.int = TRUE, 
#            cor.coef = TRUE, cor.method = "spearman",
#            xlab = xname, ylab = yname)
  
  x <- data[,xname]
  y <- data[,yname]

  spearman <- cor.test(x,y,  method = "spearman")
  spearman
  
  df[i-1,1] = xname
  df[i-1,2] = yname
  df[i-1,3] = spearman$estimate
  df[i-1,4] = spearman$p.value
  
}

write.csv(df,"C:/<WORKING DIRECTORY>/spearman.csv", row.names = FALSE)


############################################################
####### BLUPS FOR EVERY TRAIT ########
############################################################

# Load raw dataset
setwd("C:/<WORKING DIRECTORY>")
data <- read.csv("widiv_3env_raw_UMN16.csv", header=T, stringsAsFactors=F)

# Remove metadata columns
data <- data[,c(-1, -2, -3, -5, -7, -8, -9, -12, -13)] #Keep 4, 6, 10, 11, and num data
colnames(data)[1] <- "genotype"

# Setting variables as factors for genotype, env, rep, block
data[,1] <- as.factor(data[,1])
data[,2] <- as.factor(data[,2])
data[,3] <- as.factor(data[,3])
data[,4] <- as.factor(data[,4])

# Make sure identifying data is "Factor" and dataset is "num"
str(data)

# Save numerical data column names
traits <- colnames(data)[5:ncol(data)]

# Change Working Directory to where we want file output
setwd("C:/<WORKING DIRECTORY>")


for(trait in traits){
  pdf(paste0("plots_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,2))
  
  # Histogram of raw data
  hist(data[,trait],
       main=paste0("Histogram of raw\n", trait, " values"),
       xlab=trait,
       col="cadetblue")
  
  # Stripchart of raw data across environments
  stripchart(data[,trait] ~data$env,
             main=paste0("Stripchart of raw\n", trait, " values"),
             xlab="Environment",
             ylab=trait,
             vertical=TRUE, 
             method="jitter",
             bg="cadetblue",
             pch=21
  )
  
  dev.off()
  
  if(paste0("stats_", trait, ".txt") %in% list.files()){
    system(paste0("rm stats_", trait, ".txt"))}
  # Summary statistics of the trait
  summary <- summary(data[,trait], )
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Test for normality #SAMPLE SIZE MUST BE BELOW 5000
  normality <- shapiro.test(data[,trait])
  out <- capture.output(normality)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run an ANOVA (switched : in Env:Rep and Rep:Block for nesting) #ALTERED FROM ORIGINAL
  model <- lm(get(trait) ~ genotype + rep + rep/block, data=data)
  #model <- lm(get(trait) ~ genotype + env + env/rep + rep/block + genotype:env, data=data)
  

  anova <- anova(model)
  out <- capture.output(anova)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  # Run a random effects model #ALTERED FROM ORIGINAL
  model.1 <- lmer(get(trait) ~ (1|genotype) + (1|rep) + (1|rep/block), data = data, REML = TRUE)
  #model.1 <- lmer(get(trait) ~ (1|genotype) + (1|env) + (1|env/rep) + (1|rep/block) + (1|genotype:env), data = data, REML = TRUE)
  
  
    # Decreasing stopping tolerances
  strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
  if (all(model.1@optinfo$optimizer=="nloptwrap")) {
    model <- update(model.1, control=strict_tol)
  }
  summary(model, correlation=FALSE)
  random_effects <- ranef(model)
  # Write out BLUPs for Genotypes
  write.table(random_effects$genotype, paste0("blups_", trait, ".csv"), col.names=F, row.names=T, sep=",")
  # Summary of random effects
  summary <- summary(model, correlation=FALSE)
  out <- capture.output(summary)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("resids_", trait, ".csv"), col.names=F, row.names=F, sep=",")
  # Calculate heritability 
  model_variances <- as.data.frame(VarCorr(model))
  h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/5)+(model_variances$vcov[8]/10))
  out <- capture.output(h2)
  cat(out, file=paste0("stats_",trait,".txt"), sep="\n", append=TRUE)
  
  pdf(paste0("assumptions_", trait, ".pdf"), width = 15, height = 5)
  par(mfrow=c(1,3))
  
  # Model Fit with REML
  plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
  abline(h=0,col="red", lwd=1, lty=1)
  # histogram of residuals
  hist(residuals(model),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
  x=seq(-5e-15,9e-15,5e-15)
  curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
  # qq plot
  #qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
  
  dev.off()
  
}



############################################################
####### COMBINE BLUPS INTO MASTER FILES
############################################################

setwd("C:/<WORKING DIRECTORY>")

data1 <- read.csv("blups_c21hc.csv", header=F, stringsAsFactors=F)
rowlabels <- data1[,1]
blup_combined <- data.frame(matrix(ncol = 109, nrow = 448))
blup_combined[,1] <- rowlabels

collabels <- c("Taxa","c21hc","c22hc","c23hc","c24hc","c25hc","c26hc","c27hc","c28hc","c29hc",
               "c30hc","c31hc","c33hc","c35hc","c24_1_7hc","c24_1_9hc","c24_2_unkhc","c25_1_7hc",
               "c25_1_9hc","c27_1_7hc","c27_1_9hc","c27_2_unkhc","c29_1_7hc","c29_1_9hc",
               "c29_2_ahc","c29_2_bhc","c31_1_7hc","c31_1_9hc","c31_1_unkhc","c33_1_7hc","c33_1_9hc",
               "c33_1_unkhc","c35_1_7hc","c35_1_9hc","c16fatms","c17fatms","c18fatms","c20fatms","c22fatms",
               "c24fatms","c18_1_afa","c18_1_bfa","c18ohtms","c19ohtms","c25ohtms","c26ohtms","c27ohtms","c28ohtms",
               "sl","hc","hcsat","hcunsat","hceven","hcevensat","hcevenunsat","hcodd","hcoddsat",
               "hcoddunsat","fa","fasat","faunsat","faeven","faodd","oh","oheven","ohodd","7monoenes",
               "9monoenes","79alkenes","othermonoenes","hcmonoenes","dienes","hc_psl","hcsat_psl",
               "hcunsat_psl","hceven_psl","hcevensat_psl","hcevenunsat_psl","hcodd_psl","hcoddsat_psl",
               "hcoddunsat_psl","fa_psl","fasat_psl","faunsat_psl","faeven_psl","faodd_psl","oh_psl",
               "oheven_psl","ohodd_psl","7monoenes_psl","9monoenes_psl","hcsat_phc","hcunsat_phc",
               "hcevensat_phc","hcoddsat_phc","hcevenunsat_phc","hcoddunsat_phc","hcodd_phcsat",
               "hceven_phcsat","hcodd_phcunsat","hceven_phcunsat","hcmonoenes_phcunsat","dienes_phcunsat",
               "7monoenes_phcmonoenes","9monoenes_phcmonoenes","7monoenes_p79alkenes","9monoenes_p79alkenes",
               "21hc_p22fa","23hc_p24fa"
)

colnames(blup_combined) <- collabels


runlist <- c("c21hc","c22hc","c23hc","c24hc","c25hc","c26hc","c27hc","c28hc","c29hc",
             "c30hc","c31hc","c33hc","c35hc","c24_1_7hc","c24_1_9hc","c24_2_unkhc","c25_1_7hc",
             "c25_1_9hc","c27_1_7hc","c27_1_9hc","c27_2_unkhc","c29_1_7hc","c29_1_9hc",
             "c29_2_ahc","c29_2_bhc","c31_1_7hc","c31_1_9hc","c31_1_unkhc","c33_1_7hc","c33_1_9hc",
             "c33_1_unkhc","c35_1_7hc","c35_1_9hc","c16fatms","c17fatms","c18fatms","c20fatms","c22fatms",
             "c24fatms","c18_1_afa","c18_1_bfa","c18ohtms","c19ohtms","c25ohtms","c26ohtms","c27ohtms","c28ohtms",
             "sl","hc","hcsat","hcunsat","hceven","hcevensat","hcevenunsat","hcodd","hcoddsat",
             "hcoddunsat","fa","fasat","faunsat","faeven","faodd","oh","oheven","ohodd","X7monoenes",
             "X9monoenes","X79alkenes","othermonoenes","hcmonoenes","dienes","hc_psl","hcsat_psl",
             "hcunsat_psl","hceven_psl","hcevensat_psl","hcevenunsat_psl","hcodd_psl","hcoddsat_psl",
             "hcoddunsat_psl","fa_psl","fasat_psl","faunsat_psl","faeven_psl","faodd_psl","oh_psl",
             "oheven_psl","ohodd_psl","X7monoenes_psl","X9monoenes_psl","hcsat_phc","hcunsat_phc",
             "hcevensat_phc","hcoddsat_phc","hcevenunsat_phc","hcoddunsat_phc","hcodd_phcsat",
             "hceven_phcsat","hcodd_phcunsat","hceven_phcunsat","hcmonoenes_phcunsat","dienes_phcunsat",
             "X7monoenes_phcmonoenes","X9monoenes_phcmonoenes","X7monoenes_p79alkenes","X9monoenes_p79alkenes",
             "X21hc_p22fa","X23hc_p24fa"
)

j=2

for (i in runlist){

  filename1 <- ""
  filename <- ""
  filename1 <- paste("blups", i, sep = "_")
  filename <- paste(filename1, "csv", sep = ".")
  
  data <- read.csv(file=filename, header=F, stringsAsFactors=F)
  
  blup_combined[,j] <- data[,2]
  j <- j+1
  
}

setwd("C:/<WORKING DIRECTORY>")
write.csv(blup_combined,"gwas_blup_umn16.csv", row.names = FALSE)



############################################################
####### CHECK DISTRIBUTION OF TRAITS
############################################################

# Load dataset
setwd("C:/<WORKING DIRECTORY>")
data <- read.csv("widiv_3env_raw_full.csv", header=T, stringsAsFactors=F)

# Remove metadata columns
#data <- data[,c(-1)] #remove taxa info
data <- data[,c(14:121)]

# Save numerical data column names
#traits <- colnames(data)

traits <- c("sl", "hc", "hcsat", "hcunsat", "fa", "oh",
            "hc_psl", "hcsat_psl", "hcunsat_psl", "fa_psl", "oh_psl")

# Change Working Directory to where we want file output
#setwd("C:/Users/travi/Desktop/Genetics/TJH WiDiv GWAS Manuscript/traits_isu17")

#MAKE INDIVIDUAL HISTOGRAMS OF EACH TRAIT
#for(trait in traits){
#  pdf(paste0("blup_plots_", trait, ".pdf"), width = 5, height = 5)
#  par(mfrow=c(1,1))
  # Histogram of blups
#  hist(data[,trait],
#       main=paste0("Histogram of blups\n", trait, " values"),
#       xlab=trait,
#       col="cadetblue")
#  
#  dev.off()
#}

#MAKE ONE PDF OF ALL HISTOGRAMS
pdf(paste0("raw_someplots_full.pdf"), width = 55, height = 5) #25x25
par(mfrow=c(1,11)) #11,10

for(trait in traits){
  # Histogram of blups
  hist(data[,trait],
       main=paste0(trait),
       xlab=NULL,
       ylab=NULL,
       col="gray") #maroon, blue, gold, gray
}

dev.off()

