# change directory 
# setwd("/projectnb/bf528/users/group_5/project_1/project-1-group-5/program/files")

# install.packages("BiocManager")
#BiocManager::install("c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)

#2
# read file in
data <- ReadAffy(celfile.path = "/projectnb/bf528/users/group_5/project_1/project-1-group-5/program/files")

#normalize using:
norm_data <- rma(data)
# Background correcting
# Normalizing
# Calculating Expression
# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 54675 features, 134 samples 
# element names: exprs 
# protocolData
# sampleNames: GSM971958_JS_04_U133_2.CEL.gz
# GSM971959_JS_05_U133_2.CEL.gz ...
# GSM972521_VB_335T_U133_2.CEL.gz (134 total)
# varLabels: ScanDate
# varMetadata: labelDescription
# phenoData
# sampleNames: GSM971958_JS_04_U133_2.CEL.gz
# GSM971959_JS_05_U133_2.CEL.gz ...
# GSM972521_VB_335T_U133_2.CEL.gz (134 total)
# varLabels: sample
# varMetadata: labelDescription
# featureData: none
# experimentData: use 'experimentData(object)'
# Annotation: hgu133plus2

#4
# third step: compute Relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) scores
# Pset is data produced 
Pset <- fitPLM(data, normalize = TRUE, background=TRUE)

#RLE plot
RLE(Pset)

# RLE summary stats including median and IQR - only need median 
RLE_stats <- RLE(Pset,type="stats")
hist(RLE_stats["median",], xlab = "Median RLE", main = "Distribution of RLE Medians") # sample median distribution histogram 
 
#NUSE plot
NUSE(Pset)

# NUSE summary stats including median and IQR
NUSE_stats <- NUSE(Pset,type="stats")
hist(NUSE_stats["median",], xlab = "Median NUSE", main = "Distribution of NUSE Medians") # median distribution histogram 

# 5
# transform into a matrix - expression data where probesets are rows and samples are columns
matrix_data <- exprs(norm_data)

# read in the meta data 
meta_data <- read.csv("/projectnb/bf528/users/group_5/project_1/project-1-group-5/program/files/proj_metadata.csv")

# only interested in normalizationcombatbatch and normalizationcombatmod
batch = meta_data$normalizationcombatbatch
mod = model.matrix(~normalizationcombatmod, data=meta_data)  

# running ComBat to correct for batch effects while preserving features of interest like tumor and MMR status
bf_corrected <- ComBat(dat = matrix_data, batch = batch, mod = mod)

# expression data where probesets are rows and samples are columns
write.csv(bf_corrected, file = "expressed_data.csv")

# Principal Component Analysis (PCA) on normalized data 
prcomp() 