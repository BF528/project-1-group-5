## Steps 1 & 2: install and load packages (R version 3.6.0)

install.packages("BiocManager")
BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)


## Step 3: Reading in the CEL files and normalize them together
data <- ReadAffy(celfile.path = "/projectnb/bf528/users/group_5/project_1/samples/all_CELfiles")
norm_data <- rma(data)


## Step 4: Computing relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) scores of the microarray samples 
Pset <- fitPLM(data, normalize = TRUE, background=TRUE)

# RLE and NUSE stats for each sample including median and IQR - only interested in the medians 
RLE_stats <- RLE(Pset,type="stats")
NUSE_stats <- NUSE(Pset,type="stats")

# RLE and NUSE median distribution histograms
hist(RLE_stats["median",], xlab = "Median RLE", main = "Distribution of Median RLE Scores") 

hist(NUSE_stats["median",], xlab = "Median NUSE", main = "Distribution of Median NUSE Scores")


## Step 5: Correcting for batch effects 
# transform the normalized data into a matrix: rows = probesets and columns = samples
matrix_data <- exprs(norm_data)

# read in the meta data 
meta_data <- read.csv("/projectnb/bf528/users/group_5/project_1/programmer_Mary/mary_output_data/proj_metadata.csv")

# only interested in "normalizationcombatbatch" and "normalizationcombatmod"
batch = meta_data$normalizationcombatbatch
mod = model.matrix(~normalizationcombatmod, data=meta_data)  

# correcting for batch effects while preserving certain features of interest such as tumor and MMR status
bf_corrected <- ComBat(dat = matrix_data, batch = batch, mod = mod)

# write out the expression data to a CSV file
write.csv(bf_corrected, file = "expressed_data.csv")


## Step 6: Performing Principal Component Analysis (PCA) on the batch-corrected-normalized data 
# scaling within each gene rather than within each sample
transpose_bf_corrected = t(bf_corrected)
scaled_transpose <- scale(transpose_bf_corrected) 

# re-transpose to return data to its original orientation
retranspose_data <- t(scaled_transpose) 

# running PCA 
pca_output_data <- prcomp(retranspose_data, scale. = FALSE, center = FALSE)

# view values for each of the principal components per sample 
head(pca_output_data$rotation) 


## Step 7: PC1 vs PC2 with percent variability explained on each axis labels }
PC1 <- pca_output_data$rotation[,1] # first principal component
PC2 <- pca_output_data$rotation[,2] # second principal component

# SD, Proportion of Variance and cumulative proportion of each principal component 
data_summary <- summary(pca_output_data)
data_summary$importance

# plot 
plot(PC2, PC1, 
     xlab = paste0("PC2 (", round(data_summary$importance[2,2]*100, digits = 2), "%)"),
     ylab = paste0("PC1 (", round(data_summary$importance[2,1]*100, digits = 2), "%)"))


## Step 8: Outlier examination (condition: remove 3 SD from the mean)
# checking outlier for PC1
boxplot(PC1)
which(PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1)-3*sd(PC1))

# checking outlier for PC2
boxplot(PC2)
which(PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2)-3*sd(PC2)) 


# removing the outlier from the batch effects corrected data  
id <- which(!(PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2)-3*sd(PC2) |
                PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1)-3*sd(PC1)))

final_data <- (bf_corrected[,id])

# saving the final data 
write.csv(final_data, file = "final_data.csv")