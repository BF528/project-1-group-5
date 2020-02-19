# title: "BF528 - Project 1: Microarray Based Tumor Classification"
# author: "Eetu Eklund, Mary T. Yohannes, Evie Wan, Salam Al-Abdullatif"
# date: "2/19/2020"

### Project code 

### Data preprocessing & quality control (Programmer)

##### Install and load packages that needed for the project (R version 3.6.0)

# install.packages("BiocManager")
# BiocManager::install(c("affy","affyPLM","sva","AnnotationDbi","hgu133plus2.db")
library(affy)
library(affyPLM)
library(sva)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(GSEABase)

##### Reading in the CEL files and normalize them together 

data <- ReadAffy(celfile.path = "/projectnb/bf528/users/group_5/project_1/samples/all_CELfiles")
norm_data <- rma(data)

##### Computing relative Log Expression (RLE) and Normalized Unscaled Standard Error (NUSE) scores of the microarray samples 

Pset <- fitPLM(data, normalize = TRUE, background=TRUE)

# RLE and NUSE stats for each sample including median and IQR - only interested in the medians 
RLE_stats <- RLE(Pset,type="stats")
NUSE_stats <- NUSE(Pset,type="stats")

# RLE and NUSE median distribution histograms
hist(RLE_stats["median",], xlab = "Median RLE", main = "Figure 1: Distribution of Median RLE Scores") 

hist(NUSE_stats["median",], xlab = "Median NUSE", main = "Figure 2: Distribution of Median NUSE Scores")

##### Correcting for batch effects 

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

##### Performing Principal Component Analysis (PCA) on the batch-corrected-normalized data

# scaling within each gene rather than within each sample
transpose_bf_corrected = t(bf_corrected)
scaled_transpose <- scale(transpose_bf_corrected) 

# re-transpose to return data to its original orientation
retranspose_data <- t(scaled_transpose) 

# running PCA 
pca_output_data <- prcomp(retranspose_data, scale. = FALSE, center = FALSE)

# view values for each of the principal components per sample 
#head(pca_output_data$rotation) 

##### PC1 vs PC2 with percent variability explained on each axis labels 

PC1 <- pca_output_data$rotation[,1] # first principal component
PC2 <- pca_output_data$rotation[,2] # second principal component

# SD, Proportion of Variance and cumulative proportion of each principal component 
data_summary <- summary(pca_output_data)
#data_summary$importance

# plot 
plot(PC2, PC1, 
     xlab = paste0("PC2 (", round(data_summary$importance[2,2]*100, digits = 2), "%)"),
     ylab = paste0("PC1 (", round(data_summary$importance[2,1]*100, digits = 2), "%)"),
     main = "Figure 3: PC1 vs PC2 with percent variability")

##### Outlier examination (condition: remove 3 SD from the mean)

# checking outlier for PC1
#boxplot(PC1)
which(PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1)-3*sd(PC1))

# checking outlier for PC2
#boxplot(PC2)
which(PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2)-3*sd(PC2)) 


# removing the outlier from the batch effects corrected data  
id <- which(!(PC2 > mean(PC2) + 3*sd(PC2) | PC2 < mean(PC2)-3*sd(PC2) |
                PC1 > mean(PC1) + 3*sd(PC1) | PC1 < mean(PC1)-3*sd(PC1)))

final_data <- (bf_corrected[,id])

# saving the final data 
write.csv(final_data, file = "final_data.csv")

### Noise filtering & dimensionality reduction (Analyst)

##### Filter based on expression values and variances 

## filter values < log2(15)

file <- read.csv("/projectnb/bf528/users/group_5/project_1/programmer_Mary/mary_output_data/final_data.csv", row.names=1)
data0 <- rowSums( file > log2(15) ) >= 0.2*ncol(file)
data1 <- file[data0, ]

## calculate row variance
rvar <- apply(data1, 1, var)
data1$Variance <- rvar

# median variance = 0.1585906
mdvar<- median(data1$Variance)

##### Test statistics and chi square

df <- ncol(file)-1
data1$Test_Stats <- ( df*(data1$Variance))/(mdvar)

###DELETE: data1$chiupper <- qchisq((1 - 0.99)/2, 132, lower.tail = FALSE)
chiupper <- qchisq((1 - 0.99)/2, df, lower.tail = FALSE)

#filter: variance larger than upper-bound value
chi_filtered <- data1 %>% rownames_to_column('gene') %>% filter(data1$Test_Stats >= chiupper) %>% column_to_rownames('gene')

#filter using coefficient of variation
chi_filtered2 <- as.matrix(chi_filtered[,c(1:133)])
cv_filtered <- chi_filtered %>% rownames_to_column('gene') %>% filter (sqrt(chi_filtered$Variance)/mean(chi_filtered2) >= 0.186)%>% column_to_rownames('gene')

## save to .csv file
write.csv(chi_filtered, 'expression_filtered.csv')

### Hierarchical clustering & subtype discovery (Analyst)

##### Cluster analysis

cluster_data <- select(cv_filtered, -c(Variance, Test_Stats) )

clusters <- hclust(dist(t(cluster_data)))

# Samples(patients) and their cluster(1/2)
clusters2 <- cutree(clusters, k=2)
table(clusters2)

##### Heat map

annotation <- read.csv("/projectnb/bf528/users/proj_metadata.csv")

annotation_df <- data.frame(row.names=annotation$geo_accession, subtype = annotation$cit.coloncancermolecularsubtype)
annotation_df$color <- ("red")
annotation_df$color[annotation_df$subtype=="C3"] <- ("blue")
str_split<- strsplit(colnames(cluster_data), split="_")

col_split <- do.call("rbind",strsplit(colnames(cluster_data), split="_"))
annotation_color <- annotation_df[col_split[,1], 2]
heatmap(as.matrix(cluster_data), ColSideColors = annotation_color)

##### T-test

##Cluster1
group1 <- cluster_data[,clusters2==1]

##Cluster2
group2 <- cluster_data[,clusters2==2]
group_tots <- cbind(group1, group2)

all.genes.p = apply(group_tots, 1, function(x) {t.test(x[1:56], x[57:133]) $p.value})


all.genes.t = apply(group_tots, 1, function(x) { 
  t.test(x[1:56], x[57:133]) $statistic})

p.adjust <- p.adjust(all.genes.p)

df.all.genes.p <- as.data.frame(all.genes.p)
df.all.genes.t <- as.data.frame(all.genes.t)
df.p.adjust <- as.data.frame(p.adjust)
df.all.genes <- cbind(df.all.genes.p, df.all.genes.t, df.p.adjust)

### filter p<0.05
filter_p <- df.all.genes %>% rownames_to_column('gene') %>% filter(df.all.genes$p.adjust<0.05) %>% column_to_rownames('gene')
nrow(filter_p)

#sort and filter by t-value
sorted <- filter_p[order(filter_p$all.genes.t), ]
filter.all.genes <- sorted %>% rownames_to_column('gene') %>% filter(sorted$all.genes.t>abs(15)) %>% column_to_rownames('gene')

write.csv(filter.all.genes, 'final_genes.csv')

### In-depth Analysis (Biologist)

##### Mapping probe IDs and top 10 differentially expressed genes

detach(package:dplyr, unload = TRUE)
select <- select(hgu133plus2.db, keys(hgu133plus2.db), column = c("SYMBOL"))

#results <- read.csv("differential_expression_results.csv")
results <- read.csv("final_genes.csv")

for (i in 1:dim(results)[1]){
  symbols <- select$SYMBOL[which(select$PROBEID == results[i,"X"])]
  results[i, "symbol"] <- symbols[1]
}

top10 <- tail(results, 10)
# table of top 10 differentially expressed genes 
top10