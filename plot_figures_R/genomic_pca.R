library("FactoMineR")
library("factoextra")
library(tidyverse)
library("corrplot")


df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/expanded_pca_R.csv")

#make "Genes" column the index
rownames(df) <- df$Genes

#remove the extra "Genes" column and the "target column"
new <- select (df, -c(Genes, target))
head(new)

new <- log(new)

#reduce the df into 9 orthogonal dimensions using PCA
res.pca <- PCA(new, scale.unit=TRUE, ncp=9, graph = FALSE)

#get the eigenvalues
eig.val <- get_eigenvalue(res.pca)

#plot the PC variances:
fviz_eig(res.pca, addlabels = TRUE)
var <- get_pca_var(res.pca)

#Find the contribution from each feature
feature_var<- t(var$contrib)
#Find the contribution from each PC
pc_var <- res.pca$eig[,c("percentage of variance")]
#convert the PC variance vector to decimals
pc_var <- pc_var * 0.01

#multiply the PC contribution with Feature contribution
feature_pc <- feature_var* pc_var
feat_pc_t <- t(feature_pc)
feat_pc_t
#change the col and row names

colnames(feat_pc_t) <- c("PC1", "PC2", "PC3", "PC4", "PC5", 'PC6', 'PC7', 'PC8', 'PC9')
rownames(feat_pc_t) <- c("Phos Site No", 'Domain No', 'Protein Mass', 'Transcript Length', 'Gene Length', 'Exon No', 'CDS Length', 'Isoform No', 'Amino Acid Length')

new<- feat_pc_t[,c("PC1","PC2","PC3", "PC4")]

par(mar = c(7, 7, 5, 5))
#plot the normalized contribution of each feature for each PC:
corrplot(new, is.corr=FALSE, method="circle", tl.col="black")