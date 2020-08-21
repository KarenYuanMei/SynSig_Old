library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/pred_genes_above_4.9.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/nb_validated.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/nb_validated_db.csv')



ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

df <- read.csv('/Users/karenmei/ego2.csv')
ego2 <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/ego2.csv', row.names="ID")
dotplot(ego2, showCategory=15)


#========fig 3=============
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)

#df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/validate_new_exp/adult/new_adult_pred.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/new_adult_pred.csv')
ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))
dotplot(ego2, showCategory=15)+theme(text = element_text(size=20))


p1 <- cnetplot(ego2, foldChange=df$genes)
## categorySize can be scaled by 'pvalue' or 'geneNum'
#p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
#p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
p3 <- cnetplot(ego2), node_label="all") 

#======fig 4======================

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(DOSE)

#df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/validate_new_exp/fetal/fetal_only_val.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/fetal_only_val.csv')

ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))
dotplot(ego2, showCategory=15)+theme(text = element_text(size=20))


p1 <- cnetplot(ego2, foldChange=df$genes)
## categorySize can be scaled by 'pvalue' or 'geneNum'
#p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
#p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
p3 <- cnetplot(ego2), node_label="all") 

#=============Fig 6==========================
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

#df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nb_validated_db.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nonbrain_pred_genes_above_4.67.csv')




ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

#==============Fig 7========================
#=============Fig 7==========================
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

#df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nb_validated_db.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nb_proteomics_new_syndromic.csv')




ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)


====SynSig enrichment===================
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

#df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nb_validated_db.csv')
#df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/nb_proteomics_new_syndromic.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/new_adult_pred.csv')




ego2 <- enrichGO(gene         = df$genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

