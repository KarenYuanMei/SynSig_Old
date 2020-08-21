library(ggplot2)
library(hrbrthemes)

df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/adult_diff.csv")

# Change color by groups
dp <- ggplot(df, aes(x=Group, y=Adult_Expression, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Autism Genes", y = "Adult Expression (RPKM)")
dp + theme_classic()+
scale_fill_manual(values=c("pink", "coral"))+scale_y_log10()+

theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")



df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/nonbrain_big_genes_pool_pipeline/synapse_fetal_diff.csv")

# Change color by groups
dp <- ggplot(df, aes(x=Group, y=Fetal_Expression, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Autism Genes", y = "Fetal Expression (RPKM)")
dp + theme_classic()+
scale_fill_manual(values=c("pink", "coral"))+
theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
