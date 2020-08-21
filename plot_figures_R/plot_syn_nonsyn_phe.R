library(ggplot2)
library(hrbrthemes)

              


df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/syndromic_nonsyndromic_pheno_R.csv")

# Change color by groups
dp <- ggplot(df, aes(x=Group, y=Phe_No, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Autism Genes", y = "Number of Phenotypes \n (Excluding Nervous System)")
dp + theme_classic()+
scale_fill_manual(values=c("pink", "coral"))+scale_y_log10()+theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")

df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/pred_neg_nb_R_new.csv.csv")

# Change color by groups
dp <- ggplot(df, aes(x=Group, y=Phe_No, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Autism Genes", y = "Number of Phenotypes \n (Excluding Nervous System)")
dp + theme_classic()+
scale_fill_manual(values=c("pink", "coral"))+scale_y_log10()+theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")

#====June 1,2020=====

df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/syndromic_nonsyndromic_nohead_R.csv")

# Change color by groups
dp <- ggplot(df, aes(x=Group, y=Phe_No, fill=Group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="Autism Genes", y = "Number of Phenotypes \n (Excluding Nervous System, Head or Neck)")
dp + theme_classic()+
scale_fill_manual(values=c("pink", "coral"))+scale_y_log10()+theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    theme(legend.position = "none")
