#update:April 10, 2020

library(ggplot2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(hrbrthemes)

#df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/validate_new_exp/fetal/fetal_specific_trajectory.csv")

df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/fetal_specific_trajectory.csv")
df %>%
  ggplot( aes(x=Age, y=mean)) +
    geom_line( color="grey") +
    geom_point(shape=21, color="black", fill="brown3", size=6) +
    theme_ipsum() +
    theme_bw()+scale_x_log10()+
    labs(y= "Brain Expression (TPM)", x = "Age (Days)")+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# Add the regression line
ggplot(df, aes(x=Age, y=mean)) + 
  geom_point()+
  geom_smooth(color="black")+
  geom_point(shape=21, color="black", fill="brown3", size=6) +
    theme_ipsum() +
    theme_bw()+scale_x_log10()+ylim(0,65)+
    labs(y= "Brain Expression (TPM)", x = "Age (Days)")+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
new <- select(df, Age, mean)

names(new)[names(new) == "mean"] <- "Fetal"

#df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/NetworkClass/Entry_Ontology/synapse_10/validate_new_exp/fetal/adult_specific_trajectory.csv")
df <- read.csv("/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/adult_specific_trajectory.csv")


adult_mean <- select(df, mean)
final <- cbind(new, adult_mean)
names(final)[names(final) == "mean"] <- "Adult"

final.m<-melt(final, id.vars="Age", measure.vars=c("Fetal", "Adult"))

names(final.m)[names(final.m) == "variable"] <- "Group"

ggplot(final.m, aes(x=Age, y=value, color=Group)) +
  geom_point()+
  geom_smooth(aes(fill=Group))+
  geom_point(shape=21, color="black", size=6, aes(fill=Group))+
  scale_color_manual(values=c('red', 'gray'))+
  scale_fill_manual(values=c('red', 'gray'))+
  theme_ipsum() +
    theme_bw()+scale_x_log10()+ylim(0,65)+
    labs(y= "Brain Expression (TPM)", x = "Age (Days)")+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



