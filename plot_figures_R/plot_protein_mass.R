#Goal: analyze number of domains in the human genes as reported by Ensembl to determine if there are differences between synapse and non-synapse genes in this feature

#source: Ensembl

library(ggplot2)
library(dplyr)
library(hrbrthemes)

df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/Analyze_Synapse_Features/data_for_R/expanded_ENSEMBL_isoform_no_feature_histogram_R.csv")


head(df)

library(plyr)
mu <- ddply(df, "Type", summarise, grp.mean=mean(Values))
head(mu)


p <- ggplot(df, aes(x=Type, y=Values)) + 
  geom_violin()
# Function to produce summary statistics (mean and +/- sd)

# Change color by groups
dp <- ggplot(df, aes(x=Type, y=Values, fill=Type)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")
dp + theme_classic()