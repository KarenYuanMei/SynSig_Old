library(ggplot2)
theme_set(theme_bw())  

# Data Prep
#df <- read.csv('/Users/karenmei/Documents/Synapse_Ontology/cell_location/regression_coef_cell_location.csv')
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/nonbrain_big_genes_pool_pipeline/regression_coef_cell_location.csv')

df$mean_type <- ifelse(df$mean < 0, "below", "above")  # above / below avg flag

df <- df[order(df$mean), ]  # sort

df$`Locations` <- factor(df$`Locations`, levels = df$`Locations`)  # convert to factor to retain sorted order in plot.


# Diverging Barcharts
ggplot(df, aes(x=`Locations`, y=mean, label='Mean Association')) + 
  geom_bar(stat='identity', aes(fill=mean_type), width=.5)  +
  geom_errorbar(aes(ymin=mean-stderr, ymax=mean+stderr),
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9))+
  scale_fill_manual(name="Contr. to Synapse Pred.", 
                    labels = c("Positive", "Negative"), 
                    values = c("above"="brown3", "below"="gray")) + 
  #labs(subtitle="Normalised mileage from 'mtcars'", 
      # title= "Diverging Bars") + 
      labs(x='Non-Neuronal Cell Locations', y='Mean Contribution')+
      theme(text = element_text(size=20))+
#axis.text.x = element_text(angle = 45, hjust = 1))+
theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background=element_blank())+theme(legend.position="bottom")+

  coord_flip()
