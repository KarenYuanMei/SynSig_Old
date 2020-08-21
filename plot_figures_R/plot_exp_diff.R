#===========================

df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/Analyze_Synapse_Features/GTEX_tissue_expression/compare_train_synch_nonbrain.csv")

df.m<-melt(df, id.vars="Tissues", measure.vars=c("Training", "Predicted"))

# prepare data for plotting
plotting_df <-
  df.m %>% 
  # a trick!
  mutate(value = if_else(variable == "Training", -value, value))

## find the order
temp_df <-
  plotting_df %>% 
  filter(variable == "Predicted") %>% 
  arrange(value)
the_order <- temp_df$Tissues
# plot
p <- 
  plotting_df %>% 
  ggplot(aes(x = Tissues, y = value, group = variable, fill = variable)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = the_order) +
  # another trick!
  scale_y_continuous(breaks = seq(-7, 7, 1), labels = abs(seq(-7, 7, 1))) +
  #labs(x = "Features", y = "Fold Change", title = "Genomic and Protein Features") +
  labs(x = "Non-Brain Regions", y = "Fold Change") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill =  "grey90")) +
  # reverse the order of items in legend
  # guides(fill = guide_legend(reverse = TRUE)) +
  # change the default colors of bars
  scale_fill_manual(values=c("gray", "brown3"),
                    name="",
                    breaks=c("Training", "Predicted"),
                    labels=c("Training", "SynSig"))+
  theme_ipsum() +
  theme_bw()+

  theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(panel.background=element_blank())
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


print(p)


