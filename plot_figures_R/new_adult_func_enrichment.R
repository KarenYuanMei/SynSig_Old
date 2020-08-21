

#062820================
library(ggplot2)

df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/update_enrichment/New_Adult_Pred_enrichment_R.csv')

color_table <- tibble(
  Land_cover = c("Metabolism", "Transport", "Assembly", "Morphology"),
  Color = c( "coral", "navyblue","brown3", "orchid" ))

p <- ggplot(df, aes(x=reorder(term_name, -index), y=Gene_Ratio, fill=categories))+geom_bar(stat="identity")+
scale_fill_manual(values = color_table$Color)

p<- p + expand_limits( y = 0)

p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.4))


p + labs(y="Proportion of SynSig", x = "Enriched Functions")+theme_bw()+theme(text = element_text(size=25))+coord_flip()
