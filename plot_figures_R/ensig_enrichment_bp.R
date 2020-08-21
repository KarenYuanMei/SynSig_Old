#062820================
library(ggplot2)

df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/update_enrichment/New_EnSig_enrichment_R.csv')

p <- ggplot(df, aes(x=reorder(term_name, -adjusted_p_value), y=Gene_Ratio, fill=negative_log10_of_adjusted_p_value))+geom_bar(stat="identity", fill='dodgerblue')

p<- p + expand_limits( y = 0)

p <- p + scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5))


p + labs(y="Proportion of EnSig", x = "Enriched Functions")+theme_bw()+theme(text = element_text(size=25), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+coord_flip()
