#===============
library(ggplot2)
library(hrbrthemes)

df2 <- data.frame(Group=rep(c('Gen. Dev Disorders', "Non-Syndromic Autism", "Syndromic Autism"), each=1),
                Phe_No=c(1.01, 1.51, 2.87))
                
c4 = c('gray', "pink", "darkred")

# horizontal
p <- ggplot(df2) +
  geom_bar( aes(x=Group, y=Phe_No), stat="identity", fill=c4, alpha=0.5) +

  labs(y= "Fold Enrichment", x= 'Non-Brain Predicted Synapse Proteins')+
  theme_bw()+
    theme(text = element_text(size=20))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())

