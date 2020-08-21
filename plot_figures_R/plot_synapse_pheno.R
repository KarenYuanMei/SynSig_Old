library(ggplot2)

#df <- read.csv("/Users/karenmei/Documents/Synapse_Ontology/synapse_diseases/synapse_phe_no.csv")



df2 <- data.frame(Group=rep(c("Synapse", "Negative"), each=1),
                Phe_No=c(62.136, 45.281),
                Stderr=c(4.12, 2.87))

c4 = c("royalblue", "skyblue")

# horizontal
ggplot(df2) +
  geom_bar( aes(x=Group, y=Phe_No), stat="identity", fill=c4, alpha=0.5) +

  geom_errorbar( aes(x=Group, ymin=Phe_No-Stderr, ymax=Phe_No+Stderr), width=0.4, colour="black", alpha=0.9, size=1)+
  labs(y= "Number of Phenotypes", x= 'Protein Group')+
  theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
