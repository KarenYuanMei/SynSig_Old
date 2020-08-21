library(ggplot2)
library(hrbrthemes)

df2 <- data.frame(Group=rep(c("Glutamate Receptors", "Scaffolds-Core", 'Other Scaffolds-Adaptors', 'Kinases', 'Phosphatases', 'Adhesion', 'Channels', 'GAPS-GEFs'), each=1),
                Phe_No=c(8,11, 14, 54, 24, 27, 7, 37))
                

# horizontal
ggplot(df2) +
  geom_bar( aes(x=Group, y=Phe_No), stat="identity", fill="coral", alpha=0.5) +

  labs(y= "Number of Proteins", x= 'Protein Group')+
  theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    coord_flip()

# horizontal
p <-ggplot(df2) +
  geom_bar( aes(x=Group, y=Phe_No), stat="identity", fill="coral", alpha=0.5) +

  labs(y= "Number of Proteins", x= 'Protein Group')+
  theme_bw()+
    theme(text = element_text(size=25))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x=element_blank())