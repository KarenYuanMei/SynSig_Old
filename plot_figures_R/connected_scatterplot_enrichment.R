# Libraries
    
#June 13, 2020
library(ggplot2)
library(dplyr)
library(hrbrthemes)

# Load dataset from github
df <- data.frame(Number_of_Consensus_Sources=c(0, 1, 2, 3, 4, 5, 6, 7, 8), Fold_Enrichment=c(0.37641154328732745, 0.9552238805970149, 1.3944954128440368,2.0428571428571427, 3.2448979591836733, 3.888888888888889, 5.230769230769231, 6.555555555555555, 6.666666666666667), PValue=c(0.9999999999712506, 0.7822656332319116, 7.198742237365059e-06, 3.993568924425307e-18, 2.233431381273741e-46, 1.382067671628916e-51, 2.3457102490953157e-71, 7.969174521432951e-38, 1.3100118156995085e-14))

# Plot
df %>%
  ggplot( aes(x=Number_of_Consensus_Sources, y=Fold_Enrichment, size=-log(PValue))) +
    geom_line( color="grey") +
    geom_point(shape=21, color="black", fill="darkred", size=10) +
    theme_ipsum() +
    theme_bw()+
    #xlim(0, 13)+ylim(0, 10)+
    labs(y= "SynSig Fold Enrichment", x = "Number of Consensus Sources")+
    theme(text = element_text(size=25))+
    #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    

