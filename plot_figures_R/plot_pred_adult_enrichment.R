library(ggplot2)
library(hrbrthemes)

df2 <- data.frame(Group=rep(c('Cortex', "Striatum", "CTX/STR", 'SynGO'), each=1),
                Fold_Enrichment=c(2.15, 2.68, 2.9, 2.62), significance=c(247, 178, 181, 131, 203, 125))

#===================



df2 <- data.frame(Group=rep(c('Cortex', "Striatum", "Cortex and Striatum"), each=1),
                Fold_Enrichment=c(2.27, 2.77, 3.03))

#df2 <- data.frame(Group=rep(c("CTX", "STR", "CTX/STR Overlap" ), each=1),Fold_Enrichment=c(2.62, 2.15, 2.68, 2.9))

df2$Group <- factor(df2$Group, levels = df2$Group)


c4 = c("azure3", "skyblue", "navyblue")

#p <- ggplot(df, aes(x = reorder(f.name, -age), y = age))

p <- ggplot(df2, aes(x = Group, y = Fold_Enrichment))
p <- p + geom_bar(stat="identity",fill=c4)+
labs(y= "Fold Enrichment with Predicted", x= "Synapse Proteomics")+
  theme_bw()+
    theme(text = element_text(size=20))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p+coord_flip()



#=========================
df2 <- data.frame(Group=rep(c('SynSysNet', 'SynDB', 'SynGO'), each=1),
                Fold_Enrichment=c(5.41, 4.86, 4.93), significance=c(131, 203, 125))

c4 = c("darkgreen", "darkgray", 'seagreen2')

#p <- ggplot(df, aes(x = reorder(f.name, -age), y = age))

p <- ggplot(df2, aes(x = reorder(Group, Fold_Enrichment), y = Fold_Enrichment))
p <- p + geom_bar(stat="identity",fill=c4)+
labs(y= "Fold Enrichment with Predicted", x= 'Synapse Databases')+
  theme_bw()+
    theme(text = element_text(size=20))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p+coord_flip()

