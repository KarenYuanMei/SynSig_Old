library(ggplot2)

df2 <- data.frame(Groups=rep(c("GO Synapse", "SynGO", "SynSig"), each=3),
                Model=rep(c("Fetal Brain", "NGN2-Derived", "Overlap: Fetal/NGN2"),3),
                Fold_Enrichment=c(1.69, 1.21, 1.80, 2.34, 1.53, 2.53, 2.65, 1.76,2.95))
head(df2)

p <- ggplot(data=df2, aes(x=Model, y=Fold_Enrichment, fill=Groups)) +
geom_bar(stat="identity", position=position_dodge())+
scale_fill_brewer(palette="Reds", name="Synapse Model", labels = c("GO Synapse", "SynGO", "SynSig"))+
  theme_minimal()

p+labs(x="Proteomics Screens", y="Fold Enrichment")+
theme(text = element_text(size=20))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=30, hjust=1))

p + scale_y_continuous(limits=c(0, 3))