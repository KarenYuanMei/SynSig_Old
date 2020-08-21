library(ggplot2)

df2 <- data.frame(Groups=rep(c("SynGO", "SynV"), each=4),
                Model=rep(c("Fetal Human", "Fetal NGN2-Derived", "Adult Striatum", "Adult Cortex"),2),
                Fold_Enrichment=c(2.4722222222222223, 1.5360824742268042,2.4,2.1584158415841586, 2.8958333333333335, 2.0930232558139537,3.1515151515151514, 2.5597014925373136))
head(df2)

p <- ggplot(data=df2, aes(x=Model, y=Fold_Enrichment, fill=Groups)) +
geom_bar(stat="identity", position=position_dodge())+
scale_fill_brewer(palette="Reds", name="Model", labels = c("SynGO", "SynCh"))+
  theme_minimal()

p+labs(x="Experimental Model", y="Fold Enrichment")+
theme(text = element_text(size=20))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.text.x=element_text(angle=30, hjust=1))