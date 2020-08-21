
df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/final_cells_fc.csv')

obj <- ggplot(df)+geom_point(aes(x=Fold_Change, y=reorder(Cell_Line, Fold_Change), color=Group, size=5))+scale_color_manual(name='Group', labels=c('Training', 'SynSig'), values=c("Training"="#00ba38", "SynSig"="#f8766d"))+scale_x_continuous(limits=c(0,10), breaks=c(1,2,3,4, 5, 6, 7, 8, 9, 10))

obj + labs(y="Non-Neuronal Cell Lines", x = "Synapse/Non-Synapse Fold Change")+theme_bw()+labs(color = "Group")+coord_flip()


df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/SynSig_Cell_Lines_Exp.csv')

obj <- ggplot(df)+geom_point(aes(x=Fold_Change, y=reorder(Cell_Line, Fold_Change), size=8))+scale_x_continuous(limits=c(0,10), breaks=c(1,2,3,4, 5, 6, 7, 8, 9, 10))

obj + labs(y="Non-Neuronal Cell Lines", x = "SynSig Positives/Negatives")+theme_bw()


df <- read.csv('/Users/karenmei/Documents/Synapse_Paper_Code/synapse_11/brain_RNA_big_gene_pool_pipeline/Training_Cell_Lines_Exp.csv')

obj <- ggplot(df)+geom_point(aes(x=Fold_Change, y=reorder(Cell_Line, Fold_Change), color=Corr_Significance, size=8))+scale_x_continuous(limits=c(0,10), breaks=c(1,2,3,4, 5, 6, 7, 8, 9, 10))

obj + labs(y="Non-Neuronal Cell Lines", x = "Training Positives/Negatives")+theme_bw()
