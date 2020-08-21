library(ggplot2)
theme_set(theme_bw())

# Load dataset from github
df <- data.frame(x=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), y=c(317, 128, 84, 70, 92, 94, 99, 96, 84, 76, 54, 25, 20, 8, 1))

# Draw plot
ggplot(df, aes(x=x, y=y)) + 
  geom_bar(stat="identity", width=.5, fill="tomato3") + 
  labs(title="Novel Genes Not Yet Validated by Literature", 
       subtitle="Number of Genes Validated by Literature: 931 out of 1248 Predicted Genes") + 
  labs(y= "Number of Predicted Genes", x = "Number of Supporting Sources")+coord_flip()
  
df %>% 
	group_by(x) %>% 
	mutate(highlight_flag = ifelse(x == 0, T, F)) %>% 
	ggplot(aes(x =x, y = y)) +
	geom_bar(aes(fill = highlight_flag), stat = 'identity') +
	scale_fill_manual(values = c('#595959', 'red')) +guides(fill=FALSE)+
	coord_flip()+
	labs(x = 'Number of Supporting Sources',y = 'Number of Predicted Genes' ,title = 'Novel Genes Not Yet Validated by Literature')