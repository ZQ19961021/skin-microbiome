data1 = read.table("filter_data.txt",sep = "\t",check.names = F,stringsAsFactors = F,
                   header = T,row.names = 1)
library(ggplot2)
ggplot(data1, aes(x=name, y=estimate)) +
  geom_segment( aes(x=name, xend=name, y=0, yend=estimate)) +
  geom_point( aes(color = site), size=4) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("Value of Y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
