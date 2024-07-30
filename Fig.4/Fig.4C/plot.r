data1 = read.table("filter_module.txt",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T)
#data1[,1] = gsub(".*;s_","",data1[,1])
library(ggplot2)
data1 = data1[,c(1,8,9,10)]
data1$q = log10(data1$p.value)
data1$q[data1$Elderly_mean>data1$Young_mean] = -data1$q[data1$Elderly_mean<data1$Young_mean]
data1 = data1[order(data1$q),]
name = data1$species
data1$species = factor(data1$species,levels = name)
ggplot(data1,aes(data1[,1],y = q))+geom_bar(stat = 'identity', position = 'dodge')+
  coord_flip()+theme_classic()

