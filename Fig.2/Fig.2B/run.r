alpha =read.table("shannon",sep = "\t",row.names = 1,
                  check.names = F,stringsAsFactors = F)
colnames(alpha) = c("observed gene","shannon")
alpha = alpha[grep("DP",row.names(alpha)),]
mapping = read.table("mapping_file.txt",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T,row.names = 1)
mapping = mapping[row.names(mapping)%in%row.names(alpha),]
mapping = mapping[row.names(alpha),]
alpha$type = mapping$Gender

cheek = alpha[mapping$Type=="Ck",]
forehead = alpha[mapping$Type =="Fh",]
nose = alpha[mapping$Type =="AI",]

library(ggplot2)
library(reshape2)
library(ggpubr)
data1 = alpha
data1$type1 = mapping$Type
data1$type2 = mapping$age_type
shannon = data1[,2:3]
obs = data1[,c(1,3)]

my_comparisons <- list(c("old", "young"), c("old", "older"), c("young", "older"))
data1$type2 = factor(data1$type2,levels = c("young","old","older"))
ggplot(data1,aes(x= type2,y = shannon,fill=type2))+geom_boxplot()+facet_grid(. ~ type1)+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")
 


