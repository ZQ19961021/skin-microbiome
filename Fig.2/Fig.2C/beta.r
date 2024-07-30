mgs = read.csv("species0619.profile",sep = "\t",check.names = F,stringsAsFactors = F,
               header = T,row.names = 1)
mgs = mgs[grep("Bacteria",row.names(mgs)),]
library(vegan)
mapping = read.table("mapping_file.txt",sep = "\t",check.names = F,
                     stringsAsFactors = F,header = T)
mapping = mapping[grep("DP",mapping[,1]),]
row.names(mapping) = mapping[,1]
mapping = mapping[row.names(mapping)%in%colnames(mgs),]
mgs = mgs[,row.names(mapping)]
mgs_dist = as.matrix(vegdist(t(mgs)))

library(reshape2)
for(i in 1:ncol(mgs_dist)){
  for(j in i:ncol(mgs_dist)){
    mgs_dist[i,j] = 0
  }
}
mgs_dist = melt(mgs_dist,factorsAsStrings = F)
mgs_dist = data.frame(mgs_dist,stringsAsFactors = F,check.names = F)
mgs_dist = mgs_dist[mgs_dist$value!=0,]
mgs_dist[,1] = as.character(mgs_dist[,1])
mgs_dist[,2] = as.character(mgs_dist[,2])
row.names(mapping) = mapping[,1]
mgs_dist$Var1_type = mapping[mgs_dist[,1],5] 
mgs_dist$Var2_type = mapping[mgs_dist[,2],5] 
mgs_dist$type = "inter"
mgs_dist$type[mgs_dist$Var1_type==mgs_dist$Var2_type] = "intra"
mgs_dist$type1 = paste(mgs_dist$Var1_type,mgs_dist$Var2_type,mgs_dist$type,sep = "_")
mgs_dist$type1[mgs_dist$type1=="old_young_inter"] = "young_old_inter"
mgs_dist$type1[mgs_dist$type1=="older_young_inter"] = "young_older_inter"
mgs_dist$type1[mgs_dist$type1=="older_old_inter"] = "old_older_inter"

library(ggplot2)
ggplot(mgs_dist,aes(x = type1,y = value))+geom_boxplot()+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  facet_grid(. ~ type)
wilcox.test(mgs_dist[,3][mgs_dist$type=="intra"&mgs_dist$type1=="young_old_inter"],
            mgs_dist[,3][mgs_dist$type=="intra"&mgs_dist$type1=="young_older_inter"])

