#真菌的shannon指数与年龄之间的相关性
library(vegan)
load("china_fungi.RData")
load("pheno_de.Rdata")
diver <- diversity(ch_fu,index = "shannon")
rownames(ch_fu)==rownames(pheno_p)
age <- pheno_p$age
diver_age <- data.frame(diver,age)
plot(age,diver)

library(psych)
library(ggplot2)

cor.result <- corr.test(diver,age,method ="spearman")
result.r <- cor.result$r
result.p <- cor.result$p

p1 <- ggplot(diver_age,aes(x=age,y=diver))+ 
  geom_point(size=1.5,shape=21,color="black",fill='#f29325')+
  geom_smooth(method=lm,color="#F08080")+
  labs(x='Age',y="Fungal Shannon Index")+
  annotate("text", x = 28 , y = 2.8,label = "Spearman, r=0.06, p=0.09")+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

p1
#ggsave(p1,file='fungal diversity_age.png')

ggsave(p1,file='fungal diversity_age-1.pdf',width = 5,height = 4)

#不同年龄段之间真菌多样性差距
library(ggpubr)
diver_age$group<- cut(diver_age$age, breaks = c(20,35, 50, 65), labels = c("young","middle-age","aged"), right=FALSE)
diver_age$group
table(diver_age$group)

my_comparisons <- list(c("young", "middle-age"), c("young", "aged"),c("middle-age", "aged"))

p2<- ggplot(diver_age, aes(diver_age$group,diver,fill = group)) +
  geom_boxplot(width=0.5) +
  labs(
    x = NULL,
    y = "Fungal shannon index"
  ) +
  scale_fill_manual(values = c("#e77d72",'#6e9df8',"#53b74c"))+
  stat_compare_means(comparisons=my_comparisons,
                     label.y = c(2.8, 3.1, 3.4),
                     method="wilcox.test",
                     label="p.signif"
  )+
  theme_bw(base_size = 10)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2
#ggsave(p2,file="Fungal_diversity_agegroup.png",width=4.46,height = 4)
ggsave(p2,file="Fungal_diversity_agegroup.pdf",width = 4,height = 4)


