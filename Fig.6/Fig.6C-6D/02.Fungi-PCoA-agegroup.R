#真菌的PCoA
# Load package
library(vegan)
library(ggplot2)
library(ggthemes)
# Load data
load("ch_fu_sep.Rdata")
load("China_phenotype.RData")
rownames(ch_fu)==rownames(pheno)
site <- pheno$type1
age <- pheno$age
fu_sample <- data.frame(rownames(ch_fu),age,site)
colnames(fu_sample)[1] <- "ID"

#pcoa
# vegdist函数，计算距离；method参数，选择距离类型
distance <- vegdist(ch_fu, method = 'bray')
# 对加权距离进行PCoA分析
pcoa <- cmdscale(distance, k = (nrow(ch_fu) - 1), eig = TRUE)
## plot data
# 提取样本点坐标
plot_data <- data.frame({pcoa$point})[1:2]

# 提取列名，便于后面操作。
plot_data$ID <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig

#为样本点坐标添加分组信息
plot_data$ID==fu_sample$ID
group<- cut(fu_sample$age, breaks = c(20,35, 50, 65), labels = c("young","middle_age","aged"), right=FALSE)
plot_data$group <- group
fu_sample$group <- group
# 计算加权bray-curtis距离
dune_dist <- vegdist(ch_fu, method="bray", binary=F)
dune_pcoa <- cmdscale(dune_dist, k=(nrow(ch_fu) - 1), eig=T)

dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)

colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)

dune_pcoa_result <- cbind(dune_pcoa_points,group,site)
colnames(dune_pcoa_result)
head(dune_pcoa_result)
library(ggplot2)

#PERMANOVA
ch_fu.div <- adonis2(ch_fu ~ group, data = fu_sample, permutations = 999, method="bray")
ch_fu.div
ch_fu_adonis <- paste0("adonis R2= ",round(ch_fu.div$R2,2), "; P-value= ", ch_fu.div$`Pr(>F)`)

p1 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, fill=group)) +
  geom_point(shape = 21,size=3) +
  stat_ellipse(level=0.95,geom = "polygon",alpha=0.2)+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  scale_fill_manual(values = c("#e77d72",'#6e9df8',"#53b74c"))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=ch_fu_adonis)  +
  theme_bw()+
  theme(panel.grid = element_blank())
p1
ggsave(p1,file = "fungi-agegroup-PCoA.pdf",width =6,height = 4)


#绘制箱线图观察群落Beta多样性高低水平
##根据分组获得组内距离矩阵
library(ggsignif)
dis <- as.matrix(distance)
young <- subset(fu_sample,group=="young")$ID
dis_young <- dis[young,young]
dis_young <- as.vector(dis_young) #将矩阵转化为向量，以便用于作图和统计

middle_age <- subset(fu_sample, group == "middle_age")$ID
dis_middle_age <- dis[middle_age,middle_age]
dis_middle_age <- as.vector(dis_middle_age)

aged <- subset(fu_sample,group=="aged")$ID
dis_aged <- dis[aged,aged]
dis_aged <- as.vector(dis_aged)

##构建作图数据集
dat <- data.frame(
  dis = c(dis_young, dis_middle_age,dis_aged),
  group = factor(c(
    rep('young', length(dis_young)), 
    rep('middle_age', length(dis_middle_age)),
    rep('aged', length(dis_aged))),
    levels = c('young', 'middle_age',"aged"))
)


##使用 ggplot2 绘制各组内 Bray-curtis 距离指数分布的箱线图，并添加Wilcoxon 秩和检验显著线
my_comparisons <- list(c("young", "middle_age"), 
                       c("young", "aged"),
                     c("middle_age", "aged"))

library(ggplot2)

p2 <- ggplot(dat, aes(group, dis)) +
  geom_boxplot(aes(fill = group), width = 0.5) +
  scale_fill_manual(values = c("#e77d72",'#6e9df8',"#53b74c")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'), legend.position = 'none') +
  stat_compare_means(comparisons=my_comparisons,
                     method="wilcox.test",
                     label="p.signif"
  )+
  labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')+
theme_bw()+
  theme(panel.grid = element_blank())
p2
ggsave(p2,file="Beta diversity-age group.png")
ggsave(p2,file="Beta diversity-age group.pdf",p2 <- ggplot(dat, aes(group, dis)) +
  geom_boxplot(aes(fill = group), width = 0.5) +
  scale_fill_manual(values = c("#e77d72",'#6e9df8',"#53b74c")) +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = 'black'), legend.position = 'none') +
  stat_compare_means(comparisons=my_comparisons,
                     method="wilcox.test",
                     label="p.signif"
  )+
  labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')+
theme_bw()+
  theme(panel.grid = element_blank())
p2
ggsave(p2,file="Beta diversity-age group.png")
ggsave(p2,file="Beta diversity-age group.pdf",height = 4,width = 4)
    