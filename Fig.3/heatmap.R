load("ba_top50.Rdata")
load("fu_top50.Rdata")
load("pheno_v.Rdata")
pheno <- pheno[,3:15]
rownames(ba_top50)==rownames(pheno)
rownames(fu_top50)==rownames(pheno)

library(psych)
library(ggplot2)
library(reshape2)
library(dplyr)
library(rlang)

# 计算细菌丰度与表型之间的Spearman相关系数
r_result <- corr.test(ba_top50,pheno,method = "spearman", adjust = "fdr") 

r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Microbe', 'Phenotype', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Microbe', 'Phenotype', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

# 将Microbe列因子反转排序
data$Microbe <- factor(data$Microbe, levels = rev(levels(data$Microbe)))


# 创建热图
theme <- theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        text = element_text(face='bold'),
        legend.key.size = unit(15, "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(angle=45,vjust=1, hjust=1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom",  # 将图例放置在图的下方
        axis.title.x = element_blank())  # 隐藏横坐标标题

p1 <- ggplot(data,aes(x=Phenotype, y=Microbe))+
  geom_tile(aes(fill=r),color = 'white',alpha = 0.8) +
  geom_text(aes(label = P_value_sig), color = 'black', size = 3,hjust = 0.5, vjust = 0.5) +
  ggsci::scale_fill_gsea() +
  theme +
  xlab('') +
  ylab('')+
  #coord_equal()+
  scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0)+
  labs(fill = "Spearman correlation")

p1


# 计算真菌丰度与表型之间的Spearman相关系数
r_result <- corr.test(fu_top50,pheno,method = "spearman", adjust = "fdr") 

r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Microbe', 'Phenotype', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Microbe', 'Phenotype', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

# 将Microbe列因子反转排序
data$Microbe <- factor(data$Microbe, levels = rev(levels(data$Microbe)))


# 创建热图
theme <- theme_bw(base_size = 10) +
  theme(panel.grid.major = element_blank(),
        text = element_text(face='bold'),
        legend.key.size = unit(15, "pt"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x=element_text(angle=45,vjust=1, hjust=1),
        axis.text.y = element_text(face = "italic"),
        legend.position = "bottom",  # 将图例放置在图的下方
        axis.title.x = element_blank())  # 隐藏横坐标标题

p2 <- ggplot(data,aes(x=Phenotype, y=Microbe))+
  geom_tile(aes(fill=r),color = 'white',alpha = 0.8) +
  geom_text(aes(label = P_value_sig), color = 'black', size = 3,hjust = 0.5, vjust = 0.5) +
  ggsci::scale_fill_gsea() +
  theme +
  xlab('') +
  ylab('')+
  #coord_equal()+
  scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0)+
  labs(fill = "Spearman correlation")

p2



