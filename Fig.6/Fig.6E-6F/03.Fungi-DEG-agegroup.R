# 导入必要的库
library(ggplot2)
library(dplyr)
library(ggsignif)
library(tidyverse)
#library(Cairo)
library(sysfonts)
library(showtextdb)
library(showtext)

# 假设df是你的数据，其中包含了两组的所有重复
load("ch_fu_sep.Rdata")
colnames(ch_fu)
load("pheno_de.Rdata")
data <- as.data.frame(t(ch_fu))
Taxonomy <- rownames(data)
data <- data.frame(Taxonomy,data)

metadata <- data.frame(rownames(pheno_p),pheno_p$age)
colnames(metadata) <- c("Sample",'Age')
metadata$Group<- cut(metadata$Age, breaks = c(20,35, 50, 65), labels = c("young","middle-age","aged"), right=FALSE)
table(metadata$Group)
metadata <- metadata[!grepl("middle-age",metadata$Group),]
dim(metadata)
metadata <- metadata[,-2]

data_long <- data %>%
  pivot_longer(
    cols = -Taxonomy,
    names_to = "Sample",
    values_to = "Value"
  )

data_merge <- merge(data_long, metadata, by.x = "Sample", by.y = "Sample")

# 创建一个空的数据框来存储结果
result <- data.frame(Taxonomy = character(), p_value = numeric())

# 对每个物种在两个组中进行Wilcoxon秩和检验，并计算每个组的平均丰度
for (TAX in unique(data_merge$Taxonomy)) {
  df <- data_merge[data_merge$Taxonomy == TAX, ]
  p_value <- round(wilcox.test(df$Value ~ df$Group,exact = TRUE,correct = TRUE,paired = FALSE,conf.level = 0.95)$p.value, 4)
  result <- rbind(result, data.frame(Taxonomy = TAX, p_value = p_value)) #rbind()按行添加且合并表格
}


# 筛选出差异显著的物种
significant_TAX <- result %>%
  filter(p_value < 0.05)

# 添加显著性标记
significant_TAX <- significant_TAX %>%
  mutate(signif_mark = case_when(
    p_value < 0.05 & p_value >= 0.01 ~ "*",
    p_value < 0.01 & p_value >= 0.001 ~ "**",
    p_value < 0.001 ~ "***",
    TRUE ~ ""  # 对于其他情况，保持空白
  ))

#write.csv(significant_TAX, file = "significant_TAX.csv", row.names = FALSE)

# 将原数据表只留下显著的数据
data_merge_new <- merge(data_merge, significant_TAX, by = "Taxonomy", all.x = FALSE)

# 计算每个分类单元在每个组中的平均相对丰度
data_merge_new_summary <- data_merge_new %>%
  group_by(Taxonomy, Group) %>%
  summarise(mean_abundance = mean(Value))

# 创建柱状图
p1 <- ggplot(data_merge_new_summary, aes(x = Taxonomy, y =mean_abundance, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge()) + #stat = "identity"意味着条形的高度就是数据框中的数据值。position = position_dodge()意味着条形会被并排放置。
  coord_flip() +
  theme_minimal() +
  #scale_y_continuous(trans = "log10")+
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),axis.text = element_text(size = 14, face = "bold",family = 'arial')) + # 添加黑色的轴线
  labs(y = "Relative Abundance") # 改变x轴的标签


# 在柱状图上添加点，点有黑色的描边
#p <- p + geom_point(data = data_merge_new, aes(x = Taxonomy, y = Value),
                    #position = position_dodge(width=0.9), shape=21, colour="black")


# 添加p值和显著性标记
df_signif <- significant_TAX %>%
  select(Taxonomy, p_value, signif_mark) %>%
  rename(Tax = Taxonomy)

for(i in 1:nrow(df_signif)){ 
  p1 <- p1 + 
    annotate("text", 
             x = df_signif$Tax[i], 
             y = max(data_merge_new$Value)+1, 
             label = paste0("p=", df_signif$p_value[i], " ", df_signif$signif_mark[i]),
             hjust = 1) 
}

# 显示图形
print(p1)
#ggsave("属显著.pdf", device = cairo_pdf, width = 10, height = 8)

#————————————————————————————————————————————————————————————--------
#以下是按照差异真菌的Foldchange进行排序、做图

df_summary <- data_merge_new_summary

# 提取年轻组和年老组的平均表达量
young_expression <- df_summary %>% filter(Group == "young")
old_expression <- df_summary %>% filter(Group == "aged")

# 合并年轻组和年老组的表达量
merged_expression <- merge(young_expression, old_expression, by = "Taxonomy", suffixes = c("_young", "_old"))

# 计算分组间的差异倍数
merged_expression$fold_change <- log2(merged_expression$mean_abundance_young / merged_expression$mean_abundance_old)

# 根据差异倍数的正负值进行分组
merged_expression$Group <- ifelse(merged_expression$fold_change>0, "Young enriched", "Aged enriched")
table(merged_expression$fold_change>0)

# 根据差异倍数的正负值重新排序
merged_expression <- merged_expression[order(merged_expression$fold_change, decreasing = TRUE), ]

# 绘制柱状图
library(ggplot2)

p2 <- ggplot(merged_expression, aes(x = reorder(Taxonomy, fold_change), y = fold_change, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5 )) + # 添加黑色的轴线
  labs(x = "",y = "log2(Fold Change)")
p2
#ggsave(p2,file="Fungi-age-signif.png")
ggsave(p2,file="Fungi-age-signif-2.pdf",height = 11,width = 8)
###########################################################
#做差异真菌与表型的相关性
library(psych)
library(ggplot2)
library(reshape2)
#读取物种与环境因子数据
# 假设您有一个名为 microbiome_data 的数据框，其中包含微生物丰度数据，以及一个名为 phenotype_data 的数据框，其中包含表型数据
df_sig <-ch_fu[,rev(merged_expression$Taxonomy)]
dim(df_sig)


load("pheno_de.Rdata")
pheno <- pheno_p[,4:17]
#a=c("pores_diameter")
pheno <- pheno[,-12]
colnames(pheno) <- c("porphyrin","age","lentigines","telangiectasia","elasticity","TEWL","hydration","pH","L*","a*","b*","pores_area","sebum")
colnames(pheno)
pheno <- pheno[,c(2,13,1,7,6,8,5,9,10,11,3,4,12)]
rownames(df_sig)==rownames(pheno)


# 计算微生物丰度与表型之间的Spearman相关系数
r_result <- corr.test(df_sig,pheno,method = "spearman", adjust = "bonferron") 

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

p3 <- ggplot(data,aes(x=Phenotype, y=Microbe))+
  geom_tile(aes(fill=r),color = 'white',alpha = 0.8) +
  geom_text(aes(label = P_value_sig), color = 'black', size = 3,hjust = 0.5, vjust = 0.5) +
  ggsci::scale_fill_gsea() +
  theme +
  xlab('') +
  ylab('')+
  #coord_equal()+
  scale_fill_gradient2(low = "#2c7bb6", high = "#d7191c", mid = "white", midpoint = 0)+
  labs(fill = "Spearman correlation")

p3
#ggsave(p3,file="fungi-pheno-cor.png")

ggsave(p3,file="fungi-pheno-cor-3.pdf",height = 11,width = 8)

