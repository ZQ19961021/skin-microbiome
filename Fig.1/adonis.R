library(vegan)
library(ggplot2)
library(reshape2)

# 加载数据
load("china_micr.RData")
load("pheno_v.Rdata")

# 处理性别
gender <- ifelse(pheno$gender == "Male", 1, 2)
pheno$gender <- as.factor(gender)

# 选择表型数据
phenotype <- pheno[, c(1, 3:15)]

# 去除缺失值的行
complete_cases <- complete.cases(phenotype)
phenotype_complete <- phenotype[complete_cases, ]
china_mic_complete <- china_mic[rownames(phenotype_complete), ]

# 确保数据框行名匹配
stopifnot(all(rownames(china_mic_complete) %in% rownames(phenotype_complete)))
stopifnot(all(rownames(phenotype_complete) %in% rownames(china_mic_complete)))

# 计算细菌、真菌和病毒的 adonis2 结果
results <- list()

for (microbe_type in c("Bacteria", "Eukaryota", "Viruses")) {
  # 提取对应类别的数据
  microbe_data <- china_mic_complete[, grepl(microbe_type, colnames(china_mic_complete))]
  
  # 执行 adonis2 分析
  adonis_result <- adonis2(microbe_data ~ ., data = phenotype_complete)
  
  # 打印 adonis2 结果以检查输出
  print(adonis_result)
  
  # 确认 adonis_result 的结构
  aov_tab <- adonis_result$aov.tab
  print(aov_tab)
  
  # 检查 aov_tab 的行名和列名
  print(rownames(aov_tab))
  print(colnames(aov_tab))
  
  # 提取 R2 值
  # 确保行名和列名与代码匹配
  if ("Total" %in% rownames(aov_tab) && "R2" %in% colnames(aov_tab)) {
    r2_value <- aov_tab["Total", "R2"]
  } else {
    r2_value <- NA
  }
  
  # 保存结果
  results[[microbe_type]] <- r2_value
}


# 绘制堆积图
results <- read.csv("R2.csv")
results$microbe_type <- factor(results$microbe_type, levels = c("Viruses", "Eukaryota", "Bacteria"))
ggplot(results, aes(x = phenotype, y = r2, fill =microbe_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(y = "Adonis R2") +
  theme_bw() +
  scale_fill_manual(values = c("Bacteria" = "#68ac56", "Eukaryota" = "#497eb3", "Viruses" = "#d2352a"))




