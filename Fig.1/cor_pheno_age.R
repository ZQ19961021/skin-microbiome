# 加载必要的库
library(dplyr)
library(reshape2)
library(ggplot2)
library(pheatmap)

# 加载数据
load("pheno_v.Rdata")

# 创建一个函数来计算相关性和显著性
calculate_cor_pval <- function(df, var) {
  valid_data <- df %>% filter(!is.na(age) & !is.na(.data[[var]]))
  if (nrow(valid_data) > 2) {
    cor_test <- cor.test(valid_data$age, valid_data[[var]], method = "spearman")
    return(list(cor = cor_test$estimate, pval = cor_test$p.value))
  } else {
    return(list(cor = NA, pval = NA))
  }
}

# 按部位分组计算相关性和显著性
cor_results <- pheno %>%
  group_by(site) %>%
  summarise(across(c(sebum, porphyrin, hydration, TEWL, pH, elasticity, `L*`, `a*`, `b*`, lentigines, telangiectasia, pores_area), 
                   ~ calculate_cor_pval(cur_data(), cur_column()), .names = "{col}"))

# 分离相关性和p值
cor_matrix <- cor_results %>%
  summarise(across(everything(), ~ .x$cor)) %>%
  pivot_longer(cols = -site, names_to = "variable", values_to = "correlation")

pval_matrix <- cor_results %>%
  summarise(across(everything(), ~ .x$pval)) %>%
  pivot_longer(cols = -site, names_to = "variable", values_to = "p_value")

# 合并相关性和p值
cor_pval_matrix <- left_join(cor_matrix, pval_matrix, by = c("site", "variable"))

# 转换为宽格式
cor_pval_wide <- cor_pval_matrix %>%
  pivot_wider(names_from = variable, values_from = c(correlation, p_value))

# 提取相关性和显著性矩阵
cor_matrix <- cor_pval_wide %>% select(starts_with("correlation_"))
pval_matrix <- cor_pval_wide %>% select(starts_with("p_value_"))

# 转换行名
row.names(cor_matrix) <- cor_pval_wide$site
row.names(pval_matrix) <- cor_pval_wide$site

# 去除列名前缀
colnames(cor_matrix) <- gsub("correlation_", "", colnames(cor_matrix))
colnames(pval_matrix) <- gsub("p_value_", "", colnames(pval_matrix))


# 自定义颜色条
heatmap_colors <- colorRampPalette(c("#6c89b0", "white", "#c9716b"))(50)

# p值显著性标记
signif_mark <- function(pvals) {
  ifelse(pvals < 0.001, "***", 
         ifelse(pvals < 0.01, "**", 
                ifelse(pvals < 0.05, "*", "")))
}

# 生成热图并将图形旋转为竖向
pheatmap(t(cor_matrix),
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         display_numbers = signif_mark(t(pval_matrix)),
         color = heatmap_colors,
         cellwidth = 20,
         cellheight = 20)
