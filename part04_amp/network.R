setwd("/share/data1/limin/vagina/mmseq/analysis/05.amp/")
load("../03.taxo/taxonomy.known.norm.RData")
taxo <- dt_norm

##############
result <- align_dt_sample(taxo, sample_map1, ID = "Sample")
aligned_dt <- result$dt
aligned_dt1 <- aligned_dt

result <- align_dt_sample(taxo, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- aligned_dt

dt2 <- rownames_to_column(aligned_dt1,var="name")

high_abundance_features <- dt2 %>%
  # 将除了第1列（物种名）之外的所有列转换为数值型（防止文本干扰）
  mutate(across(-1, as.numeric)) %>%
  # 计算行平均值
  mutate(Mean_Abundance = rowMeans(across(-1), na.rm = TRUE)) %>%
  # 只保留物种名和平均丰度
  select(1, Mean_Abundance)

high_abundance_features1 <- high_abundance_features
#save(high_abundance_features1,file="high_abundance_features1")
high_abundance_features2 <- high_abundance_features
#save(high_abundance_features2,file="high_abundance_features2")

#####################
load("amp.known.norm.RData")
amp_profile <- amp_dt#流行率大于10%
amp <- rownames_to_column(amp_profile,var="gene")

###########
sample_map1.1 <- sample_map1[sample_map1$Country == "China",]
sample_map1.2 <- sample_map1[sample_map1$Country == "U.S.",]
result <- align_dt_sample(taxo, sample_map1.2, ID = "Sample")
aligned_dt <- result$dt
taxo_dt1 <- aligned_dt
#write.table(taxo_dt1,file = "network/taxo_China", sep = "\t", row.names = T, quote = FALSE)
#write.table(taxo_dt1,file = "network/taxo_US", sep = "\t", row.names = T, quote = FALSE)

result <- align_dt_sample(amp_profile, sample_map1.2, ID = "Sample")
aligned_dt <- result$dt
amp_dt1 <- aligned_dt
amp_dt1 <- amp_dt1[, colnames(taxo_dt1) ]  # 对齐样本顺序
#write.table(amp_dt1,file = "network/amp_China", sep = "\t", row.names = T, quote = FALSE)
#write.table(amp_dt1,file = "network/amp_US", sep = "\t", row.names = T, quote = FALSE)

###############
unique(sample_map2$Body_Site_Category)
sample_map2.1 <- sample_map2[sample_map2$Body_Site_Category == "Cervix",]
sample_map2.2 <- sample_map2[sample_map2$Body_Site_Category == "Upper vaginal segment",]
sample_map2.3 <- sample_map2[sample_map2$Body_Site_Category == "Mid vagina",]
sample_map2.4 <- sample_map2[sample_map2$Body_Site_Category == "Lower segment of the vagina",]

result <- align_dt_sample(taxo, sample_map2.3, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map
taxo_dt2 <- aligned_dt

#write.table(taxo_dt2,file = "network/taxo_Cervix", sep = "\t", row.names = T, quote = FALSE)
#write.table(taxo_dt2,file = "network/taxo_UP", sep = "\t", row.names = T, quote = FALSE)
#write.table(taxo_dt2,file = "network/taxo_MID", sep = "\t", row.names = T, quote = FALSE)
#write.table(taxo_dt2,file = "network/taxo_LOW", sep = "\t", row.names = T, quote = FALSE)

result <- align_dt_sample(amp_profile, sample_map2.3, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map
amp_dt2 <- aligned_dt
amp_dt2 <- amp_dt2[, colnames(taxo_dt2) ]  # 对齐样本顺序

#write.table(amp_dt2,file = "network/amp_Cervix", sep = "\t", row.names = T, quote = FALSE)
#write.table(amp_dt2,file = "network/amp_UP", sep = "\t", row.names = T, quote = FALSE)
#write.table(amp_dt2,file = "network/amp_MID", sep = "\t", row.names = T, quote = FALSE)
#write.table(amp_dt2,file = "network/amp_LOW", sep = "\t", row.names = T, quote = FALSE)

#############用py脚本计算spearman
setwd("/share/data1/limin/vagina/mmseq/analysis/05.amp/network/")
rcorr_result = read.table("China_amp_taxo_corr.txt",fill=T, header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

China_result <- threshold
China_result$Site <- "China"
China_result <- China_result %>% 
  mutate(group = paste(A, B, sep = "|"))
China_result <- China_result[China_result$r.label == "negative",]

rcorr_result = read.table("US_amp_taxo_corr.txt",fill=T, header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

USA_result <- threshold
USA_result$Site <- "USA"
USA_result <- USA_result %>% 
  mutate(group = paste(A, B, sep = "|"))
USA_result <- USA_result[USA_result$r.label == "negative",]

rcorr_result = read.table("Cervix_amp_taxo_corr.txt",fill=T, header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

Cervix_result <- threshold
Cervix_result$Site <- "Cervix"
Cervix_result <- Cervix_result %>% 
  mutate(group = paste(A, B, sep = "|"))
Cervix_result <- Cervix_result[Cervix_result$r.label == "negative",]

rcorr_result = read.table("UP_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

UP_result <- threshold
UP_result$Site <- "UP"
UP_result <- UP_result %>% 
  mutate(group = paste(A, B, sep = "|"))
UP_result <- UP_result[UP_result$r.label == "negative",]

rcorr_result = read.table("MID_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

MID_result <- threshold
MID_result$Site <- "MID"
MID_result <- MID_result %>% 
  mutate(group = paste(A, B, sep = "|"))
MID_result <- MID_result[MID_result$r.label == "negative",]

rcorr_result = read.table("LOW_amp_taxo_corr.txt", header = T,  sep = "\t")
rcorr_result$qvalue <- p.adjust(rcorr_result$pvalue, method = "fdr")
# 添加 r.label 和 r.value
rcorr_result$r.label  <- ifelse(rcorr_result$corr > 0, "positive", "negative")
rcorr_result$r.value  <- abs(rcorr_result$corr)  # 取绝对值

threshold <- rcorr_result %>%
  filter( qvalue < 0.05 & r.value > 0.4) 

LOW_result <- threshold
LOW_result$Site <- "LOW"
LOW_result <- LOW_result %>% 
  mutate(group = paste(A, B, sep = "|"))
LOW_result <- LOW_result[LOW_result$r.label == "negative",]

library(purrr) 
library(tidyverse) 
combind_dt1 <- list(China_result[,9:8], USA_result[,9:8]) %>%
  reduce(full_join, by = "group")
names(combind_dt1) <- c("group","China","U.S.")
#write.table(combind_dt1,file = "cervix_negative_amp", sep = "\t", row.names = F, quote = FALSE)

combind_dt2 <- list(Cervix_result[,9:8],UP_result[,9:8], MID_result[,9:8],LOW_result[,9:8]) %>%
  reduce(full_join, by = "group")
names(combind_dt2) <- c("group","Cervix","UP","MID","LOW")
#write.table(combind_dt2,file = "us_negative_amp", sep = "\t", row.names = F, quote = FALSE)

combind_dt1 <- myread("network/cervix_negative_amp")
combind_dt1 <- rownames_to_column(combind_dt1,var="group")
combind_dt2 <- myread("network/us_negative_amp")
combind_dt2 <- rownames_to_column(combind_dt2,var="group")

################
#cols <- c( "China", "U.S.")
#cols <- c( "Cervix","UP","MID","LOW")

result1 <- combind_dt1 %>% 
  separate(col = group, into = c("A", "B"), sep = "\\|", remove = TRUE)

high_abundance_features3 <- high_abundance_features1 %>%
  filter(Mean_Abundance > 0.01) 
result1 <- result1[result1$B %in% high_abundance_features3$name,]

result2 <- combind_dt2 %>% 
  separate(col = group, into = c("A", "B"), sep = "\\|", remove = TRUE)

high_abundance_features3 <- high_abundance_features2 %>%
  filter(Mean_Abundance > 0.01) 
result2 <- result2[result2$B %in% high_abundance_features3$name,]

###############
result_filtered <- result2 %>%
  filter(
    !grepl("[0-9]", B),            
    !grepl("uncultured", B, ignore.case = TRUE),
    !grepl("sp\\.", B, ignore.case = TRUE) 
  )
length(unique(result_filtered$A))
length(unique(result_filtered$B))

df_long <- result_filtered %>%
  pivot_longer(
    cols = starts_with(c("China", "U.S.")), #"China", "U.S."  "Cervix","UP","MID","LOW"
    names_to = "sample_type",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>% # 移除没有样本类型的记录
  select(-value)            # 删除不必要的列

base_colors <- c(
  "China" = "#e31a1c", 
  "U.S." = "#fdbf6f", 
  "China+U.S." = "#ff7f00"
)

base_colors <- c(
  "Cervix" = "#a6cee3", 
  "UP" = "#b2df8a", 
  "MID" = "#fb9a99", 
  "LOW" = "#cab2d6",
  
  # --- Cervix 同色系 (由浅到深) ---
  "Cervix+UP" = "#82b1d4",
  "Cervix+MID" = "#5a96c5",
  "Cervix+LOW" = "#337bb6",
  "Cervix+MID+UP" = "#1f63a0",
  "Cervix+LOW+UP" = "#0a4a8a",
  "Cervix+LOW+MID+UP" = "#003366", # 最深色
  
  # --- LOW 同色系 ---
  "LOW+UP" = "#a682b8", 
  
  # --- MID 同色系 ---
  "MID+UP" = "#e06b6b" 
)

# 2. 数据汇总与组合标签生成
df_summary <- df_long %>%
  group_by(A, B) %>%
  summarise(
    types = list(sort(unique(sample_type))), 
    .groups = 'drop'
  ) %>%
  mutate(
    n_types = lengths(types),
    combo_label = sapply(types, paste, collapse = "+")
  )

# 3. 动态处理可能出现的未知组合 (以防万一)
all_combos <- unique(df_summary$combo_label)
new_combos <- setdiff(all_combos, names(base_colors))
if (length(new_combos) > 0) {
  auto_colors <- scales::hue_pal()(length(new_combos))
  names(auto_colors) <- new_combos
  final_color_mapping <- c(base_colors, auto_colors)
} else {
  final_color_mapping <- base_colors
}

y_order <- df_summary %>%
  count(B, sort = TRUE) %>%   # 按 B 分组计数，并按频次降序排列
  pull(B)                     # 提取排好序的 B 的向量

df_summary$B <- factor(df_summary$B, levels = rev(y_order))

ggplot(df_summary, aes(x = A, y = B)) +
  # 将 fill 映射到自动识别的 combo_label 上
  geom_tile(aes(fill = combo_label), color = "white") +
  
  # 使用动态生成的颜色字典
  scale_fill_manual(
    values = final_color_mapping,
    name = "Type Combination"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    aspect.ratio = 2,
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    legend.position = "right"
  ) +
  labs(x = "AMP", y = "")

