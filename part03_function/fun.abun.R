align_dt_sample <- function(dt, sample_map, ID = NA) {
  intersect_id <- intersect(sample_map[[ID]], colnames(dt))
  if (length(intersect_id) != nrow(sample_map)) {
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有", length(intersect_id), "个样本可以匹配\033[0m")
    sample_map <- sample_map[sample_map[[ID]] %in% intersect_id, , drop = FALSE]
  }
  dt <- dt[, sample_map[[ID]], drop = FALSE] %>% dplyr::filter(rowSums(.) != 0)
  list(dt = dt, sample_map = sample_map)
}
library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)      # 用于 deframe 等（也可用 tibble::deframe）
library(ggrepel)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

load("kegg.known.norm.RData")
kegg_dt <- dt_norm
load("cazy.known.norm.RData")
cazy_dt <- dt_norm

load("vfdb.known.norm.RData")
vfdb_dt <- dt_norm
load("rgi.known.norm.RData")
rgi_dt <- dt_norm

sample_map1 = read.table("cervix.sample", header = T,  sep = "\t")
sample_map2 = read.table("usa.sample", header = T,  sep = "\t")

##############
result <- align_dt_sample(kegg_dt, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

dt2 <- aligned_dt %>%
  rownames_to_column(var = "name")

##########rgi
map <- read.table("rgi_category", header = T, fill = TRUE, sep = "\t",check.names = F)
map2 <- map[,-1]
map3 <- map2 %>%
  separate_rows('Drug Class', sep = ";\\s*")
map4 <- map3 %>% 
  distinct(`Drug Class`, ARO, .keep_all = TRUE)

merge2 <- merge(dt2,map4,by.x = "name",by.y = "ARO")
merge3 <- merge2[,-1]
names(merge3)[ncol(merge3)] <- "name"

result <- merge3 %>%
  group_by(name) %>%   # 第一步：按照 name 列进行分组
  summarise(           # 第二步：对组内数据进行汇总计算
    across(where(is.numeric), sum, na.rm = TRUE) # 对所有数值型列求和，并自动忽略 NA 值
  )

###########vfdb
map <- read.table("VFDB_setB.latest.map", header = T, fill = TRUE, sep = "\t",check.names = F)
merge <- merge(dt2,map[,c(1,6)],by.x = "name",by.y = "GeneID")
merge <- merge[,-1]
names(merge)[ncol(merge)] <- "name"

result <- merge %>%
  group_by(name) %>%   # 第一步：按照 name 列进行分组
  summarise(           # 第二步：对组内数据进行汇总计算
    across(where(is.numeric), sum, na.rm = TRUE) # 对所有数值型列求和，并自动忽略 NA 值
  )

###########cazy
dt3 <- dt2 %>%
  mutate(name = str_extract(name, "^[A-Za-z]+")) 

result <- dt2 %>%
  group_by(name) %>%   # 第一步：按照 name 列进行分组
  summarise(           # 第二步：对组内数据进行汇总计算
    across(where(is.numeric), sum, na.rm = TRUE) # 对所有数值型列求和，并自动忽略 NA 值
  )

###########kegg
map <- read.table("KO_level_A_B_C_D_Description", header = F, fill = TRUE, sep = "\t",check.names = F)
map <- map[!grepl("Brite", map$V2), ]

merge1 <- merge(dt2,map[,c(2,7)],by.x = "name",by.y = "V7")
merge2 <- merge(dt2,map[,c(4,7)],by.x = "name",by.y = "V7")
merge3 <- merge(dt2,map[,c(6,7)],by.x = "name",by.y = "V7")

merge <- merge3[,-1]
names(merge)[ncol(merge)] <- "name"

result <- merge %>%
  group_by(name) %>%   # 第一步：按照 name 列进行分组
  summarise(           # 第二步：对组内数据进行汇总计算
    across(where(is.numeric), sum, na.rm = TRUE) # 对所有数值型列求和，并自动忽略 NA 值
  )

##################
result2 <- result %>%
  mutate(
    # 创建新列 row_mean，计算除 name 列外所有数值的行平均
    row_mean = rowMeans(select(., -name), na.rm = TRUE) 
  )

dt_long <- result %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

# 合并 sample_map 获取每个样本对应的 Group
dt_grouped <- dt_long %>%
  left_join(aligned_sample, by = "Sample")

mean_abundance <- dt_grouped %>%
  group_by(name, Body_Site_Category) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

# 将结果转为宽格式，便于比较两组
mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Body_Site_Category, values_from = MeanAbundance) %>%
  replace(is.na(.), 0)  # 缺失值补零（如果某 fungi 在某组中无样本）
kegg3_mean_abundance_wide <- mean_abundance_wide %>% 
  mutate(name = gsub("\\s*\\[.*?\\]", "", name))
#write.table(mean_abundance_wide,file = "kegg.cervix_abun_sort2.txt", sep = "\t", row.names = T, quote = FALSE)
#write.table(kegg3_mean_abundance_wide,file = "kegg.us_abun_sort3.txt", sep = "\t", row.names = T, quote = FALSE)

mean_abundance_wide1 <- read.table("kegg.cervix_abun_sort3.txt", header = T, fill = TRUE, sep = "\t",check.names = F)
result1 <- mean_abundance_wide1 %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE))) %>%
  mutate(row_idx = row_number()) %>%
  group_by(group = ifelse(row_idx <= 20, name, "Others")) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  rename(name = group) %>%
  select(-row_idx) %>%
  select(name, everything()) %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))                      
data1 <- result1 %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))

mean_abundance_wide2 <- read.table("kegg.us_abun_sort3.txt", header = T, fill = TRUE, sep = "\t",check.names = F)
result2 <- mean_abundance_wide2 %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE))) %>%
  mutate(row_idx = row_number()) %>%
  group_by(group = ifelse(row_idx <= 20, name, "Others")) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE), .groups = "drop") %>%
  rename(name = group) %>%
  select(-row_idx) %>%
  select(name, everything()) %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))    
data2 <- result2 %>%
  arrange(desc(rowMeans(select(., -name), na.rm = TRUE)))

data <- data1 %>%
  full_join(data2, by = "name")   

###########
df_long <- data %>%
  pivot_longer(cols = -name, names_to = "Site", values_to = "Abundance") %>%
  arrange(factor(name, levels = unique(result$name))) 

df_long2 <- df_long %>%
  mutate(Abundance100 = Abundance)

original_levels <- unique(data$name)
x_levels <- c(original_levels[original_levels != "Others"], "Others")

y_levels <- c("name","China","U.S.","Cervix","Lower segment of the vagina","Mid vagina","Upper vaginal segment")

library(scales)
df_plot <- df_long2 %>%
  mutate(
    Abundance_raw = Abundance100,  # 保留原始值用于标注
    Abundance_log_neg = Abundance  # 用于填色
  )

p <- ggplot(df_plot, aes(x = factor(name, levels = x_levels), y = factor(Site, levels = rev(y_levels)), fill = Abundance_log_neg)) +
  geom_tile(color = "white") +
  
  # 修改此处：使用 ifelse 判断，如果是 NA (原数据为0) 则不显示文字，否则保留两位小数
  geom_text(aes(label = ifelse(is.na(Abundance_log_neg), "", sprintf("%.2f", Abundance_raw))), 
            size = 3.5, 
            color = "black", 
            na.rm = TRUE) + 
  
  scale_fill_gradientn(
    colours = c("#0077B5","#E6F5FF","#FFE6E6","#DC143C"),
    values = rescale(range(df_plot$Abundance_log_neg, na.rm = TRUE)),
    name = "-Log(Relative Abundance)",
    na.value = "grey90"  # 0值依然填充为浅灰色
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  coord_fixed()+
  labs(title = "", x = "Function", y = "Country") 

print(p)
