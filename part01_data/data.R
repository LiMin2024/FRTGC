setwd("/share/data1/limin/vagina/mmseq/analysis/00.data/curve/")
dt = read.table("roc.data.tsv", sep="\t", header=T)

dtf <- dt %>%
  group_by(nspecies, group) %>%
  summarise(value = mean(obs)) %>%
  data.frame()

dtf_filtered <- dtf %>%
  filter(group == "all" | (group == "nosingle"))

p1 <- ggplot(data = dtf_filtered, aes(x = nspecies / 1000000, y = value / 1000000, color = group)) + # 1. Y轴数值除以1000
  geom_line(size = 1) +
  #geom_point()+ 
  theme_bw() +
  theme(
    aspect.ratio = 0.6,
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.grid = element_blank() # 2. 去掉所有网格线
  ) +
  scale_x_continuous(breaks = seq(0, max(dtf_filtered$nspecies, na.rm = T), by = 100)) +
  xlab("Number of genes (million)") +
  ylab("Number of non-redundant genes (million)") 
p1

#ggsave("vagina.roc.pdf", p1, width=5, height=5)

###############pip
df <- data.frame(
  group = c("Bacteria", "Eukaryota","Virus"),
  value = c(4760147, 207921, 302945)
)

# 3. 数据预处理：计算百分比和标签位置
# 计算百分比（保留两位小数）
df$percentage <- round(df$value / sum(df$value) * 100, 2)
# 计算标签在饼图中的垂直位置（保证标签显示在每个扇形的中间）
df$ypos <- cumsum(df$percentage) - 0.5 * df$percentage

# 4. 绘制饼图
p <- ggplot(df, aes(x = "", y = percentage, fill = group)) +
  # 绘制堆叠柱状图，width=1表示柱子之间没有间隙，color="white"给扇形之间加白色边框
  geom_bar(stat = "identity", width = 1, color = "white") +
  # 核心步骤：将直角坐标系转换为极坐标系，生成饼图
  coord_polar("y", start = 0) +
  # 添加百分比标签，并设定标签颜色为白色，大小为4
  geom_text(aes(y = ypos, label = paste0(percentage, "%")), color = "white", size = 4) +
  # 使用空白主题，去掉背景、网格和坐标轴，让饼图更干净
  theme_void() +
  # 自定义填充颜色（你可以换成自己喜欢的颜色）
  scale_fill_brewer(palette = "Set2")

# 5. 展示图片
print(p)

###########map.rate
setwd("/share/data1/limin/vagina/mmseq/analysis/00.data/map/")
data <- read.table("/share/data1/limin/vagina/metagenome/metadata_new",header=TRUE,sep="\t")
names(data)[1] <- "Sample_ID"
data <- data[!duplicated(data$Sample_ID),]

map <- read.table("map.rate",header=FALSE,sep="\t")
names(map) <- c("Sample_ID","rate")
map[, 1] <- sub("\\..*", "", map[, 1])
map <- map[!duplicated(map$Sample_ID),]
map2 <- map[map$Sample_ID %in% data$Sample_ID,]

not_in_data <- map$Sample_ID[!map$Sample_ID %in% data$Sample_ID]
not_in_map <- data$Sample_ID[!data$Sample_ID %in% map$Sample_ID]

data2 <- merge(map2,data,by="Sample_ID")

#########
table(data2$Country)
data2_top10 <- data2 %>%
  group_by(Country) %>%
  slice_sample(n = 10)  %>%      # 提取前10行
  ungroup() 
#write.table(data2_top10,"select.sample",sep = "\t", row.names = FALSE, quote = FALSE)

#########
data_clean <- data2[!is.na(data2$rate), ]
data_clean <- data_clean[,c(2,3,5)]
data_clean$rate_num <- as.numeric(gsub("%", "", data_clean$rate))
names(data_clean)
library(patchwork)

p_all <- ggplot(data_clean, aes(x = "All Samples", y = rate_num)) +
  geom_boxplot(
    alpha = 0.7, 
    fill = "lightblue", 
    outlier.shape = NA
  ) +
  labs(title = "", x = "", y = "Map Rate (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    aspect.ratio = 8
  )

# 3. 绘制第二张图：按 site 分组
p_site <- ggplot(data_clean, aes(x = Body_Site_Category, y = rate_num)) +  # 去掉 aes 里的 fill = site
  geom_boxplot(
    alpha = 0.7, 
    fill = "lightblue", 
    outlier.shape = NA
  ) +
  labs(title = "", x = "", y = "") +
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    aspect.ratio = 2
    # 删掉 legend.position = "none"，因为已经没有图例了
  )

# 4. 绘制第三张图：按 country 分组
p_country <- ggplot(data_clean, aes(x = Country, y = rate_num)) + # 去掉 aes 里的 fill = country
  geom_boxplot(
    alpha = 0.7, 
    fill = "lightblue", 
    outlier.shape = NA
  ) +
  labs(title = "", x = "", y = "") +
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    aspect.ratio = 0.5
    # 删掉 legend.position = "none"，因为已经没有图例了
  )

final_plot <- p_all + p_site + p_country
print(final_plot)

################
data <- read.table("/share/data1/limin/vagina/metagenome/metadata_new",header=TRUE,sep="\t")
names(data)[1] <- "Sample_ID"
data <- data[!duplicated(data$Sample_ID),]
table(data$Body_Site_Category)
map <- read.table("map.rate",header=FALSE,sep="\t")
map <- read.table("/share/data2/limin/vagina/VIRGO2/map/02.profile/map.rate",header=FALSE,sep="\t")

names(map) <- c("Sample_ID","rate")
map[, 1] <- sub("\\..*", "", map[, 1])
map <- map[!duplicated(map$Sample_ID),]
map2 <- map[map$Sample_ID %in% data$Sample_ID,]

not_in_data <- map$Sample_ID[!map$Sample_ID %in% data$Sample_ID]
not_in_map <- data$Sample_ID[!data$Sample_ID %in% map$Sample_ID]

data2 <- merge(map2,data,by="Sample_ID")

data_clean <- data2[!is.na(data2$rate), ]
data_clean <- data_clean[,c(1:3,5)]
data_clean$rate_num <- as.numeric(gsub("%", "", data_clean$rate))

FVGC_data_clean <- data_clean
FVGC_data_clean$Source <- "FVGC"
VIRGO2_data_clean <- data_clean
VIRGO2_data_clean$Source <- "VIRGO2"

# 合并数据
combined_data <- rbind(FVGC_data_clean, VIRGO2_data_clean)

##############
mean_rates <- combined_data %>%
  group_by(Source) %>%       # 按 site 和 Source 联合分组
  summarise(
    Mean_Rate = median(rate_num, na.rm = TRUE),  # 计算平均值，忽略缺失值
    Count = n()                              # 顺便统计每组的样本量（可选）
  ) %>%
  ungroup()  

median_rates <- combined_data %>%
  group_by(Source) %>%       # 按 site 和 Source 联合分组
  summarise(
    Median_Rate = median(rate_num, na.rm = TRUE),  # 计算平均值，忽略缺失值
    Count = n()                              # 顺便统计每组的样本量（可选）
  ) %>%
  ungroup() 

median_site_table <- combined_data %>%
  group_by(Body_Site_Category,Source) %>%
  summarise(Median_Map_Rate = median(rate_num, na.rm = TRUE))

median_country_table <- combined_data %>%
  group_by(Country,Source) %>%
  summarise(Median_Map_Rate = median(rate_num, na.rm = TRUE))

df_wide <- median_country_table %>%
  pivot_wider(names_from = Source, values_from = Median_Map_Rate)
df_ratio <- df_wide %>%
  mutate(Ratio_FVGC_vs_VIRGO2 = FVGC / VIRGO2)

###############
p_all <- ggplot(combined_data, aes(x = "All Samples", y = rate_num, fill = Source)) +
  
  # 1. 箱线图：使用 position_dodge 并排错开，隐藏默认异常值
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA
  ) +
  # 3. 统一配色（与 p_site_combined 保持一致）
  scale_fill_manual(values = c("FVGC" = "lightblue", "VIRGO2" = "lightcoral")) +
  
  labs(title = "", x = "", y = "Map Rate (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 19, color = "black"),
    axis.title = element_text(size = 19),
    plot.title = element_text(size = 19, hjust = 0.5),
    aspect.ratio = 13,
    legend.position = "none"
  )

p_site_combined <- ggplot(combined_data, aes(x = Country, y = rate_num, fill = Source)) +  
  
  # 1. 绘制箱线图：使用 position_dodge() 使不同 Source 的箱子并排错开
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA  # 【建议】隐藏箱线图自带的异常值，避免和下面的散点重叠
  ) +  
  
  # 3. 自定义颜色
  scale_fill_manual(values = c("FVGC" = "lightblue", "VIRGO2" = "lightcoral")) + 
  
  # 5. 主题美化
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(), 
    aspect.ratio = 0.6,
    legend.position = "none"
  )

p_country_combined <- ggplot(combined_data, aes(x = Body_Site_Category, y = rate_num, fill = Source)) +  
  
  # 1. 绘制箱线图：使用 position_dodge() 使不同 Source 的箱子并排错开
  geom_boxplot(
    alpha = 0.7, 
    color = "black", 
    position = position_dodge(width = 0.8),
    outlier.shape = NA  # 【建议】隐藏箱线图自带的异常值，避免和下面的散点重叠
  ) +  
  
  # 3. 自定义颜色
  scale_fill_manual(values = c("FVGC" = "lightblue", "VIRGO2" = "lightcoral")) + 
  
  # 5. 主题美化
  theme_bw() +                         
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.x = element_text(size = 19, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(), 
    aspect.ratio = 3,
    legend.position = "top"
  )

final_plot <- p_all + p_country_combined + p_site_combined 
print(final_plot)
