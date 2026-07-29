sample_map <- read.table("/share/data1/limin/vagina/metagenome/metadata_new",header=TRUE,sep="\t")
names(sample_map)[1] <- "Sample"

########
sample_map <- sample_map[sample_map$Body_Site_Category != "Others",]

table(sample_map$Country,sample_map$Body_Site_Category)
sample_map1 <- sample_map[sample_map$Body_Site_Category == "Cervix",]
sample_map1 <- sample_map1[sample_map1$Country == "china"|sample_map1$Country == "usa",]

sample_map2 <- sample_map[sample_map$Country == "usa",]

###########
result <- align_dt_sample(dt_norm, sample_map1, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- t(aligned_dt)
#write.table(aligned_dt2,file = "cazy.cervix.dt", sep = "\t", row.names = T, quote = FALSE)

aligned_sample <- result$sample_map
aligned_sample <- aligned_sample %>%
  mutate(Country = recode(Country,
                          "china" = "China",
                          "usa"   = "U.S."))
#write.table(aligned_sample,"cervix.sample",sep = "\t", row.names = FALSE, quote = FALSE)

result <- align_dt_sample(dt_norm, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- t(aligned_dt)
#write.table(aligned_dt2,file = "cazy.us.dt", sep = "\t", row.names = T, quote = FALSE)

aligned_sample <- result$sample_map
#write.table(aligned_sample,"usa.sample",sep = "\t", row.names = FALSE, quote = FALSE)

###############
data <- myread("../01.abun/kegg.profiles.known.norm")

sample_map1 <- read.table("cervix.sample",header=TRUE,sep="\t")
sample_map2 <- read.table("usa.sample",header=TRUE,sep="\t")

###############
source("/share/data1/limin/Script/R_function/alpha_calculate.R")
alpha_result <- calculate_alpha_by_sample(dt_norm, sample_map2)
#write.table(alpha_result,"cazy.us.alpha_result",sep = "\t", row.names = FALSE, quote = FALSE)

alpha_result1 <- read.table("cazy.cervix.alpha_result",header=TRUE,sep="\t")
alpha_result2 <- read.table("cazy.us.alpha_result",header=TRUE,sep="\t")

################
data2 <- alpha_result1
my_colors1 <- c("#e31a1c", "#fdbf6f")
data2$Country <- factor(data2$Country, 
                        levels = c("China", "U.S."))

data3 <- alpha_result2
my_colors2 <- c("#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6")  
data3$Body_Site_Category <- factor(data3$Body_Site_Category, 
                                   levels = c("Cervix", "Upper vaginal segment", "Mid vagina","Lower segment of the vagina"))
###############
# ================== 绘制图1 (Obs Diversity Index) ==================
# 1. 计算 Obs 的两两比较结果
stat_results_obs <- data2 %>%
  wilcox_test(obs ~ Country) %>%      
  add_significance("p") %>%           
  add_xy_position(x = "Country")      

# 2. 绘制图1
p1 <- ggplot(data2, aes(x = Country, y = obs)) + 
  geom_violin(aes(fill = Country), trim = FALSE, alpha = 0.8, size = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA, fill = "white", size = 0.5) + 
  
  # 【注意】这里通过 '+' 正确连接下一个图层
  stat_pvalue_manual(
    stat_results_obs,               
    label = "p.signif",         
    tip_length = 0.01,          
    hide.ns = TRUE              
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 19),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    aspect.ratio = 2
  )  +
  labs(title = "", x = "", y = "Obs Diversity Index") +
  scale_fill_manual(values = my_colors1) + 
  scale_color_manual(values = my_colors1) + 
  guides(color = FALSE, fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(ylim = c(0, max(stat_results_obs$y.position) * 1)) 

print(p1)


# ================== 绘制图2 (Shannon Diversity Index) ==================
# 1. 计算 Shannon 的两两比较结果
stat_results_shannon <- data2 %>%
  wilcox_test(shannon ~ Country) %>%      
  add_significance("p") %>%           
  add_xy_position(x = "Country")      

# 2. 绘制图2
p2 <- ggplot(data2, aes(x = Country, y = shannon)) + 
  geom_violin(aes(fill = Country), trim = FALSE, alpha = 0.8, size = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA, fill = "white", size = 0.5) + 
  
  # 【注意】传入对应的 shannon 统计结果表
  stat_pvalue_manual(
    stat_results_shannon,               
    label = "p.signif",         
    tip_length = 0.01,          
    hide.ns = TRUE              
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 19),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    aspect.ratio = 2
  )  +
  labs(title = "", x = "", y = "Shannon Diversity Index") +
  scale_fill_manual(values = my_colors1) + 
  scale_color_manual(values = my_colors1) + 
  guides(color = FALSE, fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(ylim = c(2, max(stat_results_shannon$y.position) * 1)) 

print(p2)

# ================== 绘制图3 (Body_Site_Category - Obs Diversity Index) ==================
# 1. 计算 Body_Site_Category 分组下 obs 的两两比较结果
stat_results_site_obs <- data3 %>%
  wilcox_test(obs ~ Body_Site_Category) %>%      
  add_significance("p") %>%           
  add_xy_position(x = "Body_Site_Category")      

# 2. 绘制图3
p3 <- ggplot(data3, aes(x = Body_Site_Category, y = obs)) + 
  geom_violin(aes(fill = Body_Site_Category), trim = FALSE, alpha = 0.8, size = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA, fill = "white", size = 0.5) + 
  
  # 添加显著性标注
  stat_pvalue_manual(
    stat_results_site_obs,               
    label = "p.signif",         
    tip_length = 0.01,          
    hide.ns = TRUE              
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 19),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    aspect.ratio = 2
  )  +
  labs(title = "", x = "", y = "Obs Diversity Index") +
  scale_fill_manual(values = my_colors2) + 
  scale_color_manual(values = my_colors2) + 
  guides(color = FALSE, fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(ylim = c(0, max(stat_results_site_obs$y.position) * 1)) 

print(p3)

# ================== 绘制图4 (Site - Shannon Diversity Index) ==================
# 1. 计算 Site 分组下 shannon 的两两比较结果
stat_results_site_shannon <- data3 %>%
  wilcox_test(shannon ~ Body_Site_Category) %>%      
  add_significance("p") %>%           
  add_xy_position(x = "Body_Site_Category")      

# 2. 绘制图4
p4 <- ggplot(data3, aes(x = Body_Site_Category, y = shannon)) + 
  geom_violin(aes(fill = Body_Site_Category), trim = FALSE, alpha = 0.8, size = 0.5) +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA, fill = "white", size = 0.5) + 
  
  # 添加显著性标注
  stat_pvalue_manual(
    stat_results_site_shannon,               
    label = "p.signif",         
    tip_length = 0.01,          
    hide.ns = TRUE              
  ) +
  
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 19),
    axis.text = element_text(size = 17),
    axis.title = element_text(size = 19),
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    legend.position = "none",
    aspect.ratio = 2
  ) +
  labs(title = "", x = "", y = "Shannon Diversity Index") +
  scale_fill_manual(values = my_colors2) + 
  scale_color_manual(values = my_colors2) + 
  guides(color = FALSE, fill = guide_legend(override.aes = list(alpha = 1))) +
  coord_cartesian(ylim = c(2, max(stat_results_site_shannon$y.position) * 1)) 

print(p4)

library(cowplot)
combined_plot <- plot_grid(p1, p2, p3, p4, ncol = 4, align = "hv", rel_widths = c(1, 1, 1, 1))
print(combined_plot)

##########################
source("/share/data1/limin/Script/R_function/beta_plot2.R")
dist_matrix <- read.table("cazy.cervix.dist.txt", header = TRUE, row.names = 1,fill = TRUE)
mydist <- as.dist(dist_matrix)
sample_map <- sample_map1
ado_result <- adonis2(mydist ~ sample_map$Country, data = sample_map, permutations = 999, parallel = 40)
#save(ado_result,file="cazy.cervix_ado.Rdata")
r2 <- ado_result$R2[1]
df_model <- ado_result$Df[1] 
n <- nrow(dist_matrix)
ado_r2 <- RsquareAdj(r2, m = df_model, n = n)
ado_p <- ado_result$`Pr(>F)`[1]
p_label <- ifelse(ado_p < 0.001, "< 0.001", round(ado_p, 3))
pcoa_plot_result <- zy_pcoa(
  dt = NULL,                  
  mydist = mydist,            
  sample_map = sample_map1,
  group = "Country",          
  ID = "Sample",           
  
  sample.color = NULL,        
  pca_method = "bray",        
  title = "Oropharynx",    # 建议加上标题
  
  x = 1, y = 2,               
  star_plot = FALSE,
  ellipse_plot = TRUE,
  levels = 0.95,              
  
  ados = FALSE,               # 关键：关闭内部 adonis 计算
  
  return_dist = FALSE
)

#save(pcoa_plot_result,file="cazy.cervix_PCoA.Rdata")

##############
#load("cervix_ado.Rdata")
r2 <- ado_result$R2[1]
df_model <- ado_result$Df[1] 
n <- nrow(dist_matrix)
ado_r2 <- RsquareAdj(r2, m = df_model, n = n)
ado_p <- ado_result$`Pr(>F)`[1]
p_label <- ifelse(ado_p < 0.001, "< 0.001", round(ado_p, 3))
#load("cervix_PCoA.Rdata")
p <- pcoa_plot_result$plot+
  theme(axis.text.x = element_text(size = 17), # 设置x轴文本大小
        axis.text.y = element_text(size = 17), # 设置y轴文本大小
        axis.title.x = element_text(size = 19), # 设置x轴标题大小
        axis.title.y = element_text(size = 19),
        aspect.ratio = 1 )

target_sites <- c("China", "U.S.")
named_colors <- c(
  "China" = "#e31a1c", 
  "U.S." = "#fdbf6f"
)

filtered_data <- p$data %>% 
  dplyr::filter(Country %in% target_sites)

# 3. 整合绘图与美化操作
p_final <- p %+% filtered_data +   
  geom_point(aes(color = Country), size = 1) +    
  scale_color_manual(values = named_colors) +
  scale_fill_manual(values = named_colors) +  
  labs(
    title = "Cervix",         # 设置主标题文字
    subtitle = bquote(R^2 == .(ado_r2) ~ "," ~ italic(p) == .(p_label))) +
  theme(
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    axis.title.x = element_text(size = 19), 
    axis.title.y = element_text(size = 19),
    aspect.ratio = 1
  )+
  scale_x_continuous(position = "top")

print(p_final)
p1 <- p_final

########
dm <- filtered_data 

# 更新因子水平顺序（根据实际需要调整这6个的顺序）
dm$Country <- factor(dm$Country, 
                     levels = c("China", "U.S."))

# 使用对应的 6 色配色方案
my_colors_site <- c(
  "China" = "#e31a1c", 
  "U.S." = "#fdbf6f"
)

# 生成两两比较的组合
groups <- levels(dm$Country)
comparisons <- combn(groups, 2, simplify = FALSE)

max_y2 <- max(dm$pcoa.2, na.rm = TRUE)
y_positions_side <- seq(max_y2 * 1.2, max_y2 * 2.5, length.out = length(comparisons))

# 动态计算下方水平箱线图 (PCoA1) 的位置
max_x1 <- abs(min(dm$pcoa.1, na.rm = TRUE))
y_positions_horizontal <- seq(max_x1 * 1.2, max_x1 * 2.5, length.out = length(comparisons))

# ====== 绘制右侧竖直箱线图 (PCoA2) ======
p_side <- ggplot(dm, aes(x = Country, y = pcoa.2, fill = Country)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) + 
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_side,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    aspect.ratio = 2
  ) +
  labs(y = "PCoA2")

# ====== 绘制下方水平箱线图 (PCoA1) ======
p_bottom <- ggplot(dm, aes(x = Country, y = pcoa.1, fill = Country)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) +
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_horizontal,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    aspect.ratio = 0.5
  ) +
  labs(x = "Site") +
  coord_flip()

final_plot <- (
  (p_final + p_side) /        
    (p_bottom + plot_spacer())     
) +
  plot_layout(
    widths = c(9, 1),   
    heights = c(9, 3)   
  ) +
  theme(plot.margin = margin(10, 10, 10, 10))

print(final_plot)

####################
dist_matrix <- read.table("cazy.us.dist.txt", header = TRUE, row.names = 1,fill = TRUE)
mydist <- as.dist(dist_matrix)
sample_map <- sample_map2
ado_result <- adonis2(mydist ~ sample_map$Body_Site_Category, data = sample_map, permutations = 999, parallel = 40)
#save(ado_result,file="cazy.us_ado.Rdata")
r2 <- ado_result$R2[1]
df_model <- ado_result$Df[1] 
n <- nrow(dist_matrix)
ado_r2 <- RsquareAdj(r2, m = df_model, n = n)
ado_p <- ado_result$`Pr(>F)`[1]
p_label <- ifelse(ado_p < 0.001, "< 0.001", round(ado_p, 3))
pcoa_plot_result <- zy_pcoa(
  dt = NULL,                  
  mydist = mydist,            
  sample_map = sample_map2,
  group = "Body_Site_Category",          
  ID = "Sample",           
  
  sample.color = NULL,        
  pca_method = "bray",        
  title = "Oropharynx",    # 建议加上标题
  
  x = 1, y = 2,               
  star_plot = FALSE,
  ellipse_plot = TRUE,
  levels = 0.95,              
  
  ados = FALSE,               # 关键：关闭内部 adonis 计算
  
  return_dist = FALSE
)

#save(pcoa_plot_result,file="cazy.us_PCoA.Rdata")

##############
#load("cervix_PCoA.Rdata")
p <- pcoa_plot_result$plot+
  theme(axis.text.x = element_text(size = 17), # 设置x轴文本大小
        axis.text.y = element_text(size = 17), # 设置y轴文本大小
        axis.title.x = element_text(size = 19), # 设置x轴标题大小
        axis.title.y = element_text(size = 19),
        aspect.ratio = 1 )

target_sites <- c("Cervix", "Upper vaginal segment", "Mid vagina","Lower segment of the vagina")
named_colors <- c(
  "Cervix" = "#a6cee3", 
  "Upper vaginal segment" = "#b2df8a",
  "Mid vagina" = "#fb9a99", 
  "Lower segment of the vagina" = "#cab2d6"
)

filtered_data <- p$data %>% 
  dplyr::filter(Body_Site_Category %in% target_sites)

# 3. 整合绘图与美化操作
p_final <- p %+% filtered_data +   
  geom_point(aes(color = Body_Site_Category), size = 1) +    
  scale_color_manual(values = named_colors) +
  scale_fill_manual(values = named_colors) +  
  labs(
    title = "U.S.",         # 设置主标题文字
    subtitle = bquote(R^2 == .(ado_r2) ~ "," ~ italic(p) == .(p_label))) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 17), 
    axis.text.y = element_text(size = 17), 
    axis.title.x = element_text(size = 19), 
    axis.title.y = element_text(size = 19),
    aspect.ratio = 1
  )+
  scale_x_continuous(position = "top")

print(p_final)
p2 <- p_final

library(patchwork)
combined_plot <- p1 + p2
combined_plot

###################
dm <- filtered_data 

# 更新因子水平顺序（根据实际需要调整这6个的顺序）
dm$Body_Site_Category <- factor(dm$Body_Site_Category, 
                                levels = c("Cervix", "Upper vaginal segment", "Mid vagina","Lower segment of the vagina"))

# 使用对应的 6 色配色方案
my_colors_site <- c(
  "Cervix" = "#a6cee3", 
  "Upper vaginal segment" = "#b2df8a",
  "Mid vagina" = "#fb9a99", 
  "Lower segment of the vagina" = "#cab2d6"
)

# 生成两两比较的组合
groups <- levels(dm$Body_Site_Category)
comparisons <- combn(groups, 2, simplify = FALSE)

max_y2 <- max(dm$pcoa.2, na.rm = TRUE)
y_positions_side <- seq(max_y2 * 1.2, max_y2 * 2.5, length.out = length(comparisons))

# 动态计算下方水平箱线图 (PCoA1) 的位置
max_x1 <- abs(min(dm$pcoa.1, na.rm = TRUE))
y_positions_horizontal <- seq(max_x1 * 1.2, max_x1 * 2.5, length.out = length(comparisons))

# ====== 绘制右侧竖直箱线图 (PCoA2) ======
p_side <- ggplot(dm, aes(x = Body_Site_Category, y = pcoa.2, fill = Body_Site_Category)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) + 
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_side,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    aspect.ratio = 2
  ) +
  labs(y = "PCoA2")

# ====== 绘制下方水平箱线图 (PCoA1) ======
p_bottom <- ggplot(dm, aes(x = Body_Site_Category, y = pcoa.1, fill = Body_Site_Category)) +
  geom_boxplot(width = 0.6, outlier.alpha = 0.6) +
  scale_fill_manual(values = my_colors_site) +
  geom_signif(
    comparisons = comparisons,
    y_position = y_positions_horizontal,
    map_signif_level = TRUE,
    textsize = 3, tip_length = 0.015, linewidth = 0.5, vjust = -0.2
  ) +
  theme_bw() +
  theme(
    legend.position = "none", panel.grid = element_blank(),
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    aspect.ratio = 0.5
  ) +
  labs(x = "Site") +
  coord_flip()

final_plot <- (
  (p_final + p_side) /        
    (p_bottom + plot_spacer())     
) +
  plot_layout(
    widths = c(9, 1),   
    heights = c(9, 3)   
  ) +
  theme(plot.margin = margin(10, 10, 10, 10))

print(final_plot)
