setwd("/share/data1/limin/vagina/Geneset/analysis/02.kegg/")
myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

align_dt_sample <- function(dt, sample_map, ID=NA){
  intersect_id = intersect(sample_map[,ID],colnames(dt))
  if(length(intersect_id) != nrow(sample_map)){
    message("\033[31m警告\n\tdt和sample_map有数据不匹配\033[0m")
    message("\033[31m\t一共有",length(intersect_id),"个样本可以匹配\033[0m")
    sample_map = sample_map[sample_map[,ID] %in% intersect_id,]
  }
  dt = dt[,sample_map[,ID]] %>% filter(rowSums(.) !=0)
  list(dt=dt, sample_map=sample_map)
}

sample_map <- read.table("/share/data1/limin/vagina/metagenome/metadata_new",header=TRUE,sep="\t")
names(sample_map)[1] <- "Sample"

########
sample_map <- sample_map[sample_map$Body_Site_Category != "Others",]
table(sample_map$Country,sample_map$Body_Site_Category)
sample_map1 <- sample_map[sample_map$Body_Site_Category == "Cervix",]
sample_map1 <- sample_map1[sample_map1$Country == "china"|sample_map1$Country == "usa",]

sample_map2 <- sample_map[sample_map$Country == "usa",]

############
data <- myread("../01.abun/kegg.profiles.known.norm")
colnames(data) <- sub("\\..*", "", colnames(data))
norm_data <- function(dt){
  
  prof <- as.data.frame(apply(dt, 2, function(x) {
    col_sum <- sum(x)
    if (col_sum == 0) {
      return(rep(0, length(x)))  # 列和为0时返回全0
    } else {
      return(x / col_sum * 100)  # 否则正常归一化
    }
  }))
  rownames(prof) = rownames(dt)
  prof
}
dt_norm <- norm_data(data)
#save(dt_norm,file="kegg.known.norm.RData")

load("kegg.known.norm.RData")
dt_norm <- dt_norm[, !duplicated(colnames(dt_norm))]
result <- align_dt_sample(dt_norm, sample_map, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

df_unique <- aligned_sample[!duplicated(aligned_sample$Sample), ]
df_unique <- df_unique[df_unique$Body_Site_Category != "Others",]
table(df_unique$Country,df_unique$Body_Site_Category)

data2 <- rownames_to_column(aligned_dt,var="name")

result <- data2 %>%
  mutate(average = rowMeans(select(., -name), na.rm = TRUE)) %>% 
  select(name, average)
sum(result$average, na.rm = TRUE)

map <- read.table("/share/data1/Database/KEGG/20230401/KO_level_A_B_C_D_Description", 
                  header = F, 
                  fill = TRUE, 
                  sep = "\t", 
                  check.names = FALSE, 
                  quote = "") 
map <- map[,c(2,6,7)]
names(map) <- c("leverA","leverC","name")
merge <- merge(result,map,by="name")

merge <- merge[!grepl("Brite", merge$leverA), ]
merge2 <- merge %>%
  group_by(leverC,leverA) %>%
  summarise(total_average = sum(average, na.rm = TRUE))

df_top50 <- merge2 %>%
  ungroup() %>%              # 【关键】先取消分组，确保后续操作针对整个数据框
  arrange(desc(total_average)) %>%  # 按 total_average 降序排列
  slice_head(n = 50) 

df_top50 <- df_top50 %>% 
  mutate(leverC = gsub("\\s*\\[.*?\\]", "", leverC))

df_top50$leverC <- factor(df_top50$leverC, levels = df_top50$leverC)

# 4. 定义颜色映射 (参考原图配色)
custom_colors <- c(
  "Metabolism" = "#FFCC99",                 # 浅橙色
  "Environmental Information Processing" = "#CDA4DE", # 紫色
  "Cellular Processes" = "#99CC66",         # 绿色
  "Genetic Information Processing" = "#99CCFF", # 蓝色
  "Human Diseases" = "#FF9999"              # 粉红色
)

# 5. 绘图
p <- ggplot(df_top50, aes(x = leverC, y = total_average, fill = leverA)) +
  geom_bar(stat = "identity", width = 0.8) +   # 绘制柱状图
  scale_fill_manual(values = custom_colors) +  # 应用自定义颜色
  theme_bw(base_size = 12) +                   # 使用黑白主题作为基础
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11), # X轴文字倾斜45度
    axis.text.y = element_text(size = 11),  
    legend.position = "right",                     # 图例放在顶部
    legend.title = element_blank(),              # 去掉图例标题
    panel.grid = element_blank(), 
    panel.border = element_blank(),              # 去掉边框
    axis.line = element_line(color = "black"),    # 保留坐标轴线
    aspect.ratio = 0.3
  ) +
  labs(
    x = NULL,                                    # 去掉X轴标题
    y = "Relative abundance (%)",                # Y轴标题
    title = "KEGG Pathway Enrichment Analysis"   # 可选的标题
  )
p

##########pip
merge3 <- merge %>%
  group_by(leverA) %>%
  summarise(total_average = sum(average, na.rm = TRUE))

labels <- merge3$leverA
values <- merge3$total_average
custom_colors <- c(
  "Metabolism" = "#FFCC99",                 # 浅橙色
  "Environmental Information Processing" = "#CDA4DE", # 紫色
  "Cellular Processes" = "#99CC66",         # 绿色
  "Genetic Information Processing" = "#99CCFF", # 蓝色
  "Human Diseases" = "#FF9999",
  "Organismal Systems" = "yellow"
)

pie_colors <- custom_colors[merge3$leverA]

pie(
  x = values, 
  labels = labels, 
  main = "KEGG LeverA",
  col = pie_colors,   # 使用映射好的颜色向量
  cex = 0.7           # 缩小标签字体，防止重叠
)

###############
merge <- merge(result,map,by="name")
gene.map <- read.table("../01.abun/kegg.map", header = F, fill = TRUE, sep = "\t",check.names = F)

merge2 <- merge(gene.map,merge,by.x = "V2",by.y = "name")
merge3 <- merge2[,c(2,4)]

taxo <- read.table("/share/data1/limin/vagina/mmseq/taxon/geneset.taxonomy", fill = TRUE, header = FALSE, sep = "\t")
merge3 <- merge(merge3,taxo[,c(1,7)],by = "V1")

merge4 <- merge3[,-1]
names(merge4) <- c("Function","Genus")
merge4 <- merge4 %>% 
  mutate(Function = gsub("\\s*\\[.*?\\]", "", Function))
merge4 <- merge4[!grepl("Brite", merge4$Function), ]
merge4 <- merge4[merge4$Genus != "g__", ]

merge_dt <- merge4 %>%
  group_by(Function, Genus) %>%      # 1. 按这两列分组
  summarise(Count = n(), .groups = "drop") %>% # 2. 统计每组的数量，并解除分组
  arrange(desc(Count))   

# 2. 计算 Function 和 Genus 的排序（按总流量降序）
func_total <- merge_dt %>%
  group_by(Function) %>%
  summarise(total_count = sum(Count)) %>%
  arrange(desc(total_count))

top_n <- 10
top_functions <- func_total %>%
  slice(1:top_n) %>%
  pull(Function)

# 3. 创建映射：Function -> Group（前10保留，其余为 Others）
df_top10_plus_others <- func_total %>%
  mutate(Group = ifelse(Function %in% top_functions, Function, "Others")) %>%
  select(Function, Group)

top_functions <- df_top10_plus_others %>%
  filter(Group != "Others") %>%
  pull(Group)

func_species_count_filtered <- merge_dt %>%
  mutate(Function = ifelse(Function %in% top_functions, Function, "Others")) %>%
  group_by(Function, Genus) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  ungroup()

###############
species_total <- func_species_count_filtered %>%
  group_by(Genus) %>%
  summarise(total_count = sum(Count), .groups = 'drop') %>%
  arrange(desc(total_count))

top_species <- species_total %>%
  slice(1:10) %>%
  pull(Genus)

func_species_count_final <- func_species_count_filtered %>%
  mutate(Genus = ifelse(Genus == "Unclassified", "Others", Genus))%>%
  mutate(Genus = ifelse(Genus %in% top_species, Genus, "Others")) %>%
  group_by(Function, Genus) %>%
  summarise(Count = sum(Count), .groups = 'drop') %>%
  ungroup()

###############
node_totals <- func_species_count_final %>%
  pivot_longer(cols = c(Function, Genus), names_to = "variable", values_to = "node") %>%
  group_by(node) %>%
  summarise(total_flow = sum(Count), .groups = 'drop') %>%
  arrange(desc(total_flow))

# 分离 Others 和非 Others
non_others <- node_totals %>%
  filter(node != "Others") %>%
  pull(node)

# ✅ 构建最终排序：非 Others 降序 + Others 在最后
desired_order <- c(non_others, "Others")

library(ggsankey) 
# 筛选实际存在的节点
sankey_data_temp <- func_species_count_final %>%
  make_long(x = "Function", next_x = "Genus", value = "Count")

existing_nodes <- unique(sankey_data_temp$node)
desired_order <- desired_order[desired_order %in% existing_nodes]
desired_order <- rev(desired_order)
# =================== 5. 构建桑基图数据 ===================
sankey_data <- func_species_count_final %>%
  make_long(x = "Function", next_x = "Genus", value = "Count")

# 设置 node 因子顺序（决定垂直位置）
sankey_data$node <- factor(sankey_data$node, levels = rev(desired_order), ordered = TRUE)

# 设置 x 轴顺序
sankey_data$x <- forcats::fct_relevel(sankey_data$x, "Function", "Genus")

# =================== ✅ 6. 配色：Others 用灰色，其余用彩虹色 ===================
all_nodes <- levels(sankey_data$node)
n_colors <- length(all_nodes)

# 你的彩虹色
node_colors <- c(
  "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",
  "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
  "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#ed1299",
  "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e", "#4aef7b", "#e86502",
  "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f", "#ce2523",
  "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", "#d66551",
  "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559", "#8e3af4", "#f9a270", "#22547f",
  "#db5e92", "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f", "#7b34c1",
  "#0cf29a", "#d80fc1", "#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5", "#925bea", "#63ff4f"
)

n_rainbow_needed <- if ("Others" %in% all_nodes) {
  n_colors - 1
} else {
  n_colors
}

# 3. 截取或循环颜色 (修正了这里的逻辑)
if (n_rainbow_needed > length(node_colors)) {
  # 如果需要的颜色比色板多，则循环使用色板
  node_colors <- rep(node_colors, length.out = n_rainbow_needed)
} else {
  # 如果色板够用，直接截取
  node_colors <- node_colors[1:n_rainbow_needed]
}

# 4. 构建最终映射
# 找出非-Others的节点
non_others_nodes <- all_nodes[all_nodes != "Others"]

# 确保颜色数量与节点数量一致 (Debug 检查)
if (length(non_others_nodes) != length(node_colors)) {
  stop("颜色数量与节点数量不匹配！请检查 all_nodes 是否有重复值。")
}

# 构建向量
node_color_mapping <- c(
  setNames(node_colors, non_others_nodes),
  "Others" = "gray50"
)


sankey_data$node <- fct_rev(sankey_data$node)
# =================== 7. 绘图 ===================
p <- ggplot(sankey_data,
            aes(x = x,
                next_x = next_x,
                node = node,
                next_node = next_node,
                value = value,
                fill = node)) +
  
  geom_sankey(flow.alpha = 0.8,
              node.color = "black",
              node.linewidth = 0.5,
              show.legend = FALSE) +
  
  # --- 修改点 1: 添加 angle 参数旋转标签 ---
  geom_sankey_label(aes(label = node),
                    hjust = 0.5,
                    size = 5,
                    color = "black",
                    fill = "white",
                    alpha = 0.7,
                    show.legend = FALSE) +
  
  scale_fill_manual(values = node_color_mapping, guide = "none") +
  
  scale_x_discrete(expand = expansion(add = c(0.2, 0.6)),
                   labels = c("x" = "Function", "next_x" = "Species")) +
  
  theme_sankey(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    aspect.ratio = 1.2
  ) +
  
  labs(
    title = "Sankey Diagram: Function → Species",
    subtitle = "Nodes ordered by flow (Others at bottom, gray)",
    x = ""
  )
print(p)
