setwd("/share/data1/limin/vagina/mmseq/analysis/03.taxo/maaslin/")
library(Maaslin2)
library(permute)
library(lattice)
library(stringr)
library(tidyr)
library(dplyr)
library(tibble)

myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

result <- read.table("cervix.data", sep="\t", header=T, check.names=F)
result <- result %>% distinct(feature, .keep_all = TRUE)
result <- result[,-c(2:10)]
label <- result[,1:2]

result2 <- result[,-2]
names(result2) <- c("species","superkingdom","phylum","class","order","family","genus")

####################tree

taxonomy_df <- result2

char_cols <- sapply(taxonomy_df, is.character)

#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("data.tree")
library(data.tree)
library(ape)

# 构建路径列表
tree_nodes <- list()
for (i in 1:nrow(taxonomy_df)) {
  row <- taxonomy_df[i, ]
  
  path <- paste0(
    row$superkingdom, "/", 
    row$phylum, "/",
    row$class, "/",
    row$order, "/",
    row$family, "/",
    row$genus, "/",
    row$species
  )
  
  node <- unlist(strsplit(path, "/"))
  tree_nodes[[row$species]] <- as.list(node)
}

# 合并为路径字符串
paths <- lapply(tree_nodes, function(x) paste(x, collapse = "/"))
path_df <- data.frame(path = unlist(paths), stringsAsFactors = FALSE)

clean_path_df <- data.frame(
  path = path_df$path,
  stringsAsFactors = FALSE
)

# 提取叶子名（最后一个 / 后面的内容）
clean_path_df$id <- sub(".*/", "", clean_path_df$path)

# 定义函数：将路径添加到树中
add_path_to_tree <- function(tree, parts, leaf_name) {
  if (length(parts) == 1) {
    tree[[leaf_name]] <- list()
  } else {
    current <- parts[1]
    rest <- parts[-1]
    
    if (is.null(tree[[current]])) {
      tree[[current]] <- list()
    }
    
    # 递归处理子路径，并将结果正确赋值回当前节点
    tree[[current]] <- add_path_to_tree(tree[[current]], rest, leaf_name)
  }
  return(tree)
}

# 定义函数：将树结构转为 Newick 格式
tree_to_newick <- function(tree) {
  if (is.null(tree) || length(tree) == 0) return("")
  
  children <- names(tree)
  child_strings <- sapply(children, function(child) {
    child_tree <- tree[[child]]
    if (length(child_tree) == 0) {
      return(child)
    } else {
      return(paste0("(", tree_to_newick(child_tree), ")", child))
    }
  })
  
  if (length(child_strings) == 0) return("")
  
  return(paste0("(", paste(child_strings, collapse = ","), ")"))
}

# 初始化空树
tree <- list()

# 遍历所有行，添加路径到树中
for (i in seq_len(nrow(clean_path_df))) {
  parts <- unlist(strsplit(clean_path_df$path[i], "/"))
  leaf_name <- clean_path_df$id[i]
  tree <- add_path_to_tree(tree, parts, leaf_name)
}

# 转换为 Newick 格式
newick_str <- paste0(tree_to_newick(tree), ";")

#writeLines(newick_str, "cervix.newick")

###################分支颜色
species_to_phylum <- taxonomy_df[, c("species", "phylum")]
species_to_phylum <- as.data.frame(species_to_phylum)


# === 自定义颜色列表 ===
total_color1 = c("#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                 "#6a3d9a", "#ffff99", "#b15928","#8dd3c7",
                 "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
                 "#fdb462", "#b3de69", "#fccde5", "#bc80bd",
                 "#ccebc5", "#ffed6f", "#a6cee3", "#1f78b4",
                 "#b2df8a", "#33a02c", "#fb9a99","#ed1299",
                 "#09f9f5","#246b93","#cc8e12","#d561dd",
                 "#c93f00","#ddd53e","#4aef7b","#e86502",
                 "#9ed84e","#39ba30","#6ad157","#8249aa",
                 "#99db27","#e07233","#ff523f","#ce2523",
                 "#f7aa5d","#cebb10","#03827f","#931635",
                 "#373bbf","#a1ce4c","#ef3bb6","#d66551",
                 "#1a918f","#ff66fc","#2927c4","#7149af",
                 "#57e559","#8e3af4","#f9a270","#22547f",
                 "#db5e92","#edd05e","#6f25e8","#0dbc21",
                 "#280f7a","#6373ed","#5b910f","#7b34c1",
                 "#0cf29a","#d80fc1","#dd27ce","#07a301",
                 "#167275","#391c82","#2baeb5","#925bea","#63ff4f")

# === 提取唯一 phylum ===
unique_phyla <- unique(species_to_phylum$phylum)
n_phyla <- length(unique_phyla)
n_colors <- length(total_color1)

# 检查颜色是否足够
if (n_colors < n_phyla) {
  stop("错误：颜色数量不足！有 ", n_phyla, " 个门，但只有 ", n_colors, " 个颜色。")
}

# 使用你的自定义颜色（按顺序或随机）
assigned_colors <- total_color1[1:n_phyla]  # 按顺序取前 n_phyla 个颜色
# 或者随机打乱：assigned_colors <- sample(total_color1, n_phyla)

# 正确创建 phylum_to_color 数据框（不要覆盖！）
phylum_to_color <- data.frame(
  phylum = unique_phyla,
  color = assigned_colors,
  stringsAsFactors = FALSE
)

# === 现在可以安全合并 ===
data <- left_join(species_to_phylum, phylum_to_color, by = "phylum")

# 再与 clean_path_df 合并
merge_data <- merge(clean_path_df, data, by.x = "id", by.y = "species")
merge_data <- merge_data[3:4]

df <- cbind(merge_data[1], clade = "clade", merge_data[2:ncol(merge_data)])
df <- cbind(df[1:3], normal = "normal")
df <- cbind(df[1:4], size = "3")
length(unique(df$phylum))
table(df$phylum)
# 打开连接（写入文件）
file_conn <- file("branch_color.txt", "w")

# 写入前三行
writeLines("TREE_COLORS", file_conn)
writeLines("SEPARATOR TAB", file_conn)  
writeLines("DATA", file_conn)

# 写入数据（用制表符分隔）
write.table(df, file_conn, 
            sep = "\t",             # 使用制表符分隔
            row.names = FALSE,      # 不写行名
            col.names = FALSE,      # 不重复写列名（因为前面已有 DATA）
            quote = FALSE)          # 不给字符串加引号

# 关闭连接
close(file_conn)

############
color_df <- df
color_df <- color_df[,c(1,3)]
color_df <- color_df[!duplicated(color_df), ]

# 创建命名向量
my_colors <- setNames(color_df$color, color_df$phylum)

# 绘制纯图例
ggplot(color_df, aes(x = 1, y = phylum, color = phylum)) + 
  geom_point() + # 这里用固定的 x=1 和 V1 作为 y 轴来触发图例生成
  scale_color_manual(values = my_colors) +
  labs(color = "Phylum") +
  theme(
    axis.title = element_blank(),   # 隐藏坐标轴标题
    axis.text.y = element_blank(),  # 隐藏 Y 轴文字
    axis.ticks = element_blank(),   # 隐藏刻度线
    panel.background = element_blank() # 隐藏背景网格
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 4))) # 调大色块大小

################
tax <- label

file_name <- "us.labels.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入注释和头部信息
writeLines(
  c(
    "LABELS",
    "#use this template to change the leaf labels, or define/change the internal node names (displayed in mouseover popups)",
    "#lines starting with a hash are comments and ignored during parsing",
    "",
    "#=================================================================#",
    "#          MANDATORY SETTINGS              #",
    "#=================================================================#",
    "#select the separator which is used to delimit the data below (TAB,SPACE or COMMA).This separator must be used throughout this file (except in the SEPARATOR line, which uses space).",
    "SEPARATOR TAB",
    "#SEPARATOR SPACE",
    "#SEPARATOR COMMA",
    "",
    "#=================================================================#",
    "#    Actual data follows after the \"DATA\" keyword       #",
    "#=================================================================#",
    "DATA"
  ),
  con = file_conn
)

# 5. 写入数据部分（用制表符分隔）
write.table(tax, file_conn,
            sep = "\t",           # 使用 TAB 分隔
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

close(file_conn)

##############abun
mean_abundance_wide = read.table("us_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(result,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,c(2,9:ncol(merge_data))]

merge_data2 <- column_to_rownames(merge_data,var="feature")

row_scaled <- t(apply(merge_data2, 1, function(x) {
  (x - mean(x)) / sd(x)
}))

row_scaled_df <- rownames_to_column(as.data.frame(row_scaled), var = "name")
neg_log_df <- merge_data2 %>%
  mutate(across(everything(), ~ log10(.)))%>%
  rownames_to_column(var="name")

file_name <- "us_abun_scale.txt"

# 3. 打开连接写入文件
file_conn <- file(file_name, "w")

# 4. 写入头部信息
writeLines(
  c(
    "DATASET_HEATMAP",
    "SEPARATOR SPACE",
    "DATASET_LABEL Effect Size Heatmap",
    "# 渐变颜色设置",
    "COLOR_MIN #1f78b4",#lightblue #1f78b4
    "COLOR_MAX #e31a1c", ##1f78b4 #e31a1c
    "USE_MID_COLOR 1",
    "COLOR_MID #FFFFFF",
    "STRIP_WIDTH 25",
    "MARGIN 0",
    "SHOW_TREE 1",
    "FIELD_LABELS Cervix Lower Mid Upper",#Cervix Upper  Mid  Lower
    "DATA"
  ),
  con = file_conn
)
write.table(row_scaled_df, file_conn,
            sep = " ",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 6. 关闭连接
close(file_conn)

###############
enriched_data <- read.table("us_enriched_deficient.txt", header = T,  sep = "\t",check.names = FALSE)
names(enriched_data)[1] <- "name"
enriched_data_clean <- enriched_data %>%
  mutate(
    # 将除了大小写字母以外的所有字符替换为空字符串
    name_clean = str_replace_all(name, "[^A-Za-z]", "")
  )
us_id_clean <- us_id %>%
  mutate(
    # 将除了大小写字母以外的所有字符替换为空字符串
    name_clean = str_replace_all(name, "[^A-Za-z]", "")
  )

enriched_data2 <- merge(enriched_data_clean,us_id_clean,by="name_clean")
enriched_data3 <- enriched_data2[,c(6,3,4)]

df_long <- enriched_data3 %>%
  pivot_longer(cols = c(enriched, deficient), 
               names_to = "status", 
               values_to = "Body_Site_Category") %>%
  # 移除缺失值（如果有 NA）
  filter(!is.na(Body_Site_Category))

# 转换为宽格式：Body_Site_Category 作为列，status 作为值
df_wide <- df_long %>%
  pivot_wider(names_from = Body_Site_Category, 
              values_from = status, 
              values_fill = NA)

#write.table(df_wide,file = "us_enriched_deficient_wide.txt", sep = "\t", row.names = T, quote = FALSE)

data1 <- df_wide%>%
  select(id,Cervix)
names(data1)[2] <- "enriched"
data1$enriched <- ifelse(data1$enriched == "deficient", "1,4,#391c82,1,-1", data1$enriched)
data1$enriched <- ifelse(data1$enriched == "enriched", "1,4,#a6cee3,1,-1", data1$enriched)

data2 <- df_wide%>%
  select(id,'Lower vaginal segment')
names(data2)[2] <- "enriched"
data2$enriched <- ifelse(data2$enriched == "deficient", "1,4,#2baeb5,1,-1", data2$enriched)
data2$enriched <- ifelse(data2$enriched == "enriched", "1,4,#cab2d6,1,-1", data2$enriched)

data3 <- df_wide%>%
  select(id,'Mid vagina')
names(data3)[2] <- "enriched"
data3$enriched <- ifelse(data3$enriched == "deficient", "1,4,#925bea,1,-1", data3$enriched)
data3$enriched <- ifelse(data3$enriched == "enriched", "1,4,#fb9a99,1,-1", data3$enriched)

data4 <- df_wide%>%
  select(id,'Upper vaginal segment')
names(data4)[2] <- "enriched"
data4$enriched <- ifelse(data4$enriched == "deficient", "1,4,#63ff4f,1,-1", data4$enriched)
data4$enriched <- ifelse(data4$enriched == "enriched", "1,4,#b2df8a,1,-1", data4$enriched)


#colors <- ("#391c82","#2baeb5","#925bea","#63ff4f")
#named_colors <- c(
#  "Cervix" = "#a6cee3", 
#  "Upper vaginal segment" = "#b2df8a",
#  "Mid vagina" = "#fb9a99", 
#  "Lower segment of the vagina" = "#cab2d6"
#)

file_name <- "up_enriched.txt"

# 5. 打开连接，写入头信息 + DATA + 数据
file_conn <- file(file_name, "w")

writeLines(
  c(
    "DATASET_SYMBOL",
    "#lines starting with a hash are comments and ignored during parsing",
    "#=================================================================#",
    "#                    MANDATORY SETTINGS                           #",
    "#=================================================================#",
    "SEPARATOR COMMA",
    "",
    "# 标签（图例中显示）",
    "DATASET_LABEL,Two-color Star Symbols",
    "",
    "# 数据集整体颜色（可选，不影响单个符号）",
    "COLOR,#ffff00",
    "",
    "#=================================================================#",
    "#                    OPTIONAL SETTINGS                            #",
    "#=================================================================#",
    "",
    "# 最大符号大小",
    "MAXIMUM_SIZE,40",
    "",
    "#=================================================================#",
    "#                       LEGEND 设置                              #",
    "#=================================================================#",
    "",
    "# 图例标题",
    "LEGEND_TITLE,Data Types",
    "",
    "# 图例形状：两个都是星号（实际用三角形）",
    "LEGEND_SHAPES,1,1",
    "",
    "# 颜色：黄色代表类型A，红色代表类型B",
    "LEGEND_COLORS,#ffff00,#ff0000",
    "",
    "# 图例标签",
    "LEGEND_LABELS,Type A,Type B",
    "",
    "# 可选：调整形状大小比例（保持一致）",
    "LEGEND_SHAPE_SCALES,1,1",
    "",
    "# 图例自动定位（不固定位置）",
    "# LEGEND_POSITION_X,80",
    "# LEGEND_POSITION_Y,80",
    "",
    "#=================================================================#",
    "#                       ACTUAL DATA                              #",
    "#=================================================================#",
    "DATA"
  ),
  con = file_conn
)

# 写入数据部分（逗号分隔，不带列名）
write.table(data4, file_conn,
            sep = ",", 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 关闭连接
close(file_conn)

##############cervix
data <- read.table("cervix.data",header = T,sep = "\t")

data <- data %>%
  mutate(enriched = case_when(
    coef < 0 ~ "China",
    coef > 0 ~ "U.S.",
    TRUE ~ ""
  ))

data2 <- data%>%
  select(feature,enriched)

data2$enriched <- ifelse(data2$enriched == "China", "3,40,#DC143C,1,-1", data2$enriched)

data2$enriched <- ifelse(data2$enriched == "U.S.", "3,40,#0077B5,1,-1", data2$enriched)

file_name <- "cervix.enriched.star.txt"

# 5. 打开连接，写入头信息 + DATA + 数据
file_conn <- file(file_name, "w")

writeLines(
  c(
    "DATASET_SYMBOL",
    "#lines starting with a hash are comments and ignored during parsing",
    "#=================================================================#",
    "#                    MANDATORY SETTINGS                           #",
    "#=================================================================#",
    "SEPARATOR COMMA",
    "",
    "# 标签（图例中显示）",
    "DATASET_LABEL,Two-color Star Symbols",
    "",
    "# 数据集整体颜色（可选，不影响单个符号）",
    "COLOR,#ffff00",
    "",
    "#=================================================================#",
    "#                    OPTIONAL SETTINGS                            #",
    "#=================================================================#",
    "",
    "# 最大符号大小",
    "MAXIMUM_SIZE,40",
    "",
    "#=================================================================#",
    "#                       LEGEND 设置                              #",
    "#=================================================================#",
    "",
    "# 图例标题",
    "LEGEND_TITLE,Data Types",
    "",
    "# 图例形状：两个都是星号（实际用三角形）",
    "LEGEND_SHAPES,3,3",
    "",
    "# 颜色：黄色代表类型A，红色代表类型B",
    "LEGEND_COLORS,#ffff00,#ff0000",
    "",
    "# 图例标签",
    "LEGEND_LABELS,Type A,Type B",
    "",
    "# 可选：调整形状大小比例（保持一致）",
    "LEGEND_SHAPE_SCALES,1,1",
    "",
    "# 图例自动定位（不固定位置）",
    "# LEGEND_POSITION_X,80",
    "# LEGEND_POSITION_Y,80",
    "",
    "#=================================================================#",
    "#                       ACTUAL DATA                              #",
    "#=================================================================#",
    "DATA"
  ),
  con = file_conn
)

# 写入数据部分（逗号分隔，不带列名）
write.table(data2, file_conn,
            sep = ",", 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 关闭连接
close(file_conn)

###########
mean_abundance_wide = read.table("cervix_mean_abun.txt", sep="\t", header=T, check.names=F)

merge_data <- merge(result,mean_abundance_wide,by.x="species",by.y="name")
merge_data <- merge_data[,c(2,9:ncol(merge_data))]

merge_data2 <- column_to_rownames(merge_data,var="feature")

merge_data3 <- merge_data2 %>%
  mutate(across(everything(), ~ sqrt(.)))%>%
  rownames_to_column(var="name")

china <- merge_data3[,1:2]
us <- merge_data3[,c(1,3)]

###########
file_name <- "carvix.china_bars.txt"

file_conn <- file(file_name, "w")

# 写入头部信息
writeLines(
  c(
    "DATASET_SIMPLEBAR",
    "SEPARATOR COMMA",
    "DATASET_LABEL,Donor enriched Relative Abundance (%)",
    "DATASET_SCALE,0-Low-#FFFFFF-5-1-1,3-Medium-#A9A9A9-5-1-1,6-High-#4D4D4D-5-1-1",
    "COLOR,#0077B5",
    "WIDTH,100",
    "MARGIN,0",
    "HEIGHT_FACTOR,1",
    "BAR_SHIFT,0",
    "BAR_ZERO,hide",
    "DATA"
  ),
  con = file_conn
)
#DC143C #0077B5
# 写入数据部分（逗号分隔，不带列名）
write.table(china, file_conn,
            sep = ",", 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 关闭连接
close(file_conn)

############
file_name <- "bac_pre_enriched_bars.txt"

file_conn <- file(file_name, "w")

# 写入头部信息
writeLines(
  c(
    "DATASET_SIMPLEBAR",
    "SEPARATOR COMMA",
    "DATASET_LABEL,pre-FMT enriched Relative Abundance (%)",
    "DATASET_SCALE,0-Low-#FFFFFF-5-1-1,1.5-Medium-#A9A9A9-5-1-1,3-High-#4D4D4D-5-1-1",
    "COLOR,#DC143C",
    "WIDTH,100",
    "MARGIN,0",
    "HEIGHT_FACTOR,1",
    "BAR_SHIFT,0",
    "BAR_ZERO,hide",
    "DATA"
  ),
  con = file_conn
)

# 写入数据部分（逗号分隔，不带列名）
write.table(pre, file_conn,
            sep = ",", 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# 关闭连接
close(file_conn)
