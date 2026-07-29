####################ko
taxo <- myread("/share/data1/limin/vagina/Geneset/taxon/geneset.taxonomy")
tax <- rownames_to_column(taxo,var="name")
tax <- tax[,1:2]
names(tax) <- c("name","kin")

map <- read.table("../01.abun/kegg.map", sep="\t", header=F, check.names=F)
names(map) <- c("name","ko")

dt3 <- left_join(tax,map,by="name")
dt4 <- dt3[,-1]

dt4$kin[dt4$kin == "Archaea"] <- "Bacteria"
dt4 <- dt4 %>% filter(!is.na(ko))
dt4 <- dt4[dt4$kin !="Unknown",]

data_list <- split( dt4$ko,dt4$kin)

#options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#install.packages("ggvenn", dependencies = TRUE)
library(ggvenn)
p <- ggvenn(data = data_list,              # 传入处理好的列表
            show_elements = FALSE,         # 不显示具体元素名称
            show_percentage = FALSE,        # 显示百分比
            digits = 1,                    # 保留1位小数
            fill_color = c("#4EBBA6", "#F27D52", "#5C5C99"), # 沿用之前的配色
            stroke_size = 0.5,             # 边框粗细
            set_name_size = 5,             # 集合名称字体大小
            text_size = 4) +               # 内部数字字体大小
  ggtitle("Overlap of KOs Across Bacteria, Virus, and Fungi") +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5)
  )
p

####################veen
setwd("/share/data1/limin/vagina/mmseq/analysis/02.kegg/")
setwd("/share/data1/limin/vagina/mmseq/analysis/04.function/")
setwd("/share/data1/limin/vagina/mmseq/analysis/05.amp/")
load("amp.known.norm.RData")
dt <- dt_norm

sample_map1 <- read.table("../02.kegg/cervix.sample",header=TRUE,sep="\t")
sample_map2 <- read.table("../02.kegg/usa.sample",header=TRUE,sep="\t")

result <- align_dt_sample(dt, sample_map1, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- t(aligned_dt)

data <- as.data.frame(t(aligned_dt))
data <- rownames_to_column(data,var="Sample")
merge_data <- merge(data,sample_map1[,c(1,2)],by="Sample")

df_binary <- merge_data
df_binary[,-c(1,ncol(df_binary))] <- (merge_data[,-1] > 0) * 1
merge_data <- df_binary

merge_data1 <- merge_data[merge_data$Country == "China",]
merge_data1 <- merge_data1[,-ncol(merge_data1)]
merge_data2 <- merge_data[merge_data$Country == "U.S.",]
merge_data2 <- merge_data2[,-ncol(merge_data2)]

get_expressed_genes <- function(data, threshold = 0) {
  # 提取基因列（去掉 Sample 列）
  gene_cols <- data[, -1]  # 假设第一列是 Sample
  # 找出在任一样本中丰度 > 阈值的基因
  expressed <- apply(gene_cols, 2, function(x) any(x > threshold, na.rm = TRUE))
  return(names(gene_cols)[expressed])
}

genes_list <- list(
  "China" = get_expressed_genes(merge_data1),
  "U.S."           = get_expressed_genes(merge_data2)
)

lapply(genes_list, length)

all_genes <- unique(unlist(genes_list))

# 创建逻辑矩阵：每个基因在每个集合中是否出现
df_venn <- as.data.frame(lapply(genes_list, function(x) all_genes %in% x))
rownames(df_venn) <- all_genes

#install.packages("ggvenn")
library(ggvenn)

p <- ggvenn(df_venn,
            fill_color = c("#e31a1c", "#fdbf6f"),#"#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99"
            stroke_size = 0.5,
            show_elements = FALSE)+ 
  ggtitle("Gene Expression Overlap Across 5 Sample Types") +
  labs(subtitle = "(Expressed: Abundance > 1 in at least one sample)") +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )
p

###############
result <- align_dt_sample(dt, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- t(aligned_dt)

data <- as.data.frame(t(aligned_dt))
data <- rownames_to_column(data,var="Sample")
merge_data <- merge(data,sample_map2[,c(1,4)],by="Sample")

df_binary <- merge_data
df_binary[,-c(1,ncol(df_binary))] <- (merge_data[,-1] > 0) * 1
merge_data <- df_binary

unique(merge_data$Body_Site_Category)
merge_data1 <- merge_data[merge_data$Body_Site_Category == "Cervix",]
merge_data1 <- merge_data1[,-ncol(merge_data1)]
merge_data2 <- merge_data[merge_data$Body_Site_Category == "Upper vaginal segment",]
merge_data2 <- merge_data2[,-ncol(merge_data2)]
merge_data3 <- merge_data[merge_data$Body_Site_Category == "Mid vagina",]
merge_data3 <- merge_data3[,-ncol(merge_data3)]
merge_data4 <- merge_data[merge_data$Body_Site_Category == "Lower segment of the vagina",]
merge_data4 <- merge_data4[,-ncol(merge_data4)]

get_expressed_genes <- function(data, threshold = 0) {
  # 提取基因列（去掉 Sample 列）
  gene_cols <- data[, -1]  # 假设第一列是 Sample
  # 找出在任一样本中丰度 > 阈值的基因
  expressed <- apply(gene_cols, 2, function(x) any(x > threshold, na.rm = TRUE))
  return(names(gene_cols)[expressed])
}

genes_list <- list(
  "Cervix"                          = get_expressed_genes(merge_data1),
  "Upper vaginal segment"           = get_expressed_genes(merge_data2),
  "Mid vagina"                      = get_expressed_genes(merge_data3),
  "Lower segment of the vagina"     = get_expressed_genes(merge_data4)
)

lapply(genes_list, length)

all_genes <- unique(unlist(genes_list))

# 创建逻辑矩阵：每个基因在每个集合中是否出现
df_venn <- as.data.frame(lapply(genes_list, function(x) all_genes %in% x))
rownames(df_venn) <- all_genes

#install.packages("ggvenn")
library(ggvenn)

p <- ggvenn(df_venn,
            fill_color = c("#a6cee3", "#b2df8a", "#fb9a99", "#cab2d6") ,
            stroke_size = 0.5,
            show_elements = FALSE,
            show_percentage = FALSE)+ 
  theme(
    plot.title = element_text(size = 12, hjust = 1),
    plot.subtitle = element_text(size = 12, hjust = 1)
  )
p

library(UpSetR)
library(gridExtra)
upset(fromList(genes_list), 
      nsets = 5,           # 显示全部5个集合
      nintersects = 30,    # 显示前30个最显著的交集（可根据需要调整或设为NA显示全部）
      order.by = "freq",   # 按照交集的频率（元素数量）降序排列
      mb.ratio = c(0.5, 0.5),
      sets.bar.color = "#56B4E9", # 左侧集合大小条形图的颜色
      point.size = 3,      # 矩阵中点的大小
      line.size = 1,       # 连线的粗细
      text.scale = 1.5     # 整体文字标签的缩放比例
)
