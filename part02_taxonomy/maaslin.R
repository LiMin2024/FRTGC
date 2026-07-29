setwd("/share/data1/limin/vagina/mmseq/analysis/03.taxo/maaslin/")

library(permute)
library(lattice)
library(vegan)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)
library(tibble)
library(Maaslin2)

myread <- function(x){
  library(data.table)
  dt <- fread(x, sep="\t", header=T, check.names=F, data.table = F, nThread = 80)
  rownames(dt) = dt[,1]
  dt[,-1]
}

###############
load("../taxonomy.known.norm.RData")
dt <- dt_norm

result <- align_dt_sample(dt_norm, sample_map1, ID = "Sample")
aligned_dt <- result$dt
aligned_dt1 <- aligned_dt

result <- align_dt_sample(dt_norm, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_dt2 <- aligned_dt

#############
dt2 <- rownames_to_column(aligned_dt1,var="name")
dt2$id <- paste0("V", rownames(dt2))
dt3 <- dt2 %>%
  select(name,id)
cervix_id <- dt3

dt <- dt2 %>%
  select(-name,-id)

samp_map <- sample_map1
names(samp_map)[1] <- "SampleID"

xx = intersect(colnames(dt), samp_map$SampleID)

## 整体做检验
sampf = subset(samp_map, SampleID %in% xx & (!is.na(Country)))
dtf = t(dt[,sampf$SampleID])
rownames(sampf) = sampf$SampleID

options(contrasts = c("contr.sum", "contr.poly"))
sampf$Country = factor(sampf$Country)
library(Maaslin2)
dir.create("global", recursive = TRUE, showWarnings = FALSE)
res <- Maaslin2(input_data = dtf, input_metadata = sampf,
                output = "global",
                min_abundance = 1e-6, min_prevalence = 0.1, 
                #random_effects = c("BioProjectID", "Country"),
                fixed_effects = c("Country"),
                reference = c("Country,China"),   
                standardize = FALSE, cores=1,
                plot_heatmap = F, plot_scatter = F)

ress = res$results
sig.features <- subset(ress, qval<0.05)
save(ress, file="cervix_masslin.RData")

########################
dt2 <- rownames_to_column(aligned_dt2,var="name")
dt2$id <- paste0("V", rownames(dt2))
dt3 <- dt2 %>%
  select(name,id)
us_id <- dt3

dt <- dt2 %>%
  select(-name,-id)
rownames(dt) <- paste0("V", rownames(dt))

samp_map <- sample_map2
names(samp_map)[1] <- "SampleID"

xx = intersect(colnames(dt), samp_map$SampleID)

## 整体做检验
sampf = subset(samp_map, SampleID %in% xx & (!is.na(Body_Site_Category)))
dtf = t(dt[,sampf$SampleID])
rownames(sampf) = sampf$SampleID

options(contrasts = c("contr.sum", "contr.poly"))
sampf$Body_Site_Category = factor(sampf$Body_Site_Category)
library(Maaslin2)
dir.create("global", recursive = TRUE, showWarnings = FALSE)
res <- Maaslin2(input_data = dtf, input_metadata = sampf,
                output = "global",
                min_abundance = 1e-6, min_prevalence = 0.1, 
                #random_effects = c("BioProjectID", "Country"),
                fixed_effects = c("Body_Site_Category"),
                reference = c("Body_Site_Category,Cervix"),   
                standardize = FALSE, cores=1,
                plot_heatmap = F, plot_scatter = F)

ress <- res$results
sig.features <- subset(ress, qval < 0.05) %>%pull(feature)

# 2. 准备元数据和表达矩阵
sampf <- subset(samp_map, SampleID %in% xx & (!is.na(Body_Site_Category)))
dtf <- t(dt[, sampf$SampleID])
rownames(sampf) <- sampf$SampleID

# 3. 设置因子与获取【原始】分组名称（⚠️ 核心修改：千万不要在这里用 make.names）
options(contrasts = c("contr.treatment", "contr.poly"))
sampf$Body_Site_Category <- factor(sampf$Body_Site_Category)
sample_types <- levels(sampf$Body_Site_Category) 
all_site_res <- list()

# 4. 循环运行 Maaslin2
for (ref in sample_types) {
  # 动态生成并强制创建输出目录（⚠️ 注意：只有文件夹路径需要 make.names 防止空格报错）
  output_path <- file.path("xx", "global", make.names(ref, unique = TRUE))
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  
  message("Running Maaslin2 with reference: ", ref)
  
  # 准备表达矩阵
  dtff <- dtf
  colnames(dtff) <- make.names(colnames(dtff))
  
  # ⚠️ 核心修改：取交集，防止特征名不匹配导致 subscript out of bounds
  common_features <- intersect(sig.features, as.character(colnames(dtff)))
  
  if (length(common_features) == 0) {
    warning(paste("No matching features found for reference:", ref))
    next
  }
  
  fit <- Maaslin2(
    input_data = dtff[, common_features, drop = FALSE], # 使用交集特征
    input_metadata = sampf,
    output = output_path,
    fixed_effects = c("Body_Site_Category"), 
    reference = c(paste0("Body_Site_Category,", ref)), # ⚠️ 核心修改：直接使用原始的 ref
    min_abundance = 1e-6,
    min_prevalence = 0.1, 
    normalization = "TSS",
    transform = "LOG",
    plot_heatmap = FALSE, 
    plot_scatter = FALSE
  )
  
  # 读取结果
  res_file <- file.path(output_path, "all_results.tsv")
  if (file.exists(res_file)) {
    res_temp <- read.csv(res_file, sep = "\t") # 建议改名为 res_temp 防止覆盖外部的 res
    res_temp$Reference <- ref
    all_site_res[[ref]] <- res_temp
  } else {
    warning(paste("Result file not found for reference:", ref))
  }
}

final_res <- bind_rows(all_site_res)
save( res, all_site_res, final_res, file="us_masslin.RData")

###########################
tax <- myread("/share/data1/limin/vagina/mmseq/taxon/geneset.taxonomy")
names(tax) <- c("domain","phylum","class","order","family","genus","species")
tax2 <- tax %>% distinct()

load("cervix_masslin.RData")
cervix <- ress[ress$qval < 0.05,]

dt2 <- rownames_to_column(aligned_dt1,var="name")

high_abundance_features <- dt2 %>%
  # 将除了第1列（物种名）之外的所有列转换为数值型（防止文本干扰）
  mutate(across(-1, as.numeric)) %>%
  # 计算行平均值
  mutate(Mean_Abundance = rowMeans(across(-1), na.rm = TRUE)) %>%
  # 只保留物种名和平均丰度
  select(1, Mean_Abundance)

high_abundance_features2 <- high_abundance_features %>%
  filter(Mean_Abundance > 0.01) 

high_abundance_features3 <- merge(high_abundance_features2,cervix_id,by="name")
names(high_abundance_features3)[1] <- "species"

cervix2 <- merge(cervix,high_abundance_features3[,-2],by.x="feature",by.y="id")
cervix3 <- left_join(cervix2,tax2,by="species")
cervix3 <- cervix3[cervix3$species !="s__",]
#write.table(cervix3,file = "cervix.data", sep = "\t", row.names = T, quote = FALSE)

###########
load("us_masslin.RData")
us <- res[res$qval < 0.05,]

#names(us_id)[1] <- "species"
#us2 <- merge(us,us_id,by.x="feature",by.y="id")
#us3 <- left_join(us2,tax2,by="species")
#unique(us3$feature)

dt2 <- rownames_to_column(aligned_dt2,var="name")

high_abundance_features <- dt2 %>%
  # 将除了第1列（物种名）之外的所有列转换为数值型（防止文本干扰）
  mutate(across(-1, as.numeric)) %>%
  # 计算行平均值
  mutate(Mean_Abundance = rowMeans(across(-1), na.rm = TRUE)) %>%
  # 只保留物种名和平均丰度
  select(1, Mean_Abundance)

high_abundance_features2 <- high_abundance_features %>%
  filter(Mean_Abundance > 0.01) 

high_abundance_features3 <- merge(high_abundance_features2,us_id,by="name")
names(high_abundance_features3)[1] <- "species"

us2 <- merge(us,high_abundance_features3[,-2],by.x="feature",by.y="id")
us3 <- left_join(us2,tax2,by="species")
us3 <- us3[us3$species !="s__",]
#write.table(us3,file = "us.data", sep = "\t", row.names = T, quote = FALSE)

#########
library(tidyr)
library(dplyr)
dt2 <- rownames_to_column(aligned_dt1,var="name")
dt_long <- dt2 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

dt_grouped <- dt_long %>%
  left_join(sample_map1, by = "Sample")

mean_abundance <- dt_grouped %>%
  group_by(name, Country) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Country, values_from = MeanAbundance) %>%
  replace(is.na(.), 0) 

mean_abundance_wide2 <- mean_abundance_wide[mean_abundance_wide$name %in% cervix3$species,]
mean_abundance_wide2 <- mean_abundance_wide2[mean_abundance_wide2$name !="s__",]
#write.table(mean_abundance_wide2,file = "cervix_mean_abun.txt", sep = "\t", row.names = T, quote = FALSE)

#################
library(tidyr)
library(dplyr)
dt2 <- rownames_to_column(aligned_dt2,var="name")
dt_long <- dt2 %>%
  pivot_longer(cols = -name, names_to = "Sample", values_to = "Abundance")

dt_grouped <- dt_long %>%
  left_join(sample_map2, by = "Sample")

mean_abundance <- dt_grouped %>%
  group_by(name, Body_Site_Category) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup()

mean_abundance_wide <- mean_abundance %>%
  pivot_wider(names_from = Body_Site_Category, values_from = MeanAbundance) %>%
  replace(is.na(.), 0) 

mean_abundance_wide2 <- mean_abundance_wide[mean_abundance_wide$name %in% us3$species,]
mean_abundance_wide2 <- mean_abundance_wide2[mean_abundance_wide2$name !="s__",]
#write.table(mean_abundance_wide2,file = "us_mean_abun.txt", sep = "\t", row.names = T, quote = FALSE)

#############
mean_abundance_wide2 <- read.table("us_mean_abun.txt", fill = TRUE, header = T, sep = "\t")
cols <- c("Cervix", "Upper.vaginal.segment", "Mid.vagina","Lower.segment.of.the.vagina")

result2 <- mean_abundance_wide2 %>%
  rowwise() %>%
  mutate(
    high = {
      vals <- c_across(all_of(cols))
      if (all(is.na(vals))) NA_character_ else cols[which.max(vals)]
    },
    low = {
      vals <- c_across(all_of(cols))
      if (all(is.na(vals))) NA_character_ else cols[which.min(vals)]
    }
  ) %>%
  ungroup()

result3 <- result2

############high
Cervix <- result3[result3$high == "Cervix",]
UP <- result3[result3$high == "Upper.vaginal.segment",]
MID <- result3[result3$high == "Mid.vagina",]
LOW <- result3[result3$high == "Lower.segment.of.the.vagina",]

############low
Cervix <- result3[result3$low == "Cervix",]
UP <- result3[result3$low == "Upper.vaginal.segment",]
MID <- result3[result3$low == "Mid.vagina",]
LOW <- result3[result3$low == "Lower.segment.of.the.vagina",]

data2 <- rownames_to_column(aligned_dt2,var="name")

data3 <- data2[data2$name %in% Cervix$name,]#################################
rownames(data3) <- NULL
data3 <- column_to_rownames(data3,var="name")

result <- align_dt_sample(data3, sample_map2, ID = "Sample")
aligned_dt <- result$dt
aligned_sample <- result$sample_map

intersect_id = intersect(colnames(aligned_dt), sample_map$Sample)
sampf = unique(subset(sample_map, Sample %in% intersect_id, c("Body_Site_Category","Sample")))
rownames(sampf) = sampf$Sample
dtf = t(aligned_dt[,sampf$Sample])
dtf <- as.data.frame(dtf)
rownames(dtf) = rownames(sampf)

sampf2 <- sampf %>%
  mutate(Body_Site_Category = ifelse(Body_Site_Category == "Cervix", Body_Site_Category, "Others"))
#Cervix Upper vaginal segment Mid vagina Lower segment of the vagina

output_dir <- "Cervix_low"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 运行 MaAsLin2，Day 为固定效应，PRJ 作为协变量
fit_data <- Maaslin2(
  input_data = dtf,
  input_metadata = sampf2,
  output = output_dir,
  min_prevalence = 0,
  normalization = "NONE",
  fixed_effects = c("Body_Site_Category"),       
  random_effects = c(),           
  reference = c("Body_Site_Category,Others"), 
  plot_heatmap = FALSE,
  plot_scatter = FALSE
)

Cervix_high <- read.table("Cervix_high/all_results.tsv", sep="\t", header=T, check.names=F)
Cervix_high$high <- "Cervix"

UP_high <- read.table("UP_high/all_results.tsv", sep="\t", header=T, check.names=F)
UP_high$high <- "Upper vaginal segment"

MID_high <- read.table("MID_high/all_results.tsv", sep="\t", header=T, check.names=F)
MID_high$high <- "Mid vagina"

LOW_high <- read.table("LOW_high/all_results.tsv", sep="\t", header=T, check.names=F)
LOW_high$high <- "Lower vaginal segment"

high_data <- rbind(Cervix_high,UP_high,MID_high,LOW_high)
high_data$enriched <- ifelse(high_data$qval < 0.05, high_data$high, NA_character_)

Cervix_low <- read.table("Cervix_low/all_results.tsv", sep="\t", header=T, check.names=F)
Cervix_low$low <- "Cervix"

UP_low <- read.table("UP_low/all_results.tsv", sep="\t", header=T, check.names=F)
UP_low$low <- "Upper vaginal segment"

MID_low <- read.table("MID_low/all_results.tsv", sep="\t", header=T, check.names=F)
MID_low$low <- "Mid vagina"

LOW_low <- read.table("LOW_low/all_results.tsv", sep="\t", header=T, check.names=F)
LOW_low$low <- "Lower vaginal segment"

low_data <- rbind(Cervix_low,UP_low,MID_low,LOW_low)
low_data$deficient <- ifelse(low_data$qval < 0.05, low_data$low, NA_character_)

merge_dt <- merge(high_data[,c(1,ncol(high_data))],low_data[,c(1,ncol(low_data))],by="feature")
#write.table(merge_dt,file = "us_enriched_deficient.txt", sep = "\t", row.names = F, quote = FALSE)
