library(Seurat)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(tibble)
library(dplyr)
library(readxl)
library(glue)


# chr5q_list 복구

chr5q <- read_excel("chr5q(11~35).v2025.1.Hs.xlsx", 
                                     col_names = FALSE)

chr5q_list <- list()

for (i in c(1:13)) {
  
  chr5q_list[[i]] <- as.character(chr5q[i, c(-1, -2)])
  chr5q_list[[i]] <- chr5q_list[[i]][!is.na(chr5q_list[[i]])]
  
}

names(chr5q_list) <- chr5q[[1]]

# KEGG_metabolism 리스트 복구

KEGG_metabolism_nc <- read_excel("KEGG_metabolism_nc.xlsx", 
                                 col_names = FALSE)

KEGG_metabolism_nc <- KEGG_metabolism_nc %>%
  column_to_rownames(var = "...1")

KEGG_metabolism_nc <- KEGG_metabolism_nc[, -1]
colnames(KEGG_metabolism_nc) <- c(1:174)

KEGG_metabolism <-  apply(KEGG_metabolism_nc, MARGIN = 1, function(x) x[!is.na(x)])

KEGG_56 <- names(KEGG_metabolism)
KEGG_abc <- order(names(KEGG_metabolism))
KEGG_sort <- names(KEGG_metabolism)[KEGG_abc]
KEGG_sort <- KEGG_sort[-c(2,5,7,8,10,12,15,16,17,25,26,28,32,38,39,42,43,44,48,51,60,61,63,70,76,77,82,83,85)]

KEGG_56 <- KEGG_metabolism[KEGG_sort]

# MDS_35678_list 복구

MDS_celltypes_added <- readRDS("MDS_Celltypes_added.RDS")
MDS_35678 <- subset(MDS_celltypes_added, subset = SubType %in% c("311", "335", "511", "611", "711", "811", "832"))

MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "311", "DX_3", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "335", "Res_3", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "511", "DX_5", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "611", "DX_6", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "711", "Res_7", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "811", "DX_8", MDS_35678@meta.data$Sample_INFO)
MDS_35678@meta.data$Sample_INFO <- ifelse(MDS_35678@meta.data$Sample_INFO == "832", "Res_8", MDS_35678@meta.data$Sample_INFO)

table(MDS_35678@meta.data$Sample_INFO)

celltype_order <- c("Early Erythroid", "Late Erythroid1", "Late Erythroid2",
                    "HSC", "LMPP", "GMP", "CMP", "ERP",
                    "CD4 T cell", "CD8 T cell", "T cell", "Regulatory T cell","NK cell", "B cell", "Plasma cell",
                    "CD14 Monocyte", "CD16 Monocyte", "Macrophage", "Neutrophil",
                    "preDC", "pDC", "Platelet")

MDS_35678_list <- list()

for (celltype in celltype_order) {
  
  MDS_35678_list[[celltype]] <- subset(MDS_35678, subset = Celltypes == celltype)
  
}

# a_chr5q_KEGG : chr5q와 KEGG_56의 교집합.
# KEGG_DEG_list : KEGG_56와 deg_list_MDS(Res와 DX 간 DEG)의 교집합.

# Res와 DX 간의 DEG.

sheet_names <- getSheetNames("MDS_Res_DX_DEG.xlsx")

deg_list_restored <- lapply(sheet_names, function(sn) {
  read.xlsx(
    "MDS_Res_DX_DEG.xlsx",
    sheet = sn,
    colNames = TRUE,     # 반드시 첫 행을 컬럼 이름으로
    rowNames = FALSE,    # 절대 행 이름으로 쓰지 않음
    detectDates = FALSE  # 날짜 자동 인식 방지
  )
})
names(deg_list_restored) <- sheet_names

# p-value 0.05 미만만 선정하기.
deg_list_MDS <- list()

for (celltype in celltype_order) {
  
  deg_list_MDS[[celltype]] <- subset(deg_list_restored[[celltype]], p_val_adj < 0.05)
  
}

# celltype별, KEGG pathway별 DEG 교집합 분류.
KEGG_DEG_list <- list()

for (celltype in celltype_order) {
  for (KEGG in names(KEGG_56)) {
    
    KEGG_DEG_list[[celltype]][[KEGG]] <- intersect(deg_list_MDS[[celltype]]$gene, KEGG_56[[KEGG]])
    
  }
}


# DotPlot list 생성하기

a_KEGG_DotPlot_list <- list()

for (celltype in celltype_order) {
  for (KEGG in names(KEGG_56)) {
    
    a_KEGG_DotPlot_list[[celltype]][[KEGG]] <- 
      DotPlot(MDS_35678_list[[celltype]], features = KEGG_56[[KEGG]], 
              cols = c("skyblue", "red"), group.by = "Sample_INFO") +
      ggtitle(glue("Genes in {celltype} : {KEGG}")) +
      theme(plot.title = element_text(hjust = 0.5,size=20,face='bold')) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
  }
}


a_KEGG_DEG_DotPlot_list <- list()

for (celltype in celltype_order) {
  for (KEGG in names(KEGG_56)) {
    
    if (length(KEGG_DEG_list[[celltype]][[KEGG]]) > 0) {
      
      a_KEGG_DEG_DotPlot_list[[celltype]][[KEGG]] <- 
        DotPlot(MDS_35678_list[[celltype]], features = KEGG_DEG_list[[celltype]][[KEGG]], 
                cols = c("skyblue", "red"), group.by = "Sample_INFO") +
        ggtitle(glue("DEGs in {celltype} : {KEGG}")) +
        theme(plot.title = element_text(hjust = 0.5,size=20,face='bold')) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
      
    }
    
    else next
  }
}

a_KEGG_DEG_DotPlot_list[["ERP"]][["Glycolysis / Gluconeogenesis"]]

print(a_KEGG_DotPlot_list[["CD14 Monocyte"]][["Inositol phosphate metabolism"]])


a_chr5q_KEGG <- list()

for (KEGG in names(KEGG_56)) {
  
  for (location in names(chr5q_list)) {
    
    a_chr5q_KEGG[[KEGG]][[location]] <- intersect(chr5q_list[[location]], KEGG_56[[KEGG]])
    
  }
  
}

# DotPlot 이미지 저장

for (celltype in names(a_KEGG_DEG_DotPlot_list)) {
  for (KEGG in names(a_KEGG_DEG_DotPlot_list[[celltype]])) {
    
    plot_object <- a_KEGG_DEG_DotPlot_list[[celltype]][[KEGG]]
    
    # 유전자 목록이 비어있지 않고 실제로 ggplot 객체인 경우에만 저장
    if (inherits(plot_object, "ggplot")) { 
      
      # 저장할 파일 이름 정의
      file_name <- glue("ggsave/{celltype}_{KEGG}_DotPlot.png")
      
      ggsave(
        filename = file_name, 
        plot = plot_object, 
        width = 12, 
        height = 8
      )
    }
  }
}

DotPlot(MDS_35678_list[["CD14 Monocyte"]], features = c( "COX7C", "NDUFA13", "COX6C",  "NDUFC2", "ATP5PF", "UQCRH", "COX6A1", 
                                                         "NDUFA11", "ATP5F1C", "NDUFB11", "NDUFA12", "NDUFB10", "SDHD", "ATP6V1D", "NDUFAB1", "ATP6V0E2"), 
        cols = c( "skyblue", "red"), group.by = "Sample_INFO") +
  ggtitle(glue("DEGs in CD14 Monocyte : Oxidative phosphorylation")) +
  theme(plot.title = element_text(hjust = 0.5,size=20,face='bold')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

