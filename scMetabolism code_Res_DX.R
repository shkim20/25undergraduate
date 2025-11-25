library(scMetabolism)
library(ggplot2)
library(rsvd)
library(tidyr)
library(tibble)

# scMetabolism으로 분석, K-S 정규성 검정, Wilcoxson Rank-Sum test 분석.
# wilcoxson p-value와 metabolic score를 구간별로 나누어 pheatmap으로 그리는 과정을 포함한 코드.

# 1. 데이터 준비(LNLD_Res, DX를 함께 묶어 scMetabolism으로 분석해야 한다.)

MDS_Res_DX <- subset(`MDS_Celltypes added`, subset = Sample_State %in% c("LNLD_Res", "DX"))

MDS_Res_DX_scMeta<-sc.metabolism.Seurat(obj = MDS_Res_DX, method = "AUCell", 
                                        imputation = F, ncores = 2, metabolism.type = "KEGG")

MDS_Res_DX_3 <- subset(MDS_Res_DX, subset = SampleID %in% c("MDS-311", "MDS-335"))
MDS_Res_DX_567 <- subset(MDS_Res_DX, subset = SampleID %in% c("MDS-511", "MDS-611", "MDS-711"))
MDS_Res_DX_8 <- subset(MDS_Res_DX, subset = SampleID %in% c("MDS-811", "MDS-832"))

MDS_3_scMeta <- sc.metabolism.Seurat(obj = MDS_Res_DX_3, method = "AUCell", 
                                          imputation = F, ncores = 2, metabolism.type = 
                                            "KEGG")

MDS_567_scMeta <- sc.metabolism.Seurat(obj = MDS_Res_DX_567, method = "AUCell", 
                                       imputation = F, ncores = 2, metabolism.type = 
                                         "KEGG")

MDS_8_scMeta <- sc.metabolism.Seurat(obj = MDS_Res_DX_8, method = "AUCell", 
                                     imputation = F, ncores = 2, metabolism.type = 
                                       "KEGG")

scMeta_score_3 <-  MDS_3_scMeta@assays[["METABOLISM"]][["score"]]
scMeta_score_567 <-  MDS_567_scMeta@assays[["METABOLISM"]][["score"]]
scMeta_score_8 <-  MDS_8_scMeta@assays[["METABOLISM"]][["score"]]

# 만약 모든 마침표(.)를 하이픈(-)으로 바꾸어도 문제가 없다면:
colnames(scMeta_score_8) <- gsub("\\.", "-", colnames(scMeta_score_8))

Res_metadata <- MDS_Res_DX_8@meta.data[MDS_Res_DX_8@meta.data$Sample_State == "LNLD_Res", ]
DX_metadata <- MDS_Res_DX_8@meta.data[MDS_Res_DX_8@meta.data$Sample_State == "DX", ]

scMeta_score_Res <- scMeta_score_8[, colnames(scMeta_score_8) %in% rownames(Res_metadata) ]
scMeta_score_DX <- scMeta_score_8[, colnames(scMeta_score_8) %in% rownames(DX_metadata) ]


# scMetabolism 결과들을 다시 celltype에 따라 분류해야 한다...

# meta.data에서 각 셀의 celltype 정보 가져오기 (Res & DX)
celltype_info_Res <- Res_metadata[colnames(scMeta_score_Res), "Celltypes"]
names(celltype_info_Res) <- rownames(Res_metadata)

celltype_info_DX <- DX_metadata[colnames(scMeta_score_DX), "Celltypes"]
names(celltype_info_DX) <- rownames(DX_metadata)

# celltypes 순서대로 scMeta_score_8 열을 묶어 list로 저장

celltypes <- unique(MDS_Res_DX@meta.data$Celltypes)

scMeta_split_Res <- lapply(celltypes, function(ct) {
  cols <- names(celltype_info_Res)[celltype_info_Res == ct]
  scMeta_score_Res[, cols, drop = FALSE]   # 각 celltype별 subset
})

scMeta_split_DX <- lapply(celltypes, function(ct) {
  cols <- names(celltype_info_DX)[celltype_info_DX == ct]
  scMeta_score_DX[, cols, drop = FALSE]   # 각 celltype별 subset
})

# celltypes 이름을 list에 부여
names(scMeta_split_Res) <- celltypes
names(scMeta_split_DX) <- celltypes

# 2. Kolmogorov-Smirnov Test(KS test, KS 검정) 실시. 정규성 검정.

ks_pnorm_results_Res <- list()
ks_pnorm_results_DX <- list()

for (pathway in rownames(scMeta_score_8)) {
  
  # 하위 리스트를 비어있는 리스트로 명시적으로 초기화
  ks_pnorm_results_Res[[pathway]] <- list()
  ks_pnorm_results_DX[[pathway]] <- list()
  
  for (celltype in celltypes) {
    
    Res_vector <- as.numeric(as.matrix(scMeta_split_Res[[celltype]][pathway, ]))
    DX_vector <- as.numeric(as.matrix(scMeta_split_DX[[celltype]][pathway, ]))
    
    ks_pnorm_results_Res[[pathway]][[celltype]] <- ks.test(Res_vector, "pnorm")
    ks_pnorm_results_DX[[pathway]][[celltype]] <- ks.test(DX_vector, "pnorm")
    
  }
  
}

for (pathway in rownames(scMeta_score_8)) {
  
  for (celltype in celltypes) {
    
    if (ks_pnorm_results_Res[[pathway]][[celltype]][["p.value"]] >= 0.05) {
      
      cat( pathway, " - " , celltype, " : ", 
          ks_pnorm_results_Res[[pathway]][[celltype]][["p.value"]], "\n")
      
    }
    
  }
  
}

# 출력문이 없음 : 모든 pathway, celltype에 대해 p.value < 0.05 -> 정규성이 없음.
# 그러므로 독립표본 t-검정 대신 Wilcoxson Rank-Sum Test을 실시.



# 2-2. Shapiro-Wilk test 실시. 정규성 검정.

sw_pnorm_results_Res <- list()
sw_pnorm_results_DX <- list()

for (pathway in rownames(scMeta_score_8)) {
  
  # 하위 리스트를 비어있는 리스트로 명시적으로 초기화
  sw_pnorm_results_Res[[pathway]] <- list()
  sw_pnorm_results_DX[[pathway]] <- list()
  
  for (celltype in celltypes) {
    
    Res_vector <- as.numeric(as.matrix(scMeta_split_Res[[celltype]][pathway, ]))
    DX_vector <- as.numeric(as.matrix(scMeta_split_DX[[celltype]][pathway, ]))
    
    sw_pnorm_results_Res[[pathway]][[celltype]] <- shapiro.test(Res_vector)
    sw_pnorm_results_DX[[pathway]][[celltype]] <- shapiro.test(DX_vector)
    
  }
  
}

for (pathway in rownames(scMeta_score_8)) {
  
  for (celltype in celltypes) {
    
    if (sw_pnorm_results_Res[[pathway]][[celltype]][["p.value"]] >= 0.05) {
      
      cat( pathway, " - " , celltype, " : ", 
           sw_pnorm_results_Res[[pathway]][[celltype]][["p.value"]], "\n")
      
    }
    
  }
  
}



# 근데 3번 환자의 경우는 pDC, CD16 Monocyte에서 p.value > 0.05, cell 수 적음.
# 5,6,7번 환자는 Res group의 CD16 Monocyte, pDC, Macrophage cell 수가 0임.


# 3. Wilcoxson Rank-Sum Test 실시.

wilcox_scMeta <- list()

for (pathway in rownames(scMeta_score_8)) {
  
  wilcox_scMeta[[pathway]] <- list()
  
  for (celltype in celltypes) {
    
    Res_vector <- as.numeric(as.matrix(scMeta_split_Res[[celltype]][pathway, ]))
    DX_vector <- as.numeric(as.matrix(scMeta_split_DX[[celltype]][pathway, ]))
    
    wilcox_scMeta[[pathway]][[celltype]] <- wilcox.test(Res_vector, DX_vector, 
                                                        alternative = "two.sided", paired = FALSE)
    
  }
  
}

for (pathway in rownames(scMeta_score_8)) {
  
  for (celltype in celltypes)
    
    if (wilcox_scMeta[[pathway]][[celltype]][["p.value"]] >= 0.05) {
      
      cat(pathway, " - ", celltype, " : ", 
          wilcox_scMeta[[pathway]][[celltype]][["p.value"]], "\n")
    }
  
}

wilcox_scMeta_pvalue <- data.frame(NULL)
wilcox_scMeta_pvalue <- as.data.frame(
  t(sapply(wilcox_scMeta, function(x) {
    sapply(x, function(y) y[["p.value"]])
  }))
)



# for-loop 버전 :

# # pathway × celltype 크기 지정
# pval_matrix <- matrix(
#   NA, 
#   nrow = length(rownames(scMeta_score_8)), 
#   ncol = length(celltypes),
#   dimnames = list(rownames(scMeta_score_8), celltypes)
# )
# 
# for (pathway in rownames(scMeta_score_8)) {
#   for (celltype in celltypes) {
#     pval_matrix[pathway, celltype] <- wilcox_scMeta[[pathway]][[celltype]][["p.value"]]
#   }
# }
# 
# # data.frame으로 변환
# wilcox_scMeta_pvalue <- as.data.frame(pval_matrix)

wilcox_scMeta_pvalue$pvalue <- apply(wilcox_scMeta_pvalue, 1, function(row) {
  sum(row < 0.05, na.rm = TRUE)  # 0.05 미만 값 개수
})
wilcox_scMeta_pvalue <- rbind(
  wilcox_scMeta_pvalue,
  pvalue_celltype = apply(wilcox_scMeta_pvalue, 2, function(x) sum(x < 0.05, na.rm = TRUE))
)

wilcox_scMeta_heatmap <- wilcox_scMeta_pvalue[1:84, 1:22]





# scMetabolism score를 LNLD_Res, DX에 따라 구분해두기.
names(scMeta_split_Res) <- paste0(names(scMeta_split_Res), "_Res")
names(scMeta_split_DX) <- paste0(names(scMeta_split_DX), "_DX")
scMeta_split <- c(scMeta_split_Res, scMeta_split_DX)

MDS_Res_DX_scMeta_save <- MDS_Res_DX_scMeta

MDS_Res_DX_scMeta@meta.data$Celltypes <- ifelse(
  MDS_Res_DX_scMeta@meta.data$Sample_State == "LNLD_Res",
  paste0(MDS_Res_DX_scMeta@meta.data$Celltypes, "_Res"),
  paste0(MDS_Res_DX_scMeta@meta.data$Celltypes, "_DX")
)
MDS_Res_DX_scMeta@meta.data$Celltypes <- as.character(
  MDS_Res_DX_scMeta@meta.data$Celltypes
)



plot_obj_8 <- DotPlot.metabolism(obj = MDS_8_scMeta, 
                               pathway = rownames(MDS_8_scMeta@assays[["METABOLISM"]][["score"]]), 
                               phenotype = "Celltypes", 
                               norm = "y")

print(plot_obj_8)


# 2. DotPlot 내부 데이터 (색상 매핑용) 추출

# 사용자 정의 CellType 순서 벡터를 다시 정의합니다.
celltype_order <- c("Early Erythroid", "Late Erythroid1", "Late Erythroid2",
                    "HSC", "LMPP", "GMP", "CMP", "ERP",
                    "CD4 T cell", "CD8 T cell", "T cell", "Regulatory T cell","NK cell", "B cell", "Plasma cell",
                    "CD14 Monocyte", "CD16 Monocyte", "Macrophage", "Neutrophil",
                    "preDC", "pDC", "Platelet")

# Res/DX 그룹을 추가하여 새로운 벡터 생성
celltype_order_Res_DX <- c()

for (celltype in celltype_order) {
  celltype_order_Res_DX <- c(celltype_order_Res_DX, 
                             paste0(celltype, "_Res"), 
                             paste0(celltype, "_DX"))
}

plot_data <- plot_obj_8$data
colnames(plot_data) <- c("CellType", "Pathway", "NormValue") # 열 이름 명시

# DotPlot.metabolism에 쓰인 celltype별 정규화된 중앙값을 data.frame으로 정리.

norm_score_matrix <- plot_data %>%
  # CellType을 열 이름으로, NormValue를 값으로 피벗합니다.
  pivot_wider(names_from = CellType, 
              values_from = NormValue) %>%
  
  # Pathway 열을 행 이름으로 설정하고 Pathway 열 자체는 제거합니다.
  column_to_rownames(var = "Pathway")

norm_8_dotplot <- norm_score_matrix[rownames(nscore3), celltype_order]
# 77x22 data.frame인지 확인.
norm_8_sig <- trunc(10 * norm_8_dotplot)

norm_567_sig <- ifelse(norm_567_sig >= 0.9, 1, 0)

wilcox_77 <- wilcox_scMeta_pvalue[rownames(nscore3), celltype_order]
wilcox_77_8 <- wilcox_77

wilcox_77 <- ifelse(wilcox_77 < 0.001, 1, 0)

wilcox_77[wilcox_77 < 0.001] <- 3
wilcox_77[wilcox_77 < 0.01] <- 2
wilcox_77[wilcox_77 < 0.05] <- 1.01
wilcox_77[wilcox_77 <= 1] <- 0
wilcox_77[wilcox_77 == 1.01] <- 1

wilcox_norm_8 <- wilcox_77 * norm_8_sig

pheatmap(wilcox_77,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c(
           "white", "palegoldenrod", "orange"))(100),
         main = "Meaningful_***(pvalue) (Patient8)")

pheatmap(norm_3_sig,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c(
           "white", "palegoldenrod", "orange"))(100),
         main = "Significant(metabolic score >= 0.9) (Patient3)")

pheatmap(norm_8_sig,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c(
           "white", "palegoldenrod", "orange", "red"))(100),
         main = "Significant(metabolic score, 0.1 scale) (Patient8)")

pheatmap(wilcox_norm_8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c(
           "white", "palegoldenrod", "orange", "red"))(100),
         main = "Meaningful x Significant (metabolic score x p-value) (Patient8)")

wilcox_sqrt_all <- wilcox_heatmap
rm(wilcox_heatmap)
