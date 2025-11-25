library(scMetabolism)
library(ggplot2)
library(rsvd)

# 1. 데이터 준비(LNLD_Res, DX를 함께 묶어 scMetabolism으로 분석해야 한다.)

MDS_Res_DX_scMeta<-sc.metabolism.Seurat(obj = MDS_Res_DX, method = "AUCell", 
                                        imputation = F, ncores = 2, metabolism.type = "KEGG")

MDS_Res_DX_Reactome<-sc.metabolism.Seurat(obj = MDS_Res_DX, method = "AUCell", 
                                          imputation = F, ncores = 2, metabolism.type = 
                                            "REACTOME")

scMeta_score = MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]]
Reactome_score = MDS_Res_DX_Reactome@assays[["METABOLISM"]][["score"]]

# 만약 모든 마침표(.)를 하이픈(-)으로 바꾸어도 문제가 없다면:
colnames(scMeta_score) <- gsub("\\.", "-", colnames(scMeta_score))

scMeta_score_Res <- scMeta_score[, colnames(scMeta_score) %in% rownames(MDS_Res@meta.data) ]
scMeta_score_DX <- scMeta_score[, colnames(scMeta_score) %in% rownames(MDS_DX@meta.data) ]

colnames(Reactome_score) <- gsub("\\.", "-", colnames(Reactome_score))

Reactome_score_Res <- Reactome_score[, colnames(Reactome_score) %in% rownames(MDS_Res@meta.data) ]
Reactome_score_DX <- Reactome_score[, colnames(Reactome_score) %in% rownames(MDS_DX@meta.data) ]

# scMetabolism 결과들을 다시 celltype에 따라 분류해야 한다...

# meta.data에서 각 셀의 celltype 정보 가져오기 (Res & DX)
celltype_info_Res <- MDS_Res@meta.data[colnames(scMeta_score_Res), "Celltypes"]
names(celltype_info_Res) <- rownames(MDS_Res@meta.data)

celltype_info_DX <- MDS_DX@meta.data[colnames(scMeta_score_DX), "Celltypes"]
names(celltype_info_DX) <- rownames(MDS_DX@meta.data)

# celltypes 순서대로 scMeta_score 열을 묶어 list로 저장

celltypes <- unique(MDS_Res@meta.data$Celltypes)

scMeta_split_Res <- lapply(celltypes, function(ct) {
  cols <- names(celltype_info_Res)[celltype_info_Res == ct]
  scMeta_score_Res[, cols, drop = FALSE]   # 각 celltype별 subset
})

# celltypes 이름을 list에 부여
names(scMeta_split_Res) <- celltypes

names(Reactome_split_Res) <- celltypes

# 2. Kolmogorov-Smirnov Test(KS test, KS 검정) 실시. 정규성 검정.

ks_pnorm_results_Res <- list()
ks_pnorm_results_DX <- list()

for (pathway in rownames(scMeta_score)) {
  
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

for (pathway in rownames(scMeta_score)) {
  
  for (celltype in celltypes) {
    
    if (ks_pnorm_results_Res[[pathway]][[celltype]][["p.value"]] >= 0.05) {
      
      cat(pathway, " - " , celltype, " : ", 
          ks_pnorm_results_Res[[pathway]][[celltype]][["p.value"]], "\n")
      
    }
    
  }
  
}

# 출력문이 없음 : 모든 pathway, celltype에 대해 p.value < 0.05 -> 정규성이 없음.
# 그러므로 독립표본 t-검정 대신 Wilcoxson Rank-Sum Test을 실시.


# 3. Wilcoxson Rank-Sum Test 실시.

wilcox_scMeta <- list()

for (pathway in rownames(scMeta_score)) {
  
  wilcox_scMeta[[pathway]] <- list()
  
  for (celltype in celltypes) {
    
    Res_vector <- as.numeric(as.matrix(scMeta_split_Res[[celltype]][pathway, ]))
    DX_vector <- as.numeric(as.matrix(scMeta_split_DX[[celltype]][pathway, ]))
    
    wilcox_scMeta[[pathway]][[celltype]] <- wilcox.test(Res_vector, DX_vector, 
                                                         alternative = "two.sided", paired = FALSE)
    
  }
  
}

for (pathway in rownames(scMeta_score)) {
  
  for (celltype in celltypes)
    
    if (wilcox_scMeta[[pathway]][[celltype]][["p.value"]] >= 0.05) {
      
      cat(pathway, " - ", celltype, " : ", 
          wilcox_scMeta[[pathway]][[celltype]][["p.value"]], "\n")
    }
  
}

wilcox_scMeta_pvalue <- as.data.frame(
  t(sapply(wilcox_scMeta, function(x) {
    sapply(x, function(y) y[["p.value"]])
  }))
)

# for-loop 버전 :

# # pathway × celltype 크기 지정
# pval_matrix <- matrix(
#   NA, 
#   nrow = length(rownames(scMeta_score)), 
#   ncol = length(celltypes),
#   dimnames = list(rownames(scMeta_score), celltypes)
# )
# 
# for (pathway in rownames(scMeta_score)) {
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

# pheatmap으로 heatmap 그리기

# 1. -(p value)로 변환
wilcox_scMeta_heatmap <- -(wilcox_scMeta_heatmap)
pheatmap(wilcox_scMeta_heatmap)

# 1-1. -(p value)/4로 변환
wilcox_scMeta_heatmap <- -(wilcox_scMeta_heatmap)/4
pheatmap(wilcox_scMeta_heatmap)

# p value >= 0.05인 값들은 중요하지 않으니 일괄적으로 회색으로 처리하라 하셨다.
# 이에 값을 0.05 미만은 모두 빨간색, 0.05 이상은 회색으로 표시할 것이다.

# 만약 0.01로 기준을 높여본다면?

# 1) 이진화: 0.01 미만 → 1, 이상 → 0
wilcox_scMeta_heatmap <- ifelse(wilcox_scMeta_heatmap < 0.01, 1, 0)

# 2) 색상 정의 (0 = 회색, 1 = 빨강)
annotation_colors <- colorRampPalette(c("grey80", "red"))(2)

# 3) pheatmap 실행
pheatmap(
  wilcox_scMeta_heatmap,
  color = annotation_colors,
  cluster_rows = FALSE,
  legend_breaks = c(0, 1),
  legend_labels = c("p >= 0.01", "p < 0.01"),
  main = "Pathways displaying difference (Wilcoxon p-values < 0.01)"
)


# 2. -log2(p value)로 변환. 망
?pheatmap
pheatmap(wilcox_scMeta_heatmap)
wilcox_scMeta_heatmap <- -log(wilcox_scMeta_heatmap, base = 2)

# log10을 씌워도 값이 0.1에서 200까지도 튄다. 상한값/하한값 조정이 필요.

wilcox_scMeta_heatmap[wilcox_scMeta_heatmap > 50] <- 50

# 그래도 안 된다.

# 3. -sqrt(log2(p value))로 변환.
wilcox_scMeta_heatmap <- -(log(wilcox_scMeta_heatmap, base = 2))
wilcox_scMeta_heatmap <- sqrt(wilcox_scMeta_heatmap)

pheatmap(wilcox_scMeta_heatmap)

# 상한값 조정 : 25를 10까지
wilcox_scMeta_heatmap[wilcox_scMeta_heatmap > 10] <- 10
pheatmap(wilcox_scMeta_heatmap)

pheatmap_sqrt <- pheatmap(wilcox_scMeta_heatmap)

wilcox_scMeta_heatmap <- wilcox_scMeta_heatmap %>%
  arrange(desc(n))


# pheatmap에서 붉은색으로 나타난 pathway들.
rownames(wilcox_scMeta_heatmap)[c(1:11, 78:84)]
# 또는 rownames(wilcox_scMeta_heatmap)[c(1:11, 36:38, 40, 42:53, 78:84)]
rowseq <- pheatmap_sqrt[["tree_row"]][["order"]][c(1:11, 78:84)]
pathway_sqrt <- rownames(wilcox_scMeta_heatmap)[rowseq]
pathway_sqrt


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

DotPlot.metabolism(obj = MDS_Res_DX_scMeta, pathway = 
                     rownames(MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]]), 
                   phenotype = "Celltypes", norm = "y") + 
  scale_color_gradientn(colors = c( "white", "red"))

# celltype 너무 많아서 실패. 정제 필요.

pvalue_celltype <- wilcox_scMeta_pvalue[85, 1:22]
pvalue_pathway <- wilcox_scMeta_pvalue[1:84, 23, drop = F]

col_order <- order(as.numeric(pvalue_celltype[1, ]), decreasing = TRUE)
pvalue_celltype <- pvalue_celltype[, col_order]
pvalue_celltype <- as.data.frame(t(pvalue_celltype))
colnames(pvalue_celltype) <-  "pvalue"

pvalue_pathway <- pvalue_pathway %>%
  arrange(desc(pvalue))
rownames(pvalue_celltype)

pheatmap(wilcox_scMeta_heatmap[rownames(pvalue_pathway), rownames(pvalue_celltype)],
         cluster_rows = FALSE, cluster_cols = FALSE)

print(rownames(pvalue_pathway)[1:42])
print(rownames(pvalue_celltype)[1:12])
phenotype_12 <- c("Early Erythroid_Res", "Early Erythroid_DX", "ERP_Res", "ERP_DX",
                  "CD14 Monocyte_Res", "CD14 Monocyte_DX", "CD4 T cell_Res", "CD4 T cell_DX",
                  "Late Erythroid1_Res", "Late Erythroid1_DX", "HSC_Res", "HSC_DX",
                  "NK cell_Res", "NK cell_DX", "CD8 T cell_Res", "CD8 T cell_DX",
                  "T cell_Res", "T cell_DX", "LMPP_Res", "LMPP_DX", "Late Erythroid2_Res",
                  "Late Erythroid2_DX", "GMP_Res", "GMP_DX")

# phenotype_12 <- factor(phenotype_12, levels = phenotype_12)

MDS_Res_DX_scMeta <- MDS_Res_DX_scMeta_save
Metabolism <- MDS_Res_DX_scMeta@assays[["METABOLISM"]]

colnames(MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]])[1:10] # cell 이름

rownames(MDS_Res_DX_scMeta@meta.data)[1:10] # 얘도 cell 이름

MDS_Res_DX_scMeta_24 <- subset(
  MDS_Res_DX_scMeta, 
  subset = Celltypes %in% phenotype_12 # object$column_name 대신 subset 인자를 사용
)

# [["METABOLISM"]]에서 cell 이름의 모든 마침표(.)를 하이픈(-)으로 바꾸기:
colnames(MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]]) <- 
  gsub("\\.", "-", colnames(MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]]))

rownames(MDS_Res_DX_scMeta_24@meta.data)[1:10]

# 1. list 구조에서 score 매트릭스를 직접 추출
metabolism_score_matrix <- MDS_Res_DX_scMeta@assays[["METABOLISM"]][["score"]]

# 2. 필터링할 cell ID 리스트 정의 (subset된 객체의 cell names)
cells_to_keep <- colnames(MDS_Res_DX_scMeta_24)

# 3. score 행렬에서 해당 cell ID에 해당하는 열만 추출
Metabolism_filtered_matrix <- metabolism_score_matrix[, cells_to_keep, drop = FALSE]

# 4. 필터링된 행렬을 원래의 "score" key를 가진 list로 복원
new_metabolism_list <- list(score = Metabolism_filtered_matrix)

# 5. MDS_Res_DX_scMeta_24 객체의 assays 슬롯에 복원된 list 할당
MDS_Res_DX_scMeta_24@assays[["METABOLISM"]] <- new_metabolism_list

MDS_Res_DX_scMeta_24@meta.data$Celltypes <- factor(
  MDS_Res_DX_scMeta_24@meta.data$Celltypes,
  levels = phenotype_12
)

# Early Erythroid ~ GMP, 12개 celltype & 40 + 2개 pathway 기준
DotPlot.metabolism(obj = MDS_Res_DX_scMeta_24, pathway = 
                     rownames(pvalue_pathway)[1:42], 
                   phenotype = "Celltypes", norm = "y") + 
  scale_color_gradientn(colors = c("white", "red")) +
  labs(color = "Metabolic score", size = "Metabolic score")

table(MDS_Res_DX_scMeta@meta.data$Celltypes)

trace(DotPlot.metabolism)
