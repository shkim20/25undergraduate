library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(openxlsx)
library(readxl)

install.packages("Seurat")

.libPaths()

# 0. Seurat 파일 준비 및 list로 변환

batch_list = SplitObject(`MDS_Celltypes added`, split.by = "SampleID")
MDS_Res_DX <- subset(`MDS_Celltypes added`, subset = Sample_State %in% c("DX", "LNLD_Res"))
saveRDS(MDS_Res_DX, file = "MDS_Res_DX.RDS")

MDS_Res <- subset(`MDS_Celltypes added`, subset = Sample_State == "LNLD_Res")
saveRDS(MDS_Res, file = "MDS_Res.RDS")

MDS_DX <- subset(`MDS_Celltypes added`, subset = Sample_State == "DX")
saveRDS(MDS_DX, file = "MDS_DX.RDS")

MDS_Res_DX_patient3 <- subset(MDS_Res_DX, subset = SampleID %in% c("MDS-311", "MDS-335"))
saveRDS(MDS_Res_DX_patient3, file = "MDS_Res_DX_patient3.RDS")

MDS_Res_DX_patient8 <- subset(MDS_Res_DX, subset = SampleID %in% c("MDS-811", "MDS-832"))
saveRDS(MDS_Res_DX_patient8, file = "MDS_Res_DX_patient8.RDS")

# 리스트에 저장
MDS_Res_DX_list <- list()

Celltypes <- unique(`MDS_Celltypes added`@meta.data[["Celltypes"]])
for (ct in Celltypes) {
  MDS_Res_DX_list[[ct]] <- subset(MDS_Res_DX, subset = Celltypes == ct)
}

# Seurat으로 저장

'''
for (ct in celltypes) {
  assign(
    paste0("MDS_Res_DX_", gsub(" ", "_", ct)),  # 변수명 공백을 _로 변환
    subset(MDS_Res_DX, subset = Celltype == ct)
  )
}
'''

MDS_Res_DX_patient3_list <- list()

Celltypes <- unique(`MDS_Celltypes added`@meta.data[["Celltypes"]])
for (ct in Celltypes) {
  MDS_Res_DX_patient3_list[[ct]] <- subset(MDS_Res_DX_patient3, subset = Celltypes == ct)
}


# DEG 계산 코드


MDS_Res <- subset(`MDS_Celltypes added`, subset = Sample_State == "LNLD_Res")
# 1. Celltype별 DEG 계산
deg_list <- lapply(names(MDS_Res_DX_list), function(ct) {
  obj <- MDS_Res_DX_list[[ct]]
  
  # ident 설정: Sample_State 기준
  Idents(obj) <- obj@meta.data$Sample_State
  
  # LNLD_Res vs DX DEG 계산 (기존 설정 반영)
  FindMarkers(
    obj,
    group.by = "Sample_State",     # group.by는 meta.data 컬럼
    ident.1 = "LNLD_Res",
    ident.2 = "DX",
    logfc.threshold = 0.5,
    only.pos = FALSE,
    test.use = "bimod"
  )
})

# 2. 리스트 이름: Celltype 이름 공백 제거
# names(deg_list) <- gsub(" ", "_", names(MDS_Res_DX_list))
names(deg_list) <- names(MDS_Res_DX_list)

# 3. Excel 파일 생성
wb <- createWorkbook()

for (ct in names(deg_list)) {
  # 엑셀 시트에 쓰기 전에 데이터프레임 수정
  deg_df <- deg_list[[ct]]
  deg_df$gene <- rownames(deg_df)
  
  # 원하는 열 순서로 재정렬 (gene 열을 맨 앞으로)
  deg_df <- deg_df[, c("gene", setdiff(colnames(deg_df), "gene"))]
  
  addWorksheet(wb, ct)
  writeData(wb, sheet = ct, deg_df)
}

# 4. 저장
saveWorkbook(wb, file = "MDS_Res_DX_DEG.xlsx", overwrite = TRUE)

# deg_list 객체를 RDS 파일로 저장
saveRDS(deg_list, file = "MDS_Res_DX_DEG.rds")
