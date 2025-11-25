library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(tidyverse)
library(writexl)
library(readxl)

MDS <- readRDS("/home/shared/Intern_LYS/MDS_SC.RDS")
MDS_save <-  MDS
DimPlot(MDS, reduction = "umap")
DimPlot(MDS, reduction = "umap") + NoLegend()
?RunFastMNN

batch_list = SplitObject(MDS, split.by = "SampleID") # 데이터 가독성을 위해

# 일단은 Seurat tutorial대로
top10 <- head(MDS@assays[["mnn.reconstructed"]]@var.features, 10)

plot1 <- VariableFeaturePlot(MDS@assays[["mnn.reconstructed"]]@meta.features)
#FindVariableFeatures를 해달랜다. 패스
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# all.genes <- rownames(MDS)
# MDS <- ScaleData(MDS, features = all.genes)
# MDS <- RunPCA(MDS, features = VariableFeatures(object = MDS))

# MDS@reductions와 MDS@commands를 보니 FastMNN, UMAP, FindNeighbors, FindClusters
# 가 실행되었다. PCA가 안 보이는데, UMAP를 실행하면 사라지나?
head(Idents(MDS), 5)
clusternames = levels(Idents(MDS))

# Finding DE

?FindMarkers
clusterCD8.markers <- FindMarkers(MDS, ident.1 = "CD8 T cell") # cluster 0~38 중 6, 2, 25
head(clusterCD8.markers, n = 5)
clusterCD8_2.markers <- FindMarkers(MDS, ident.1 = "CD8 T cell", ident.2 = c("CD4 T cell",
                                                                            "T cell"))
# ident.2를 넣었을 때의 결과 비교용
head(clusterCD8_2.markers, n = 5)
type(clusterCD8.markers)

MDS.markers <- FindAllMarkers(MDS, only.pos = TRUE)
MDS.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & pct.2 > 0.01) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)

# MDS@reductions에 PCA 정보가 없다.
MDS <- RunPCA(MDS, features = VariableFeatures(object = MDS))
# 할려고 했는데 ScaleData가 RunFastMNN 이후에 값이 날라갔는지
# Error in PrepDR(object = object, features = features, verbose = verbose) : 
# Data has not been scaled. Please run ScaleData and retry 라는 오류 출력.
# 선배님 안내: UMAP > t-SNE > PCA 순으로 상위 차원 축소 방법이라,
# UMAP를 실행했다면 PCA를 할 필요는 없다.

VlnPlot(MDS, features = c("GZMA", "CST7"), pt.size = 0) # slot = "counts" 추가 가능.

DimPlot(MDS, label=T, group.by = 'RNA_snn_res.1.5',raster=T)
DimPlot(MDS, label=T, raster=T)

View(MDS)

colnames(MDS@meta.data)[10]

?FindMarkers
#for(i in 1:NROW(ctypeGDM)){
aa1=FindMarkers(MDS, group.by="RNA_snn_res.1.5", 
                  ident.1 = '1', logfc.threshold = 0.5, only.pos=T)
head(aa1)
aa25=FindMarkers(MDS, group.by="RNA_snn_res.1.5", 
                ident.1 = '25', logfc.threshold = 0.5, only.pos=T)
head(aa21)
View(aa00)

aa0=aa0[order(aa0$avg_log2FC),]
aa0=aa0[aa0$p_val_adj<0.05,]
aa00=aa0[aa0$pct.1>0.75,]
aa00=aa00[aa00$pct.2<0.5,]
write.csv(aa00, "/home/welcome1/cluster0.csv")
View(aa00)


NROW(MDSmeta[MDSmeta$RNA_snn_res.1.5 %in% '0',])/NROW(MDSmeta)



DimPlot(MDS, label=T, raster=T)


# 같은 NK cell type 안인 '1', '21'도 marker가 조금씩 다른데, cell type을 어떻게 구별할까?
# 졸업논문 3.1에 따르면 < 10 DEG면 먼저 cluster를 합치고, cell type specific marker를
# 활용했다고 한다. 강의 자료에서 Marker 종류는 연구논문들 참고하면 나와 있다고 한다.
# 사이트 예시: https://dice-database.org/ , https://www.proteinatlas.org/ 등등

# Marker를 사이트에서 찾아 비교 중... 근데 가끔씩 cell type이 안 맞는 경우도 있음 ㅋ

# aa1[,6]=ctypeGDM[i,1]
# aa1[,7]=rownames(aa1)
# colnames(aa1)[6]="Celltype"
# colnames(aa1)[7]="Gene_Name"
# aa2=rbind(aa1, aa2)
#}

# print(MDS[["umap"]], dims = 1:5, nfeatures = 5)
# VizDimLoadings(MDS, dims = 1:2, reduction = "umap") 실행이 안 돼서 일시정지

?DimHeatmap
DimHeatmap(MDS, reduction = "mnn", dims = 1, cells = 500, balanced = TRUE)
# 얘도 pca 대신 umap이나 mnn 넣으면 오류 나온다.

FeaturePlot(MDS, features = rownames(clusterCD8.markers)[1:5])
FeaturePlot(MDS, features = rownames(aa1)[1:5])
VlnPlot(MDS, features = rownames(clusterCD8.markers)[1:5], pt.size = 0)
VlnPlot(MDS, features = rownames(aa25)[1:5], pt.size = 0, raster = F)

for (i in c(0:6)) {
  
  assign(paste0("aa", i), FindMarkers(MDS, group.by="RNA_snn_res.1.5", 
                                      ident.1 = as.character(i), logfc.threshold = 0.5, only.pos=T))
  
} 
# 앞에 적어놓은 cell type, marker에 대한 의문 다시. cluster 2-6 사이.
# 그리고 cluster 2, 6, 25 사이에서 regulatory T cell는 어떻게 분리해 냈지?

# Tabula Sapiens - Immune 또는 Azimuth 사용.
# Web에서 Azimuth를 사용하려면 H5 Seurat 형식으로 파일을 변환해야 한다.

# library(remotes)
# remotes::install_github("satijalab/azimuth") 에러 뜨고 안 된다...

Sys.which("make")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("TFBSTools", "DirichletMultinomial"))


aa1=aa1[order(desc(aa1$avg_log2FC)),]
aa1=aa1[order(desc(aa1$pct.1)),]
aa1=aa1[aa1$p_val_adj<0.05,]

# cluster1의 sample 1
aa11=aa1[aa1$pct.1 - aa1$pct.2 > 0.2, ]
aa11=aa11[aa11$pct.2 < 0.5, ]
aa11 = aa11[order(dplyr::desc(aa11$pct.1 - aa11$pct.2)), ] # 맨 위 30개 선별

# cluster1의 sample 2
aa12 = aa1[aa1$pct.1 / aa1$pct.2 > 2.5, ]
aa12 = aa12[aa12$pct.2 < 0.5, ]
aa12 = aa12[order(dplyr::desc(aa12$pct.1 - aa12$pct.2)), ] # 역시 30개 선별

write.csv(aa11, "/home/welcome1/cluster1-1.csv")
write.csv(aa12, "/home/welcome1/cluster1-2.csv")

for (i in c(35:38)) {
  
  aa <- FindMarkers(MDS, group.by="RNA_snn_res.1.5", ident.1 = c("5", "33", "34"), 
                    logfc.threshold = 0.5, only.pos=T)
  
  #aa1 <- read_excel(paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i, "-1.xlsx"))
  #aa2 <- read_excel(paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i, "-2.xlsx"))
  
  aa1=aa[aa$p_val_adj<0.05, ]
  aa2=aa[aa$p_val_adj<0.05, ]
  
  # 해당 cluster의 sample 1
  aa11 = aa1[aa1$pct.1 - aa1$pct.2 > 0.2, ]
  aa11 = aa11[aa11$pct.2 < 0.5, ]
  aa11 = aa11[order(dplyr::desc(aa11$pct.1 - aa11$pct.2)), ] 
  # xlsx에서 쓸 때는 맨 위 30개 선별
  
  # 해당 cluster의 sample 2
  aa12 = aa2[aa2$pct.1 / aa2$pct.2 > 2.5, ]
  aa12 = aa12[aa12$pct.2 < 0.5, ]
  aa12 = aa12[order(dplyr::desc(aa12$pct.1 - aa12$pct.2)), ] 
  # 역시 30개 선별
  
  write_xlsx(aa11, paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i ,"-1.xlsx"))
  write_xlsx(aa12, paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i, "-2.xlsx"))
}

colnames(MDSmeta)[7]
sub=subset(MDS, subset = RNA_snn_res.1.5 %in% c("1"))
saveRDS(sub, file = "/home/welcome1/cluster1.RDS")

for (i in c(0:38)) {
  
  sub = subset(MDS, subset = RNA_snn_res.1.5 %in% c(i))
  saveRDS(sub, file = paste0("/home/welcome1/cluster0~38 RDS files/cluster", 
                                        as.character(i),".RDS"))
}

View(sub@meta.data)

View(aa11)

for (i in c(1:30)) {
  
  if (rownames(aa11)[i] != rownames(aa12)[i]) {
    print(paste0(rownames(aa11)[i], " and ", rownames(aa12)[i]))
  }
}

?install.packages
install.packages("/home/welcome1/TFBSTools_1.44.0.tar.gz", repos = NULL, type = "source")

colnames(MDSmeta)
DimPlot(MDS, label=T, raster=T, group.by = "RNA_snn_res.3")

length(levels(MDSmeta$RNA_snn_res.1.5))
table(MDS$RNA_snn_res.1.5)
# cluster3의 gene 개수 aa1 - 3개, aa2 - 6개. cluster4는 aa1 - 1개, aa2 - 0개.
# 나머지는 모두 30개 이상.

for (i in c(0:33, 35:38)) {
  
  aa1 <- read_excel(paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i, "-1.xlsx"))
  aa2 <- read_excel(paste0("/home/welcome1/cluster0~38 xlsx files/cluster", i, "-2.xlsx"))
  
  if (nrow(aa1) < 30) print(paste(i, "aa1", nrow(aa1), nrow(aa2)) )
  if (nrow(aa2) < 30) print(paste(i, "aa2", nrow(aa2), nrow(aa1)) )
}

aa3 <- FindMarkers(MDS, group.by="RNA_snn_res.1.5", ident.1 = "3", 
                  logfc.threshold = 0.5, only.pos=T)
aa33=aa3[aa3$p_val_adj<0.05, ]

aa331 = aa33 %>%
  filter(pct.1 > 0.75 & pct.2 < 0.5) %>%
  arrange(desc(pct.1 - pct.2))

# 해당 cluster의 sample 2
aa332 = aa33 %>%
  filter(pct.1 / pct.2 > 2.5 & pct.2 < 0.5) %>%
  arrange(desc(pct.1 - pct.2))

i <- 5
