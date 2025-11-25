file_list<-list.files("/home/welcome1/MDSfile/")
file_list
#sampleNames = mapply(function(x) x[length(x)], strsplit(file_list, split = '_'))
#sampleNames
meta.info = data.table(file_Path = file_list)
meta.info

meta.info[, sampleID:=mapply(function(x) x[length(x)], strsplit(file_list, split = 'x'))]
meta.info[, SubType:=mapply(function(x) x[length(x)], strsplit(sampleID, '-'))]
DT::datatable(meta.info)
#VlnPlot(dataMDSF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size=0, group.by = 'MDS_CON', raster=F)
#VlnPlot(dataMDSF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MT_cut = 25
MT_cut
batch_list = list()
#getwd()
setwd("/home/welcome1/MDSfile/")
getwd()

for(i in 1:nrow(meta.info)){
  dir_of_10X = meta.info$file_Path[i]
  sampleID = meta.info$sampleID[i]
  subtype = meta.info$SubType[i]
  
  meta_mat = Read10X(dir_of_10X)
  sc_aml_f <- CreateSeuratObject(counts = meta_mat, project = sampleID, min.cells = 3, min.features = 200)
  sc_aml_f[["percent.mt"]] <- PercentageFeatureSet(sc_aml_f, pattern = "^MT-")
  
  FeatureScatter(sc_aml_f, feature1 = "nCount_RNA", feature2 = "percent.mt")
  FeatureScatter(sc_aml_f, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  sc_aml_f <- subset(sc_aml_f, subset = nFeature_RNA > 200 & percent.mt <= MT_cut)
  
  sc_aml_f <- NormalizeData(sc_aml_f, verbose = F)
  sc_aml_f <- FindVariableFeatures(sc_aml_f, selection.method = "vst", 
                                   nfeatures = 2000, verbose = FALSE)
  sc_aml_f <- ScaleData(sc_aml_f, verbose = F)
  
  sc_aml_f <- RenameCells(object = sc_aml_f, add.cell.id = sampleID)  # add sample name as prefix
  
  ### add information
  sc_aml_f@meta.data[,'SubType'] = subtype
  sc_aml_f@meta.data[,'SampleID'] = sampleID
  
  ### add batch data into a list to merge
  batch_list <- append(batch_list, sc_aml_f)
  print(paste("Finishing processing", sampleID, 'at', Sys.time()))
}

#아래실패시 여기부터다시
MDS <- Reduce(function(x,y){merge(x,y,merge.data = TRUE)}, batch_list) 

MDS <- FindVariableFeatures(MDS, selection.method = "vst", 
                                    nfeatures = 2000, verbose = FALSE)

#MDS <- ScaleData(MDS, verbose = F)
#t <- merge(x=batch_list[1], y=batch_list[2], add.cell.ids = c("4K", "8K"), project = "PBMC12K")

#setwd("/home/yesung2023/")
#data_aml <- ScaleData(data_aml, vars.to.regress = c('nCount_RNA', 'percent.mt'),  verbose = F)

#DimPlot(data.combined[,data.combined$SubType==c('811', '832')], reduction='pca', raster=F)
#batch fastmnn 실험적 요인 제거 // 여기부터 
MDS <- RunFastMNN(object.list = SplitObject(MDS, split.by = "SampleID"))
##갑자기 RunFastMNN이 안되서 다음에 할 땐 RunHarmony로 바꿔야 할 듯

saveRDS(MDS_LNLD, "/home/yesung2023/MDS_LNLD.RDS")
table(MDS_LNLD3$SampleID)

#data_MDS@assays$RNA@scale.data<-as.matrix(data_MDS@assays$RNA@data)

?RunPCA
data.combined <- RunPCA(MDS, npcs=50, verbose=F)
data.combined
data.combined <- RunPCA(data.combined, npcs=50, verbose=F)
ElbowPlot(data.combined,ndims = 50)
JackStrawPlot(data.combined, dims = 1:50)

num_PC = 30
#20221019 여기부터 HSC resolution 높혀서 다시 시도
data_MDSF <- FindNeighbors(data_MDSF, reduction = 'mnn', dims=1:num_PC)
data_MDSF <- FindClusters(data_MDSF, resolution = 1.5)
# resolution 조절 
data_MDSF <- RunUMAP(data_MDSF, reduction = 'mnn', dims=1:num_PC)
DimPlot(dataMDSF, label = TRUE, raster = F, group.by = "RNA_snn_res.3")