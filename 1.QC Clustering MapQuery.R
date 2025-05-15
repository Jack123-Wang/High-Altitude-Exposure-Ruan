
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(dplyr)
library(viridis)
library(hrbrthemes)
library(Seurat)
library(ggforce)
library(magrittr)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(ggpubr)
# step1:导入数据 --------------------------------------------------------------
expression_matrix_1<-Read10X('B-17-1')
expression_matrix_3<-Read10X('B-39-1')
expression_matrix_5<-Read10X('B-4-1')
expression_matrix_7<-Read10X('B-5-1')
expression_matrix_9<-Read10X('QX-39-3')
expression_matrix_11<-Read10X('QX-4-3')
expression_matrix_13<-Read10X('QX-5-3')

cds_1 <- CreateSeuratObject(counts = expression_matrix_1,min.features = 150,min.cells = 10)
cds_3 <- CreateSeuratObject(counts = expression_matrix_3,min.features = 150,min.cells = 10)
cds_5 <- CreateSeuratObject(counts = expression_matrix_5,min.features = 150,min.cells = 10)
cds_7 <- CreateSeuratObject(counts = expression_matrix_7,min.features = 150,min.cells = 10)
cds_9 <- CreateSeuratObject(counts = expression_matrix_9,min.features = 150,min.cells = 10)
cds_11 <- CreateSeuratObject(counts = expression_matrix_11,min.features = 150,min.cells = 10)
cds_13 <- CreateSeuratObject(counts = expression_matrix_13,min.features = 150,min.cells = 10)

cds_1@meta.data$percent.mt<-PercentageFeatureSet(cds_1, pattern = "^MT-")
cds_3@meta.data$percent.mt<-PercentageFeatureSet(cds_3, pattern = "^MT-")
cds_5@meta.data$percent.mt<-PercentageFeatureSet(cds_5, pattern = "^MT-")
cds_7@meta.data$percent.mt<-PercentageFeatureSet(cds_7, pattern = "^MT-")
cds_9@meta.data$percent.mt<-PercentageFeatureSet(cds_9, pattern = "^MT-")
cds_11@meta.data$percent.mt<-PercentageFeatureSet(cds_11, pattern = "^MT-")
cds_13@meta.data$percent.mt<-PercentageFeatureSet(cds_13, pattern = "^MT-")

dim(cds_1@assays$RNA)
dim(cds_3@assays$RNA)
dim(cds_5@assays$RNA)
dim(cds_7@assays$RNA)
dim(cds_9@assays$RNA)
dim(cds_11@assays$RNA)
dim(cds_13@assays$RNA)

cds_1 <- subset(cds_1, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 10)
cds_3 <- subset(cds_3, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 10)
cds_5 <- subset(cds_5, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 10)
cds_7 <- subset(cds_7, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 10)
cds_9 <- subset(cds_9, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
cds_11 <- subset(cds_11, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)
cds_13 <- subset(cds_13, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 10)

# step2:合并 ----------------------------------------------------------------

cds.combined <- merge(cds_1, y = c(cds_3,cds_5,cds_7,
                                   cds_9,cds_11,cds_13), 
                      add.cell.ids = c("1-B","2-B","3-B","4-B",
                                       "2-HB","3-HB","4-HB"))

cds.combined@meta.data$type<-unlist(lapply(strsplit(rownames(cds.combined@meta.data),'_'),function(x) x[1]))
cds.combined@meta.data$type
# step2.1:质控 --------------------------------------------------------------
head(cds.combined@meta.data)

# 可视化每个样本的细胞计数
library(ggplot2)

cds.combined@meta.data %>% 
  ggplot(aes(x=type, fill=type)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# 可视化每个细胞的UMI/转录本数目
cds.combined@meta.data %>% 
  ggplot(aes(color=type, x=nCount_RNA, fill= type)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept =300)


# 通过频数图可视化每个细胞检测出的基因数分布
cds.combined@meta.data %>% 
  ggplot(aes(color=type, x=nFeature_RNA, fill= type)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = c(300,3000))

# 可视化检测到的基因数和UMI数之间的关系，并且观察是否存在大量低数目的基因数/UMI数的细胞
cds.combined@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = c(500,10000) )+
  geom_hline(yintercept = c(300,3000) )+
  facet_wrap(~type)

# step3:整合 ----------------------------------------------------------------

cds.combined$mitoRatio <- PercentageFeatureSet(cds.combined, pattern = "^MT-")
VlnPlot(object = cds.combined,features = "mitoRatio",pt.size = 0,raster=FALSE)

split_seurat <- SplitObject(cds.combined, split.by = "type")

split_seurat <- split_seurat[c("1-B","2-B","3-B","4-B","2-HB","3-HB","4-HB")]

#现在，对所有样品进行细胞周期评分和sctransform。这可能需要一些时间（〜10分钟）。
cc.genes.updated.2019
options(future.globals.maxSize = 1000 * 1024^2)
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio","S.Score","G2M.Score"))
}

##identify variable features for each dataset independently
split_seurat<-lapply(X=split_seurat, FUN=function(x){
  x<-FindVariableFeatures(x,verbose=FALSE)
})

# 选择变异最大的基因进行聚合
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 

#RPCA要求在整合前在每个单独的数据集上单独运行PCA
split_seurat<-lapply(X=split_seurat, FUN=function(x){
  x<-ScaleData(x,features=integ_features,verbose=FALSE)
  x<-RunPCA(x,features=integ_features,verbose=FALSE)
})

split_seurat<-PrepSCTIntegration(split_seurat,anchor.features = integ_features)

# 寻找最佳伙伴
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        reduction="rpca",dim=1:50)
# 跨条件聚合
filter_genes <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, features = integ_features)
  return(seurat_obj)
}

filtered_seurat_list <- lapply(split_seurat, filter_genes)

integ_anchors <- FindIntegrationAnchors(object.list = filtered_seurat_list, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features,
                                        reduction="rpca",dim=1:50)

seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
seurat_integrated <- ScaleData(seurat_integrated,verbose=FALSE)
seurat_integrated <- FindVariableFeatures(seurat_integrated, 
                                          selection.method = "vst", nfeatures = 3000) 
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")
PCAPlot(seurat_integrated,
        split.by = "type")  

# 绘制UMAP,观察批次效应是否缓解                             
DimPlot(seurat_integrated,
        split.by = "type")  
DimPlot(seurat_integrated,
        split.by = "people") 

# step4：聚类 ----------------------------------------------------------------

#聚类细胞
# 确定k-近邻图
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:40)

# 确定聚类的不同分辨率                               
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.5,0.8,1.5,2))

# 探索分辨率

table(seurat_integrated@meta.data$integrated_snn_res.0.8)

# 分配类群的身份
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"


# 绘制UMAP图
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size=3,raster=FALSE,group.by ="integrated_snn_res.0.8" )
DefaultAssay(seurat_integrated) <- "RNA"


# step5:映射 ----------------------------------------------------------------

reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
options(future.globals.maxSize = 5*1000 * 1024^2)

seurat_integrated<- PrepSCTFindMarkers(seurat_integrated)
seurat_integrated <- SCTransform(seurat_integrated, verbose = FALSE)


anchors <- FindTransferAnchors(
  reference = reference,
  query = seurat_integrated,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50)


seurat_integrated <- MapQuery(
  anchorset = anchors,
  query = seurat_integrated,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

seurat_integrated <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference,
  query = seurat_integrated, 
  new.reduction.name = "ref.spca"
)

seurat_integrated <- ProjectUMAP(
  query = seurat_integrated, 
  query.reduction = "ref.spca", 
  reference = reference, 
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)
p1 = DimPlot(seurat_integrated, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
p2 = DimPlot(seurat_integrated, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend()
p1 + p2



library(stringr)

metadata <- seurat_integrated@meta.data
names(metadata)
unique(metadata$type)
metadata$origin<-metadata$type
metadata$origin<-str_replace(metadata$origin,"1-B","plain_before")
metadata$origin<-str_replace(metadata$origin,"2-B","plain_before")
metadata$origin<-str_replace(metadata$origin,"3-B","plain_before")
metadata$origin<-str_replace(metadata$origin,"4-B","plain_before")
metadata$origin<-str_replace(metadata$origin,"2-HB","plateau_before")
metadata$origin<-str_replace(metadata$origin,"3-HB","plateau_before")
metadata$origin<-str_replace(metadata$origin,"4-HB","plateau_before")

metadata$people<-metadata$type
table(metadata$people)
unique(metadata$people)[2]

metadata$people<-str_replace(metadata$people,"2-B","2")
metadata$people<-str_replace(metadata$people,"1-B","1")
metadata$people<-str_replace(metadata$people,"3-B","3")
metadata$people<-str_replace(metadata$people,"4-B","4")
metadata$people<-str_replace(metadata$people,"2-HB","2")
metadata$people<-str_replace(metadata$people,"3-HB","3")
metadata$people<-str_replace(metadata$people,"4-HB","4")

metadata$celltype<-metadata$integrated_snn_res.0.8
metadata$celltype<-str_replace(metadata$celltype,"15","Neutrophils")
metadata$celltype<-str_replace(metadata$celltype,"14","Neutrophils")
metadata$celltype<-str_replace(metadata$celltype,"12","Neutrophils")
metadata$celltype<-str_replace(metadata$celltype,"0","Neutrophils")
metadata$celltype<-str_replace(metadata$celltype,"1","Neutrophils")
metadata$celltype<-str_replace(metadata$celltype,"5","Neutrophils")

metadata$predicted.celltype.l1[metadata$celltype=='Neutrophils'] <- 'Neutrophils'
metadata$predicted.celltype.l2[metadata$celltype=='Neutrophils'] <- 'Neutrophils'

metadata$celltype2<-metadata$integrated_snn_res.0.8
metadata$celltype2<-str_replace(metadata$celltype,"24","Platelet")

metadata$predicted.celltype.l1[metadata$celltype2=='Platelet'] <- 'Platelet'
metadata$predicted.celltype.l2[metadata$celltype2=='Platelet'] <- 'Platelet'


table(metadata$predicted.celltype.l1)

seurat_integrated@meta.data <-metadata

seurat_integrated1<-subset(x=seurat_integrated, predicted.celltype.l1=="other",invert=T)
seurat_integrated1<-subset(x=seurat_integrated1, predicted.celltype.l1=="other T",invert=T)
saveRDS(seurat_integrated1,'RDS/11_20.rds')

# step6：可视化 ---------------------------------------------------------------

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by="predicted.celltype.l1",
        label = T,split.by="origin",
        raster=F,
        alpha=0.3)

library(ggsci)
library(scales)

pal= pal_npg("nrc")(10)
cols =pal_jco("default")(8)
cols =pal_jco("default")(4)
cols =pal_lancet("lanonc")(9)	
cols =pal_d3("category10" )(9)
show_col(cols)
colors <- c("#ba642f", '#f17679','#e5422e','#3e9e47','#3f8baa','#ebad30','#70bf47','#b3976f')
cols <- c('#7AA6DCFF','#A73030FF')

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by="predicted.celltype.l1",cols = cols,
        label = T,
        raster=F,
        alpha=0.4)

DimPlot(seurat_integrated,
        reduction = "umap",
        group.by="people",cols = cols,
        label = F,
        raster=F,
        alpha=0.4,shuffle = T)

DimPlot(seurat_integrated,
        reduction = "umap",cols = cols,
        group.by="origin",
        label = F,
        raster=F,
        alpha=0.4,shuffle = T)

DimPlot(seurat_integrated1,
        reduction = "umap",
        group.by="integrated_snn_res.0.8",
        label = T,
        raster=F,
        alpha=0.3)


