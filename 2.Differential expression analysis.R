
# step7:差异分析 --------------------------------------------------------------
seurat_integrated1=readRDS("seurat_integrated1.rds")

NK_single <- subset(seurat_integrated1, predicted.celltype.l1=="NK")
B_single <- subset(seurat_integrated1, predicted.celltype.l1=="B")
CD4T_single <- subset(seurat_integrated1, predicted.celltype.l1=="CD4 T")
CD8T_single <- subset(seurat_integrated1, predicted.celltype.l1=="CD8 T")
Mono_single<- subset(seurat_integrated1, predicted.celltype.l1=="Mono")
Neutrophils_single<- subset(seurat_integrated1, predicted.celltype.l1=="Neutrophils")
platelet_single<- subset(seurat_integrated1, predicted.celltype.l1=="Platelet")
DC_single <- subset(seurat_integrated1, predicted.celltype.l1=="DC")

nrow(NK_single@meta.data)
nrow(B_single@meta.data)
nrow(CD4T_single@meta.data)
nrow(CD8T_single@meta.data)
nrow(Mono_single@meta.data)
nrow(Neutrophils_single@meta.data)
nrow(platelet_single@meta.data)

table(seurat_integrated1@meta.data$predicted.celltype.l1)
DefaultAssay(NK_single) <- "RNA"
DefaultAssay(B_single) <- "RNA"
DefaultAssay(CD4T_single) <- "RNA"
DefaultAssay(CD8T_single) <- "RNA"
DefaultAssay(Mono_single) <- "RNA"
DefaultAssay(Neutrophils_single) <- "RNA"
DefaultAssay(platelet_single) <- "RNA"
DefaultAssay(DC_single) <- "RNA"


NK_single <- JoinLayers(NK_single) #seurat v5需要先跑joinlayer，结构有改变。
B_single <- JoinLayers(B_single) #seurat v5需要先跑joinlayer，结构有改变。
CD4T_single <- JoinLayers(CD4T_single) #seurat v5需要先跑joinlayer，结构有改变。
CD8T_single <- JoinLayers(CD8T_single) #seurat v5需要先跑joinlayer，结构有改变。
Mono_single <- JoinLayers(Mono_single)
Neutrophils_single <- JoinLayers(Neutrophils_single)
platelet_single <- JoinLayers(platelet_single)
DC_single <- JoinLayers(DC_single)

diff_DC <- FindMarkers(DC_single, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "origin",
                       ident.1 ="plateau_before",
                       ident.2="plain_before")
diff_DC$type<-"DC"
diff_DC$Symbol<-rownames(diff_DC)

diff_platelet <- FindMarkers(platelet_single, min.pct = 0.1, 
                             logfc.threshold = 0.1,
                             group.by = "origin",
                             ident.1 ="plateau_before",
                             ident.2="plain_before")
diff_platelet$type<-"platelet"
diff_platelet$Symbol<-rownames(diff_platelet)

diff_Neutrophils <- FindMarkers(Neutrophils_single, min.pct = 0.1, 
                                logfc.threshold = 0.1,
                                group.by = "origin",
                                ident.1 ="plateau_before",
                                ident.2="plain_before")
diff_Neutrophils$type<-"Neutrophils"
diff_Neutrophils$Symbol<-rownames(diff_Neutrophils)

diff_NK <- FindMarkers(NK_single, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "origin",
                       ident.1 ="plateau_before",
                       ident.2="plain_before")
diff_NK$type<-"NK"
diff_NK$Symbol<-rownames(diff_NK)



diff_B <- FindMarkers(B_single, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "origin",
                      ident.1 ="plateau_before",
                      ident.2="plain_before")
diff_B$type<-"B"
diff_B$Symbol<-rownames(diff_B)



diff_CD4T <- FindMarkers(CD4T_single, min.pct = 0.1, 
                         logfc.threshold = 0.1,
                         group.by = "origin",
                         ident.1 ="plateau_before",
                         ident.2="plain_before")
diff_CD4T$type<-"CD4T"
diff_CD4T$Symbol<-rownames(diff_CD4T)


diff_CD8T <- FindMarkers(CD8T_single, min.pct = 0.1, 
                         logfc.threshold = 0.1,
                         group.by = "origin",
                         ident.1 ="plateau_before",
                         ident.2="plain_before")
diff_CD8T$type<-"CD8T"
diff_CD8T$Symbol<-rownames(diff_CD8T)

diff_Mono <- FindMarkers(Mono_single, min.pct = 0.1, 
                         logfc.threshold = 0.1,
                         group.by = "origin",
                         ident.1 ="plateau_before",
                         ident.2="plain_before")
diff_Mono$type<-"Mono"
diff_Mono$Symbol<-rownames(diff_Mono)


total<-rbind(diff_Neutrophils,diff_NK,
             diff_B,
             diff_CD4T,
             diff_CD8T,diff_Mono,diff_platelet,diff_DC)

total<-filter(total,p_val<0.05)

table(total$type)
total_Up<-filter(total,avg_log2FC>0)
total_Down<-filter(total,avg_log2FC<0)
total_Up$Change<-"Up"
total_Down$Change<-"Down"
total<-rbind(total_Up,total_Down)
table(total_Up$type)
table(total_Down$type)


# step8：富集分析 --------------------------------------------------------------

library(org.Hs.eg.db)
library(clusterProfiler)
library(openxlsx)
library(DOSE)
library(ggplot2)
go_all_total<-NA
table(need$type)
table(total$type)
?enrichGO
for(i in 1:8){
  a<-filter(total,type==unique(total$type)[i])
  a_Up<-filter(a,Change =='Up')
  a_Down<-filter(a,Change =='Down')
  
  geneDEG_Up <- a_Up$Symbol
  geneDEG_Down <- a_Down$Symbol
  
  res_go_Up <- enrichGO(gene = geneDEG_Up,  #待富集的基因列表
                        OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                        keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                        qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                        readable = FALSE)
  
  res_go_Down <- enrichGO(gene = geneDEG_Down,  #待富集的基因列表
                          OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                          keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                          ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                          pAdjustMethod = 'fdr',  #指定 p 值校正方法
                          pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                          qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                          readable = FALSE)
  go_results_Up <- res_go_Up@result
  go_results_Down <- res_go_Down@result
  
  
  go_results_Up$cluster<-unique(total$type)[i]
  go_results_Down$cluster<-unique(total$type)[i]
  go_results_Up$type<-"Up"
  go_results_Down$type<-"Down"
  go_all_total<-rbind(go_all_total,go_results_Up,go_results_Down)
}
go_results_Up<-go_results_Up[c(1,2,5,9)]
go_results_Down<-go_results_Down[c(1,2,5,9)]
