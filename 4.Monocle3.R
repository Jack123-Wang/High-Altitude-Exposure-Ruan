# step1：Monocle3导入数据 ------------------------------------------------------
library(monocle3)
library(Seurat)
library(ggplot2)

seurat_integrated_Neutrophils<-readRDS("seurat_integrated1_Neutrophils.rds") #monocle3导入数据只需要Seurat导出的这个rds文件

data<-GetAssayData(seurat_integrated_Neutrophils,assay='RNA',slot='counts')
cell_metadata<-seurat_integrated_Neutrophils@meta.data
gene_annotation<-data.frame(gene_short_name=rownames(data))
rownames(gene_annotation)<-rownames(data)
cds<-new_cell_data_set(data,
                       cell_metadata=cell_metadata,
                       gene_metadata=gene_annotation)


# step2：Monocle3降维，整合Seurat数据 ---------------------------------------------
library(ggplot2)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")

colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP",
                 color_cells_by="RNA_snn_res.0.15",
                 show_trajectory_graph = F) + ggtitle('cds.umap')
p1
#plot_cells(cds,gene=c('TSPO','S100A11','CSF3R'),cell_size =0.5,graph_label_size = 1)
#plot_cells(cds,gene=c('NAMPT','FCGR3B','CXCR2'),cell_size =0.5,graph_label_size = 2)
#plot_cells(cds,gene=c('PLAC8','HEXA','PRSS57'),cell_size =0.5,graph_label_size = 2)
#plot_cells(cds,gene=c('CTSC','PRDX1','ITGA4'),cell_size =0.5,graph_label_size = 2)

##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat_integrated_Neutrophils, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP",
                 color_cells_by="subtype",
                 show_trajectory_graph = FALSE) + ggtitle('int.umap')
p2
p1|p2 


# step3：聚类，识别轨迹 -----------------------------------------------------------

## Monocle3聚类分区
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p1|p2

## 识别轨迹
cds <- learn_graph(cds)
#cds <- learn_graph(cds,verbose = T,learn_graph_control = list(minimal_branch_len=80))

p=plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = T, 
             label_branch_points = T)


##选起点
cds <- order_cells(cds)

p1=plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
              label_leaves = FALSE,  label_branch_points = FALSE,
              trajectory_graph_color = "white",
              cell_size=0.3) +
  theme_classic()+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))+
  scale_color_manual(values  = c("red","blue"))

pdata<-p1$data

ggplot(pdata,aes(x=data_dim_1,data_dim_2,color=cell_color))+
  geom_point(size=0.5)+
  labs(x="umap_1",y="umap_2",color="pseudotime",title="pseudotime")+
  scale_color_gradientn(colours  = c("#96b6e2","#dbe5f5","#f6c7c0","#eda69d","#ce564f"))+
  theme_classic()+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"),
        plot.title = element_text(hjust = 0.5)  # This centers the title
  )

ggplot(pdata,aes(x=subtype,y=cell_color,fill=subtype))+
  geom_violin(scale=T,width=1)+
  scale_fill_manual(values =c("#7aa6dc","#efc000","#cd534c","#8e7700"))+
  theme_classic()+
  labs(x="",y="pseudotime")+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))


#计算每个origin下的pse组成
pdata <- pdata %>%
  mutate(cell_color_category = case_when(
    cell_color >= 0 & cell_color < 10 ~ "pse 0-10",
    cell_color >= 10 & cell_color < 20 ~ "pse 10-20",
    cell_color >= 20 & cell_color <= 27 ~ "pse 20-30",
    TRUE ~ NA_character_  # 对于不在上述范围内的值，标记为 NA
  ))


# 计算每个 origin 分组下的 subtype 组成百分比
percent_data <- pdata %>%
  group_by(origin, cell_color_category) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

percent_data$origin<-str_replace_all(percent_data$origin,"plain_before","LA")
percent_data$origin<-str_replace_all(percent_data$origin,"plateau_before","HA")
percent_data$origin=factor(percent_data$origin,levels=c("LA","HA"))

# 使用 ggplot2 绘制柱形图
p1=ggplot(percent_data, aes(x = origin, y = percentage, fill = cell_color_category)) +
  geom_bar(stat = "identity", position = "stack",color="black") +
  labs(title = "          Pseudotime",
       x = "",
       y = "Percentage",
       fill = "Subtype") +
  scale_fill_manual(values =c("#7aa6dc","#efc000","#cd534c","#8e7700"))+
  theme_classic()+theme(axis.text.x = element_text(size=10,color="black"),
                        axis.text.y = element_text(size=10,color="black"))


##统计两种环境下的G0-3的组成
metadata<-seurat_integrated1_Neutrophils@meta.data
metadata<-metadata[c("subtype","origin")]


# 计算每个 origin 分组下的 subtype 组成百分比
percent_data <- metadata %>%
  group_by(origin, subtype) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100)

percent_data$origin<-str_replace_all(percent_data$origin,"plain_before","LA")
percent_data$origin<-str_replace_all(percent_data$origin,"plateau_before","HA")
percent_data$origin=factor(percent_data$origin,levels=c("LA","HA"))
# 使用 ggplot2 绘制柱形图
p2=ggplot(percent_data, aes(x = origin, y = percentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack",color="black") +
  labs(title = "          Subtype",
       x = "",
       y = "Percentage",
       fill = "Subtype") +
  scale_fill_manual(values =c("#7aa6dc","#efc000","#cd534c","#8e7700"))+
  theme_classic()+theme(axis.text.x = element_text(size=10,color="black"),
                        axis.text.y = element_text(size=10,color="black"))
p1|p2



