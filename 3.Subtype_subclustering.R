library(stringr)
seurat_integrated1<-readRDS("seurat_integrated1.rds")
metadata<-seurat_integrated1@meta.data
library(Seurat)
DimPlot(seurat_integrated1,
        reduction = "umap",
        group.by="predicted.celltype.l1",
        label = TRUE,
        label.size = 4)

# step1:提取Neutrophils细胞 -----------------------------------------------------------
Idents(seurat_integrated1)<-"predicted.celltype.l1"
Idents(seurat_integrated1)
table(seurat_integrated1@meta.data$predicted.celltype.l1)
table(seurat_integrated1@meta.data$origin)

table(seurat_integrated1$predicted.celltype.l1)

seurat_integrated1_Neutrophils<-subset(seurat_integrated1, 
                                       cells = WhichCells(object = seurat_integrated1, idents = c("Neutrophils")))

table(seurat_integrated1_Neutrophils@meta.data$predicted.celltype.l2)

metadata<-seurat_integrated1_Neutrophils@meta.data
table(metadata$subtype)

#核对两组是否有差异
nrow(filter(metadata,predicted.celltype.l1=="NK"))
nrow(filter(seurat_integrated1_Neutrophils@meta.data,predicted.celltype.l1=="NK"))

# step2：分群 ----------------------------------------------------------------
seurat_integrated1_Neutrophils <- NormalizeData(seurat_integrated1_Neutrophils, normalization.method = "LogNormalize", scale.factor = 1e4) 
seurat_integrated1_Neutrophils <- FindVariableFeatures(seurat_integrated1_Neutrophils, selection.method = 'vst', nfeatures = 2000)
seurat_integrated1_Neutrophils <- ScaleData(seurat_integrated1_Neutrophils, vars.to.regress = "percent.mt")
seurat_integrated1_Neutrophils <- RunPCA(seurat_integrated1_Neutrophils, features = VariableFeatures(object = seurat_integrated1_Neutrophils)) 
seurat_integrated1_Neutrophils <- FindNeighbors(seurat_integrated1_Neutrophils, dims = 1:10)
seurat_integrated1_Neutrophils <- FindClusters(seurat_integrated1_Neutrophils, resolution = 0.15 )

metadata<-seurat_integrated1_Neutrophils@meta.data

table(seurat_integrated1_Neutrophils$seurat_clusters) 
seurat_integrated1_Neutrophils <- RunUMAP(seurat_integrated1_Neutrophils, dims = 1:10)

DimPlot(seurat_integrated1_Neutrophils, 
        reduction = 'umap',
        label = TRUE)


##Figure3A
library(ggsci)
cols =pal_jco("default")(9)
cols

DimPlot(seurat_integrated1_Neutrophils, 
        group.by="seurat_clusters",
        reduction = 'umap',
        label = TRUE,
        cols=c("#7aa6dc","#efc000","#cd534c","#8e7700","#003c67"),
        shuffle = T,seed=1)


# step3:鉴定每一群的marker基因 ----------------------------------------------------------------
seurat_integrated1_Neutrophils<-readRDS("seurat_integrated1_Neutrophils.rds")

diff_0 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "seurat_clusters",
                      ident.1 ="0")

diff_0$type<-"diff_0"
diff_0$Symbol<-rownames(diff_0)

diff_1 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "seurat_clusters",
                      ident.1 ="1")

diff_1$type<-"diff_1"
diff_1$Symbol<-rownames(diff_1)

diff_2 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "seurat_clusters",
                      ident.1 ="2")

diff_2$type<-"diff_2"
diff_2$Symbol<-rownames(diff_2)


diff_3 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "seurat_clusters",
                      ident.1 ="3")

diff_3$type<-"diff_3"
diff_3$Symbol<-rownames(diff_3)


diff_4 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                      logfc.threshold = 0.1,
                      group.by = "seurat_clusters",
                      ident.1 ="4")

diff_4$type<-"diff_4"
diff_4$Symbol<-rownames(diff_4)

#Neu共有marker
DotPlot(object = seurat_integrated1_Neutrophils, 
        features = c("S100A2","CXCR2","G0S2"),
        group.by = "seurat_clusters")

seurat_integrated1_Neutrophils@meta.data$seurat_clusters<-factor(seurat_integrated1_Neutrophils@meta.data$seurat_clusters,levels=rev(c("3","4","0","1","2")))

DotPlot(seurat_integrated1_Neutrophils, 
        features = c("RPL12","ITGA4","PLAC8","PRDX1","SOX4", #G0 共有
                     "MMP9","EGR1","PADI4","CATIP", #0-G1 
                     "RSAD2","IFIT1","ISG15", #1-G2
                     "CFD","FCER1G","S100P",#2-G3
                     "CXCR2","FCGR3B","CSF3R","NAMPT"),#G123共有
        group.by = "subtype") + 
  scale_size_continuous(range=c(0.5,6))+
  scale_color_gradientn(colors =c("#7aa6dc","#7aa6dc","white","#cd534c","#cd534c"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

DimPlot(seurat_integrated1_Neutrophils, 
        group.by="subtype",
        reduction = 'umap',
        label = TRUE,
        cols=c("#7aa6dc","#efc000","#cd534c","#8e7700","#003c67"),
        shuffle = T,seed=1)

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

# 使用 ggplot2 绘制柱形图
ggplot(percent_data, aes(x = origin, y = percentage, fill = subtype)) +
  geom_bar(stat = "identity", position = "stack",color="black") +
  labs(title = "          Subtype",
       x = "",
       y = "Percentage",
       fill = "Subtype") +
  scale_fill_manual(values =c("#7aa6dc","#efc000","#cd534c","#8e7700"))+
  theme_classic()+theme(axis.text.x = element_text(size=10,color="black"),
                        axis.text.y = element_text(size=10,color="black"))

# step4：每一群的差异分析及富集分析 ------------------------------------------------------

DefaultAssay(seurat_integrated1_Neutrophils) <- "RNA"
diff_G0 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "subtype",
                       ident.1 ="G0")

diff_G0$type<-"diff_G0"
diff_G0$Symbol<-rownames(diff_G0)



diff_G1 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "subtype",
                       ident.1 ="G1")

diff_G1$type<-"diff_G1"
diff_G1$Symbol<-rownames(diff_G1)


diff_G2 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "subtype",
                       ident.1 ="G2")

diff_G2$type<-"diff_G2"
diff_G2$Symbol<-rownames(diff_G2)


diff_G3 <- FindMarkers(seurat_integrated1_Neutrophils, min.pct = 0.1, 
                       logfc.threshold = 0.1,
                       group.by = "subtype",
                       ident.1 ="G3")

diff_G3$type<-"diff_G3"
diff_G3$Symbol<-rownames(diff_G3)

total<-rbind(diff_G0,diff_G1,diff_G2,diff_G3)
total<-filter(total,avg_log2FC>0&p_val<0.05)
table(total$type)

go_all_total<-NA
for(i in 1:4){
  a<-filter(total,type==unique(total$type)[i])
  a_Up<-filter(a,avg_log2FC>0)
  
  geneDEG_Up <- a_Up$Symbol
  
  res_go_Up <- enrichGO(gene = geneDEG_Up,  #待富集的基因列表
                        OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                        keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'BP',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                        qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                        readable = FALSE)
  
  go_results_Up <- res_go_Up@result
  
  go_results_Up<-go_results_Up[c(1,2,5,9)]
  go_results_Up$cluster<-unique(total$type)[i]
  go_results_Up$type<-"Up"
  go_all_total<-rbind(go_all_total,go_results_Up)
}

go_all_total<-as.data.frame(go_all_total)
go_all_total<-filter(go_all_total,pvalue<0.05)
table(go_all_total$cluster)
go_all_total$logP<- -log2(go_all_total$pvalue)


# step5:AddmoduleScore ----------------------------------------------------

genes<-read_excel("geneset.xlsx",sheet="Sheet1")
genes[is.na(genes)] <- "XXXXXX"

genes1 <- lapply(genes, function(column) {
  column[column != "" & column != "XXXXXX"]
})

seurat_integrated1_Neutrophils<-AddModuleScore(seurat_integrated1_Neutrophils,
                                               features = genes1)
#seurat_integrated1_Neutrophils@meta.data<-seurat_integrated1_Neutrophils@meta.data[-c(25:52)]
colnames(seurat_integrated1_Neutrophils@meta.data)
colnames(seurat_integrated1_Neutrophils@meta.data)[25:39]<-colnames(genes)

DotPlot(seurat_integrated1_Neutrophils, 
        features =rev(c("Cell Proliferation","NET formation","Azurophilic granule","ISG","Neutrophil Maturation","Specific granuly","Gelatinase granule","Secretory vesicle","Phagocytosis","Chemotaxis")),
        group.by = "subtype") + 
  scale_size_continuous(range=c(0.5,6))+
  scale_color_gradientn(colors =c("#7aa6dc","#7aa6dc","white","#cd534c","#cd534c"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+coord_flip()


# step4:细胞密度图 -------------------------------------------------------------
library(viridis)
library(ggpointdensity)
umap<-cbind(Embeddings(object=seurat_integrated1_Neutrophils@reductions[["umap"]]),FetchData(seurat_integrated1_Neutrophils,"subtype"),FetchData(seurat_integrated1_Neutrophils,"origin"))
umap_plain_before<-filter(umap,origin=="plain_before")
umap_plateau_before<-filter(umap,origin=="plateau_before")

p1<-ggplot(umap_plain_before,aes(x=umap_1,y=umap_2))+
  geom_pointdensity(size=0.5)+
  scale_color_gradientn(colors = c("#7aa6dc","#fafbb0","#cd534c"))+
  labs(title="Plain_before")+
  theme_classic()+
  theme(legend.position = 'none')+theme(axis.text.x = element_text(size=10,color="black"),
                                        axis.text.y = element_text(size=10,color="black"))
p3<-ggplot(umap_plateau_before,aes(x=umap_1,y=umap_2))+
  geom_pointdensity(size=0.5)+
  scale_color_gradientn(colors = c("#7aa6dc","#fafbb0","#cd534c"))+
  labs(title="Plateau_before")+
  theme_classic()+
  theme(legend.position = 'none')+theme(axis.text.x = element_text(size=10,color="black"),
                                        axis.text.y = element_text(size=10,color="black"))

library(patchwork)
p1+p3+plot_layout(ncol = 2)

saveRDS(seurat_integrated1_Neutrophils,"seurat_integrated1_Neutrophils.rds")

# step5：每个亚群在两种环境下的比较(富集分析) -----------------------------------------------------
seurat_integrated1_Neutrophils<-readRDS("seurat_integrated1_Neutrophils.rds")

Idents(seurat_integrated1_Neutrophils)<-"subtype"
Seurat_G0<-subset(seurat_integrated1_Neutrophils,subtype=="G0")
Seurat_G1<-subset(seurat_integrated1_Neutrophils,subtype=="G1")
Seurat_G2<-subset(seurat_integrated1_Neutrophils,subtype=="G2")
Seurat_G3<-subset(seurat_integrated1_Neutrophils,subtype=="G3")

table(seurat_integrated1_Neutrophils@meta.data$subtype)
table(Seurat_G0@meta.data$subtype)
table(Seurat_G1@meta.data$subtype)
table(Seurat_G2@meta.data$subtype)
table(Seurat_G3@meta.data$subtype)

DefaultAssay(Seurat_G0) <- "RNA"
DefaultAssay(Seurat_G1) <- "RNA"
DefaultAssay(Seurat_G2) <- "RNA"
DefaultAssay(Seurat_G3) <- "RNA"

Seurat_G0 <- JoinLayers(Seurat_G0) #seurat v5需要先跑joinlayer，结构有改变。
Seurat_G1 <- JoinLayers(Seurat_G1) #seurat v5需要先跑joinlayer，结构有改变。
Seurat_G2 <- JoinLayers(Seurat_G2) #seurat v5需要先跑joinlayer，结构有改变。
Seurat_G3 <- JoinLayers(Seurat_G3) #seurat v5需要先跑joinlayer，结构有改变。

cluster<-c("G0","G1","G2","G3")
test<-list(Seurat_G0,Seurat_G1,Seurat_G2,Seurat_G3)


total<-NA
for(i in 1:4){
  diff.plain<-FindMarkers(test[[i]], min.pct = 0.1, 
                          logfc.threshold = 0.1,
                          group.by = "origin",
                          ident.1 ="plateau_before",
                          ident.2="plain_before")
  diff.plain$type<-cluster[i]
  
  diff.plain$Symbol<-rownames(diff.plain)
  
  total<-rbind(total,diff.plain)
  print(cluster[i])
}

total<-filter(total,p_val<0.05)

write_xlsx(total,"Neu_subtype_diff.xlsx")

#统计上下调分别变化的基因
table(total$type)
total_Up<-filter(total,avg_log2FC>0)
total_Down<-filter(total,avg_log2FC<0)
total_Up$Change<-"Up"
total_Down$Change<-"Down"
total<-rbind(total_Up,total_Down)
table(total_Up$type)
table(total_Down$type)


library(clusterProfiler)
set.seed(12345)
library(org.Hs.eg.db)
#用enrichGO富集所有GO BP条目，挑选后呈现
go_all_total<-NA
for(i in 1:4){
  a<-filter(total,type==unique(total$type)[i])
  a_Up<-filter(a,avg_log2FC>0)
  a_Down<-filter(a,avg_log2FC<0)
  
  geneDEG_Up <- a_Up$Symbol
  geneDEG_Down <- a_Down$Symbol
  
  res_go_Up <- enrichGO(gene = geneDEG_Up,  #待富集的基因列表
                        OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                        keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                        ont = 'BP',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                        pAdjustMethod = 'fdr',  #指定 p 值校正方法
                        pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                        qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                        readable = FALSE)
  
  res_go_Down <- enrichGO(gene = geneDEG_Down,  #待富集的基因列表
                          OrgDb = org.Hs.eg.db,  #指定物种的基因数据库，示例物种是绵羊（sheep）
                          keyType = 'SYMBOL',  #指定给定的基因名称类型，例如这里以 entrze id 为例
                          ont = 'BP',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                          pAdjustMethod = 'fdr',  #指定 p 值校正方法
                          pvalueCutoff = 0.05,  #指定 p 值阈值（可指定 1 以输出全部）
                          qvalueCutoff = 1,  #指定 q 值阈值（可指定 1 以输出全部）
                          readable = FALSE)
  go_results_Up <- res_go_Up@result
  go_results_Down <- res_go_Down@result
  
  go_results_Up<-go_results_Up[c(1,2,5,9)]
  go_results_Down<-go_results_Down[c(1,2,5,9)]
  go_results_Up$cluster<-unique(total$type)[i]
  go_results_Down$cluster<-unique(total$type)[i]
  go_results_Up$type<-"Up"
  go_results_Down$type<-"Down"
  go_all_total<-rbind(go_all_total,go_results_Up,go_results_Down)
  print(i)
  print(unique(total$type)[i])
}

go_all_total<-as.data.frame(go_all_total)
go_all_total<-filter(go_all_total,pvalue<0.05)
table(go_all_total$cluster)
go_all_total$logP<- -log2(go_all_total$pvalue)
go_all_total$Cluster_type<-paste(go_all_total$cluster,go_all_total$type,sep="_")
library(reshape2)
wide_df <- dcast(go_all_total, ID + Description ~ Cluster_type, value.var = "logP")
wide_df_count <- dcast(go_all_total, ID + Description ~ Cluster_type, value.var = "Count")

head(wide_df)
wide_df_select<-wide_df[c(839,122,1214,1298,2549,390,
                          979,941,1446,
                          2273,234,151,
                          501,291,1252,492,27,102,1224),]
wide_df_select<-wide_df_select[c(2,4,6,8,10,3,5,7,9)]
wide_df_select[is.na(wide_df_select)]<-0

wide_df_count_select<-wide_df_count[c(839,122,1214,1298,2549,390,
                                      979,941,1446,
                                      2273,234,151,
                                      501,291,1252,492,27,102,1224),]
wide_df_count_select<-wide_df_count_select[c(2,4,6,8,10,3,5,7,9)]
wide_df_count_select[is.na(wide_df_count_select)]<-0

#下调变为负值
wide_df_select$G0_Down<- -wide_df_select$G0_Down
wide_df_select$G1_Down<- -wide_df_select$G1_Down
wide_df_select$G2_Down<- -wide_df_select$G2_Down
wide_df_select$G3_Down<- -wide_df_select$G3_Down
wide_df_select$Description<-str_to_title(wide_df_select$Description)


long_df_select<-melt(wide_df_select)
long_df_count_select<-melt(wide_df_count_select)

long_df_select$variable<-factor(long_df_select$variable,levels=colnames(wide_df_select)[2:9])
long_df_select$Description<-factor(long_df_select$Description,levels=rev(wide_df_select$Description))

long_df_select$Count<-long_df_count_select$value
colnames(long_df_select)[3]<-"LogP"


bk <- c(seq(-40,-0.1,by=0.1),seq(0,10,by=0.1))
bk

Up<-max(long_df_select$LogP)/(max(long_df_select$LogP)-min(long_df_select$LogP))
Down<- -min(long_df_select$LogP)/(max(long_df_select$LogP)-min(long_df_select$LogP))

long_df_select<-filter(long_df_select,Count>0)

ggplot(long_df_select,aes(x=variable,y=Description,size=Count,color=LogP))+
  geom_point()+
  scale_color_gradientn( colours= c(colorRampPalette(colors = c("#0073c2","#0073c2","#7aa6dc","white"))(length(bk)*Down),
                                    colorRampPalette(colors = c("white","#cd534c","#cd534c","#cd534c"))(length(bk)*Up)))+
  theme_bw()+
  labs(x="",y="")+  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1,size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))+
  scale_size (range=c (1,7)) +
  theme(panel.border = element_blank())+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

# step6：多组火山图 -------------------------------------------------------------------
test<-list(Seurat_G0,Seurat_G1,Seurat_G2,Seurat_G3)

total<-NA
for(i in 1:4){
  diff.plain<-FindMarkers(test[[i]], min.pct = 0.1, 
                          logfc.threshold = 0.1,
                          group.by = "origin",
                          ident.1 ="plateau_before",
                          ident.2="plain_before")
  diff.plain$type<-cluster[i]
  
  diff.plain$Symbol<-rownames(diff.plain)
  
  total<-rbind(total,diff.plain)
  print(cluster[i])
}

total<-na.omit(total)

#total
total$type<-factor(total$type,levels=unique(total$type))

total_Up<-filter(total,avg_log2FC>0,p_val<0.05)
total_Down<-filter(total,avg_log2FC<0,p_val<0.05)
total_Nonsig<-filter(total,p_val>0.05)
total_Up$Type<-"Up"
total_Down$Type<-"Down"
total_Nonsig$Type<-"Nonsig"

mycol <- c("#1bd66c","#e4d554","#d1362d","#8c559f")
dfcol<-data.frame(x=c(1:4),
                  y=0,
                  label=unique(total$type))

total_Nonsig$type<-factor(total_Nonsig$type,levels=c("G3","G2","G1","G0"))
total_Up$type<-factor(total_Up$type,levels=c("G3","G2","G1","G0"))
total_Down$type<-factor(total_Down$type,levels=c("G3","G2","G1","G0"))

p1=ggplot()+
  geom_jitter(data = total_Nonsig,
              aes(x = type, y = avg_log2FC),
              color = "grey",
              size = 0.3,
              width =0.4)+
  geom_jitter(data = total_Up,
              aes(x = type, y = avg_log2FC),
              color = "#cd534c",
              size = 0.3,
              width =0.4)+
  geom_jitter(data = total_Down,
              aes(x = type, y = avg_log2FC),
              color = "#0073c2",
              size = 0.3,
              width =0.4)+
  ylim(-2,2)+
  geom_tile(data = dfcol,
            aes(x=x,y=y),
            height=0.2,
            #color = "black",
            fill = mycol,
            alpha = 0.6,
            show.legend = F)+
  theme_classic()+
  labs(x="",y="Average log2FoldChange",title="")+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))+
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_blank(),  # 隐藏 Y 轴标签
        axis.ticks.y = element_blank(),  # 隐藏 Y 轴刻度
        axis.line.y = element_blank()) +  # 隐藏 Y 轴线
  coord_flip()


#统计上下调分别变化的基因
total<-filter(total,p_val<0.05)
table(total$type)
total_Up<-filter(total,avg_log2FC>0)
total_Down<-filter(total,avg_log2FC<0)
total_Up$Change<-"Up"
total_Down$Change<-"Down"
total<-rbind(total_Up,total_Down)
Num_Up<-as.data.frame(table(total_Up$type))
Num_Up$Type="Up"
Num_Down<-as.data.frame(table(total_Down$type))
Num_Down$Type="Down"
Num_Down
Num<-rbind(Num_Up,Num_Down)
Num$Var1<-factor(Num$Var1,levels=c("G3","G2","G1","G0"))

p2=ggplot(Num, aes(x = Freq, y = Var1, fill = Type)) +
  geom_bar(stat = "identity", position = "identity", aes(x = ifelse(Type == "Up", Freq, -Freq)),color="black") +
  labs(x = "Number of Differential genes", y = "Subtype", title = "") +
  theme_classic()+
  scale_fill_manual(values = c("#0073c2","#cd534c"))+
  theme(axis.text.x = element_text(size=10,color="black"),
        axis.text.y = element_text(size=10,color="black"))

p1|p2  

