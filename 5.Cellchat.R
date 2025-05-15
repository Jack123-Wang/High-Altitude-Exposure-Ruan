options(timeout=1000)
library(CellChat)
library(patchwork)
library(Seurat)
# step1：读取数据 ----------------------------------------------------

getwd()
cellChat.Plateau<-readRDS("cellChat.plateau_b.rds")
cellChat.Plain<-readRDS("cellChat.plain_b.rds")

df.plateau <- subsetCommunication(cellChat.Plateau)
df.plain <- subsetCommunication(cellChat.Plain)

df.netp.plateau <- subsetCommunication(cellChat.Plateau, slot.name = "netP")
df.netp.plain <- subsetCommunication(cellChat.Plain, slot.name = "netP")

object.list<-list(LA=cellChat.Plain,HA=cellChat.Plateau)

cellChat<-mergeCellChat(object.list,add.names = names(object.list))

levels(cellChat@idents) # show factor levels of the cell labels
cellChat@idents

# Part1:预测细胞间通讯的一般原理 ------------------------------------------------------

#比较交互总数和交互强度
gg1 <- compareInteractions(cellChat, show.legend = F, group = c(1))
gg2 <- compareInteractions(cellChat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

p=compareInteractions(cellChat, show.legend = F, group = c(2), measure = "weight")
bar=p[["data"]]

# 绘制竖向柱形图
library(ggplot2)
ggplot(bar, aes(x = dataset, y = count, fill = dataset)) +
  geom_col(color="black") +  # 使用geom_col()绘制竖向柱形图
  labs(title = "", x = "", y = "Interaction strength") +  # 添加标题和坐标轴标签
  theme_classic2() +  # 使用简约主题
  geom_text(aes(label = round(count, 3)), vjust = -0.5, color = "black", size = 3.5)+
  scale_fill_manual(values = c("LA" = "#7aa6dc", "HA" = "#cd534c"))+  # 自定义颜色
  theme(
    axis.text.x = element_text(color = "black", size = 10),  # 设置x轴文本为黑色，大小为10
    axis.text.y = element_text(color = "black", size = 10),  # 设置y轴文本为黑色，大小为10
    axis.title.y = element_text(color = "black", size = 12),  # 设置y轴标题为黑色，大小为10
    legend.position = "none"  # 去除图例
  )

#不同细胞群之间的相互作用数量或强度的差异
#两个数据集之间的细胞-细胞通信网络中相互作用的差异数量或差异作用强度可以使用圆图来可视化
#其中红色的红色的（或者蓝色的) 彩色边缘代表与第一个数据集相比第二个数据集中的信号增加（或者减少) 。
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellChat, weight.scale = T, measure = "weight",
                          color.use = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                        "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"))

#分别呈现两组的内部的交流情况
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2],
                   color.use = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                 "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"),
                   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#unloadNamespace("leiden")
#library(CellChat)
#用热图在更大的细节中显示交互的差异数或交互强度
#顶部彩色条形图表示热图（传入信号）中显示的列值的总和。
#右边的彩色条形图表示一行值（传出信号）的总和。
#在色条中红色或蓝色表示第二个数据集中与第一个数据集相比增加或[减少]信号。

netVisual_heatmap(cellChat, measure = "weight")

#上图色标有问题，自己e定义颜色断点和对应的颜色
breaks <- c(-0.1, 0, 0.1, 0.2, 0.3)
colors <- c("#2166ac", "#f7f7f7", "#e9ada8", "#c95754", "#b2182b")

# 创建颜色渐变函数
color_palette <- colorRampPalette(colors)

# 生成更细致的颜色梯度（例如100个颜色）
n_colors <- 100
gradient_colors <- color_palette(n_colors)

# 创建色标图
library(ggplot2)

df <- data.frame(x = seq(-0.1, 0.3, length.out = n_colors),
                 y = 1,
                 color = gradient_colors)

ggplot(df, aes(x = x, y = y, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  scale_x_continuous(breaks = breaks) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Value", title = "Color Gradient from -0.1 to 0.3")



cellChat<-updateCellChat(cellChat)

## 检查每种细胞发出的信号 （每种细胞和其他细胞的互作情况）
mat <- cellChat@net$Plain$count
mat <- cellChat@net$Plateau$count

Idents(cellChat) <- 'predicted.celltype.l1'
cellChat@
  table(cellChat@idents)
par(mfrow =c(4,3),xpd=T)
groupSize <- as.numeric(table(cellChat@idents))

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   arrow.width = 0.2,
                   arrow.size = 0.1, 
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}


#比较二维空间中的主要源和目标
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()

for (i in 1:length(object.list)) {
  object.list[[i]]<-netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i],
                                               weight.MinMax = weight.MinMax,
  )+ylim(0, 2) +xlim(0,2.5)
}
patchwork::wrap_plots(plots = gg)
p=gg
LA=p[[1]][["data"]]
HA=p[[2]][["data"]]
LA$Type="LA"
HA$Type="HA"
LAHA<-rbind(LA,HA)

ggplot()+
  geom_point(LA,mapping=aes(x=x,y=y,color=labels,size=Count),alpha=0.3)+
  geom_point(HA,mapping=aes(x=x,y=y,color=labels,size=Count),alpha=1)+
  geom_text(HA,mapping=aes(x=x,y=y,label=labels), vjust = -0.7,alpha=1)+
  scale_color_manual(values = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"))+  # 自定义颜色
  theme_classic()+labs(x="Outgoing intersection strength",y="Incoming interaction strength")+
  geom_segment(aes(x=0.13730, y=0.2205, xend=0.22062, yend=0.3565), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=0.21239, y=0.2235, xend=0.43837, yend=0.2515), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=0.72033, y=0.3284, xend=1.28425, yend=0.4384), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=1.52491, y=1.0132, xend=1.71406, yend=1.3773), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=0.19348, y=0.9918, xend=0.47611, yend=1.6795), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=0.43267, y=0.2698, xend=0.49866, yend=0.3244), linetype="dashed", color="black") +  # 添加虚线
  geom_segment(aes(x=0.04306, y=0.2170, xend=0.04693, yend=0.2513), linetype="dashed", color="black") +  # 添加虚线
  
  theme(
    axis.text.x = element_text(color = "black", size = 10),  # 设置x轴文本为黑色，大小为10
    axis.text.y = element_text(color = "black", size = 10))  # 设置y轴文本为黑色，大小为10

#x-axis and y-axis are respectively the total outgoing or incoming communication probability associated with each cell group. 
#Dot size is proportional to the number of inferred links (both outgoing and incoming) associated with each cell group.

# Part2:识别保守和特定环境的信号通路 ----------------------------------------------------
#根据功能相似性识别信号组
cellChat <- computeNetSimilarityPairwise(cellChat, type = "functional",comparison = NULL)
?computeNetSimilarityPairwise

library(uwot)
cellChat <- netEmbedding(cellChat, type = "functional",umap.method="uwot")
cellChat <- netClustering(cellChat, type = "functional",do.parallel = F)
netVisual_embeddingPairwise(cellChat, type = "functional", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellChat, type = "functional", nCol = 2)

#根据结构相似性识别信号组
cellChat <- computeNetSimilarityPairwise(cellChat, type = "structural")
cellChat <- netEmbedding(cellChat, type = "structural",umap.method="uwot")
cellChat <- netClustering(cellChat, type = "structural",do.parallel = F)
netVisual_embeddingPairwise(cellChat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellChat, type = "structural", nCol = 2)


#计算并可视化学习的联合流形中的路径距离
rankSimilarity(cellChat, type = "functional")
rankSimilarity(cellChat, type = "structural",comparison2 = c(1,2))

#比较每个信号通路的整体信号流
gg1 <- rankNet(cellChat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellChat, mode = "comparison", stacked = F, do.stat = TRUE)

gg1+gg2
gg1
gg2

gg1_data=gg1[["data"]]

percent_data <- gg1_data %>%
  group_by(name) %>%
  mutate(
    total = sum(contribution.scaled),
    LA_percent = ifelse(group == "LA", (contribution.scaled / total) * 100, NA),
    HA_percent = ifelse(group == "HA", (contribution.scaled / total) * 100, NA)
  ) %>%
  ungroup()

percent_data <- percent_data %>%
  mutate(percent = ifelse(group == "LA", LA_percent, HA_percent)) 

percent_data <- percent_data %>%
  group_by(name) %>%
  mutate(
    significant = ifelse(pvalues < 0.05, TRUE, FALSE),
    color = ifelse(significant, 
                   ifelse(contribution.scaled[group == "HA"] > contribution.scaled[group == "LA"], "#cd534c", "#7aa6dc"), 
                   "black")
  ) %>%
  ungroup()

gg1=ggplot(percent_data, aes(x = name, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "stack",color="black") +  # 绘制堆叠条形图
  coord_flip() +  # 转换坐标轴，使条形图水平显示
  labs(x = "", y = "Relative information flow") +  # 添加坐标轴标签
  theme_classic() +  # 使用简约主题
  theme(
    axis.text.y = element_text(size = 10,color = percent_data$color),
    axis.text.x = element_text(size = 10,color='black'), # 设置y轴文本大小
    legend.position = "none"  # 设置图例位置
  ) +
  scale_color_manual(values = c( "#7aa6dc", "#cd534c"))+  # 自定义颜色
  scale_fill_manual(values = c("LA" = "#7aa6dc", "HA" = "#cd534c"))  # 自定义颜色

gg1

gg2=ggplot(percent_data, aes(x = name, y = contribution.scaled, fill = group)) +
  geom_bar(stat = "identity", position = "dodge",color="black") +  # 绘制堆叠条形图
  coord_flip() +  # 转换坐标轴，使条形图水平显示
  labs(x = "", y = "Relative information flow") +  # 添加坐标轴标签
  theme_classic() +  # 使用简约主题
  theme(
    axis.text.y = element_text(size = 10,color = percent_data$color),
    axis.text.x = element_text(size = 10,color='black'), # 设置y轴文本大小
    legend.position = "right"  # 设置图例位置
  ) +
  scale_fill_manual(values = c("LA" = "#7aa6dc", "HA" = "#cd534c"))  # 自定义颜色

gg1+gg2

# Part 4:使用层次图、圈图或者弦图直观地比较细胞间通信 -------------------------------------------

# 观察指定信号通路下，KO、WT的细胞传出或传入信号是否发生了改变？
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = c(1,2,3,4) # a numeric vector. 
par(mfrow = c(1,2), xpd=TRUE)

netVisual_aggregate(object.list[[1]], signaling = "CXCL", layout = "circle",
                    color.use = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                  "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"))

netVisual_aggregate(object.list[[2]], signaling = "CXCL", layout = "circle",
                    color.use = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                  "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"))

netVisual_aggregate(object.list[[1]], signaling = "BTLA", layout = "circle") #没有
netVisual_aggregate(object.list[[2]], signaling = "BTLA", layout = "circle",
                    color.use = c("B" = "#0073c2", "CD4 T" = "#efc000","CD8 T"="#868686","DC"="#cd534c",
                                  "Mono"="#7aa6dc","Neutrophils"="#003c67","NK"="#8f7700","Platelet"="#3b3b3b"))

netAnalysis_contribution(object.list[[1]], signaling = "CXCL")
netAnalysis_contribution(object.list[[2]], signaling = "CXCL")
netAnalysis_contribution(object.list[[1]], signaling = "BTLA")
netAnalysis_contribution(object.list[[2]], signaling = "BTLA")


#观察指定信号通路下，每个细胞类型在KO、WT的身份（Sender、Receiver、Mediator、Influencer）
# Compute the network centrality scores
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]],
                                                  slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]],
                                                  slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(object.list[[1]],
                                  signaling = "CXCL", 
                                  width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list[[2]],
                                  signaling = "CXCL", 
                                  width = 8, height = 2.5, font.size = 10)

netAnalysis_signalingRole_network(object.list[[1]],
                                  signaling = "BTLA", 
                                  width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(object.list[[2]],
                                  signaling = "BTLA", 
                                  width = 8, height = 2.5, font.size = 10)


# Part5:比较不同数据集之间地信号基因表达分布 ------------------------------------------------

cellChat@meta$datasets = factor(cellChat@meta$datasets, 
                                levels = c("LA", "HA")) # set factor level



p=plotGeneExpression(cellChat,
                     signaling = "CXCL",
                     split.by = "datasets", 
                     colors.ggplot = T)

CXCL=p[[1]][[1]][["data"]]
CXCL$CXCR1=p[[2]][[1]][["data"]]$CXCR1
CXCL$CXCR2=p[[3]][[1]][["data"]]$CXCR2

average_values= CXCL %>%
  group_by(ident,split) %>%
  summarise(avg_CXCL8 = median(CXCL8))

p1=ggplot(CXCL,aes(x=split,y=CXCL8,fill=split))+
  geom_line(average_values,mapping=aes(x=split,y=avg_CXCL8,group=ident))+
  geom_violin(scale = 'width')+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5)+
  stat_summary(mapping=aes(group=split),                      ##设置分组列
               fun="median",                                   ##箱线图添加均值使用的公式mean
               geom="point",size=3,fill="white",    ##均值图形的设置
               position=position_dodge(0.8))+
  stat_compare_means(
    aes(label = ..p.signif..),  # 显示 * 号（而不是 p 值）
    method = "wilcox.test",     # 使用 Wilcoxon 秩和检验
    comparisons = list(c("LA", "HA")),  # 替换成你的分组名称
    size = 4,
    vjust = 0.5
  ) +
  theme_classic()+facet_wrap(ident ~ .,scales = "free",nrow = 1)+scale_fill_manual(values = c('#7aa6dc',"#cd534c"))+
  theme(axis.text.x = element_text(colour = 'black',size=10),axis.text.y = element_text(colour = 'black',size=10),
        legend.title = element_blank(),legend.position = 'none')+labs(x="")

average_values= CXCL %>%
  group_by(ident,split) %>%
  summarise(avg_CXCR1 = median(CXCR1))

p2=ggplot(CXCL,aes(x=split,y=CXCR1,fill=split))+
  geom_line(average_values,mapping=aes(x=split,y=avg_CXCR1,group=ident))+
  geom_violin(scale = 'width')+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5)+
  stat_summary(mapping=aes(group=split),                      ##设置分组列
               fun="median",                                   ##箱线图添加均值使用的公式mean
               geom="point",size=3,fill="white",    ##均值图形的设置
               position=position_dodge(0.8))+
  stat_compare_means(
    aes(label = ..p.signif..),  # 显示 * 号（而不是 p 值）
    method = "wilcox.test",     # 使用 Wilcoxon 秩和检验
    comparisons = list(c("LA", "HA")),  # 替换成你的分组名称
    size = 4,
    vjust = 0.5
  ) +
  theme_classic()+facet_wrap(ident ~ .,scales = "free",nrow = 1)+scale_fill_manual(values = c('#7aa6dc',"#cd534c"))+
  theme(axis.text.x = element_text(colour = 'black',size=10),axis.text.y = element_text(colour = 'black',size=10),
        legend.title = element_blank(),legend.position = 'none')+labs(x="")

average_values= CXCL %>%
  group_by(ident,split) %>%
  summarise(avg_CXCR2 = median(CXCR2))

p3=ggplot(CXCL,aes(x=split,y=CXCR2,fill=split))+
  geom_line(average_values,mapping=aes(x=split,y=avg_CXCR2,group=ident))+
  geom_violin(scale = 'width')+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5)+
  stat_summary(mapping=aes(group=split),                      ##设置分组列
               fun="median",                                   ##箱线图添加均值使用的公式mean
               geom="point",size=3,fill="white",    ##均值图形的设置
               position=position_dodge(0.8))+
  stat_compare_means(
    aes(label = ..p.signif..),  # 显示 * 号（而不是 p 值）
    method = "wilcox.test",     # 使用 Wilcoxon 秩和检验
    comparisons = list(c("LA", "HA")),  # 替换成你的分组名称
    size = 4,
    vjust = 0.5
  ) +
  theme_classic()+facet_wrap(ident ~ .,scales = "free",nrow = 1)+scale_fill_manual(values = c('#7aa6dc',"#cd534c"))+
  theme(axis.text.x = element_text(colour = 'black',size=10),axis.text.y = element_text(colour = 'black',size=10),
        legend.title = element_blank(),legend.position = 'none')+labs(x="")

p1/p2/p3

p=plotGeneExpression(cellChat,
                     signaling = "BTLA",
                     split.by = "datasets", 
                     colors.ggplot = T)


BTLA=p[[1]][[1]][["data"]]
BTLA$TNFRSF14=p[[2]][[1]][["data"]]$TNFRSF14

average_values= BTLA %>%
  group_by(ident,split) %>%
  summarise(avg_BTLA = median(BTLA))

p1=ggplot(BTLA,aes(x=split,y=BTLA,fill=split))+
  geom_line(average_values,mapping=aes(x=split,y=avg_BTLA,group=ident))+
  geom_violin(scale = 'width')+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5)+
  stat_summary(mapping=aes(group=split),                      ##设置分组列
               fun="median",                                   ##箱线图添加均值使用的公式mean
               geom="point",size=3,fill="white",    ##均值图形的设置
               position=position_dodge(0.8))+
  stat_compare_means(
    aes(label = ..p.signif..),  # 显示 * 号（而不是 p 值）
    method = "wilcox.test",     # 使用 Wilcoxon 秩和检验
    comparisons = list(c("LA", "HA")),  # 替换成你的分组名称
    size = 4,
    vjust = 0.5
  ) +
  theme_classic()+facet_wrap(ident ~ .,scales = "free",nrow = 1)+scale_fill_manual(values = c('#7aa6dc',"#cd534c"))+
  theme(axis.text.x = element_text(colour = 'black',size=10),axis.text.y = element_text(colour = 'black',size=10),
        legend.title = element_blank(),legend.position = 'none')+labs(x="")

average_values= BTLA %>%
  group_by(ident,split) %>%
  summarise(avg_TNFRSF14 = median(TNFRSF14))

p2=ggplot(BTLA,aes(x=split,y=TNFRSF14,fill=split))+
  geom_line(average_values,mapping=aes(x=split,y=avg_TNFRSF14,group=ident))+
  geom_violin(scale = 'width')+
  stat_boxplot(geom="errorbar",width=0.1,size=0.5)+
  stat_summary(mapping=aes(group=split),                      ##设置分组列
               fun="median",                                   ##箱线图添加均值使用的公式mean
               geom="point",size=3,fill="white",    ##均值图形的设置
               position=position_dodge(0.8))+
  stat_compare_means(
    aes(label = ..p.signif..),  # 显示 * 号（而不是 p 值）
    method = "wilcox.test",     # 使用 Wilcoxon 秩和检验
    comparisons = list(c("LA", "HA")),  # 替换成你的分组名称
    size = 4,
    vjust = 0.5
  ) +
  theme_classic()+facet_wrap(ident ~ .,scales = "free",nrow = 1)+scale_fill_manual(values = c('#7aa6dc',"#cd534c"))+
  theme(axis.text.x = element_text(colour = 'black',size=10),axis.text.y = element_text(colour = 'black',size=10),
        legend.title = element_blank(),legend.position = 'none')+labs(x="")

p1/p2

