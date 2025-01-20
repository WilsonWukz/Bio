Cell.integrated <- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.1stAnnotated.rds")

# Step 1: Subset
cluster_mix <- WhichCells(Cell.integrated, idents = "Chondrocytes/Osteoblasts/Fibroblasts")
Cell.integrated.mix <- subset(Cell.integrated, cells = cluster_mix)

Cell.integrated.mix@active.ident<- Cell.integrated.mix$seurat_clusters
levels(Cell.integrated.mix@active.ident)
# Step 2: PCA and UMAP
Cell.integrated.mix <- RunPCA(Cell.integrated.mix, npcs = 30, verbose = TRUE)
Cell.integrated.mix <- RunUMAP(Cell.integrated.mix, reduction = "pca", dims = 1:30, verbose = TRUE)

# Step 3: Re-clustering
Cell.integrated.mix <- FindNeighbors(Cell.integrated.mix, reduction = "pca", dims = 1:30, verbose = TRUE)
Cell.integrated.mix <- FindClusters(Cell.integrated.mix, resolution = 0.9, verbose = 0)  # 提高分辨率以细分第27簇

# Step 4: Visualization
DimPlot(Cell.integrated.mix, reduction = "umap", label = TRUE, pt.size = 1)

cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")


#Fibroblasts: ..................................................................
#"COL1A1", "COL3A1", "THY1", "DCN", "LUM"
genes.Fibroblast <- c("COL1A1", "COL3A1", "THY1", "DCN", "LUM")
tidy.VlnPlot(Cell.integrated.mix, features = genes.Fibroblast, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.Fibroblast, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.Fibroblast, cols = c("blue","red"), dot.scale = 8) + RotatedAxis()

#Chondrocytes: .................................................................0 1 2 8
#"COL2A1","SOX9","ACAN","AGC1","Aggrecan","COMP","MIA"
get_marker(spc = "Human", cellname = "Chondrocyte")
genes.Chondrocytes <- c("COL2A1","SOX9","ACAN","COMP","MIA")
tidy.VlnPlot(Cell.integrated.mix, features = genes.Chondrocytes, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.Chondrocytes, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.Chondrocytes, cols = c("blue","red"), dot.scale = 8) + RotatedAxis()

#Proliferative Chondrocytes: ...................................................1 8
#"MKI67","PCNA","SOX9"
get_marker(spc = "Human", cellname = "proliferative chondrocyte cell")
genes.PC <- c("MKI67","PCNA","SOX9")
tidy.VlnPlot(Cell.integrated.mix, features = genes.PC, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.PC, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.PC, cols = c("blue","red"), dot.scale = 8) + RotatedAxis()


#Fibrochondrocytes :............................................................
#"COL1A1","COL3A1","COL6A1"
get_marker(spc = "Human", cellname = "Fibrochondrocyte")
genes.FC <- c("COL1A1","COL3A1","COL6A1")
tidy.VlnPlot(Cell.integrated.mix, features = genes.FC, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.FC, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.FC, cols = c("blue","red"), dot.scale = 8) + RotatedAxis()

#Proliferate fibrochondrocytes: ................................................2
#"COL1A1","FGF7","CTGF"
get_marker(spc = "Human", cellname = "proliferate fibrochondrocyte cell")
genes.ProFC <- c("COL1A1","FGF7","CTGF")
tidy.VlnPlot(Cell.integrated.mix, features = genes.ProFC, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.ProFC, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.ProFC, cols = c("blue","red"), dot.scale = 8) + RotatedAxis()

#Fibrochondrocyte progenitor cells:.............................................10
#"CD90","RUNX2","ACAN","GABA","B220"," CD105"," CD166","CD44","CD45"," CD45.1","CD51","CD73","CDH11","DLX5","EBF2"
get_marker(spc = "Human", tsuClass = "Bone", cellname = "Fibrochondrocyte progenitor cell")
genes.FCP <- c("MYLK","MCAM","COL3A1","COL1A1")
tidy.VlnPlot(Cell.integrated.mix, features = genes.FCP, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.FCP, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.FCP, cols = c("blue","red"), dot.scale = 8)

#Fibrocartilage chondrocytes:....................................................4
#"IGFBP5","LMCD1","MYL9","S100A6","SH3BGRL3"
get_marker(spc = "Human", tsuClass = "Bone", cellname = "Fibrocartilage chondrocyte")
genes.FCC <- c("IGFBP5","LMCD1","MYL9","S100A6","SH3BGRL3")
tidy.VlnPlot(Cell.integrated.mix, features = genes.FCC, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.FCC, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.FCC, cols = c("blue","red"), dot.scale = 8)

#Degradative Chondrocytes.......................................................12
#"COL10A1","RUNX2","MMP13","CDCP1"
get_marker(spc = "Human", cellname = "degradative chondrocyte cell")
genes.DegChondrocyte <- c("COL10A1","RUNX2","MMP13")
tidy.VlnPlot(Cell.integrated.mix, features = genes.DegChondrocyte, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.DegChondrocyte, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.DegChondrocyte, cols = c("blue","red"), dot.scale = 8)

#Hypertrophic Chondrocytes:
#"COL10A1","SPP1","MMP13"
get_marker(spc = "Human", cellname = "hypertrophic chondrocyte cell")
genes.HPChondrocyte <- c("COL10A1","SPP1","MMP13")
tidy.VlnPlot(Cell.integrated.mix, features = genes.HPChondrocyte, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.HPChondrocyte, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.HPChondrocyte, cols = c("blue","red"), dot.scale = 8)

#Homeostatic chondrocytes:.......................................................2 10
#"TMIP1","GDF15","IFITM3","TXNIP"
get_marker(spc = "Human", tsuClass = "Bone", cellname = "Homeostatic chondrocyte")
genes.HomeostaticChondrocyte <- c("GDF15","IFITM3","TXNIP")
tidy.VlnPlot(Cell.integrated.mix, features = genes.HomeostaticChondrocyte, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.HomeostaticChondrocyte, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.HomeostaticChondrocyte, cols = c("blue","red"), dot.scale = 8)


#Regulatory chondrocyte:........................................................
#"BMP2","FOSL1"
genes.RegC <- c("BMP2","FOSL1")
tidy.VlnPlot(Cell.integrated.mix, features = genes.RegC, pt.size = 0)
tidy.FeaturePlot(Cell.integrated.mix, features = genes.RegC, 
                 cols = c("grey", "red"))
DotPlot(Cell.integrated.mix, features = genes.RegC, cols = c("blue","red"), dot.scale = 8)

###################################################################
current.cluster.ids <- c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9)

new.cluster.ids <- c("Fibroblasts","Chondrocytes","Chondrocytes","Chondrocytes","Fibrocartilage chondrocytes","Homeostatic chondrocytes","Chondrocytes","Chondrocytes","Fibroblasts","Osteoblasts")  

Cell.integrated.mix@active.ident <- plyr::mapvalues(x = Cell.integrated.mix@active.ident, from = current.cluster.ids, to = new.cluster.ids)

head(Cell.integrated.mix@meta.data)

levels(Cell.integrated.mix@active.ident)

DimPlot(Cell.integrated.mix, reduction = "umap", raster = F, label = T)
# 
# Cell.integrated.mix$refined_clusters <- ifelse(
#   Cell.integrated.mix@active.ident== 4,
#   "27", #Pod
#   "28" #MES
# )

saveRDS(Cell.integrated.mix,"F:/R/Osteoarthritis/3. RDS Files/integrated/4.1 Seurat.integrated.reclustered.rds")

Cell.integrated.mix <- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4.1 Seurat.integrated.reclustered.rds")

Cell.integrated <- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.1stAnnotated.rds")



# 

# 初始化 refined_clusters 列为原始标签
Cell.integrated@meta.data$refined_clusters <- as.character(Idents(Cell.integrated))

# 将 Cell.integrated.mix 的细分标签映射回 Cell.integrated
# Cell.integrated.mix@active.ident
refined_labels <- as.character(Cell.integrated.mix@active.ident)
names(refined_labels) <- rownames(Cell.integrated.mix@meta.data)
# 
# # 用细分标签覆盖第27簇细胞在 refined_clusters 列的值
Cell.integrated@meta.data[cluster_mix, "refined_clusters"] <- refined_labels

table(Cell.integrated@meta.data$refined_clusters)
# 
DimPlot(Cell.integrated, reduction = "umap", group.by = "refined_clusters",
        pt.size = 0.5, raster = FALSE, cols = cb_palette,
        label = TRUE, label.size = 4, label.box = FALSE)


DimPlot(Cell.integrated, reduction = "umap", pt.size = 0.5, raster = F, cols = cb_palette,
        label = T, label.size = 4, label.box = F)

# 确保 cluster_mix 和 refined_labels 的名称匹配
if (!all(names(refined_labels) %in% cluster_mix)) {
  stop("Mismatch between refined_labels names and cluster_mix indices")
}

# 确保 refined_clusters 列存在，并赋值 refined_labels
Cell.integrated@meta.data$refined_clusters <- as.character(Cell.integrated@meta.data$refined_clusters)
Cell.integrated@meta.data[cluster_mix, "refined_clusters"] <- refined_labels

# 检查 refined_clusters 的值
table(Cell.integrated@meta.data$refined_clusters)

# 确保 refined_clusters 是因子并有正确的级别
valid_levels <- unique(Cell.integrated@meta.data$refined_clusters)
Cell.integrated@meta.data$refined_clusters <- factor(Cell.integrated@meta.data$refined_clusters, levels = valid_levels)
Idents(Cell.integrated) <- "refined_clusters"

# 可视化
if (!exists("cb_palette") || length(cb_palette) < length(valid_levels)) {
  stop("cb_palette is not defined or lacks sufficient colors for all clusters")
}

DimPlot(Cell.integrated, reduction = "umap", group.by = "refined_clusters", 
        pt.size = 0.5, raster = FALSE, cols = cb_palette,
        label = TRUE, label.size = 4, label.box = FALSE)

DimPlot(Cell.integrated, reduction = "umap",
        pt.size = 0.5, raster = FALSE, cols = cb_palette,
        label = TRUE, label.size = 4, label.box = FALSE)

levels(Cell.integrated@active.ident)

saveRDS(Cell.integrated, "F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.2ndAnnotated.rds")

Cell.integrated <- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.2ndAnnotated.rds")


plot3 <- DimPlot(Cell.integrated, reduction = "umap", split.by = "condition", pt.size = 0.5, label = T)

plot3 & theme(legend.position = "bottom") & guides(color = guide_legend(nrow = 2, byrow = TRUE, 
                                                                        override.aes = list(size = 3)))

#############################################################################
Cell.integrated.2ndAnnotated <- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.2ndAnnotated.rds")
##########################################################################
Cell.integrated<- readRDS("F:/R/Osteoarthritis/3. RDS Files/integrated/4. Seurat.integrated.2ndAnnotated.rds")

conserved.gene <- scCEGs(Cell.integrated, grouping.var = "condition")
View(conserved.gene)

############################################################################################ 
####################################   Draw stacked bar plot  ####################################   
####Calculate cell numbers for each cluster, based on clusters and groups/conditions
## extract meta data
md <- Cell.integrated.2ndAnnotated@meta.data %>% as.data.table
#将数据框或其他类型的数据对象转换为 data.table 对象。data.table 是R中一个高效的数据存储和操作框架。

# count the number of cells per unique combinations of "Sample" and "seurat_clusters"
Cell.integrated.2ndAnnotated_cluster.cell.numbers_1 <- md[, .N, by = c("condition", "refined_clusters")]
#理解不同条件下不同聚类中细胞的分布情况非常有用
view(Cell.integrated.2ndAnnotated_cluster.cell.numbers_1)

write.csv(Cell.integrated.2ndAnnotated_cluster.cell.numbers_1, 
          "F:/R/Osteoarthritis/4. CSV files/3.1 Cell.integrated.2ndAnnotated_cluster.cell.numbers_1.csv", quote = F)


# with additional casting after the counting
Cell.integrated.2ndAnnotated_cluster.cell.numbers_2 <- md[, .N, by = c("condition", "refined_clusters")] %>%
  dcast(., condition ~ refined_clusters, value.var = "N")
#使用 dcast() 函数将数据从长格式转换为宽格式

view(Cell.integrated.2ndAnnotated_cluster.cell.numbers_2)

write.csv(Cell.integrated.2ndAnnotated_cluster.cell.numbers_2, 
          "F:/R/Osteoarthritis/4. CSV files/3.1 Cell.integrated.2ndAnnotated_cluster.cell.numbers_2.csv", quote = F)
#########################################################################################################

#####stacked barplot
#Healthypare .csv file, containing 3 columns: Group (condition/Condition), sub-group(Cluster names/PodocyteSubset), Number (Cell numbers) 

data <- Cell.integrated.2ndAnnotated_cluster.cell.numbers_1
view(data)

#draw a stacked barplot, vertical
P_stackedbarplot <- ggplot(data, aes(fill = refined_clusters, y = N , x = condition)) + 
  geom_bar(position="fill", stat="identity") 
#创建一个堆叠条形图（stacked bar plot），用于展示不同条件下各个聚类的细胞计数或比例。
#data：传递给ggplot()函数的数据集，包含至少三列：seurat_clusters（聚类标识）、N（计数或比例数值）、condition（条件标识）。
#geom_bar()：这是ggplot2中的一个几何对象（geom），用于绘制条形图。
#position = "fill"：指定条形的堆叠方式，这里是按填充颜色堆叠，即每个 condition 下的 seurat_clusters 将堆叠在一起。
#stat = "identity"：指定条形图的统计方法，"identity" 使用数据中直接给定的y值，而不是计算它们（如求和或平均）。
P_stackedbarplot

P_stackedbarplot + theme(panel.background = element_rect(fill = NA)) + 
  theme(panel.grid.major = element_line(colour = NA))



#draw a stacked barplot, lateral
P_stackedbarplot_2 <- ggplot(data, aes(fill = refined_clusters, y = condition, x= N)) + 
  geom_bar(position="fill", stat="identity")
#y = N：指定y轴的数值基于 N 列，这代表每个条形的高度或计数。
#x = condition：指定x轴的类别基于 condition 列，不同的条件将形成不同的条形。
P_stackedbarplot_2 + theme(panel.grid.major = element_line(colour = "lemonchiffon1"))

P_stackedbarplot_2 + theme(panel.background = element_rect(fill = NA)) + 
  theme(panel.grid.major = element_line(colour = NA))+labs(title = "Frequencies", 
                                                           x = "Frequency", y = NULL, fill = "Seurat_clusters", 
                                                           size = 15)

############################################################################################
############################################################################################ 