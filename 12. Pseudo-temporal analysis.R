library(Seurat)

Cell.integrated <- readRDS("F:/R/Osteoarthritis/3. RDS Files/Integrated/4. Seurat.integrated.Annotated.rds")

cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", 
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")

DimPlot(Cell.integrated, reduction = "umap",  pt.size = 0.5, label = T, raster = F,  
        cols = cb_palette)

################################################################################
################################################################################
library(monocle3)

# Change Seurat object to Monocle's cell_data_set object
cds <- as.cell_data_set(Cell.integrated)

# View gene expression data
head(cds@assays@data$counts)
head(colData(cds))
head(rowData(cds))

cds$cell_type <- Cell.integrated@active.ident

# If we haven't done the clustering by UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Learn trajectory
cds <- learn_graph(cds)

# Visualize the trajectory and color them by original cluster labels
plot_cells(cds, color_cells_by = "cell_type")

# Calculate the pseudo time
cds <- order_cells(cds)

# Draw the complete map and color them by cell types
plot_cells(cds,
           color_cells_by = "pseudotime",       
           label_groups_by_cluster = FALSE,   
           label_leaves = TRUE,              # Display（cell fates）
           label_branch_points = TRUE,       # Display branches
           graph_label_size = 4)             

ciliated_genes <- c("RUNX2")

plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


saveRDS(cds, file = "F:/R/Osteoarthritis/3. RDS Files/Integrated/Pseudotime.order.cds.rds")

################################################################################
################################################################################
library(monocle3)

# 从Seurat对象转换到Monocle的cell_data_set对象
cds <- as.cell_data_set(Cell.integrated)

# 查看Seurat对象中的condition属性
table(Cell.integrated$condition)  # 查看 HC 和 OA 的样本数量

HC <- Cell.integrated[, Cell.integrated$condition == "Healthy"]
OA <- Cell.integrated[, Cell.integrated$condition == "Osteoarthritis"]

DimPlot(HC, reduction = "umap",  pt.size = 0.5, label = T, raster = F,  
        cols = cb_palette)

DimPlot(OA, reduction = "umap",  pt.size = 0.5, label = T, raster = F,  
        cols = cb_palette)

# 根据condition属性将数据分开为HC和OA的子集
cds_HC <- cds[, cds$condition == "Healthy"]   # 选取condition为HC的样本
cds_OA <- cds[, cds$condition == "Osteoarthritis"]   # 选取condition为OA的样本

# 查看分割后的cell_data_set对象
head(colData(cds_HC))
head(colData(cds_OA))

# 对HC子集进行伪时间追踪
cds_HC <- cluster_cells(cds_HC, reduction_method = "UMAP")
cds_HC <- learn_graph(cds_HC)
cds_HC <- order_cells(cds_HC)

# 对OA子集进行伪时间追踪
cds_OA <- cluster_cells(cds_OA, reduction_method = "UMAP")
cds_OA <- learn_graph(cds_OA)
cds_OA <- order_cells(cds_OA)

# 可视化HC子集的伪时间轨迹
plot_cells(cds_HC,
           color_cells_by = "pseudotime",       # 根据伪时间着色
           label_groups_by_cluster = FALSE,     # 不显示聚类标签
           label_leaves = TRUE,                 # 显示轨迹终点
           label_branch_points = TRUE,          # 显示分支点
           graph_label_size = 4)                # 调整标签字体大小

# 可视化OA子集的伪时间轨迹
plot_cells(cds_OA,
           color_cells_by = "pseudotime",       # 根据伪时间着色
           label_groups_by_cluster = FALSE,     # 不显示聚类标签
           label_leaves = TRUE,                 # 显示轨迹终点
           label_branch_points = TRUE,          # 显示分支点
           graph_label_size = 4)                # 调整标签字体大小
