
# This file was used to store figures and examples

## For example of single-cell multiomics data

This dataset contains 400 cells, with 200 monocytes and 200 CD8+T cells, used for example learning.

```r
## codes for visualizing example single-cell data
#load data on 10x pbmc example data
#Load scMultiomic data
pbmc_10x <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_PBMC/10X_PBMC.rds")
single_cell_10x <- pbmc_10x
Idents(single_cell_10x) <- single_cell_10x$cell_type

cell_type2<- as.character(single_cell_10x$cell_type)
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes","FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell_10x$cell_type2 <- cell_type2


#subset cells of two cell types
pbmc_10x_real_data <- subset(single_cell_10x,idents=c("Monocytes","CD8+ T cells"))
#downsample cell list
cell.list <- WhichCells(pbmc_10x_real_data,idents = c("Monocytes","CD8+ T cells"),downsample = 1000)
pbmc_10x_real_data_downsampled_2000 <- pbmc_10x_real_data[,cell.list]
table(pbmc_10x_real_data_downsampled_2000$cell_type2)
pbmc_10x_real_data_downsampled_2000 #this real dataset used for assessing ctDRTF performance


#DimPlot(pbmc_10x_real_data_downsampled_2000)

#re-clustering the UMAP
pbmc_10x_real_data_downsampled_2000 <- NormalizeData(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- FindVariableFeatures(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- ScaleData(pbmc_10x_real_data_downsampled_2000)
pbmc_10x_real_data_downsampled_2000 <- RunPCA(pbmc_10x_real_data_downsampled_2000,reduction.name = "pca")
pbmc_10x_real_data_downsampled_2000 <- FindNeighbors(pbmc_10x_real_data_downsampled_2000,
                                                     dims = 1:30,reduction = "pca")
pbmc_10x_real_data_downsampled_2000 <- FindClusters(pbmc_10x_real_data_downsampled_2000,
                                                    resolution = 0.5)
#Run UMAP
pbmc_10x_real_data_downsampled_2000 <- RunUMAP(pbmc_10x_real_data_downsampled_2000, dims = 1:30,
                                               reduction = "pca",reduction.name = "umap.new")

#Visualization
Idents(pbmc_10x_real_data_downsampled_2000) <- pbmc_10x_real_data_downsampled_2000$cell_type2

DimPlot(pbmc_10x_real_data_downsampled_2000, reduction = "umap.new",cols = c("#67ADB7","#E77A77"))

single_cell<-pbmc_10x_real_data_downsampled_2000




```
