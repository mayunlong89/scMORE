
# This file was used to store figures and examples

## For example of single-cell multiomics dataset (10X_pbmc_400cells.rds)

This dataset contains 400 cells, with 200 monocytes and 200 CD8+T cells, used for example learning.

```r
# Load 10x PBMC scMultiomic data
single_cell_10x <- readRDS("10X_PBMC.rds")

# Update cell type identities
Idents(single_cell_10x) <- single_cell_10x$cell_type
cell_type2 <- as.character(single_cell_10x$cell_type)
# Combine "CD14+ Monocytes" and "FCGR3A+ Monocytes" into "Monocytes"
cell_type2[which(cell_type2 %in% c("CD14+ Monocytes", "FCGR3A+ Monocytes"))] <- "Monocytes"
single_cell_10x$cell_type2 <- cell_type2

# Subset to include only Monocytes and CD8+ T cells
pbmc_10x_real_data <- subset(single_cell_10x, idents = c("Monocytes", "CD8+ T cells"))

# Downsample to 200 cells per cell type
cell.list <- WhichCells(pbmc_10x_real_data, idents = c("Monocytes", "CD8+ T cells"), downsample = 200)
single_cell <- pbmc_10x_real_data[, cell.list]

# Re-clustering and dimensionality reduction for UMAP
# Normalize data
single_cell <- NormalizeData(single_cell)
# Identify highly variable features
single_cell <- FindVariableFeatures(single_cell)
# Scale the data
single_cell <- ScaleData(single_cell)
# Perform PCA for dimensionality reduction
single_cell <- RunPCA(single_cell, reduction.name = "pca")
# Find neighbors based on PCA dimensions
single_cell <- FindNeighbors(single_cell, dims = 1:30, reduction = "pca")
# Find clusters at resolution 0.5
single_cell <- FindClusters(single_cell, resolution = 0.5)
# Run UMAP on PCA-reduced data
single_cell <- RunUMAP(single_cell, dims = 1:30, reduction = "pca", reduction.name = "umap.new")

# Visualization
# Set identities to the updated cell type labels
Idents(single_cell) <- single_cell$cell_type2
# Generate UMAP plot with customized colors
DimPlot(single_cell, reduction = "umap.new", cols = c("#67ADB7", "#E77A77"))

# Save the processed single-cell data to an RDS file
saveRDS(single_cell, file = "10X_pbmc_400cells.rds")


```

## Single_cell object
```r
> single_cell
An object of class Seurat 
134520 features across 400 samples within 2 assays 
Active assay: RNA (22815 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 1 other assay present: peaks
 7 dimensional reductions calculated: pca, umap.rna, lsi, umap.atac, umap.biMod, tsne.biMod, umap.new

```


