## 1. Data Acquisition

```r
# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(GEOquery)
library(BiocManager)
library(readr)
library(tibble)

# Downloading dataset from GEO
geo_id <- "GSE123456"  # Replace with actual GEO ID

# Fetch the dataset
geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
expr_data <- exprs(geo_data[[1]])

# Convert expression data to a Seurat object
sc_data <- CreateSeuratObject(counts = expr_data, project = "SingleCellProject")

# Save raw data for backup
saveRDS(sc_data, "raw_sc_data.rds")
```

## 2. Preprocessing

```r
# Quality Control
sc_data <- subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
sc_data <- NormalizeData(sc_data)

# Feature Selection
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)

# Scaling the data
sc_data <- ScaleData(sc_data)
```

## 3. Clustering and Dimensionality Reduction

```r
# PCA
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))

# Visualizing PCA
VizDimLoadings(sc_data, dims = 1:2, reduction = "pca")
DimPlot(sc_data, reduction = "pca")

# UMAP
sc_data <- RunUMAP(sc_data, dims = 1:10)
DimPlot(sc_data, reduction = "umap", label = TRUE, pt.size = 0.5)

# t-SNE (Optional)
sc_data <- RunTSNE(sc_data, dims = 1:10)
DimPlot(sc_data, reduction = "tsne", label = TRUE, pt.size = 0.5)

# Clustering
sc_data <- FindNeighbors(sc_data, dims = 1:10)
sc_data <- FindClusters(sc_data, resolution = 0.5)
DimPlot(sc_data, reduction = "umap", label = TRUE)
```

## 4. Cell Type Annotation

```r
# Annotating Cell Types
# Use known markers for annotation (replace with relevant markers)
FeaturePlot(sc_data, features = c("CD3D", "CD14", "LYZ", "MS4A1", "PPBP"))

# Assigning cell types
sc_data$celltype <- Idents(sc_data)

# Visualization of annotated cell types
DimPlot(sc_data, reduction = "umap", label = TRUE, group.by = "celltype")
```

## 5. Differential Expression Analysis

```r
# Finding markers for each cluster
cluster_markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Displaying top markers
head(cluster_markers)

# Saving results to file
write.csv(cluster_markers, "cluster_markers.csv")
```

## 6. Visualization

```r
# Visualizing differentially expressed genes
FeaturePlot(sc_data, features = c("CD3D", "LYZ", "MS4A1"), min.cutoff = "q9")

# Heatmap of top markers
DoHeatmap(sc_data, features = head(cluster_markers$gene, 20)) + NoLegend()

# Visualizing annotated clusters
DimPlot(sc_data, reduction = "umap", label = TRUE, group.by = "celltype")
