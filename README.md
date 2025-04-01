# Single-Cell-Data-Analysis-Tutorial

Overview
This project provides a complete pipeline for Single-Cell RNA-Seq Data Analysis, including preprocessing, clustering, dimensionality reduction, cell type annotation, differential expression analysis, and visualization. This tutorial is structured to provide a reproducible workflow that can be deployed to GitHub.

Objectives
- Perform quality control and preprocessing of scRNA-Seq data.
- Apply clustering and dimensionality reduction techniques.
- Annotate cell types using known markers.
- Perform differential expression analysis.
- Visualize the results through various plots.


📚 Tools and Packages
R Packages: Seurat, dplyr, ggplot2, patchwork, GEOquery, BiocManager, readr, tibble


📊 Dataset
A publicly available single-cell RNA-Seq dataset is used, downloaded from the Gene Expression Omnibus (GEO) database.
Replace the geo_id in the script with your desired dataset.

Workflow
1. Data Acquisition
library(Seurat)
library(GEOquery)

geo_id <- "GSE123456"  # Replace with your actual GEO ID
geo_data <- getGEO(geo_id, GSEMatrix = TRUE)
expr_data <- exprs(geo_data[[1]])

sc_data <- CreateSeuratObject(counts = expr_data, project = "SingleCellProject")
saveRDS(sc_data, "raw_sc_data.rds")

2. Preprocessing
sc_data <- subset(sc_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
sc_data <- NormalizeData(sc_data)
sc_data <- FindVariableFeatures(sc_data, selection.method = "vst", nfeatures = 2000)
sc_data <- ScaleData(sc_data)

3. Clustering and Dimensionality Reduction
sc_data <- RunPCA(sc_data, features = VariableFeatures(object = sc_data))
sc_data <- RunUMAP(sc_data, dims = 1:10)
sc_data <- FindNeighbors(sc_data, dims = 1:10)
sc_data <- FindClusters(sc_data, resolution = 0.5)

4. Cell Type Annotation
FeaturePlot(sc_data, features = c("CD3D", "CD14", "LYZ", "MS4A1", "PPBP"))
sc_data$celltype <- Idents(sc_data)
DimPlot(sc_data, reduction = "umap", label = TRUE, group.by = "celltype")

6. Differential Expression Analysis
cluster_markers <- FindAllMarkers(sc_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, "cluster_markers.csv")

8. Visualization
FeaturePlot(sc_data, features = c("CD3D", "LYZ", "MS4A1"), min.cutoff = "q9")
DoHeatmap(sc_data, features = head(cluster_markers$gene, 20)) + NoLegend()
DimPlot(sc_data, reduction = "umap", label = TRUE, group.by = "celltype")
