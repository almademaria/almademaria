# Seurat Tutorial Following
# install.packages("dplyr")
# install.packages("Seurat")
# install.packages("patchwork")

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
setwd("/home/alma/Desktop/Practicas/tutorial-files-Seurat")
pbmc.data <- Read10X(data.dir = "/home/alma/Desktop/Practicas/tutorial-files-Seurat/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Exploring the columns inside the pbmc
pbmc$orig.ident

pbmc$nCount_RNA

pbmc$nFeature_RNA

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# The '.' values in the matrix represent 0s (no molecules detected). 
# Since most values in an scRNA-seq matrix are 0, Seurat uses a sparse-matrix 
# representation whenever possible. This results in significant memory and 
# speed savings for Drop-seq/inDrop/10x data.
# The dense size:
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
# The sparse size
sparse.size <- object.size(pbmc.data)
sparse.size
# Difference between sizes
dense.size/sparse.size

## QC control

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc[["percent.mt"]]

# Question: Where are QC metrics stored in Seurat?
# Show QC metrics for the first 5 cells
meta_values <- pbmc@meta.data
head(meta_values, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, 
        features = c("nFeature_RNA", 
                     "nCount_RNA", 
                     "percent.mt"), 
        ncol = 3)

# Feature Scatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores, etc
plot1 <- FeatureScatter(pbmc,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, 
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1 + plot2

# Normalizing the data: after removing the unwanted cells we need to normalize the values. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# In Seurat Normalized values are stored in pbmc[["RNA]]$data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc[["RNA"]]$data

pbmc <- NormalizeData(pbmc)

# Identification of the HV features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes (HVGs)
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 # Red --> + informative genes

# Now we need to scale the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc[["RNA"]]$scale.data


# Question: How can I remove unwanted sources of variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

## Perform linear dimensional reduction 

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") + NoLegend()

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(pbmc, dims = 1:3, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset
ElbowPlot(pbmc)

## Cluster the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# Run non-linear dimensional reduction (UMAP/tSNE) 
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")

# You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
saveRDS(pbmc, file = "pbmc_tutorial.rds")

# Finding differentially expressed features (cluster biomarkers) 
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), layer = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

DoHeatmap(pbmc, features = top10$gene) + NoLegend()

## Assigning cell type identity to clusters 

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", 
                label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "pbmc3k_final.rds")

# Get the Session information:
sessionInfo()


