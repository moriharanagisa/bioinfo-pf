# PBMC scRNA-seq Analysis Workflow (Seurat)

## Overview

This workflow covers:

* Cell Ranger processing
* Quality control
* Data integration
* Clustering and UMAP visualization
* ScType cell type annotation
* Marker gene identification
* Differential expression analysis

---

# 1. Data Import

Load Cell Ranger filtered matrices.

```r
donor1.data <- Read10X(...)
donor2.data <- Read10X(...)

donor1 <- CreateSeuratObject(...)
donor2 <- CreateSeuratObject(...)
```

---

# 2. Quality Control

## Calculate Mitochondrial Percentage

```r
donor1[["percent.mt"]] <- PercentageFeatureSet(donor1, pattern = "^MT-")
donor2[["percent.mt"]] <- PercentageFeatureSet(donor2, pattern = "^MT-")
```

## Generate QC Plots

* Violin plots
* FeatureScatter plots

Inspect:

* nCount_RNA
* nFeature_RNA
* percent.mt

Determine filtering thresholds based on distributions.

---

# 3. Cell Filtering

Recommended PBMC thresholds:

```r
nCount_RNA < 50000
nFeature_RNA > 500
nFeature_RNA < 7000
percent.mt < 10
```

Apply identical filters to all samples.

---

# 4. Data Integration

## Normalize and Identify Variable Features

```r
obj.list <- list(donor1, donor2)

obj.list <- lapply(obj.list, function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(
    obj,
    nfeatures = 2000,
    verbose = FALSE
  )
  obj
})
```

## Select Integration Features

```r
features <- SelectIntegrationFeatures(
  object.list = obj.list
)
```

## Find Integration Anchors

```r
anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  verbose = FALSE
)
```

## Integrate Datasets

```r
combined <- IntegrateData(
  anchorset = anchors,
  verbose = FALSE
)
```

---

# 5. PCA, UMAP and Clustering

Use the integrated assay.

```r
DefaultAssay(combined) <- "integrated"
```

## Scale and Run PCA

```r
combined <- ScaleData(combined)
combined <- RunPCA(combined)
```

## Determine Number of PCs

```r
ElbowPlot(combined)
```

Example:

```r
dims = 1:8
```

## Run UMAP

```r
combined <- RunUMAP(
  combined,
  reduction = "pca",
  dims = 1:8
)
```

## Clustering

```r
combined <- FindNeighbors(
  combined,
  reduction = "pca",
  dims = 1:8
)

combined <- FindClusters(
  combined,
  resolution = 0.5
)
```

---

# 6. ScType Cell Type Annotation

## Load ScType Scripts

```r
source("gene_sets_prepare.R")
source("sctype_score_.R")
```

## Load Immune Reference

```r
gs_list <- gene_sets_prepare(
  ".../ScTypeDB_full.xlsx",
  "Immune system"
)
```

## Prepare RNA Assay

```r
DefaultAssay(combined) <- "RNA"

combined[["RNA"]] <- JoinLayers(
  combined[["RNA"]]
)

combined <- ScaleData(
  combined,
  verbose = FALSE
)
```

## Calculate ScType Scores

```r
es.max <- sctype_score(
  scRNAseqData = LayerData(
    combined,
    assay = "RNA",
    layer = "scale.data"
  ),
  scaled = TRUE,
  gs = gs_list$gs_positive,
  gs2 = gs_list$gs_negative
)
```

Assign the highest-scoring cell type to each cluster.

---

# 7. UMAP Visualization

## Cluster UMAP

```r
DimPlot(
  combined,
  group.by = "seurat_clusters",
  label = TRUE
)
```

## Sample UMAP

```r
DimPlot(
  combined,
  group.by = "orig.ident"
)
```

## ScType UMAP

```r
DimPlot(
  combined,
  group.by = "sctype_classification",
  label = TRUE
)
```

## Split by Sample

```r
DimPlot(
  combined,
  group.by = "sctype_classification",
  split.by = "orig.ident"
)
```

---

# 8. Marker Gene Identification

Switch back to the RNA assay.

```r
DefaultAssay(combined) <- "RNA"

combined[["RNA"]] <- JoinLayers(
  combined[["RNA"]]
)

Idents(combined) <- "sctype_classification"
```

## Find Marker Genes

```r
all_markers <- FindAllMarkers(
  combined,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
```

Save results:

```r
write.csv(
  all_markers,
  "all_markers.csv"
)
```

---

# 9. Marker Visualization

## Dot Plot

```r
DotPlot(
  combined,
  features = unique(top5$gene)
)
```

## Violin Plot

Representative PBMC markers:

* CD3D
* CD4
* CD8A
* CD19
* MS4A1
* CD14
* LYZ
* GNLY
* NKG7
* FCER1A
* CST3

```r
VlnPlot(...)
```

## Feature Plot

```r
FeaturePlot(...)
```

---

# 10. Differential Expression Analysis

## donor1 vs donor2 Within Each Cell Type

```r
Idents(combined) <- "sctype_classification"
```

Loop through all annotated cell types:

```r
for (ct in cell_types) {

  obj_ct <- subset(
    combined,
    idents = ct
  )

  Idents(obj_ct) <- "orig.ident"

  deg <- FindMarkers(
    obj_ct,
    ident.1 = "donor1",
    ident.2 = "donor2",
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

  write.csv(...)
}
```

Output:

```text
DEG_<celltype>_donor1_vs_donor2.csv
```

---

# 11. Save Workspace

```r
save(
  combined,
  all_markers,
  top5,
  file = "integrate.RData"
)
```

## Final Outputs

* integrate.RData
* all_markers.csv
* UMAP figures
* DotPlot
* Violin plots
* Feature plots
* Cell type-specific DEG tables
