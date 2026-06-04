# Seurat Integration Tutorial

## 1. Install Micromamba

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)

# When prompted:
# Micromamba binary folder? [~/.local/bin]
# /DATA/cfam000*/bin

# Init shell (bash)? [Y/n]
# Y

# Configure conda-forge? [Y/n]
# Y

# Prefix location? [~/micromamba]
# /DATA/cfam000*/micromamba

# re-login
```
```bash
micromamba self-update
```

## 2. Create Seurat Environment

```bash
micromamba create -n seurat -c conda-forge -c bioconda \
  r-base=4.3 r-seurat r-seuratobject r-tidyverse r-patchwork

# Type Y and press Enter
```
```bash
micromamba activate seurat
```
```bash
R
```

---

## 3. Install Required R Packages

```r
install.packages("HGNChelper")  # First time only
install.packages("openxlsx")    # First time only

# Type 48 and press Enter

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(HGNChelper)
library(openxlsx)
```

---

## 4. Load Data

```r
set.seed(1234)
donor1.data <- Read10X(data.dir = "./donor1/outs/filtered_feature_bc_matrix/")
donor2.data <- Read10X(data.dir = "./donor2/outs/filtered_feature_bc_matrix/")

donor1 <- CreateSeuratObject( counts = donor1.data, project = "donor1", min.cells = 3, min.features = 200)
donor2 <- CreateSeuratObject( counts = donor2.data, project = "donor2", min.cells = 3, min.features = 200)
```
---

## 5. Quality Control (QC)

```r
donor1[["percent.mt"]] <- PercentageFeatureSet(donor1, pattern = "^MT-")
donor2[["percent.mt"]] <- PercentageFeatureSet(donor2, pattern = "^MT-")
```

### Violin Plots

```r
vln1 <- VlnPlot( donor1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vln2 <- VlnPlot( donor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

ggsave("VlnPlot_QC_donor1.pdf", vln1, width = 16, height = 10)
ggsave("VlnPlot_QC_donor2.pdf", vln2, width = 16, height = 10)
```
<img width="1142" height="708" alt="image" src="https://github.com/user-attachments/assets/f6bbdb50-caf6-4cce-a3e7-4f4f2b41bb4f" />


### Scatter Plots

```r
scatter1 <- FeatureScatter( donor1, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter( donor1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("FeatureScatter_QC_donor1.pdf", scatter1, width = 16, height = 10)

scatter2 <- FeatureScatter( donor2, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  FeatureScatter( donor2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("FeatureScatter_QC_donor2.pdf", scatter2, width = 16, height = 10)
```
<img width="1136" height="706" alt="image" src="https://github.com/user-attachments/assets/f6b14bc6-fbe3-4019-9126-23c48929df19" />

### Download QC PDFs

Open a new terminal:

```bash
scp cfam000*@133.41.125.54:/DATA/cfam000*/*.pdf .
```
---

## 6. Cell Filtering

```r
donor1 <- subset( donor1, subset = nCount_RNA < 50000 & nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)
donor2 <- subset( donor2, subset = nCount_RNA < 50000 & nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 10)
```
---

## 7. Dataset Integration
```r
donor1_sct <- SCTransform(donor1)
donor2_sct <- SCTransform(donor2)
```
```r
list_sct <- c(donor1_sct, donor2_sct)
features <- SelectIntegrationFeatures( object.list = list_sct)
list_sct <- PrepSCTIntegration( object.list = list_sct, anchor.features = features)
```
```r
anchors_sct <- FindIntegrationAnchors( object.list = list_sct, anchor.features = features, normalization.method = "SCT")
```
```r
combined <- IntegrateData( anchorset = anchors_sct, normalization.method = "SCT")
```
### Save R Data (optional)
```r
save(combined, file = "combined.RData")
```
---

## 8. PCA and UMAP
```r
combined <- RunPCA(combined, verbose = FALSE)
```
```r
p <- ElbowPlot(combined, ndims = 15)
ggsave("ElbowPlot.pdf", p)
```
Open a new terminal:
```bash
scp cfam0001@133.41.125.54:/DATA/cfam000*/ElbowPlot.pdf .
```
<img width="751" height="740" alt="image" src="https://github.com/user-attachments/assets/813cf94a-548d-4198-9ef4-3db45988deb0" />

```r
combined <- RunUMAP( combined, reduction = "pca", dims = 1:8)
combined <- FindNeighbors( combined, reduction = "pca", dims = 1:8)
combined <- FindClusters( combined, graph.name = "integrated_snn", resolution = 0.5)
```
---

## 9. Cell Type Annotation with ScType
```r
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined)
combined <- ScaleData(combined)
```
### Load ScType Functions

```r
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

### Prepare Marker Database

```r
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx", "Immune system")
```

### Calculate Scores

```r
es.max <- sctype_score(scRNAseqData = LayerData(combined, assay = "RNA", layer = "scale.data"), scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
```

### Assign Cell Types

```r
cL_results <- do.call(rbind,
  lapply(unique(combined@meta.data$seurat_clusters), function(cl) {
    scores <- sort(rowSums(es.max[, rownames(combined@meta.data[combined$seurat_clusters == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(scores), scores = scores), 10)}))

sctype_scores <- cL_results %>%
  group_by(cluster) %>%
  slice_max(scores, n = 1)
```

### Add Annotation

```r
combined@meta.data$sctype_classification <- ""
for (j in unique(sctype_scores$cluster)) {
  cl_type <- sctype_scores[
    sctype_scores$cluster == j,]

  combined@meta.data$sctype_classification[
    combined@meta.data$seurat_clusters == j
  ] <- as.character(cl_type$type[1])
}
```

### Visualize

```r
p4 <- DimPlot(combined, reduction = "umap", group.by = c("sctype_classification", "orig.ident"), label = TRUE, repel = TRUE)
p5 <- DimPlot( combined, reduction = "umap", group.by = "sctype_classification", split.by = "orig.ident", label = TRUE, repel = TRUE)
ggsave("UMAP_sctype.pdf", p4, width = 16, height = 10)
ggsave("UMAP_sctype_split.pdf", p5, width = 16, height = 10)
```
<img width="1140" height="712" alt="image" src="https://github.com/user-attachments/assets/d61d51be-3849-4ef6-88e3-ac359e7fbf84" />
<img width="1138" height="713" alt="image" src="https://github.com/user-attachments/assets/45544663-4a88-485c-a1e2-414eb912a40d" />

---

## 10. Identify conserved cell type markers

```r
combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
Idents(combined) <- "sctype_classification"
```

### Marker Detection

```r
all_markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, "all_markers.csv", row.names = FALSE)
```

### Top 5 Markers

```r
top5 <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
```

---

## 11. Visualization

```r
p_dot <- DotPlot(combined, features = unique(top5$gene)) + RotatedAxis()
ggsave("DotPlot_top5.pdf", p_dot, width = 16, height = 10)
```
<img width="1148" height="708" alt="image" src="https://github.com/user-attachments/assets/1aea1962-b15a-4066-b486-c3352314b8cf" />

```r
pbmc_markers <- c("CD3D", "CD4", "CD8A", "CD19", "MS4A1", "CD14", "LYZ", "GNLY", "NKG7", "FCER1A", "CST3")
p_vln <- VlnPlot(combined, features = pbmc_markers, group.by = "sctype_classification", pt.size = 0, ncol = 4)
ggsave("VlnPlot_markers.pdf", p_vln, width = 16, height = 10)
```
<img width="1141" height="710" alt="image" src="https://github.com/user-attachments/assets/b3fee6ad-d74d-46e9-9fda-73b708781dc3" />

```r
p_feat <- FeaturePlot(combined, features = pbmc_markers, ncol = 4)
ggsave("FeaturePlot_markers.pdf", p_feat, width = 16, height = 10)
```
<img width="1138" height="708" alt="image" src="https://github.com/user-attachments/assets/11b1eaae-e465-4e19-a3f6-e1fa138a80dd" />

---

## 12. DEG Analysis

```r
Idents(combined) <- "sctype_classification"
cell_types <- unique(combined@meta.data$sctype_classification)

for (ct in cell_types) {

  obj_ct <- subset(combined, idents = ct)
  Idents(obj_ct) <- "orig.ident"
  tryCatch({deg <- FindMarkers(obj_ct, ident.1 = "donor1", ident.2 = "donor2", min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
    write.csv( deg, paste0("DEG_", gsub(" ", "_", ct), "_donor1_vs_donor2.csv"))})}
```
---

## 13. Visualization

```r
combined$cluster_sample <- paste0("cluster", combined$seurat_clusters, "-", combined$orig.ident)
agg <- AggregateExpression(combined, group.by = "cluster_sample", assays = "RNA", return.seurat = TRUE)
DefaultAssay(agg) <- "RNA"

agg_mat <- GetAssayData(agg, layer = "data")

expressed_genes <- rownames(agg_mat)[
  agg_mat[, "cluster0-donor1"] > 2 &   
  agg_mat[, "cluster0-donor2"] > 2     
]

top5 <- all_markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(pct.1 > 0.25) %>%
  filter(avg_log2FC > 0.5) %>%
  filter(gene %in% expressed_genes) %>%   
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

genes.to.label <- unique(top5$gene)
genes.to.label2 <- intersect(genes.to.label, rownames(agg))
```

```r
p1 <- CellScatter(agg, cell1 = "cluster0-donor1", cell2 = "cluster0-donor2", highlight = genes.to.label2)
p2 <- LabelPoints(plot = p1, points = genes.to.label2, repel = TRUE)

ggsave("cluster0_scatter_top5.pdf", p2, width = 8, height = 6)
```
<img width="950" height="710" alt="image" src="https://github.com/user-attachments/assets/ab3ada11-c8ac-410f-af34-4e27fdcd3408" />

```r
genes.manual <- c("S100A9", "ETV6", "CD74", "HLA-B")

# FeaturePlot
p_feat <- FeaturePlot(combined, features = genes.manual, split.by = "orig.ident", max.cutoff = 3, reduction = "umap")
ggsave("cluster_FeaturePlot_manual.pdf", p_feat, width = 10, height = 12)

p_vln <- VlnPlot(combined, features = genes.manual, split.by = "orig.ident", pt.size  = 0, ncol     = 2)
ggsave("cluster_VlnPlot_manual.pdf", p_vln, width = 10, height = 12)
```
<img width="472" height="762" alt="image" src="https://github.com/user-attachments/assets/e5bbf9c5-6826-4977-894e-bb60670ba4b6" />
<img width="531" height="639" alt="image" src="https://github.com/user-attachments/assets/e8df5085-d043-46cb-a8d3-d9b23d20ed64" />

## 14. Save Final Results
```r
save(combined, all_markers, top5, file = "final.RData")
```
Open a new terminal:
```bash
scp cfam0001@133.41.125.54:/DATA/cfam000*/*.pdf .
scp cfam0001@133.41.125.54:/DATA/cfam000*/*.csv .
```
