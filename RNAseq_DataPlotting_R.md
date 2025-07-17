# RNA-seq Data Visualization using R

This workflow demonstrates how to analyze RNA-seq expression data using **DESeq2**, <br>and visualize the results through a **volcano plot**, **PCA**, and **heatmaps**.  
The input file should contain raw read counts along with gene annotations.

---

##  Input Format (example: `all_expression.txt`)

| Gene ID | Control-1(raw) | ... | Treat-3(raw) | Gene Symbol |
|---------|----------------|-----|--------------|--------------|
| 497097  | 8.5            | ... | 9.0          | Xkr4         |
| ...     | ...            |     | ...          | ...          |

---

##  Differential Expression Analysis using DESeq2

```r
library(DESeq2)
library(readr)
library(dplyr)

df <- read.table("all_expression.txt", header = TRUE, sep = "\t", check.names = FALSE)
count_data <- df[, 2:7] # columns 2 to 7 contain read counts
rownames(count_data) <- df$`Gene ID`
sample_info <- data.frame(condition = factor(c(rep("Control", 3), rep("Treat", 3))))
rownames(sample_info) <- colnames(count_data)

dds <- DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = sample_info,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ] # the threshold for low-expression filtering depends on the dataset
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = "DESeq2_results.csv")
```
##  Volcano Plot

```r
library(ggplot2)
library(ggrepel)

res_df <- transform(as.data.frame(res), GeneID = row.names(res))
gene_annotation <- setNames(df[, c("Gene ID", "Gene Symbol")], c("GeneID", "GeneSymbol"))
res_df <- merge(res_df, gene_annotation, by = "GeneID")
res_df <- na.omit(res_df)
res_df$threshold <- ifelse(res_df$padj < 1E-50 & abs(res_df$log2FoldChange) > 2, "Significant", "Not Significant")
label_df <- subset(res_df, threshold == "Significant")

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  geom_text_repel(data = label_df, aes(label = GeneSymbol), size = 3, max.overlaps = 20, show.legend = FALSE) +
  theme_bw() +
  labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10))

png("volcano_plot.png", width = 6, height = 5, units = "in", res = 300)
print(p)
dev.off()
```
<img width="1800" height="1500" alt="volcano_plot" src="https://github.com/user-attachments/assets/201017b1-6dee-480d-a88d-cdd70688d407" />

##  PCA Plot

```r
library(ggfortify)

vsd <- vst(dds, blind = FALSE)
pca_data <- t(assay(vsd))
sample_info <- as.data.frame(colData(vsd))
pca_result <- prcomp(pca_data, scale. = TRUE)

png("PCA_plot.png", width = 1000, height = 800, res = 150)
autoplot(pca_result,
         data = sample_info,
         colour = "condition",
         size = 4,
         label = FALSE,
         label.size = 5) +
  theme_bw() +
  theme(text = element_text(size = 16))
dev.off()
```
<img width="1000" height="800" alt="PCA_plot" src="https://github.com/user-attachments/assets/c9d784f7-3fa4-4c75-a9d6-35120bbde7c7" />

##  Heatmaps for Significant Genes

```r
library(pheatmap)

vsd_sig_matrix <- assay(vsd)[intersect(label_df$GeneID, rownames(assay(vsd))), , drop = FALSE]
rownames(vsd_sig_matrix) <- setNames(gene_annotation$GeneSymbol, gene_annotation$GeneID)[rownames(vsd_sig_matrix)]

# Unscaled heatmap
png("SignificantGenes_Heatmap.png", width = 1000, height = 800)
pheatmap(vsd_sig_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
dev.off()

# Row-scaled heatmap
png("SignificantGenes_Heatmap_scaled.png", width = 1000, height = 800)
pheatmap(vsd_sig_matrix,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")
dev.off()
```
<img width="1000" height="800" alt="SignificantGenes_Heatmap" src="https://github.com/user-attachments/assets/a9f2d67e-56c4-46d2-b43b-8b3c74c0051a" />
<img width="1000" height="800" alt="SignificantGenes_Heatmap_scaled" src="https://github.com/user-attachments/assets/b239b86e-6703-4f7c-b345-3e434b73cfb3" />

