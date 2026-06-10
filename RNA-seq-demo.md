# RNA-seq Analysis Pipeline
---

## 1. Preprocessing: TrimGalore

```bash
# Place fastq files in /DATA/centos/RNAseq/fastq
cd /DATA/centos/RNAseq/fastq
sudo cp /mnt/add2/*/*.fq .

# docker attach trimgalore-0.6.10
for r1 in *_1.fq; do
  r2="${r1/_1.fq/_2.fq}"
  trim_galore --paired "$r1" "$r2"
done
```

---

## 2. Alignment: STAR

```bash
cd ../STAR

for fq1 in ../fastq/*_1_val_1.fq; do
  sample=$(basename "$fq1" | sed 's/_1_val_1\.fq//')
  echo "mapping: ${sample}"
  ../STAR-2.7.0a/bin/Linux_x86_64/STAR \
    --runMode alignReads \
    --genomeDir ../ref/STAR_reference_mouse \
    --readFilesIn ../fastq/${sample}_1_val_1.fq \
                  ../fastq/${sample}_2_val_2.fq \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 128 \
    --quantMode TranscriptomeSAM \
    --outFileNamePrefix ${sample}.
done
```

---

## 3. Quantification: RSEM

```bash
for fq1 in ../fastq/*_1_val_1.fq; do
  sample=$(basename "$fq1" | sed 's/_1_val_1\.fq//')
  ../RSEM-1.3.3/rsem-calculate-expression \
    --num-threads 128 \
    --paired-end \
    --bam "${sample}.Aligned.toTranscriptome.out.bam" \
    ../ref/RSEM_reference_mouse/RSEM_reference_mouse \
    "${sample}"
done
```

---

## 4. DESeq2

> Run in the STAR directory. Prepare `sample2condition.txt` (tab-separated) and `target2gene.txt` beforehand.

```r
library(DESeq2)
library(tximport)
library(ggplot2)
library(ggrepel)

# Load metadata and gene ID mapping table
s2c <- read.table("sample2condition.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
s2c$group <- gsub(" ", "_", s2c$group)
t2g <- read.table("target2gene.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Import RSEM results (gene-level, expected_count)
files <- setNames(s2c$path, s2c$sample)
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
txi$length[txi$length == 0] <- 1  # replace length=0 to avoid division error

# Build DESeq2 object and run
sampleTable <- data.frame(condition = s2c$group, row.names = colnames(txi$counts))
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
dds <- DESeq(dds)
```

---

## 5. DEG Analysis + Volcano Plot

```r
# Define comparison pairs
conditions <- list(c("LK", "Ctrl"))

for (pair in conditions) {
  cond1 <- pair[1]
  cond2 <- pair[2]

  # Differential expression
  res <- results(dds, contrast = c("condition", cond1, cond2))
  res_df <- merge(
    data.frame(ensembl_gene_id = rownames(na.omit(res)), na.omit(res)),
    unique(t2g[, 2:3]),
    by = "ensembl_gene_id"
  )
  res_sort <- res_df[order(res_df$padj), ]
  write.table(res_sort, paste0(cond1, "-vs-", cond2, ".result.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # Volcano plot
  res_sort$log10padj <- -log10(res_sort$padj)
  res_sort$threshold <- with(res_sort,
    ifelse(padj < 0.05 & log2FoldChange >=  log2(2), "Upregulated",
    ifelse(padj < 0.05 & log2FoldChange <= -log2(2), "Downregulated", "Not significant"))
  )
  top_labels <- head(res_sort[res_sort$threshold != "Not significant", ][
    order(res_sort[res_sort$threshold != "Not significant", ]$padj), ], 10)

  p <- ggplot(res_sort, aes(x = log2FoldChange, y = log10padj, color = threshold)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(Upregulated = "red", Downregulated = "blue",
                                   "Not significant" = "gray")) +
    geom_text_repel(data = top_labels, aes(label = external_gene_name),
                    size = 3, max.overlaps = 15, show.legend = FALSE) +
    theme_classic() +
    labs(title = paste0("Volcano Plot: ", cond1, " vs ", cond2),
         x = "log2(Fold Change)", y = "-log10(padj)") +
    theme(legend.position = "bottom",
          plot.background  = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))

  ggsave(paste0(cond1, "-vs-", cond2, "_volcano_plot.png"), p, width = 7, height = 6, dpi = 300)
}
```

---

## 6. Filter DEGs (|log2FC| ≥ 1, padj < 0.05)

```r
library(dplyr)

for (input_file in list.files(pattern = "-vs-.*\\.result\\.txt$")) {
  read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
    filter(abs(log2FoldChange) >= log2(2), padj < 0.05) %>%
    write.table(paste0("filtered_", input_file), sep = "\t", quote = FALSE, row.names = FALSE)
}
```

---

## 7. Sample Clustering Heatmap & PCA

```r
library(RColorBrewer)
library(pheatmap)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Clustering heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(sampleTable$condition, rownames(sampleTable), sep = "_")
colnames(sampleDistMatrix) <- NULL

png("clustering.png")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

# PCA plot
png("PCA_plot.png")
plotPCA(vsd, intgroup = "condition") + theme_bw()
dev.off()
```

---

## 8. Top20 DEG Heatmap (per sample & per group)

```r
ens_to_ext <- setNames(t2g$external_gene_name, t2g$ensembl_gene_id)
s2c$SampleName <- s2c$sample

for (file in list.files(pattern = "^filtered_.*\\.result\\.txt$")) {
  comparison      <- gsub("^filtered_(.*)\\.result\\.txt$", "\\1", file)
  comparison_safe <- gsub("-", "_", comparison)

  filtered_results <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  present_genes    <- head(filtered_results$ensembl_gene_id, 20)
  present_genes    <- present_genes[present_genes %in% rownames(assay(vsd))]

  vsd_mat <- assay(vsd)[present_genes, rownames(sampleTable)]
  rownames(vsd_mat) <- ens_to_ext[present_genes]

  # Per-sample heatmaps
  for (scale in c("none", "row")) {
    suffix <- ifelse(scale == "none", "", "_scale-row")
    png(paste0(comparison_safe, "_sample_top20_heatmap", suffix, ".png"))
    pheatmap(vsd_mat, scale = scale,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean")
    dev.off()
  }

  # Per-group mean heatmap
  colnames(vsd_mat) <- colnames(dds)
  group_means <- sapply(unique(s2c$group), function(g) {
    samples <- s2c$SampleName[s2c$group == g]
    samples <- samples[samples %in% colnames(vsd_mat)]
    rowMeans(vsd_mat[, samples, drop = FALSE], na.rm = TRUE)
  })
  colnames(group_means) <- paste0("mean_", unique(s2c$group))

  for (scale in c("none", "row")) {
    suffix <- ifelse(scale == "none", "", "_scale-row")
    png(paste0(comparison_safe, "_group_top20_heatmap", suffix, ".png"))
    pheatmap(group_means, scale = scale,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean")
    dev.off()
  }
}
```

---

## 9. topGO (BP / CC / MF)

```r
library(topGO)
library(org.Mm.eg.db)
library(ggplot2)
library(stringr)
library(Rgraphviz)

all_genes <- rownames(dds)

for (file in list.files(pattern = "^filtered_.*\\.result\\.txt$")) {
  comparison      <- gsub("^filtered_(.*)\\.result\\.txt$", "\\1", file)
  comparison_safe <- gsub("-", "_", comparison)
  cat(sprintf("\n=== topGO: %s ===\n", comparison))

  deg_results <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(deg_results) == 0) { cat("No significant genes. Skipping.\n"); next }

  gene_list <- factor(as.integer(all_genes %in% deg_results$ensembl_gene_id))
  names(gene_list) <- all_genes

  for (ont in c("BP", "CC", "MF")) {
    cat(sprintf("\n--- %s ---\n", ont))

    GOdata <- new("topGOdata", ontology = ont, allGenes = gene_list,
                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "ENSEMBL")

    # Run 3 algorithms
    results_list <- list(
      classic  = runTest(GOdata, algorithm = "classic",  statistic = "fisher"),
      weight01 = runTest(GOdata, algorithm = "weight01", statistic = "fisher"),
      elim     = runTest(GOdata, algorithm = "elim",     statistic = "fisher")
    )

    allRes <- GenTable(GOdata,
                       classic  = results_list$classic,
                       weight01 = results_list$weight01,
                       elim     = results_list$elim,
                       orderBy = "weight01", ranksOf = "classic",
                       topNodes = 100, numChar = 1000)

    write.table(allRes,
                paste0(comparison_safe, "_topGO_", ont, "_all_algorithms.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # Convert p-values and compute plotting metrics
    for (col in c("classic", "weight01", "elim")) {
      allRes[[paste0(col, "_p")]] <- as.numeric(gsub("< ", "", allRes[[col]]))
      allRes[[paste0(col, "_p")]][is.na(allRes[[paste0(col, "_p")]])] <- 1e-50
    }
    top_res <- head(allRes, min(15, nrow(allRes)))
    top_res$GeneRatio <- top_res$Significant / top_res$Annotated
    top_res$Term_wrapped <- str_wrap(top_res$Term, width = 60)

    # Dot plots for each algorithm
    for (alg in c("classic", "weight01", "elim")) {
      p_col <- paste0(alg, "_p")
      p <- ggplot(top_res, aes(x = GeneRatio,
                                y = reorder(Term_wrapped, -log10(get(p_col))))) +
        geom_point(aes(size = Significant, color = -log10(get(p_col))), alpha = 0.8) +
        scale_color_gradient(low = "blue", high = "red", name = "-log10(p)") +
        scale_size_continuous(name = "Count", range = c(3, 8)) +
        labs(title = alg, x = "Gene Ratio", y = "") +
        theme_classic() +
        theme(axis.text.y = element_text(size = 12),
              plot.title  = element_text(hjust = 0.5, face = "bold", size = 12))

      ggsave(paste0(comparison_safe, "_topGO_", ont, "_", alg, ".png"),
             p, width = 8, height = 10, dpi = 300)

      # Network diagram (top 5 significant nodes)
      pdf(paste0(comparison_safe, "_topGO_", ont, "_", alg, "_network_top5.pdf"),
          width = 10, height = 8)
      showSigOfNodes(GOdata, score(results_list[[alg]]),
                     firstSigNodes = 5, useInfo = "all")
      dev.off()
    }

    cat(sprintf("Significant GO terms (p<0.05) — Classic: %d, Weight01: %d, Elim: %d\n",
                sum(allRes$classic_p  < 0.05),
                sum(allRes$weight01_p < 0.05),
                sum(allRes$elim_p     < 0.05)))
  }
}
```

---

## 10. clusterProfiler: GO + KEGG + GSEA

```r
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

for (file in list.files(pattern = "^filtered_.*\\.result\\.txt$")) {
  comparison      <- gsub("^filtered_(.*)\\.result\\.txt$", "\\1", file)
  comparison_safe <- gsub("-", "_", comparison)
  message("Processing: ", comparison)

  deg <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(deg) == 0) { message("No significant genes. Skipping."); next }

  # Ensembl → Entrez ID conversion
  gene_ids <- bitr(deg$ensembl_gene_id, fromType = "ENSEMBL",
                   toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  if (nrow(gene_ids) == 0) { message("ID conversion failed. Skipping."); next }
  message(nrow(gene_ids), " genes converted")

  # Helper: save result table + dotplot (+barplot if requested)
  save_enrich <- function(obj, prefix, barplot = FALSE) {
    if (is.null(obj) || nrow(obj@result) == 0) return(invisible(NULL))
    write.table(obj@result, paste0(prefix, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(paste0(prefix, "_dotplot.png"),
           dotplot(obj, showCategory = 20) + ggtitle(basename(prefix)),
           width = 10, height = 8, dpi = 300)
    if (barplot)
      ggsave(paste0(prefix, "_barplot.png"),
             barplot(obj, showCategory = 20) + ggtitle(basename(prefix)),
             width = 10, height = 8, dpi = 300)
    message("  Done: ", basename(prefix), " (", nrow(obj@result), " terms)")
  }

  # GO enrichment (BP / MF / CC)
  for (ont in c("BP", "MF", "CC")) {
    ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Mm.eg.db, ont = ont,
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                    readable = TRUE)
    save_enrich(ego, paste0(comparison_safe, "_GO_", ont), barplot = (ont == "BP"))
  }

  # KEGG enrichment
  ekegg <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "mmu",
                      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2)
  if (!is.null(ekegg) && nrow(ekegg@result) > 0)
    ekegg <- setReadable(ekegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  save_enrich(ekegg, paste0(comparison_safe, "_KEGG"), barplot = TRUE)

  # GSEA: rank genes by log2FC
  fc_vec <- setNames(deg$log2FoldChange, deg$ensembl_gene_id)
  fc_conv <- bitr(names(sort(fc_vec, decreasing = TRUE)),
                  fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  fc_list <- sort(setNames(fc_vec[fc_conv$ENSEMBL], fc_conv$ENTREZID), decreasing = TRUE)

  gsea_go <- gseGO(geneList = fc_list, OrgDb = org.Mm.eg.db, ont = "BP",
                   pvalueCutoff = 0.05, pAdjustMethod = "BH")
  save_enrich(gsea_go, paste0(comparison_safe, "_GSEA_GO"))

  gsea_kegg <- gseKEGG(geneList = fc_list, organism = "mmu",
                       pvalueCutoff = 0.05, pAdjustMethod = "BH")
  save_enrich(gsea_kegg, paste0(comparison_safe, "_GSEA_KEGG"))

  message("Finished: ", comparison)
}
```
