############################################################
# RNA-seq Exploratory Data Analysis (EDA)
# Description:
#   This script performs standard EDA for RNA-seq count data,
#   including normalization (VST), PCA, distance heatmaps,
#   correlation heatmaps, boxplots, and hierarchical clustering.
#   The script is designed for inclusion in a GitHub repository.
############################################################

# ============================================================
# 1. Load required packages
# ============================================================
library(DESeq2)        # Differential expression framework
library(ggplot2)        # Plotting
library(pheatmap)       # Heatmaps
library(genefilter)     # Filtering and distance calculations
library(vsn)            # Variance stabilizing transformation
library(RColorBrewer)   # Color palettes
library(ggrepel)        # Improved label placement in PCA plots

# ============================================================
# 2. Load raw count matrix
# ============================================================
# The counts file should contain:
# - Gene IDs as row names
# - Sample IDs as column names
counts <- read.table(
  "counts_clean.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# Print basic information about the count matrix
dim(counts)             # Number of genes and samples

# ============================================================
# 3. Define sample annotation (metadata)
# ============================================================


sample <- c(
  "SRR7821918","SRR7821919","SRR7821920",
  "SRR7821921","SRR7821922","SRR7821923",
  "SRR7821924","SRR7821925","SRR7821927",
  "SRR7821937","SRR7821938","SRR7821939",
  "SRR7821940","SRR7821941","SRR7821942"
)

# Genotype of each sample
genotype <- c(
  "WT","WT","WT",
  "WT","WT","WT",
  "DoubleKO","DoubleKO","DoubleKO",
  "WT","WT","WT",
  "DoubleKO","DoubleKO","DoubleKO"
)

# Experimental condition of each sample
condition <- c(
  "Disease","Disease","Disease",
  "Disease","Disease","Disease",
  "Disease","Disease","Disease",
  "Control","Control","Control",
  "Control","Control","Control"
)

# Construct metadata table
coldata <- data.frame(
  row.names = sample,
  genotype = factor(genotype),
  condition = factor(condition)
)

# Ensure sample names in counts match metadata
stopifnot(all(colnames(counts) == row.names(coldata)))

# ============================================================
# 4. Build DESeq2 dataset object
# ============================================================
# The design formula models the effect of genotype and condition.
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ genotype + condition
)

# ============================================================
# 5. Run normalization with DESeq2
# ============================================================
dds <- DESeq(dds)

# ============================================================
# 6. Variance Stabilizing Transformation (VST)
# ============================================================
# VST reduces heteroscedasticity in RNA-seq data,
# making samples more comparable for PCA and clustering.
vsd <- vst(dds, blind = TRUE)

# ============================================================
# 7. PCA Analysis
# ============================================================
# DESeq2 has built-in PCA, but here we use custom PCA for flexibility
pca <- prcomp(t(assay(vsd)))

# Prepare PCA data for plotting
pcaData <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  genotype  = coldata$genotype,
  condition = coldata$condition,
  sample    = rownames(coldata)
)

# Plot PCA
ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = condition, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  theme_bw() +
  labs(
    title = "PCA of VST-transformed Counts",
    x = paste0("PC1: ", round(summary(pca)$importance[2,1] * 100), "% variance"),
    y = paste0("PC2: ", round(summary(pca)$importance[2,2] * 100), "% variance")
  )

# ============================================================
# 8. Sample-to-sample distance heatmap
# ============================================================
# Computes Euclidean distance between all samples.
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDist,
  clustering_distance_cols = sampleDist,
  annotation_col = coldata,
  main = "Sample-to-sample Distance Heatmap (VST)"
)

# ============================================================
# 9. Pearson correlation heatmap
# ============================================================
# High correlation between replicates indicates good data quality.
cor_matrix <- cor(assay(vsd))

pheatmap(
  cor_matrix,
  annotation_col = coldata,
  main = "Pearson Correlation Heatmap (VST)"
)

# ============================================================
# 10. Boxplot of VST-normalized counts
# ============================================================
# Shows distribution of transformed gene expression values.
boxplot(
  assay(vsd),
  outline = FALSE,
  las = 2,
  main = "Boxplot of VST-transformed Counts"
)

# ============================================================
# 11. Hierarchical clustering
# ============================================================
# Clusters samples based on global expression similarity.
hc <- hclust(dist(t(assay(vsd))), method = "complete")
plot(hc, main = "Hierarchical Clustering of Samples")

############################################################
# End of Script
############################################################


