############################################################
# 1) Load required packages
#    - DESeq2: differential expression analysis
#    - tidyverse: data handling and plotting
#    - EnhancedVolcano: volcano plot visualization
#    - org.Mm.eg.db / biomaRt: gene annotation
############################################################

if(!require(DESeq2)) install.packages("DESeq2")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(EnhancedVolcano)) install.packages("EnhancedVolcano")
if(!require(org.Mm.eg.db)) BiocManager::install("org.Mm.eg.db")
if(!require(biomaRt)) install.packages("biomaRt")

library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(biomaRt)


############################################################
# 2) Read the count matrix
#    - Rows: genes
#    - Columns: sample IDs (currently BAM file paths)
#    - remove long file paths → keep only sample IDs
############################################################

counts <- read.table(
  "counts_clean.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Clean column names to match sample IDs
colnames(counts) <- gsub(".*/", "", colnames(counts))          # remove folder paths
colnames(counts) <- gsub(".sorted.bam", "", colnames(counts))  # remove BAM suffix

dim(counts)
head(counts)[,1:5]


############################################################
# 3) Construct sample metadata (coldata)
#    - genotype: WT vs DoubleKO
#    - condition: Control vs Disease
#    - group = genotype + condition (used in DESeq2 design)
############################################################

sampleID  <- c("SRR7821918","SRR7821919","SRR7821920",
               "SRR7821921","SRR7821922","SRR7821923",
               "SRR7821924","SRR7821925","SRR7821927",
               "SRR7821937","SRR7821938","SRR7821939",
               "SRR7821940","SRR7821941","SRR7821942")

# genotype assignment based on experiment layout
genotype  <- c(rep("WT", 6),
               rep("DoubleKO", 3),
               rep("WT", 3),
               rep("DoubleKO", 3))

# condition assignment
condition <- c(rep("Disease", 9),
               rep("Control", 6))

coldata <- data.frame(
  row.names = sampleID,
  sampleID  = sampleID,
  genotype  = factor(genotype, levels = c("WT", "DoubleKO")),
  condition = factor(condition, levels = c("Control","Disease"))
)

# Combined factor for DESeq2
coldata$group <- factor(
  paste0(coldata$genotype, "_", coldata$condition),
  levels = c("WT_Control","WT_Disease","DoubleKO_Control","DoubleKO_Disease")
)

coldata

# Ensure that count matrix columns match metadata rows — CRITICAL step
all(colnames(counts) == rownames(coldata))


############################################################
# 4) Create the DESeq2 dataset
#    - design = ~ group allows all pairwise contrasts
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ group
)

# Filter out low-expressed genes (<10 total counts)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
dds


############################################################
# 5) Run the DESeq2 pipeline
#    - normalization
#    - dispersion estimation
#    - Wald tests
############################################################

dds <- DESeq(dds)


############################################################
# 6) Extract pairwise contrasts
#    These correspond to the four biologically relevant comparisons.
############################################################

res_WT_DvsC      <- results(dds, contrast = c("group","WT_Disease","WT_Control"))
res_KO_DvsC      <- results(dds, contrast = c("group","DoubleKO_Disease","DoubleKO_Control"))
res_KO_DvsWT_D   <- results(dds, contrast = c("group","DoubleKO_Disease","WT_Disease"))  # ⭐ main contrast
res_KO_CvsWT_C   <- results(dds, contrast = c("group","DoubleKO_Control","WT_Control"))


############################################################
# 7) Helper function to count differentially expressed genes
#    - padj < 0.05
#    - count up-regulated and down-regulated genes
############################################################

count_DE <- function(res) {
  total <- sum(res$padj < 0.05, na.rm=TRUE)
  up    <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm=TRUE)
  down  <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm=TRUE)
  return(c(total=total, up=up, down=down))
}

count_DE(res_WT_DvsC)
count_DE(res_KO_DvsC)
count_DE(res_KO_DvsWT_D)   # ⭐ biologically most important contrast
count_DE(res_KO_CvsWT_C)


############################################################
# 8) Annotate ENSEMBL IDs → gene symbols
#    - needed for volcano plot labels
############################################################

mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

annotations <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(res_KO_DvsWT_D),
  mart = mart
)

# Remove duplicates
annot <- annotations %>% distinct(ensembl_gene_id, .keep_all = TRUE)

# Add gene symbols to results table
res_KO_DvsWT_D$symbol <- annot$mgi_symbol[match(rownames(res_KO_DvsWT_D), annot$ensembl_gene_id)]


############################################################
# 9) Volcano plot for the MAIN contrast
#    - Highlights Irf7, Nos2, Lcn2
#    - Shows log2FC vs adjusted p-values
############################################################

EnhancedVolcano(
  res_KO_DvsWT_D,
  lab = res_KO_DvsWT_D$symbol,
  x = "log2FoldChange",
  y = "padj",
  selectLab = c("Irf7", "Nos2", "Lcn2"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  title = "DoubleKO_Disease vs WT_Disease",
  subtitle = "Differential expression analysis",
  
  # Label styling
  labSize = 5,
  labFace = "bold",
  boxedLabels = TRUE,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf,
  
  pointSize = 2.5,
  colAlpha = 0.7
)


############################################################
# 10) Extract normalized counts for gene-level visualization
############################################################

norm_counts <- counts(dds, normalized=TRUE)


############################################################
# 11) Boxplots for Irf7, Nos2, Lcn2
#     - log scale used to stabilize variance
#     - shows expression across the four groups
############################################################

### IRF7
irf7_id <- "ENSMUSG00000025498"

df_irf7 <- data.frame(
  count = as.numeric(norm_counts[irf7_id, ]),
  group = coldata$group
)

ggplot(df_irf7, aes(x=group, y=count)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5) +
  scale_y_log10() +
  theme_bw() +
  ggtitle("Irf7 expression across groups")


### NOS2
nos2_id <- "ENSMUSG00000000056"

df_nos2 <- data.frame(
  count = as.numeric(norm_counts[nos2_id, ]),
  group = coldata$group
)

ggplot(df_nos2, aes(x=group, y=count)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5) +
  scale_y_log10() +
  theme_bw() +
  ggtitle("Nos2 expression across groups")


### LCN2
lcn2_id <- "ENSMUSG00000030704"   
df_lcn2 <- data.frame(
  count = as.numeric(norm_counts[lcn2_id, ]),
  group = coldata$group
)

ggplot(df_lcn2, aes(x=group, y=count)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, alpha=0.5) +
  scale_y_log10() +
  theme_bw() +
  ggtitle("Lcn2 expression across groups")
