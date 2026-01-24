###############################################################################
# RNA-seq Differential Expression (DESeq2) + Interaction Term + Visualization
#
#
# Goal:
#   - Fit a DESeq2 model with an interaction term: genotype + condition + genotype:condition
#   - Extract the genotype × condition interaction results
#   - Summarize significant DE genes (padj < 0.05)
#   - Visualize: per-gene normalized counts, volcano plot(s), heatmap(s)
#
# Input files expected (tab-delimited):
#   1) counts_matrix.tsv : first column = gene IDs (e.g., ENSEMBL), remaining columns = samples
#   2) samples.tsv       : must include columns: sample_id, genotype, condition
#
# Notes:
#   - This script uses mouse annotation org.Mm.eg.db (Mus musculus).
#   - Counts are rounded to integers (DESeq2 requires integer counts).
###############################################################################

# ----------------------------- 0) Packages ----------------------------------

# List of packages used in this workflow.
# DESeq2: differential expression
# tidyverse/tibble: data manipulation and tidying
# AnnotationDbi/org.Mm.eg.db: mapping between gene identifiers (SYMBOL/ENSEMBL)
# clusterProfiler/enrichplot: enrichment (loaded here; ORA steps not shown in plots below)
pkgs <- c(
  "BiocManager",
  "DESeq2",
  "tidyverse",
  "tibble",
  "AnnotationDbi",
  "org.Mm.eg.db",
  "clusterProfiler",
  "enrichplot",
  "EnhancedVolcano",
  "pheatmap"
)

# Install missing packages. Bioconductor packages are installed via BiocManager.
# (BiocManager itself must be installed via CRAN if missing.)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  # BiocManager::install handles both CRAN and Bioconductor packages.
  BiocManager::install(to_install, update = FALSE, ask = FALSE)
}

# Load packages quietly (suppresses startup messages for cleaner logs).
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(tibble)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(EnhancedVolcano)
  library(pheatmap)
})

# ----------------------------- 1) Inputs ------------------------------------

# Paths to input files.
counts_path  <- "counts_matrix.tsv"
samples_path <- "samples.tsv"

# ----------------------------- 2) Read counts --------------------------------
# The counts matrix should be:
#   - tab-delimited
#   - first column: gene IDs (often Ensembl IDs)
#   - remaining columns: sample columns (raw read counts)
counts <- read.delim(counts_path, check.names = FALSE)

# Basic sanity check: at least 1 gene column + 1 sample column.
stopifnot(ncol(counts) >= 2)

# Use first column as rownames (gene IDs), then drop it from the data frame.
rownames(counts) <- counts[[1]]
counts[[1]] <- NULL

# If Ensembl IDs contain version suffixes like "ENSMUSG00000000001.1",
# remove everything from the first dot onward to standardize IDs.
rownames(counts) <- sub("\\..*$", "", rownames(counts))

# ----------------------------- 3) Read metadata ------------------------------
# Sample metadata must contain:
#   - sample_id: matches counts column names
#   - genotype: e.g., WT, DoubleKO
#   - condition: e.g., Control, Infection
coldata <- read.delim(samples_path)

stopifnot(all(c("sample_id", "genotype", "condition") %in% colnames(coldata)))

# Align metadata rows to count matrix columns.
# This ensures coldata is in the exact same order as counts columns.
coldata <- coldata[match(colnames(counts), coldata$sample_id), ]

# Confirm alignment succeeded (no NA rows and exact match).
stopifnot(all(coldata$sample_id == colnames(counts)))

# Convert to factors and set reference levels.
# Reference levels define the baseline for interpreting coefficients.
coldata$genotype  <- relevel(factor(coldata$genotype), ref = "WT")
coldata$condition <- relevel(factor(coldata$condition), ref = "Control")

# Create a combined group factor (useful for plotting).
# Example levels: "WT_Control", "WT_Infection", "DoubleKO_Control", "DoubleKO_Infection"
coldata$group <- factor(paste(coldata$genotype, coldata$condition, sep = "_"))

# ----------------------------- 4) Build DESeq2 object ------------------------
# DESeqDataSetFromMatrix constructs the main DESeq2 container:
#   - countData: integer matrix of raw counts
#   - colData: metadata with rownames matching sample columns
#   - design: model formula
#
# Design explained:
#   ~ genotype + condition + genotype:condition
#   - genotype: main effect of genotype (DoubleKO vs WT) at baseline condition
#   - condition: main effect of infection (Infection vs Control) at baseline genotype
#   - genotype:condition: interaction: "difference of differences"
#       * captures whether infection response differs by genotype
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),                 # enforce integers
  colData   = as.data.frame(coldata) %>% column_to_rownames("sample_id"),
  design    = ~ genotype + condition + genotype:condition
)

# ----------------------------- 5) Filter low-count genes ---------------------
# Common filter: remove genes with extremely low counts across all samples.
# Rationale:
#   - genes with nearly no counts have little statistical power
#   - filtering improves multiple testing burden and stability
dds <- dds[rowSums(counts(dds)) >= 10, ]

# ----------------------------- 6) Fit the model ------------------------------
# DESeq performs:
#   - size factor estimation (library size normalization)
#   - dispersion estimation
#   - negative binomial GLM fitting
#   - statistical tests for coefficients
dds <- DESeq(dds)

# ----------------------------- 7) Extract interaction results ----------------
# resultsNames gives the available coefficient names in the fitted model.
# For interaction models, DESeq2 encodes terms like:
#   "genotype_DoubleKO_vs_WT"
#   "condition_Infection_vs_Control"
#   "genotypeDoubleKO.conditionInfection" (or similar)
coef_names <- resultsNames(dds)
print(coef_names)

# Identify the interaction term name programmatically.
# This searches for a coefficient containing both genotype and condition.
int_name <- coef_names[grepl("genotype.*condition", coef_names)]
stopifnot(length(int_name) == 1)
int_name <- int_name[[1]]
cat("Using interaction coefficient:", int_name, "\n")

# Extract DESeqResults for the interaction term.
# Interpretation:
#   log2FoldChange > 0 => infection response is stronger in the non-reference genotype (e.g., DoubleKO)
#   log2FoldChange < 0 => infection response is stronger in the reference genotype (WT)
res_int <- results(dds, name = int_name)

# ----------------------------- 8) Significant genes summary ------------------
# Filter to significant genes using adjusted p-value (FDR).
res_int_sig <- res_int[!is.na(res_int$padj) & res_int$padj < 0.05, ]

n_DE <- nrow(res_int_sig)
up_interaction   <- sum(res_int_sig$log2FoldChange > 0)
down_interaction <- sum(res_int_sig$log2FoldChange < 0)

cat(
  "Interaction term:", int_name, "\n",
  "DE genes (padj<0.05):", n_DE, "\n",
  "Stronger infection response in non-reference genotype (log2FC>0):", up_interaction, "\n",
  "Stronger infection response in reference genotype (log2FC<0):", down_interaction, "\n"
)

# ----------------------------- 9) Gene ID mapping (SYMBOL <-> ENSEMBL) -------
# Many plots and enrichment tools require ENTREZ IDs or gene symbols.
# Here we demonstrate SYMBOL -> ENSEMBL (and later ENSEMBL -> SYMBOL for labeling).
key_genes_symbol <- c("Isg15", "Ifit1", "Gbp5")

# Map SYMBOL -> ENSEMBL using org.Mm.eg.db.
# AnnotationDbi::select returns a data.frame with requested columns.
map_symbol_to_ens <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = key_genes_symbol,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
) %>%
  distinct()

print(map_symbol_to_ens)

# ----------------------------- 10) Normalized counts + per-gene plot ---------
# counts(dds, normalized=TRUE) returns size-factor-normalized counts.
# These are useful for expression plots but not used for model fitting.
norm_counts <- counts(dds, normalized = TRUE)
rownames(norm_counts) <- sub("\\..*$", "", rownames(norm_counts))

# Pick the mapped ENSEMBL IDs that exist in the normalized count matrix.
ens_ids <- unique(na.omit(map_symbol_to_ens$ENSEMBL))
ens_ids <- ens_ids[ens_ids %in% rownames(norm_counts)]

# If mapping returns multiple ENSEMBL per SYMBOL, you might not get exactly 3.
# Here we enforce 3 only if mapping is 1-to-1 in your dataset.
stopifnot(length(ens_ids) == 3)

# Subset normalized counts for the chosen genes (rows=genes, cols=samples).
norm_3 <- norm_counts[ens_ids, , drop = FALSE]

# Create a tidy long-format data frame for ggplot.
# Step A: transpose to get rows as samples and columns as genes.
df_norm <- as.data.frame(t(norm_3)) %>%
  tibble::rownames_to_column("sample_id")

# Step B: metadata to data.frame; ensure sample_id exists as a column.
df_meta <- as.data.frame(coldata)
if (!"sample_id" %in% names(df_meta)) {
  df_meta <- df_meta %>% tibble::rownames_to_column("sample_id")
}
df_meta <- df_meta[, !duplicated(names(df_meta)), drop = FALSE]

# Step C: join and pivot longer.
# IMPORTANT:
#   - norm_3 rownames are ENSEMBL IDs, so the transposed columns are ENSEMBL, not SYMBOL.
#   - To facet by SYMBOL, we rename the columns from ENSEMBL -> SYMBOL using the mapping.
ens_to_symbol <- map_symbol_to_ens %>%
  filter(!is.na(ENSEMBL)) %>%
  distinct(ENSEMBL, SYMBOL)

# Rename columns in df_norm from ENSEMBL IDs to SYMBOLs (for those genes).
for (i in seq_len(nrow(ens_to_symbol))) {
  ens <- ens_to_symbol$ENSEMBL[i]
  sym <- ens_to_symbol$SYMBOL[i]
  if (ens %in% names(df_norm)) {
    names(df_norm)[names(df_norm) == ens] <- sym
  }
}

plot_df <- df_norm %>%
  left_join(df_meta, by = "sample_id") %>%
  pivot_longer(
    cols      = all_of(key_genes_symbol),
    names_to  = "SYMBOL",
    values_to = "norm_count"
  )

stopifnot(nrow(plot_df) > 0)
print(table(plot_df$group))
print(head(plot_df))

# Boxplot + jitter per group, faceted by gene symbol.
ggplot(plot_df, aes(x = group, y = norm_count)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Group", y = "Normalized counts")

# ----------------------------- 11) Volcano plot (basic) ----------------------
# EnhancedVolcano expects a data.frame-like object with columns for x and y.
# Using res_int directly labels points by rownames (ENSEMBL IDs).
EnhancedVolcano(
  res_int,
  lab      = rownames(res_int),
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  title    = "Genotype × Infection interaction"
)

# ----------------------------- 12) Volcano plot with SYMBOL labels -----------
# Convert DESeqResults to data.frame and add ENSEMBL as a column.
res_int_df <- as.data.frame(res_int) %>%
  tibble::rownames_to_column("ENSEMBL")

# Map ENSEMBL -> SYMBOL for labeling.
# multiVals="first" chooses one symbol if multiple map (simplifies plotting).
res_int_df$SYMBOL <- mapIds(
  org.Mm.eg.db,
  keys      = res_int_df$ENSEMBL,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Create a LABEL column: prefer SYMBOL, otherwise fall back to ENSEMBL.
res_int_df$LABEL <- ifelse(is.na(res_int_df$SYMBOL), res_int_df$ENSEMBL, res_int_df$SYMBOL)

EnhancedVolcano(
  res_int_df,
  lab      = res_int_df$LABEL,
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  title    = "Genotype × Infection interaction",
  subtitle = "Interaction effect: differential infection response by genotype",
  labSize  = 3
)

# Volcano plot highlighting only selected key genes.
EnhancedVolcano(
  res_int_df,
  lab = ifelse(res_int_df$SYMBOL %in% key_genes_symbol, res_int_df$SYMBOL, ""),
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  title    = "Genotype × Infection interaction (highlighted genes)"
)

# Same but with fixed x-axis limits.
EnhancedVolcano(
  res_int_df,
  lab = ifelse(res_int_df$SYMBOL %in% key_genes_symbol, res_int_df$SYMBOL, ""),
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  xlim     = c(-8, 8),
  title    = "Genotype × Infection interaction (highlighted genes)"
)

# ----------------------------- 13) Heatmap (significant genes) ---------------
# vst() computes variance-stabilized expression values:
#   - makes counts more homoscedastic (variance ~ constant)
#   - appropriate for clustering/heatmaps (not for DE testing)
vsd <- vst(dds, blind = TRUE)

# Significant genes by interaction term.
sig_genes <- rownames(res_int_sig)

# Expression matrix for those genes (rows=genes, cols=samples).
mat_sig <- assay(vsd)[sig_genes, , drop = FALSE]

# Column annotation (sample metadata) for heatmap.
# Use colData(dds) to guarantee alignment to assay columns.
ann <- as.data.frame(colData(dds)) %>%
  select(genotype, condition, group) %>%
  mutate(across(everything(), as.factor))

# Ensure sample order matches between annotation and matrix columns.
stopifnot(all(rownames(ann) == colnames(mat_sig)))

# Heatmap of all significant genes (can be large).
# scale="row" centers/scales each gene so patterns across samples are comparable.
pheatmap(
  mat_sig,
  annotation_col = ann,
  show_rownames  = FALSE,
  scale          = "row"
)

# ----------------------------- 14) Heatmap (top N by padj) -------------------
# To keep heatmaps readable, it is common to plot the top N most significant genes.
topN <- 100
top_genes <- rownames(res_int_sig)[order(res_int_sig$padj)][1:min(topN, nrow(res_int_sig))]

mat_top <- assay(vsd)[top_genes, , drop = FALSE]

pheatmap(
  mat_top,
  annotation_col = ann,
  show_rownames  = FALSE,
  scale          = "row"
)

###############################################################################
# End of script
###############################################################################
