###############################################################################
# GO Biological Process Enrichment Analysis
# (Interaction DE genes from DESeq2: genotype × condition)
#
# Goal:
#   - Perform over-representation analysis (ORA) using clusterProfiler
#   - Use ENSEMBL IDs from the interaction term
#   - Define a proper gene universe
#   - Visualize enriched GO Biological Processes
###############################################################################

# ----------------------------- 0) Packages ----------------------------------
# clusterProfiler : GO enrichment analysis
# org.Mm.eg.db    : Mouse gene annotation database
# dplyr / tibble  : Data manipulation
# stringr         : Text wrapping for plot labels

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(tibble)
library(stringr)

# ----------------------------- 1) Significant interaction genes --------------
# res_int is a DESeqResults object from the interaction term
# Filter for statistically significant genes using adjusted p-value (FDR)

res_int_sig <- res_int[!is.na(res_int$padj) & res_int$padj < 0.05, ]

# Extract ENSEMBL gene IDs (rownames of DESeqResults)
genes_DE <- rownames(res_int_sig)

# Sanity check: number of DE genes
length(genes_DE)
# Expected: ~570 (should match DESeq2 summary above)

# ----------------------------- 2) Gene universe ------------------------------
# The gene universe should contain *all genes tested* in DESeq2,
# NOT only the significant ones.
# This avoids bias in enrichment analysis.

genes_universe <- rownames(res_int)

length(genes_universe)
# Expected: ~30,000 (mouse transcriptome scale)

# ----------------------------- 3) GO enrichment (Biological Process) ---------
# enrichGO performs over-representation analysis (ORA):
#   - gene: significant genes
#   - universe: background genes
#   - keyType: ENSEMBL IDs
#   - ont: BP = Biological Process
#   - pAdjustMethod: Benjamini–Hochberg FDR correction
#   - readable=TRUE: converts ENSEMBL IDs to gene SYMBOLs in results

ego_BP <- enrichGO(
  gene          = genes_DE,
  universe      = genes_universe,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Quick look at enrichment object
head(ego_BP)

# Convert enrichment result to a data.frame for inspection/export
ego_BP_df <- as.data.frame(ego_BP)

# Inspect key columns:
#   ID          : GO term ID
#   Description : GO term name
#   GeneRatio   : proportion of DE genes in this GO term
#   p.adjust    : FDR-adjusted p-value
head(
  ego_BP_df[, c("ID", "Description", "GeneRatio", "p.adjust")],
  10
)

# ----------------------------- 4) Dotplot visualization ----------------------
# dotplot shows:
#   - Y-axis: GO terms
#   - X-axis: GeneRatio
#   - Color: adjusted p-value
#   - Size: number of genes

dotplot(ego_BP, showCategory = 15) +
  ggplot2::ggtitle(
    "GO Biological Processes enriched in genotype × infection DE genes"
  )

# ----------------------------- 5) Improved readability -----------------------
# Long GO term names often overlap.
# We wrap text to improve readability of the y-axis labels.

p <- dotplot(ego_BP, showCategory = 15) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  ggtitle(
    "GO Biological Processes enriched in genotype × infection DE genes"
  )

p

###############################################################################
# End of GO enrichment analysis
###############################################################################
