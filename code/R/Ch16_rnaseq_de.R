# =====================================================================
# Chapter 16 - Gene Expression (RNA-seq): differential expression (R)
# The standard DESeq2 workflow on the simulated count data.
#
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#   source("rnaseq_de.R")
#
# Data: counts.csv (2000 genes x 8 samples), metadata.csv, truth.csv.
# =====================================================================
suppressMessages(library(DESeq2))

counts <- as.matrix(read.csv("counts.csv", row.names = 1))
meta   <- read.csv("metadata.csv", row.names = 1)
meta$condition <- factor(meta$condition, levels = c("control", "treated"))

# 1. build the DESeq2 dataset (genes x samples, design = ~ condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta,
                              design = ~ condition)

# 2. (optional but recommended) drop almost-empty genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 3. normalization, dispersion shrinkage, and the per-gene test -- one call
dds <- DESeq(dds)
sizeFactors(dds)                              # per-sample size factors

# 4. results: treated vs control
res <- results(dds, contrast = c("condition", "treated", "control"))
summary(res)                                  # up / down at FDR 5%
res <- res[order(res$padj), ]
head(as.data.frame(res))

# 5. shrink log2 fold changes for ranking / plotting
resLFC <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")

# 6. call DE genes and check against the known truth
sig <- which(!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1)
cat("DE genes (padj<0.05 & |LFC|>1):", length(sig), "\n")
truth <- read.csv("truth.csv", row.names = 1)
tp <- sum(truth[rownames(res)[sig], "true_DE"] == 1)
cat("true positives:", tp, "  precision:", round(tp/length(sig), 2), "\n")

# 7. plots
plotMA(res, ylim = c(-4, 4))                  # MA plot
plotPCA(vst(dds), intgroup = "condition")     # sample PCA (needs vst)
# volcano (base R):
with(as.data.frame(res),
     plot(log2FoldChange, -log10(pvalue), pch = 20, col = "grey",
          xlab = "log2 fold change", ylab = "-log10 p"))
# enrichment next:  clusterProfiler::enrichGO(rownames(res)[sig], ...)

cat("\nDone.\n")
