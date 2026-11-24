# Differential expression analysis with DESeq2 at 24 h: OHT vs Control and DPN vs Control.
# The control sample of patient 4 is removed because it was an outlier (see README).
# Usage: run it from the Apartado2 folder with: Rscript dea_deseq2.R
# Outputs: DEA/DEG_tables, DEA/volcano_plots and DEA/volcano_tables

library("DESeq2")
library("EnhancedVolcano")
library("mygene")

dir.create("./DEA/DEG_tables", recursive = TRUE, showWarnings = FALSE)
dir.create("./DEA/volcano_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("./DEA/volcano_tables", recursive = TRUE, showWarnings = FALSE)


# --- 1. Load the data ---
counts <- read.table("./input/rawcounts.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("./input/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")


# --- 2. Keep the 24 h samples and remove the control of patient 4 ---
metadata_24h <- metadata[metadata$time == "24h", ]
metadata_24h <- metadata_24h[!(metadata_24h$patient == 4 & metadata_24h$agent == "Control"), ]
counts_24h <- counts[, rownames(metadata_24h)]

metadata_24h$patient <- factor(metadata_24h$patient)
metadata_24h$agent <- factor(metadata_24h$agent)


# --- 3. DESeq2 object (patient as blocking factor) and prefiltering ---
dds <- DESeqDataSetFromMatrix(
  countData = counts_24h,
  colData = metadata_24h,
  design = ~ patient + agent
)
keep <- rowSums(counts(dds)) >= 10   # keep genes with at least 10 counts in total
dds <- dds[keep, ]


# --- 4. Fit the model (Wald test) and save the dispersion plot ---
dds <- DESeq(dds, test = "Wald")
pdf("./DEA/dispersion_estimates.pdf")
plotDispEsts(dds)
dev.off()


# --- 5. Results: OHT vs Control and DPN vs Control (alpha 0.05, Benjamini-Hochberg) ---
res_OHT <- results(object = dds,
                   contrast = c("agent", "OHT", "Control"),
                   alpha = 0.05,
                   pAdjustMethod = "BH",
                   tidy = TRUE)

res_DPN <- results(object = dds,
                   contrast = c("agent", "DPN", "Control"),
                   alpha = 0.05,
                   pAdjustMethod = "BH",
                   tidy = TRUE)


# --- 6. Replace Ensembl IDs by gene symbols (online query with mygene) ---
# If no symbol is found, the Ensembl ID is kept.
add_symbols <- function(res) {
  ids <- queryMany(res$row, scopes = "ensembl.gene", fields = "symbol", species = "human")
  ids <- ids[!duplicated(ids$query), ]
  res$row <- ifelse(is.na(ids$symbol), ids$query, ids$symbol)
  res
}
res_OHT <- add_symbols(res_OHT)
res_DPN <- add_symbols(res_DPN)


# --- 7. Significant genes (adjusted p-value < 0.05) ---
res_OHT_sig <- res_OHT[which(res_OHT$padj < 0.05), ]
res_DPN_sig <- res_DPN[which(res_DPN$padj < 0.05), ]

write.csv(res_OHT_sig[order(res_OHT_sig$padj), ],
          file = "./DEA/DEG_tables/OHT_vs_Control_24h_significant.csv",
          row.names = FALSE)
write.csv(res_DPN_sig[order(res_DPN_sig$padj), ],
          file = "./DEA/DEG_tables/DPN_vs_Control_24h_significant.csv",
          row.names = FALSE)


# --- 8. Volcano plots (genes without adjusted p-value are removed) ---
res_OHT_clean <- res_OHT[!is.na(res_OHT$padj), ]
res_DPN_clean <- res_DPN[!is.na(res_DPN$padj), ]

pdf("./DEA/volcano_plots/Volcano_OHT_24h.pdf", width = 8, height = 8)
print(EnhancedVolcano(res_OHT_clean,
                      lab = res_OHT_clean$row,
                      x = "log2FoldChange",
                      y = "padj",
                      title = "OHT vs Control (24h)",
                      subtitle = NULL,
                      FCcutoff = 1,
                      pCutoff = 0.05,
                      boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      labSize = 6.0))
dev.off()

pdf("./DEA/volcano_plots/Volcano_DPN_24h.pdf", width = 8, height = 8)
print(EnhancedVolcano(res_DPN_clean,
                      lab = res_DPN_clean$row,
                      x = "log2FoldChange",
                      y = "padj",
                      title = "DPN vs Control (24h)",
                      subtitle = NULL,
                      FCcutoff = 1,
                      pCutoff = 0.05,
                      boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      labSize = 6.0,
                      xlim = c(-3, 3)))
dev.off()


# --- 9. Full tables behind the volcano plots, sorted by adjusted p-value ---
volcano_table_OHT <- res_OHT[, c("row", "log2FoldChange", "padj")]
volcano_table_OHT <- volcano_table_OHT[order(volcano_table_OHT$padj), ]
write.csv(volcano_table_OHT,
          file = "./DEA/volcano_tables/Volcano_OHT_24h_table.csv",
          row.names = FALSE)

volcano_table_DPN <- res_DPN[, c("row", "log2FoldChange", "padj")]
volcano_table_DPN <- volcano_table_DPN[order(volcano_table_DPN$padj), ]
write.csv(volcano_table_DPN,
          file = "./DEA/volcano_tables/Volcano_DPN_24h_table.csv",
          row.names = FALSE)
