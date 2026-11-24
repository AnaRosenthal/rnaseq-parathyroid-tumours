# Exploratory analysis and outlier detection at 24 h, using all samples:
# PCA, sample-to-sample distances and Cook's distance per sample.
# The control sample of patient 4 was identified as an outlier here (see README)
# and it is removed in dea_deseq2.R.
# Usage: run it from the Apartado2 folder with: Rscript qc_pca_cooks.R
# Outputs: QC/

library("DESeq2")
library("pheatmap")
library("RColorBrewer")

dir.create("./QC", showWarnings = FALSE)


# --- 1. Load the data ---
counts <- read.table("./input/rawcounts.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("./input/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
counts <- counts[, rownames(metadata)]   # same sample order in both tables


# --- 2. Keep the 24 h samples ---
metadata_24h <- metadata[metadata$time == "24h", ]
counts_24h <- counts[, rownames(metadata_24h)]

metadata_24h$patient <- factor(metadata_24h$patient)
metadata_24h$agent <- factor(metadata_24h$agent)


# --- 3. DESeq2 object and prefiltering ---
dds <- DESeqDataSetFromMatrix(
  countData = counts_24h,
  colData = metadata_24h,
  design = ~ patient + agent
)
keep <- rowSums(counts(dds)) >= 10   # keep genes with at least 10 counts in total
dds <- dds[keep, ]


# --- 4. PCA on variance-stabilised counts ---
vsd <- vst(dds, blind = TRUE)

pdf("./QC/PCA_agent_24h.pdf")
print(plotPCA(vsd, intgroup = "agent"))
dev.off()

pdf("./QC/PCA_patient_24h.pdf")
print(plotPCA(vsd, intgroup = "patient"))
dev.off()


# --- 5. Sample-to-sample distances ---
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$patient, vsd$agent, sep = "-")
colnames(sampleDistMatrix) <- paste(vsd$patient, vsd$agent, sep = "-")

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         filename = "./QC/sample_distances_heatmap.pdf")


# --- 6. Cook's distance per sample ---
dds <- DESeq(dds, test = "Wald")
cooks <- assays(dds)$cooks
sample_labels <- paste0("P", dds$patient, "-", dds$agent)   # patient and treatment, used as labels

# Median Cook's distance of each sample; outlier if above Q3 + 1.5 * IQR
cooks_median <- apply(cooks, 2, median, na.rm = TRUE)
Q1 <- quantile(cooks_median, 0.25)
Q3 <- quantile(cooks_median, 0.75)
threshold <- Q3 + 1.5 * (Q3 - Q1)
outliers <- cooks_median > threshold

outlier_table <- data.frame(
  sample = names(cooks_median),
  label = sample_labels,
  cooks_median = cooks_median,
  outlier = outliers
)
write.csv(outlier_table, "./QC/cooks_outliers_24h.csv", row.names = FALSE)

colnames(cooks) <- sample_labels
pdf("./QC/Cooks_distance_boxplot_24h.pdf")
par(mar = c(7, 4, 4, 2))
boxplot(cooks,
        las = 2,
        main = "Cook's distance per sample",
        outline = FALSE,
        ylim = c(0, 1),
        ylab = "Cook's distance")
abline(h = threshold, col = "red", lty = 2)
dev.off()
