# Create the ranked gene list (.rnk) for GSEA Preranked: DPN vs Control at 24 h.
# The control sample of patient 4 is removed because it was an outlier (see README).
# Usage: run it from the Apartado2 folder with: Rscript gsea_create_rnk.R
# Output: GSEA/DPN_24h.rnk

library(DESeq2)

dir.create("./GSEA", showWarnings = FALSE)


# --- 1. Load the data ---
metadata <- read.delim("./input/metadata.tsv", row.names = 1)
counts <- read.delim("./input/rawcounts.tsv", row.names = 1)

metadata$agent <- factor(metadata$agent)
metadata$time <- factor(metadata$time)
metadata$patient <- factor(metadata$patient)


# --- 2. Keep 24 h samples, DPN and Control only, and remove the control of patient 4 ---
metadata_filt <- metadata[metadata$time == "24h" &
                          metadata$agent %in% c("DPN", "Control") &
                          !(metadata$patient == "4" & metadata$agent == "Control"), ]
metadata_filt <- droplevels(metadata_filt)
counts_filt <- counts[, rownames(metadata_filt)]


# --- 3. DESeq2: DPN vs Control ---
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = metadata_filt,
  design = ~ agent
)
dds <- DESeq(dds)

res <- results(dds, alpha = 0.05, contrast = c("agent", "DPN", "Control"))
summary(res)


# --- 4. Ranked list: Wald statistic of each gene ---
rnk <- data.frame(Feature = rownames(res), Stat = res$stat)
rnk$Feature <- sub("\\..*$", "", rnk$Feature)   # remove the Ensembl version suffix
rnk <- rnk[!is.na(rnk$Stat), ]                  # remove genes without statistic
rnk <- rnk[order(rnk$Stat, decreasing = TRUE), ]

# GSEA expects no header and tab-separated columns
write.table(rnk,
            file = "./GSEA/DPN_24h.rnk",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)
