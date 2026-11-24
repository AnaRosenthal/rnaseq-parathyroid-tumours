# parathyroid-tumours-rnaseq

RNA-seq analysis of primary cultures of parathyroid tumours treated with DPN (an estrogen receptor agonist) or 4-hydroxytamoxifen (OHT, a competitive estrogen receptor inhibitor) for 24 or 48 hours. Part 1 covers read quality control, trimming, alignment and quantification; Part 2 covers sample quality control, differential expression and gene set enrichment analysis (GSEA).

This is my solution to the final exercise of the transcriptomics module (Master's programme). The repository is a fork of [bioinfo-lessons/transcriptomic-final-exercise](https://github.com/bioinfo-lessons/transcriptomic-final-exercise), which provides the assignment ([docs/ASSIGNMENT.md](docs/ASSIGNMENT.md)) and the input data. My full report (in Spanish) is in [docs/report.pdf](docs/report.pdf).

## Data

The experiment consists of 24 primary cultures from ER-alpha-negative parathyroid tumours (4 patients), treated with DPN or OHT for 24 h or 48 h, plus controls (original study: Haglund et al., 2012; sequencing data in SRA).

- **Part 1:** two paired-end samples (SRR479052 and SRR479054) restricted to human chromosome 21 (GRCh38, Ensembl release 109).
- **Part 2:** raw count matrix of the 24 cultures, sample metadata, and a GMT file with the DPN 48 h gene sets (`DPN_perturbed` and `DPN_unperturbed`).

## Repository structure

    Apartado1/
      input/                     FASTQ files, genome FASTA and GTF (provided)
      qc_align_quantify.sh       QC, trimming, alignment and counting
      results/                   logs, QC reports, alignment summaries and counts
    Apartado2/
      input/                     raw counts, metadata and GMT (provided)
      qc_pca_cooks.R             PCA, sample distances and Cook's distance (outlier detection)
      dea_deseq2.R               differential expression at 24 h (DESeq2)
      gsea_create_rnk.R          ranked gene list for GSEA Preranked
      QC/  DEA/  GSEA/           outputs of the scripts above (and the GSEA run)
    docs/                        assignment and report
    environment_apartado1.yml    conda environment for Part 1
    environment_apartado2.yml    conda environment for Part 2

## Methods summary

**Part 1.** Quality control of the raw reads (FastQC, MultiQC); adapter and quality trimming (Cutadapt, `-q 20 -m 30`, Illumina adapters); indexing of chromosome 21 and paired-end alignment (HISAT2); sorting and indexing (SAMtools); read counting per gene (HTSeq, `--stranded=no`, `--mode=intersection-nonempty`, `--minaqual=10`).

**Part 2.** Analysis restricted to the 24 h samples. A PCA and the Cook's distance of each sample (median per sample above Q3 + 1.5 x IQR) identified the control sample of patient 4 as an outlier (see `Apartado2/QC/`). Only that sample was removed, leaving control n = 3, DPN n = 4 and OHT n = 4. Differential expression with DESeq2 (design `~ patient + agent`, genes with fewer than 10 total counts removed, Wald test, Benjamini-Hochberg adjusted p-value < 0.05). For GSEA, genes were ranked by the Wald statistic of DPN vs control at 24 h and tested with GSEA Preranked against the DPN 48 h gene sets.

## Main results

- **Alignment:** overall alignment rate of 91.3% and 90.8% for the two samples. About 48% of the aligned reads were assigned to genes by HTSeq, and about 35% were discarded for low alignment quality.
- **OHT vs control (24 h):** no genes with adjusted p-value < 0.05.
- **DPN vs control (24 h):** 35 genes with adjusted p-value < 0.05, none with |log2 fold change| > 1.
- **GSEA (DPN, 24 h):** `DPN_perturbed` is positively enriched (NES = 2.57, FDR q-value = 0.0) and `DPN_unperturbed` is negatively enriched (NES = -2.99, FDR q-value = 0.0). So, although no individual gene changes strongly at 24 h, DPN already shifts the expression of the genes that define its 48 h response.

## How to run

Part 1:

    conda env create -f environment_apartado1.yml
    conda activate rnaseq-apartado1
    cd Apartado1
    bash qc_align_quantify.sh

Part 2 (from the `Apartado2` folder):

    conda env create -f environment_apartado2.yml
    conda activate rnaseq-apartado2
    cd Apartado2
    Rscript qc_pca_cooks.R
    Rscript dea_deseq2.R
    Rscript gsea_create_rnk.R

The last script writes `GSEA/DPN_24h.rnk`. The enrichment itself was run with the GSEA desktop application (Preranked tool) using this file and `input/DPN_response.gmt`; the run reported in the report is in `Apartado2/GSEA/DPN_24h_preranked/` (its parameters are in the `.rpt` file).

## Software

| Software | Version |
|---|---|
| FastQC | 0.12.1 |
| MultiQC | 1.33 |
| Cutadapt | 5.2 |
| HISAT2 | 2.2.2 |
| SAMtools | 1.23.1 |
| HTSeq | 2.0.9 |
| R | 4.5.3 |
| DESeq2 | 1.50.2 |
| EnhancedVolcano | 1.28.2 |
| mygene | 1.46.0 |
| pheatmap | 1.0.13 |
| GSEA (desktop) | 4.4.0 |

## Notes

- Gene symbols in the differential expression tables are obtained with an online query (`mygene`), so they may change slightly over time.
- The ranked list for GSEA was computed from a model that only includes the treatment (`~ agent`); adding the patient as a blocking factor, as in the differential expression step, would be a natural refinement. The ranked list used in the report was generated with R 4.4.3 and DESeq2 1.46.0.
- The assignment and the input files belong to the course; the code, results and report are my own work.
