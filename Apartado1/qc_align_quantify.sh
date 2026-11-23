#!/bin/bash
# RNA-seq processing of two paired-end samples (SRR479052 and SRR479054)
# restricted to human chromosome 21 (GRCh38, Ensembl release 109).
# Steps: QC -> adapter trimming -> HISAT2 index and alignment -> samtools -> HTSeq counts -> MultiQC.
# Usage: run it from the Apartado1 folder with: bash qc_align_quantify.sh
# Required tools: fastqc, multiqc, cutadapt, hisat2, samtools, htseq-count

DATA=input
REF=input
OUT=results

mkdir -p $OUT/fastqc/raw $OUT/fastqc/trimmed $OUT/trimmed $OUT/cutadapt_logs $OUT/hisat2 $OUT/htseq $OUT/multiqc


# 1. Quality control of the raw reads
fastqc $DATA/*.fastq -o $OUT/fastqc/raw
multiqc $OUT/fastqc/raw -o $OUT/multiqc/raw


# 2. Adapter and quality trimming
# -q 20: trim low-quality 3' ends (Phred < 20); -m 30: discard reads shorter than 30 bp
# -a / -A: Illumina adapters of read 1 / read 2
for sample in SRR479052 SRR479054
do
    cutadapt \
    -q 20 \
    -m 30 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o $OUT/trimmed/${sample}_1.trimmed.fastq \
    -p $OUT/trimmed/${sample}_2.trimmed.fastq \
    $DATA/${sample}.chr21_1.fastq \
    $DATA/${sample}.chr21_2.fastq \
    > $OUT/cutadapt_logs/${sample}.cutadapt.log
done

# QC of the trimmed reads
fastqc $OUT/trimmed/*.fastq -o $OUT/fastqc/trimmed
multiqc $OUT/fastqc/trimmed -o $OUT/multiqc/trimmed


# 3. Genome index (chromosome 21)
hisat2-build --seed 123 -p 2 \
$REF/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
$OUT/hisat2/chr21_index


# 4. Alignment with HISAT2
# -k 1: report at most one alignment per read; -p 2: two threads
for sample in SRR479052 SRR479054
do
    hisat2 \
    --new-summary \
    --summary-file $OUT/hisat2/${sample}.summary \
    --seed 123 \
    --phred33 -p 2 -k 1 \
    -x $OUT/hisat2/chr21_index \
    -1 $OUT/trimmed/${sample}_1.trimmed.fastq \
    -2 $OUT/trimmed/${sample}_2.trimmed.fastq \
    -S $OUT/hisat2/${sample}.sam
done

# QC of the alignment
multiqc $OUT/hisat2 -o $OUT/multiqc/alignment


# 5. SAM -> sorted and indexed BAM
for sample in SRR479052 SRR479054
do
    samtools view -bS $OUT/hisat2/${sample}.sam > $OUT/hisat2/${sample}.bam
    samtools sort $OUT/hisat2/${sample}.bam -o $OUT/hisat2/${sample}.sorted.bam
    samtools index $OUT/hisat2/${sample}.sorted.bam
done


# 6. Read counting with HTSeq
# --stranded=no: the library is not strand-specific
# --mode=intersection-nonempty and --minaqual=10: skip reads with alignment quality < 10
# --type=exon and --idattr=gene_id: count reads over exons, summed per gene
for sample in SRR479052 SRR479054
do
    htseq-count \
    --format=bam \
    --stranded=no \
    --mode=intersection-nonempty \
    --minaqual=10 \
    --type=exon \
    --idattr=gene_id \
    --additional-attr=gene_name \
    $OUT/hisat2/${sample}.sorted.bam \
    $REF/Homo_sapiens.GRCh38.109.chr21.gtf \
    > $OUT/htseq/${sample}.counts.txt
done

# QC of the counts
multiqc $OUT/htseq -o $OUT/multiqc/counts


# 7. Final combined QC report
multiqc $OUT/fastqc/raw $OUT/fastqc/trimmed $OUT/hisat2 $OUT/htseq -o $OUT/multiqc/full_report
