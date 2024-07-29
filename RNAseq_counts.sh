#!/bin/bash

# Load necessary modules
module load apps/samtools
module load apps/subread
module load apps/R

# Define directories and files
GTF_FILE="path/to/file/Homo_sapiens.GRCh38.104.gtf"
OUTPUT_DIR="path/to/file/counts"
BAM_DIR="path/to/file/bam_files"  # Directory containing the BAM files

# Counting with featureCounts for each BAM file
for bamfile in ${BAM_DIR}/*_sorted.bam; do
    base=$(basename $bamfile _sorted.bam)
    echo "Processing: $base"
    
    featureCounts -a ${GTF_FILE} -o ${OUTPUT_DIR}/counts/${base}_counts.txt -t exon -g gene_id -s 2 -p -B -C ${bamfile}
done

echo "All samples have been processed."