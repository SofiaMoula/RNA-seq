#!/bin/bash

module load apps/trimmomatic
module load apps/multiqc
module load apps/anaconda3

# Ensure the Conda environment is correctly activated.
source activate rnaseq_env

# Directory containing the raw reads
RAW_DIR="path/to/file/Rawfiles"

# Output directories
PRE_FASTQC_DIR="${RAW_DIR}/preFASTQC"
POST_FASTQC_DIR="${RAW_DIR}/postFASTQC"
TRIMMOMATIC_DIR="${RAW_DIR}/trimmomatic"

# Ensure output directories exist
mkdir -p ${PRE_FASTQC_DIR} ${POST_FASTQC_DIR} ${TRIMMOMATIC_DIR}

# Log file
# Create log_QC working directory
mkdir log_QC
# Create analysis_log.txt file with touch command inside log_QC directory
LOG_FILE="${RAW_DIR}/log_QC/analysis_log.txt"

# Adaptor sequences
ADAPT1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPT2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Create a temporary adapters file
ADAPTERS_FILE="${TRIMMOMATIC_DIR}/adapters.fa"
echo ">Adapter1" > ${ADAPTERS_FILE}
echo "${ADAPT1}" >> ${ADAPTERS_FILE}
echo ">Adapter2" >> ${ADAPTERS_FILE}
echo "${ADAPT2}" >> ${ADAPTERS_FILE}

# Initialize log file
echo "Analysis Log" > ${LOG_FILE}
echo "============" >> ${LOG_FILE}

# Process each pair of FASTQ files
for file in ${RAW_DIR}/*_R1_001.fastq.gz; do
    # Construct file names
    base=$(basename ${file} _R1_001.fastq.gz)
    file1=${file}
    file2=${RAW_DIR}/${base}_R2_001.fastq.gz
    
    paired1=${TRIMMOMATIC_DIR}/${base}_R1_001_paired.fastq.gz
    unpaired1=${TRIMMOMATIC_DIR}/${base}_R1_001_unpaired.fastq.gz
    paired2=${TRIMMOMATIC_DIR}/${base}_R2_001_paired.fastq.gz
    unpaired2=${TRIMMOMATIC_DIR}/${base}_R2_001_unpaired.fastq.gz
    
    echo "Processing sample ${base}" >> ${LOG_FILE}
    
    # FASTQC before trimming
    fastqc -o ${PRE_FASTQC_DIR} ${file1} ${file2}
    
    # Trimmomatic for quality and adapter trimming
    trimmomatic PE -phred33 ${file1} ${file2} ${paired1} ${unpaired1} ${paired2} ${unpaired2} ILLUMINACLIP:${ADAPTERS_FILE}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    # FASTQC after Trimmomatic
    fastqc -o ${POST_FASTQC_DIR} ${paired1} ${paired2}
done

# MultiQC analysis
# Generate MultiQC report for preFASTQC 
multiqc -o ./multiqc_report_preFASTQC ./preFASTQC

# Generate MultiQC report for postFASTQC
multiqc -o ./multiqc_report_postFASTQC ./postFASTQC

echo "MultiQC report generated." >> ${LOG_FILE}

echo "QC Analysis completed." >> ${LOG_FILE}
