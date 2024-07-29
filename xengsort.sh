#!/bin/bash


index="path/to/file/index_xen"
read1="path/to/file/R1_001_paired.fastq.gz"
read2="path/to/file/R2_001_paired.fastq.gz"
prefix="path/to/file/results"

xengsort classify --index ${index} --fastq ${read1} --pairs ${read2} --prefix ${prefix} --mode count

