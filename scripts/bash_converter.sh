#!/usr/bin/env bash
# Convert all FASTQ files in the current directory to FASTA using seqtk.

for fq in *.fastq; do
    sampleName=$(basename "$fq" .fastq)
    echo "$sampleName"
    seqtk seq -A "$fq" > "$sampleName.fasta"
done
