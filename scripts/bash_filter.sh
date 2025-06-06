#!/usr/bin/env bash
# Filter sequences based on blastn results for each FASTQ file.

for fq in *.fastq; do
    sampleName=$(basename "$fq" .fastq)
    echo "$sampleName"
    seqtk subseq "$fq" "$sampleName.tsv" > "${sampleName}_nardo.fastq"
done
