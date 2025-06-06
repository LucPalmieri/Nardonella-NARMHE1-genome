#!/usr/bin/env bash
# Run blastn against the Nardonella reference database for each FASTA file.

for F in *.fasta; do
    sampleName=$(basename "$F" .fasta)
    echo "$sampleName"
    blastn -db nardonella_ALLgenome.fa -query "$F" -outfmt 6 -out "$sampleName.tsv"
done
