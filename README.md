# Nardonella NARMHE1 Genome Analysis

This repository documents the steps used to assemble and analyse the draft genome of *Nardonella dryophthoridicola* NARMHE1, the bacterial endosymbiont of *Metamasius hemipterus*. The workflow follows the protocol described in the Microbiology Resources Announcement (Vol. 11, Issue 11, e00738‑22).

The repository contains:

- `RAxML/` – alignment and tree files used for phylogenetic analysis.
- `report_FAN30463_20211108_2203_21193522.pdf` and `report_FAL84985_20211214_2120_1bccc318.pdf` – Oxford Nanopore sequencing reports.
- `scripts/` – helper bash scripts for sequence conversion, BLAST searches and filtering.

The instructions below outline how to replicate the analysis.

## 1. Sorting bacterial reads from host sequences

1. **Convert FASTQ to FASTA** using [seqtk](https://github.com/lh3/seqtk):

   ```bash
   seqtk seq -a in.fastq.gz > out.fasta
   ```
   For many FASTQ files, run `scripts/bash_converter.sh`.

2. **Create a local BLAST database** from reference *Nardonella* genomes:

   ```bash
   makeblastdb -in nardonella_ALLgenome.fa -dbtype nucl -parse_seqids
   ```

3. **Search for matches** with `blastn` (single file example):

   ```bash
   blastn -db nardonella_ALLgenome.fa -query sample.fasta -outfmt 6 -out sample.tsv
   ```
   To process multiple FASTA files, use `scripts/bash_finder.sh`.

4. **Concatenate results** and filter alignments with an E‑value ≤ 1e‑6:

   ```bash
   cat *.tsv > blastn_nardonella_sequences.tsv
   awk '{ if($11 <= 1e-6) print }' blastn_nardonella_sequences.tsv > blastn_nardonella_sequences_filtered.tsv
   ```

5. **Extract matching reads** with `seqtk` (single file example):

   ```bash
   seqtk subseq sample.fastq name.lst > sample_nardo.fastq
   ```
   For batches of FASTQ files, run `scripts/bash_filter.sh` and then merge the filtered reads:

   ```bash
   cat *_nardo.fastq > ONT_merged.fastq
   ```

## 2. File correction and quality control

- **Repair paired-end Illumina reads** (if needed) using [BBMap](https://github.com/BioInfoTools/BBMap):

  ```bash
  repair.sh in1=illumina_uncorrected_R1_nardo.fastq.gz \
           in2=illumina_uncorrected_R2_nardo.fastq.gz \
           out1=illumina_corrected_R1_nardo.fastq.gz \
           out2=illumina_corrected_R2_nardo.fastq.gz \
           outsingle=unpaired.fastq.gz
  ```

- **Remove duplicates** from Nanopore reads with `dedupe.sh`:

  ```bash
  dedupe.sh in=ONT_merged.fastq out=deduplicated.fastq ac=f -da
  ```

- **Trim short reads** and perform light QC with [Filtlong](https://github.com/rrwick/Filtlong):

  ```bash
  filtlong --min_length 500 --keep_percent 95 deduplicated.fastq > deduplicated.fastq
  ```

## 3. Genome assembly

### Canu

```bash
canu -p nardonella -d canu_assembly genomeSize=0.2m \
     correctedErrorRate=0.155 minReadLength=500 \
     minOverlapLength=250 -nanopore ONT_deduplicated.fastq
```

### Flye

```bash
conda activate medaka
flye --nano-hq ONT_deduplicated.fastq -g 0.2m -t 4 --out-dir ./flye_assembly
```

### Polishing

1. **Medaka** (long reads):

   ```bash
   medaka_consensus -i ONT_merged.fastq -d canu_contigs.fasta \
                   -o long_canu_corrected -m r103_hac_g507 -t 4
   medaka_consensus -i ONT_merged.fastq -d flye_contigs.fasta \
                   -o long_flye_corrected -m r103_hac_g507 -t 4
   ```

2. **Short read polishing** with [Polypolish](https://github.com/rrwick/Polypolish):

   ```bash
   bwa index medaka_consensus.fasta
   bwa mem -t 4 -a medaka_consensus.fasta illumina_corrected_R1_nardo.fastq.gz > alignments_R1.sam
   bwa mem -t 4 -a medaka_consensus.fasta illumina_corrected_R2_nardo.fastq.gz > alignments_R2.sam
   polypolish medaka_consensus.fasta alignments_R1.sam alignments_R2.sam > polished.fasta
   ```

3. **Additional polishing** with [POLCA](https://github.com/alekseyzimin/masurca):

   ```bash
   polca.sh -a consensus_polished.fasta -r "illumina_corrected_R1_nardo.fastq.gz illumina_corrected_R2_nardo.fastq.gz" -t 4
   ```

### Assembly summaries

| Assembly | Substitution errors | Indel errors | Size (bp) | Consensus quality |
|---------|--------------------|-------------|----------|-------------------|
| **Canu** | 89 | 9 | 205,527 | 99.9523 |
| **Flye** | 37 | 38 | 167,594 | 99.9552 |

### Combining contigs

Merge the polished assemblies with [Quickmerge](https://github.com/mahulchak/quickmerge):

```bash
merge_wrapper.py canu_assembly.fasta flye_assembly.fasta
```

The merge yielded seven contigs which were polished again. Final statistics:

| Substitution errors | Indel errors | Size (bp) | Consensus quality |
|--------------------|-------------|----------|-------------------|
| 62 | 46 | 206,311 | 99.9477 |

## 4. Scaffolding with reference genomes

Use [RagTag](https://github.com/malonge/RagTag) to correct and scaffold contigs:

```bash
ragtag.py correct nardonella_ALLgenome.fa merged_canu_flye_contigs.fasta -u -o ./corrected -t 4 -R ONT_deduplicated.fastq -T ont
ragtag.py scaffold Nardonella_RFE.fasta corrected/ragtag.correct.fasta -C -t 4 -o ./scaffold_RFEprotein
```

## 5. Gene annotation

The final assembly (`Nardonella_MHE.fasta`) was submitted to GenBank (BioProject PRJNA798699). Accession numbers are listed below:

| Contig | Length | Accession |
|-------|-------|-----------|
| contig01 | 12,016 bp | JAKMAI010000001 |
| contig02 | 28,356 bp | JAKMAI010000002 |
| contig03 | 16,370 bp | JAKMAI010000003 |
| contig04 | 6,418 bp | JAKMAI010000004 |
| contig05 | 12,133 bp | JAKMAI010000005 |
| contig06 | 13,069 bp | JAKMAI010000006 |
| contig07 | 91,945 bp | JAKMAI010000007 |

## 6. Phylogenetic analysis

The `RAxML` directory contains the alignment (`Mauve_All_Nardonella_for_RAxML.phy`) and tree files generated using RAxML.

---

These notes and scripts should allow reproduction of the assembly and analysis of the *Nardonella* NARMHE1 genome.
