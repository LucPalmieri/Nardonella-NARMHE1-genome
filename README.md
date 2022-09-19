## Supplementary information - Draft Genome Sequence of *Nardonella dryophthoridicola* NARMHE1, Endosymbiont of *Metamasius hemipterus*.
You will find below the codes, the list of programs, and parameters I used to generate and analyze the draft genome and annotation of *Nardonella* strain NARMHE1. In this repository you can also find Oxford Nanopore run reports and a folder containing the alignement used to run the phylogenetic analysis. For more information see the [publication, under review]

## 1 - Sorting bacteria from host sequences.

First, convert target fastq.gz to fasta files using 
[Seqtk](https://github.com/lh3/seqtk)
**only fasta files can be used on blastn**

For Illumina reads:

	seqtk seq -a in.fastq.gz > out.fasta (generic command)
	seqtk seq -a MTMO002_S2_L002_R1_001.fastq.gz > ovaryF.fasta
	seqtk seq -a MTMO002_S2_L002_R2_001.fastq.gz > ovaryR.fasta

For multiple Nanopore fastq files do a batch conversion using the script below.

	nano bash_converter.sh

Copy and paste code on the bash_converter.sh file.

	#!/usr/bin/env bash
		
	#running another program with the shell: seqtk
	#This command will convert multiple fastq files to the fasta format
	
	for fq in *.fastq
	do
		sampleName=$(echo $fq | cut -f 1 -d '.')
		echo $sampleName
		seqtk seq -A $fq > $sampleName\.fasta

	done

Create a local blast database using Nardonella reference genomes (downloaded previously from NCBI).

[Blastn](https://github.com/ncbi/blast_plus_docs)

	makeblastdb -in nardonella_ALLgenome.fa -dbtype nucl -parse_seqids

Search for genome matches on target fasta files using blastn.

For Illumina reads:

	blastn -db nardonella_database -query in.fasta -out results.out (generic command)
	blastn -db nardonella_genome.fasta -query ovaryF.fasta -outfmt 6 -out blastn_ovary_forward.tsv
	blastn -db nardonella_genome.fasta -query ovaryR.fasta -outfmt 6 -out blastn_ovary_reverse.tsv

Batch blast search on multiple Nanopore fasta files.

	nano bash_finder.sh

Copy and paste code on the bash_finder.sh file.

	#!/usr/bin/env bash
	
	#running another program with the shell: blastn
	#This command will search multiple fasta files for database query matches
	
	for F in *.fasta
	do
		sampleName=$(echo $F | cut -f 1 -d '.')
		echo $sampleName
		blastn -db nardonella_ALLgenome.fa -query $F -outfmt 6 -out $sampleName\.tsv
	
	done

Concatenate all Nanopore .tsv files on one single file.

	cat *.tsv > blastn_nardonella_sequences.tsv

Filter out the lines with E-value greater than 10<sup>-6</sup>.
This should be done separately on Illumina and Nanopore reads.

	awk '{ if($11 <= 0.000001 ) {print}}' blastn_nardonella_sequences.tsv -> blastn_nardonella sequences_filtered.tsv

Filter out blast identified Nardonella sequences from original fastq.gz files.
This will finally sort out bacteria sequences from host beetle sequences.

Extract sequences with names in file name.lst, one sequence name per line:

For Illumina reads:

	seqtk subseq in.fq name.lst > out.fq (generic command)
	seqtk subseq MTMO002_S2_L002_R1_001.fastq.gz ovaryF_nardo_list.csv > MTMO002_R1_nardo_only.fastq.gz
	seqtk subseq MTMO002_S2_L002_R2_001.fastq.gz ovaryR_nardo_list.csv > MTMO002_R2_nardo_only.fastq.gz


Batch extraction for multiple Nanopore files.

	nano bash_filter.sh

Copy and paste code on the bash_fiter.sh file.

	#!/usr/bin/env bash
	
	#running another program with the shell: seqtk
	#This command will find the sequences on the query .tsv file and filter them out on a new file
	
	for fq in *.fastq
	do
		sampleName=$(echo $fq | cut -f 1 -d '.')
		echo $sampleName
		seqtk subseq $fq $sampleName.tsv > $sampleName\_nardo.fastq
		
	done

Concatenate all the sorted Nanopore fastq files on one single file.

	cat *_nardo.fastq > ONT_merged.fastq

Sequence sorting is complete.
We now have two files (forward and reverse) for Illumina and/or one for Nanopore each containing only endosymbiont reads.

## 2 -  File correction and quality control.

Following sorting the reads need some attention.
Illumina forward an reverse files might not have the same number of reads, and reads might be in random order.
We need to re-pair the short sequences and remove unpaired with [BBmap](https://github.com/BioInfoTools/BBMap).

For paired-end reads in two files:

	repair.sh in1=illumina_uncorrected_R1_nardo.fastq.gz.fq in2=illumina_uncorrected_R2_nardo.fastq.gz out1=illumina_corrected_R1_nardo.fastq.gz out2=illumina_corrected_R2_nardo.fastq.gz outsingle=unpaired.fastq.gz

The Nanopore file might have very short and duplicated sequences.
To remove duplicated sequences using [BBmap] tool depupe:

	dedupe.sh in=ONT_merged.fastq out=deduplicated.fastq ac=f -da

To remove short reads and do a light QC on Nanopore reads using [Filtlong](https://github.com/rrwick/Filtlong):

	filtlong --min_length 500 --keep_percent 95 deduplicated.fastq > deduplicated.fastq

## 3 - Assembling the Nardonella genome.

### Canu assembly.
[Canu](https://github.com/marbl/canu)

	canu -p nardonella -d canu_assembly genomeSize=0.2m correctedErrorRate=0.155 minReadLength=500 minOverlapLength=250 -nanopore ONT_deduplicated.fastq

### Flye assembly.
[Flye](https://github.com/fenderglass/Flye)

	conda activate medaka

	flye --nano-hq ONT_deduplicated.fastq -g 0.2m -t 4 --out-dir ./flye_assembly

### Polishing contigs.

Correcting with long reads only using [Medaka](https://github.com/nanoporetech/medaka)

	conda activate medaka
	medaka_consensus -i ONT_merged.fastq -d canu_contigs.fasta -o long_canu_corrected -m r103_hac_g507 -t 4
	medaka_consensus -i ONT_merged.fastq -d flye_contigs.fasta -o long_flye_corrected -m r103_hac_g507 -t 4

### Additional polishing with short reads.

First short read polishing with [polypolish](https://github.com/rrwick/Polypolish)

	bwa index medaka_consensus.fasta
	bwa mem -t 4 -a medaka_consensus.fasta illumina_corrected_R1_nardo.fastq.gz > alignments_R1.sam
	bwa mem -t 4 -a medaka_consensus.fasta illumina_corrected_R2_nardo.fastq.gz > alignments_R2.sam
	polypolish medaka_consensus.fasta alignments_R1.sam alignments_R2.sam > polished.fasta

Extra polishing with [POLCA](https://github.com/alekseyzimin/masurca)

	polca.sh -a consensus_polished.fasta -r 'illumina_corrected_R1_nardo.fastq.gz illumina_corrected_R2_nardo.fastq.gz' -t 4

|**Canu assembly summary after POLCA**|
|-----------------|
|Substitution Errors: 89|
|Insertion/Deletion Errors: 9|
|Assembly Size: 205527|
|Consensus Quality: 99.9523|

|**Flye assembly summary after POLCA**|
|----------------|
|Substitution Errors: 37|
|Insertion/Deletion Errors: 38|
|Assembly Size: 167594|
|Consensus Quality: 99.9552|

### Combining the contigs of the two assemblies.

Both Canu and Flye failed to recover a single contig.
The contigs were combined using [Quickmerge](https://github.com/mahulchak/quickmerge).

	merge_wrapper.py canu_assembly.fasta flye_assembly.fasta

|Polished sequence (to stdout):|
|-----------|
|  contig_20_polypolish (16,342 bp)|
| contig_19_polypolish (14,571 bp)|
|  contig_17_polypolish (40,101 bp)|
|  seq000_polypolish (95,525 bp)|
|  contig_16_polypolish (18,125 bp)|
|  contig_18_polypolish (15,230 bp)|
|  seq006_polypolish (6,417 bp)|

The merging process resulted in seven contigs.
I passed the combined file through one extra round of polishing with long and short reads.

|**Combined assembly summary after POLCA**|
|---------|
|Substitution Errors: 62|
|Insertion/Deletion Errors: 46|
|Assembly Size: 206311|
|Consensus Quality: 99.9477|

## 4 - Scaffolding using Nardonella reference genomes.

Using homology between the sequences and reference genomes to identify and correct potential misassemblies with [RagTag](https://github.com/malonge/RagTag).

	conda activate python3

	ragtag.py correct nardonella_ALLgenome.fa merged_canu_flye_contigs.fasta -u -o ./corrected -t 4 -R ONT_deduplicated.fastq -T ont

### Scafolding the corrected file with reference nardonella proteins.

	ragtag.py scaffold Nardonella_RFE.fasta corrected/ragtag.correct.fasta -C -t 4 -o ./scaffold_RFEprotein

## 5 - Gene annotation with NCBI [prokaryot genome annotation pipeline ](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/).

The final corrected, scaffolded, and curated assemble file (Nardonella_MHE.fasta) was submitted to GenBank for official annotation.
The contigs have been deposited in GenBank accession number [JAKMAI010000000](https://www.ncbi.nlm.nih.gov/nuccore/JAKMAI000000000.1/) under the BioProject accession number PRJNA798699.

|**ID**|**Lenght**|**Acession number**|
|----------|--------|---------|
|contig01|12,016 bp|JAKMAI010000001| 
|contig02|28,356 bp|JAKMAI010000002| 
|contig03|16,370 bp|JAKMAI010000003| 
|contig04|6,418 bp|JAKMAI010000004| 
|contig05|12,133 bp|JAKMAI010000005| 
|contig06|13,069 bp|JAKMAI010000006| 
|contig00|91,945 bp|JAKMAI010000007|
