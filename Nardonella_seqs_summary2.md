### Nardonella genome assembly strain NardMH *Metamasius hemipterus*
Following the methodology of Chouaia *et al*. **Complete Genome Sequence of Rhynchophorus ferrugineus Endocytobiont “*Candidatus Nardonella dryophthoridicola*” Strain NardRF**. Microbiol Resour Announc. 2021 Jul; 10(26): e00355-21.

Illumina output summary:

|Gut|Larva|Ovary|
|:---:|:---:|:---:|
forward = 55.1 M seqs|forward = 67.1 M seqs|forward = 70.0 M seqs|
|reverse = 55.1 M seqs|reverse = 67.1 M seqs|reverse = 70.0 M seqs|
|TOTAL = 110.2 M seqs|TOTAL = 134.2 M seqs|TOTAL = 140.0 M seqs|

I used BLASTn (E value cutoff, 10<sup>-8</sup>) to identify Nardonella sequences present on raw illumina reads (all four Nardonella genomes used as reference database).

* Convert target fastq.gz to fasta files (only fasta files can be used on blastn).


    seqtk seq -a MTMG001_S2_L002_R1_001.fastq.gz > gutF.fasta
    seqtk seq -a MTMG001_S2_L002_R2_001.fastq.gz > gutR.fasta
    seqtk seq -a MTMO002_S4_L002_R1_001.fastq.gz > ovaryF.fasta
    seqtk seq -a MTMO002_S4_L002_R2_001.fastq.gz > ovaryR.fasta
    seqtk seq -a MTML003_S3_L002_R1_001.fastq.gz > larvaF.fasta
    seqtk seq -a MTML003_S3_L002_R2_001.fastq.gz > larvaR.fasta

* Create local Nardonella blast database to serach for sequences of interest.


    makeblastdb -in nardonella_genome.fasta -dbtype nucl -parse_seqids

* Search for genome matches on target fasta file using blastn.


    blastn -db nardonella_genome.fasta -query gutF.fasta -outfmt 6 -out blastn_gut_forward.tsv
    blastn -db nardonella_genome.fasta -query gutR.fasta -outfmt 6 -out blastn_gut_reverse.tsv
    blastn -db nardonella_genome.fasta -query ovaryF.fasta -outfmt 6 -out blastn_ovary_forward.tsv
    blastn -db nardonella_genome.fasta -query ovaryR.fasta -outfmt 6 -out blastn_ovary_reverse.tsv
    blastn -db nardonella_genome.fasta -query larvaF.fasta -outfmt 6 -out blastn_larva_forward.tsv
    blastn -db nardonella_genome.fasta -query larvaR.fasta -outfmt 6 -out blastn_larva_reverse.tsv


* Filter out blast identified Nardonella sequences from original fastq.gz files (to reduce the size of the files and speed up the process).



    seqtk subseq MTMG001_S2_L002_R1_001.fastq.gz gutF_nardo_list.csv > MTMG001_R1_nardo_only.fastq.gz
    seqtk subseq MTMG001_S2_L002_R2_001.fastq.gz gutR_nardo_list.csv > MTMG001_R2_nardo_only.fastq.gz
    seqtk subseq MTML003_S3_L002_R1_001.fastq.gz larvaF_nardo_list.csv > MTML003_R1_nardo_only.fastq.gz
    seqtk subseq MTML003_S3_L002_R2_001.fastq.gz larvaR_nardo_list.csv > MTML003_R2_nardo_only.fastq.gz
    seqtk subseq MTMO002_S4_L002_R1_001.fastq.gz ovaryF_nardo_list.csv > MTMO002_R1_nardo_only.fastq.gz
    seqtk subseq MTMO002_S4_L002_R2_001.fastq.gz ovaryR_nardo_list.csv > MTMO002_R2_nardo_only.fastq.gz


Blastn output summary:

 |Gut|Larva|Ovary|
 |:---:|:---:|:---:|
 |forward seqs = 9943|forward seqs = 23974|forward seqs = 126718|
 |reverse seqs = 9498|reverse seqs = 23711|reverse seqs = 124247|
 |total = **19441**|total = **47685**|total = **250965**|
|% of TOTAL reads = **0.018**|% of TOTAL reads = **0.036**|% of TOTAL reads = **0.18**|


**Based on the results above, and literature information, the ovaries are the selected tissue to proceed with the pilot project, they presented 10 times more Nardonella reads than the guts, and 5 times more reads than the larvae.**


---
### Estimation of future sequencing using the same strategy as Chouaia *et al.* ###

In the paper they save this and see if picture appears!!!


![image info](./genome_coverage.png)

Taking the ovary estimate as best case scenario, 140.0 M total sequences yield 250965 Nardonella sequences which represents roughly 5% o reference genome coverage. For 100% coverage we will need:

    Total seqs needed for 100% coverage = 100%
    250965  Nardonella sequences on ovary sample = 5%
    Total seqs needed for 100% coverage = 250965 x 100 / 5 = 5019300
    Total seqs needed for 100% coverage = 5019300 ~ 5 M seqs

To reach ~ 5 M Nardonella sequences with the current sequencing strategy we will need:

    Total ovary seqs needed = 100%
    5 M Nardo seqs = 0.18%
    Total ovary seqs needed = 2788500000 ~ 2.8 B pair-end seqs
    Total ovary seqs needed = 1.4 B reads

If this estimate is correct, at the price of US$ 6.35 per million reads we will have:

    1.4 B reads x $ 6.35 = $ 8890.00
