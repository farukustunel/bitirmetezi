====
Data
====
*Terje Stanium's* plasmids which are called:

1. 1-F5-96_S1_L001_R1_001.fastq.gz  > Library **F5**, *Forward*
2. 1-F5-96_S1_L001_R2_001.fastq.gz  > Library **F5**, *Reverse*
3. 4-F20-96_S2_L001_R1_001.fastq.gz > Library **F20**, *Forward*
4. 4-F20-96_S2_L001_R2_001.fastq.gz > Library **F20**, *Reverse*

 Used in this analysis.

Reference plasmids can be obtain from 'Refseq Database <ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/>'_.

 The name of the files are:

1. plasmid.1.1.genomic.fna.gz
2. plasmid.2.1.genomic.fna.gz
3. plasmid.3.1.genomic.fna.gz
4. plasmid.4.1.genomic.fna.gz
5. plasmid.5.1.genomic.fna.gz


=========
Protocols
=========

We compare the each libraries with references by using Burrows-Wheeler Aligner command line tool.This tool contains different
aligners with different algorithms.Two of them are Bwa-mem and Bwa-aln.

-------
Bwa-mem
-------

^^^^^^^^^^^
F5-Plasmid1
^^^^^^^^^^^

**Indexing reference, aligning Library F5 and plasmid.1.1.genomic.fna.gz**

.. code:: bash

   bwa index -p plasmid1 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz
   bwa mem plasmid1 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid1.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F5-plasmid1.bam F5-plasmid1.sam
   samtools index F5-plasmid1.bam
   samtools idxstats F5-plasmid1.bam |awk '$3>1000'|sort -k3n |grep > F5-plasmid1.stats

^^^^^^^^^^^^
F20-Plasmid1
^^^^^^^^^^^^

**Aligning Library F20 and plasmid.1.1.genomic.fna.gz**

.. code:: bash

   bwa mem plasmid1 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid1.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F20-plasmid1.bam F20-plasmid1.sam
   samtools index F20-plasmid1.bam
   samtools idxstats F20-plasmid1.bam |awk '$3>1000'|sort -k3n > F20-plasmid1.stats 

^^^^^^^^^^^
F5-Plasmid2
^^^^^^^^^^^

**Indexing reference, aligning Library F5 and plasmid.2.1.genomic.fna.gz**

.. code:: bash
   
   bwa index -p plasmid2 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.2.1.genomic.fna.gz
   bwa mem plasmid2 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid2.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F5-plasmid2.bam F5-plasmid2.sam
   samtools index F5-plasmid2.bam
   samtools idxstats F5-plasmid2.bam |awk '$3>1000'|sort -k3n > F5-plasmid2.stats

^^^^^^^^^^^^
F20-Plasmid2
^^^^^^^^^^^^

**Aligning Library F20 and plasmid.2.1.genomic.fna.gz**

.. code:: bash

   bwa mem plasmid2 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid2.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F20-plasmid2.bam F20-plasmid2.sam
   samtools index F20-plasmid2.bam
   samtools idxstats F20-plasmid2.bam |awk '$3>1000'|sort -k3n > F20-plasmid2.stats

^^^^^^^^^^^
F5-Plasmid3
^^^^^^^^^^^

**Indexing reference, aligning Library F5 and plasmid.3.1.genomic.fna.gz**

.. code:: bash

   bwa index -p plasmid3 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.3.1.genomic.fna.gz
   bwa mem plasmid3 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid3.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F5-plasmid3.bam F5-plasmid3.sam
   samtools index F5-plasmid3.bam
   samtools idxstats F5-plasmid3.bam |awk '$3>1000'|sort -k3n > F5-plasmid3.stats

^^^^^^^^^^^^
F20-Plasmid3
^^^^^^^^^^^^

**Aligning Library F20 and plasmid.3.1.genomic.fna.gz**

.. code:: bash

   bwa mem plasmid3 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid3.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F20-plasmid3.bam F20-plasmid3.sam
   samtools index F20-plasmid3.bam
   samtools idxstats F20-plasmid3.bam |awk '$3>1000'|sort -k3n > F20-plasmid3.stats

^^^^^^^^^^^
F5-Plasmid4
^^^^^^^^^^^

**Indexing reference, aligning Library F5 and plasmid.4.1.genomic.fna.gz**

.. code:: bash

   bwa index -p plasmid4 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.4.1.genomic.fna.gz
   bwa mem plasmid4 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid4.sam


**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F5-plasmid4.bam F5-plasmid4.sam
   samtools index F5-plasmid4.bam
   samtools idxstats F5-plasmid4.bam |awk '$3>1000'|sort -k3n > F5-plasmid4.stats

^^^^^^^^^^^^
F20-Plasmid4
^^^^^^^^^^^^

**Aligning Library F20 and plasmid.4.1.genomic.fna.gz**

.. code:: bash
   
   bwa mem plasmid4 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid4.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F20-plasmid4.bam F20-plasmid4.sam
   samtools index F20-plasmid4.bam
   samtools idxstats F20-plasmid4.bam |awk '$3>1000'|sort -k3n > F20-plasmid4.stats

^^^^^^^^^^^
F5-Plasmid5
^^^^^^^^^^^

**Indexing reference, aligning Library F5 and plasmid.5.1.genomic.fna.gz**

.. code:: bash

   bwa index -p plasmid5 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.5.1.genomic.fna.gz
   bwa mem plasmid5 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid5.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F5-plasmid5.bam F5-plasmid5.sam
   samtools index F5-plasmid5.bam
   samtools idxstats F5-plasmid5.bam |awk '$3>1000'|sort -k3n > F5-plasmid5.stats

^^^^^^^^^^^^
F20-Plasmid5
^^^^^^^^^^^^

**Aligning Library F20 and plasmid.5.1.genomic.fna.gz**

.. code:: bash
   
   bwa mem plasmid5 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid5.sam

**Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F20-plasmid5.bam F20-plasmid5.sam
   samtools index F20-plasmid5.bam
   samtools idxstats F20-plasmid5.bam |awk '$3>1000'|sort -k3n > F20-plasmid5.stats

-------
Bwa-aln
-------

**Indexing reference plasmids**

.. code:: bash

   bwa index -p plasmid1 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz
   bwa index -p plasmid2 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.2.1.genomic.fna.gz
   bwa index -p plasmid3 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.3.1.genomic.fna.gz
   bwa index -p plasmid4 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.4.1.genomic.fna.gz
   bwa index -p plasmid5 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.5.1.genomic.fna.gz


^^^^^^^^^^^
F5-Plasmid1
^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid1 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid1.sai
   bwa aln plasmid1 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid1.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid1 F5-R1-plasmid1.sai F5-R2-plasmid1.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz >F5-plasmid1.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid1.sam > F5-plasmid1.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid1sorted.bam  F5-plasmid1.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid1sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid1sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid1.stats


^^^^^^^^^^^^
F20-Plasmid1
^^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid1 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid1.sai
   bwa aln plasmid1 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid1.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid1 F20-R1-plasmid1.sai F20-R2-plasmid1.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid1.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid1.sam > F20-plasmid1.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid1sorted.bam  F20-plasmid1.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid1sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid1sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid1.stats


^^^^^^^^^^^
F5-Plasmid2
^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid2 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid2.sai
   bwa aln plasmid2 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid2.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid2 F5-R1-plasmid2.sai F5-R2-plasmid2.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-plasmid2.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid2.sam > F5-plasmid2.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid2sorted.bam  F5-plasmid2.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid2sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid2sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid2.stats


^^^^^^^^^^^^
F20-Plasmid2
^^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid2 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid2.sai
   bwa aln plasmid2 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid2.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid2 F20-R1-plasmid2.sai F20-R2-plasmid2.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid2.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid2.sam > F20-plasmid2.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid2sorted.bam  F20-plasmid2.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid2sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid2sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid2.stats


^^^^^^^^^^^
F5-Plasmid3
^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid3 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid3.sai
   bwa aln plasmid3 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid3.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid3 F5-R1-plasmid3.sai F5-R2-plasmid3.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-plasmid3.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid3.sam > F5-plasmid3.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid3sorted.bam  F5-plasmid3.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid3sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid3sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid3.stats


^^^^^^^^^^^^
F20-Plasmid3
^^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid3 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid3.sai
   bwa aln plasmid3 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid3.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid3 F20-R1-plasmid3.sai F20-R2-plasmid3.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid3.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid3.sam > F20-plasmid3.bam

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid3sorted.bam  F20-plasmid3.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid3sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid3sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid3.stats


^^^^^^^^^^^
F5-Plasmid4
^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid4 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid4.sai
   bwa aln plasmid4 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid4.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid4 F5-R1-plasmid4.sai F5-R2-plasmid4.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-plasmid4.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid4.sam > F5-plasmid4.bam  

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid4sorted.bam  F5-plasmid4.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid4sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid4sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid4.stats


^^^^^^^^^^^^
F20-Plasmid4
^^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid4 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid4.sai
   bwa aln plasmid4 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid4.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid4 F20-R1-plasmid4.sai F20-R2-plasmid4.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid4.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid4.sam > F20-plasmid4.bam 

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid4sorted.bam  F20-plasmid4.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid4sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid4sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid4.stats


^^^^^^^^^^^
F5-Plasmid5
^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid5 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid5.sai
   bwa aln plasmid5 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid5.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid5 F5-R1-plasmid5.sai F5-R2-plasmid5.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-plasmid5.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid5.sam > F5-plasmid5.bam  

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid5sorted.bam  F5-plasmid5.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid5sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid5sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid5.stats


^^^^^^^^^^^^
F20-Plasmid5
^^^^^^^^^^^^

**Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid5 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid5.sai
   bwa aln plasmid5 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid5.sai

**Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid5 F20-R1-plasmid5.sai F20-R2-plasmid5.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid5.sam

**Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid5.sam > F20-plasmid5.bam

**Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid5sorted.bam  F20-plasmid5.bam

**Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid5sorted.bam

**Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid5sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid5.stats


----------------
Finding Coverage
----------------

**Getting reference id and mapped reads for calculating genome coverage**

.. code:: bash

   samtools idxstats F5-plasmid1.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid1.txt
   samtools idxstats F5-plasmid2.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid2.txt
   samtools idxstats F5-plasmid3.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid3.txt
   samtools idxstats F5-plasmid4.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid4.txt
   samtools idxstats F5-plasmid5.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid5.txt
   samtools idxstats F20-plasmid1.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid1.txt
   samtools idxstats F20-plasmid2.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid2.txt
   samtools idxstats F20-plasmid3.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid3.txt
   samtools idxstats F20-plasmid4.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid4.txt
   samtools idxstats F20-plasmid5.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid5.txt


**Calculating Depth and Coverage from a BAM file**

.. code:: bash
   
   genomeCoverageBed -ibam F5-plasmid1.bam -g F5-plasmid1.txt > F5-plasmid1coverage.txt
   genomeCoverageBed -ibam F5-plasmid2.bam -g F5-plasmid2.txt > F5-plasmid2coverage.txt
   genomeCoverageBed -ibam F5-plasmid3.bam -g F5-plasmid3.txt > F5-plasmid3coverage.txt
   genomeCoverageBed -ibam F5-plasmid4.bam -g F5-plasmid4.txt > F5-plasmid4coverage.txt
   genomeCoverageBed -ibam F20-plasmid1.bam -g F20-plasmid1.txt > F20-plasmid1coverage.txt
   genomeCoverageBed -ibam F20-plasmid2.bam -g F20-plasmid2.txt > F20-plasmid2coverage.txt
   genomeCoverageBed -ibam F20-plasmid3.bam -g F20-plasmid3.txt > F20-plasmid3coverage.txt
   genomeCoverageBed -ibam F20-plasmid4.bam -g F20-plasmid4.txt > F20-plasmid4coverage.txt
   genomeCoverageBed -ibam F20-plasmid5.bam -g F20-plasmid5.txt > F20-plasmid5coverage.txt


**Calculating Breadth of Coverage and choosing best candidates for assembly**

.. code:: bash
   
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid1coverage.txt > F5-plasmid1coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid2coverage.txt > F5-plasmid2coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid3coverage.txt > F5-plasmid3coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid4coverage.txt > F5-plasmid4coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid5coverage.txt > F5-plasmid5coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid1coverage.txt > F20-plasmid1coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid2coverage.txt > F20-plasmid2coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid3coverage.txt > F20-plasmid3coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid4coverage.txt > F20-plasmid4coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid5coverage.txt > F20-plasmid5coveragesummary.txt

**Combining Library F5 results**

.. code:: bash
   
   cat F5-plasmid1coveragesummary.txt F5-plasmid2coveragesummary.txt F5-plasmid3coveragesummary.txt F5-plasmid4coveragesummary.txt F5-plasmid5coveragesummary.txt > F5-plasmidcoverageallsummary.txt
   
   sort -k2nr F5-plasmidcoverageallsummary.txt > F5-plasmidcoverageallsummarysorted.txt
  
**Combining Library F20 results**

.. code:: bash
   
   cat F20-plasmid1coveragesummary.txt F20-plasmid2coveragesummary.txt F20-plasmid3coveragesummary.txt F20-plasmid4coveragesummary.txt F20-plasmid5coveragesummary.txt > F20-plasmidcoverageallsummary.txt
   
   sort -k2nr F20-plasmidcoverageallsummary.txt > F20-plasmidcoverageallsummarysorted.txt

