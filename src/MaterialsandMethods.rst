====
Data
====
• *Terje Stanium's* plasmids which are called:

1. 1-F5-96_S1_L001_R1_001.fastq.gz  > Library **F5**, *Forward*
2. 1-F5-96_S1_L001_R2_001.fastq.gz  > Library **F5**, *Reverse*
3. 4-F20-96_S2_L001_R1_001.fastq.gz > Library **F20**, *Forward*
4. 4-F20-96_S2_L001_R2_001.fastq.gz > Library **F20**, *Reverse*

• Used in this analysis.

• Reference plasmids can be obtain from `Refseq Database`_.

.. _Refseq Database: https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid


• The name of the files are:

1. plasmid.1.1.genomic.fna.gz
2. plasmid.2.1.genomic.fna.gz
3. plasmid.3.1.genomic.fna.gz
4. plasmid.4.1.genomic.fna.gz
5. plasmid.5.1.genomic.fna.gz


=========
Protocols
=========

• We compare the each libraries with references by using Burrows-Wheeler Aligner command line tool.This tool contains different aligners with different algorithms.Two of them are Bwa-mem and Bwa-aln.

-------
Bwa-mem
-------

^^^^^^^^^^^
F5-Plasmid1
^^^^^^^^^^^

• **Indexing reference, aligning Library F5 and plasmid.1.1.genomic.fna.gz**

.. code:: bash

   bwa index -p plasmid1 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz
   bwa mem plasmid1 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz -o F5-plasmid1.sam

• **Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F5-plasmid1.bam F5-plasmid1.sam
   samtools index F5-plasmid1.bam
   samtools idxstats F5-plasmid1.bam |awk '$3>1000'|sort -k3n |grep > F5-plasmid1.stats


^^^^^^^^^^^^^^^^^
F5-Other Plasmids
^^^^^^^^^^^^^^^^^

• **Indexing reference, aligning Library F5 and other plasmids.**

.. code:: bash
   
   bwa index -p plasmid{number} {file_path}
   bwa mem plasmid{number} {Forward fastq file} {Reverse fastq file} -o F5-plasmid{number}.sam

• **Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F5-plasmid{number}.bam F5-plasmid{number}.sam
   samtools index F5-plasmid{number}.bam
   samtools idxstats F5-plasmid{number}.bam |awk '$3>1000'|sort -k3n > F5-plasmid{number}.stats


^^^^^^^^^^^^
F20-Plasmid1
^^^^^^^^^^^^

• **Aligning Library F20 and plasmid.1.1.genomic.fna.gz**

.. code:: bash

   bwa mem plasmid1 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz -o F20-plasmid1.sam

• **Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools sort -O BAM -o F20-plasmid1.bam F20-plasmid1.sam
   samtools index F20-plasmid1.bam
   samtools idxstats F20-plasmid1.bam |awk '$3>1000'|sort -k3n > F20-plasmid1.stats 


^^^^^^^^^^^^^^^^^^
F20-Other Plasmids
^^^^^^^^^^^^^^^^^^

• **Indexing reference, aligning Library F20 and other plasmids.**

.. code:: bash
   
   bwa index -p plasmid{number} {file_path}
   bwa mem plasmid{number} {Forward fastq file} {Reverse fastq file} -o F20-plasmid{number}.sam

• **Converting SAM file to BAM file, sorting BAM file and indexing,getting read ids with mapping more than 1000 reads**

.. code:: bash

   samtools sort -O BAM -o F20-plasmid{number}.bam F20-plasmid{number}.sam
   samtools index F20-plasmid{number}.bam
   samtools idxstats F20-plasmid{number}.bam |awk '$3>1000'|sort -k3n > F20-plasmid{number}.stats


-------
Bwa-aln
-------

• **Indexing reference plasmids**

.. code:: bash

   bwa index -p plasmid1 /home/db/Mirror/ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.1.genomic.fna.gz
   bwa index -p plasmid{number} {file_path}

^^^^^^^^^^^
F5-Plasmid1
^^^^^^^^^^^

• **Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid1 1-F5-96_S1_L001_R1_001.fastq.gz > F5-R1-plasmid1.sai
   bwa aln plasmid1 1-F5-96_S1_L001_R2_001.fastq.gz > F5-R2-plasmid1.sai

• **Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid1 F5-R1-plasmid1.sai F5-R2-plasmid1.sai 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz >F5-plasmid1.sam

• **Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid1.sam > F5-plasmid1.bam 

• **Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid1sorted.bam  F5-plasmid1.bam

• **Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid1sorted.bam

• **Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid1sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid1.stats


^^^^^^^^^^^^^^^^^
F5-Other Plasmids
^^^^^^^^^^^^^^^^^

• **Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid{number} {Forward fastq file} > F5-R1-plasmid{number}.sai
   bwa aln plasmid{number} {Reverse fastq file} > F5-R2-plasmid{number}.sai

• **Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid{number} F5-R1-plasmid{number}.sai F5-R2-plasmid{number}.sai {Forward fastq file} {Reverse fastq file} > F5-plasmid{number}.sam

• **Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F5-plasmid{number}.sam > F5-plasmid{number}.bam 

• **Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F5-plasmid{number}sorted.bam  F5-plasmid{number}.bam

• **Indexing sorted BAM file**

.. code:: bash
   
   samtools index F5-plasmid{number}sorted.bam

• **Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F5-plasmid{number}sorted.bam |awk '$3>1000'|sort -k3n > F5-plasmid{number}.stats


^^^^^^^^^^^^
F20-Plasmid1
^^^^^^^^^^^^

• **Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid1 4-F20-96_S2_L001_R1_001.fastq.gz > F20-R1-plasmid1.sai
   bwa aln plasmid1 4-F20-96_S2_L001_R2_001.fastq.gz > F20-R2-plasmid1.sai

• **Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid1 F20-R1-plasmid1.sai F20-R2-plasmid1.sai 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-plasmid1.sam

• **Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid1.sam > F20-plasmid1.bam 

• **Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid1sorted.bam  F20-plasmid1.bam

• **Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid1sorted.bam

• **Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid1sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid1.stats


^^^^^^^^^^^^^^^^^^
F20-Other Plasmids
^^^^^^^^^^^^^^^^^^

• **Aligning forward and reverse reads with references**

.. code:: bash
   
   bwa aln plasmid{number} {Forward fastq file} > F20-R1-plasmid{number}.sai
   bwa aln plasmid{number} {Reverse fastq file} > F20-R2-plasmid{number}.sai

• **Combining outputs in one SAM file**

.. code:: bash
   
   bwa sampe plasmid{number} F20-R1-plasmid{number}.sai F20-R2-plasmid{number}.sai {Forward fastq file} {Reverse fastq file} > F20-plasmid{number}.sam

• **Converting SAM file to BAM file**

.. code:: bash
   
   samtools view -bS F20-plasmid{number}.sam > F20-plasmid{number}.bam 

• **Sorting BAM file**

.. code:: bash
   
   samtools sort -O bam -o F20-plasmid{number}sorted.bam  F20-plasmid{number}.bam

• **Indexing sorted BAM file**

.. code:: bash
   
   samtools index F20-plasmid{number}sorted.bam

• **Getting read ids with mapping more than 1000 reads**

.. code:: bash
   
   samtools idxstats F20-plasmid{number}sorted.bam |awk '$3>1000'|sort -k3n > F20-plasmid{number}.stats


----------------
Finding Coverage
----------------

• **Getting reference id and mapped reads for calculating genome coverage**

.. code:: bash

   samtools idxstats F5-plasmid1.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid1.txt
   samtools idxstats F5-plasmid{number}.bam |awk -v OFS='\t' '{print $1, $2}' > F5-plasmid{number}.txt
   samtools idxstats F20-plasmid1.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid1.txt
   samtools idxstats F20-plasmid{number}.bam |awk -v OFS='\t' '{print $1, $2}' > F20-plasmid{number}.txt



**Calculating Depth and Coverage from a BAM file**

• The first column of output file is reference id. The following columns are respectively called depth number, number of bases
with given depth, reference length and percentage of coverage with given depth.

.. code:: bash
   
   genomeCoverageBed -ibam F5-plasmid1.bam -g F5-plasmid1.txt > F5-plasmid1coverage.txt
   genomeCoverageBed -ibam F5-plasmid{number}.bam -g F5-plasmid{number}.txt > F5-plasmid{number}coverage.txt
   genomeCoverageBed -ibam F20-plasmid1.bam -g F20-plasmid1.txt > F20-plasmid1coverage.txt
   genomeCoverageBed -ibam F20-plasmid{number}.bam -g F20-plasmid{number}.txt > F20-plasmid{number}coverage.txt



**Calculating Breadth of Coverage and choosing best candidates for assembly**

• If second column equals number zero, it means that 0 depth or unmapped regions occur for reference. We can easily calcutate
the breadth of coverage with 1-$5 condition which means that summary of percentage of all depth but 0. Then we can determine
the best candidates for each plasmids.  

.. code:: bash
   
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid1coverage.txt > F5-plasmid1coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F5-plasmid{number}coverage.txt > F5-plasmid{number}coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid1coverage.txt > F20-plasmid1coveragesummary.txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' F20-plasmid{number}coverage.txt > F20-plasmid{number}coveragesummary.txt


• **Combining Library F5 results**

.. code:: bash
   
   cat F5-plasmid1coveragesummary.txt F5-plasmid2coveragesummary.txt F5-plasmid3coveragesummary.txt F5-plasmid4coveragesummary.txt F5-plasmid5coveragesummary.txt > F5-plasmidcoverageallsummary.txt
   
   sort -k2nr F5-plasmidcoverageallsummary.txt > F5-plasmidcoverageallsummarysorted.txt
  
• **Combining Library F20 results**

.. code:: bash
   
   cat F20-plasmid1coveragesummary.txt F20-plasmid2coveragesummary.txt F20-plasmid3coveragesummary.txt F20-plasmid4coveragesummary.txt F20-plasmid5coveragesummary.txt > F20-plasmidcoverageallsummary.txt
   
   sort -k2nr F20-plasmidcoverageallsummary.txt > F20-plasmidcoverageallsummarysorted.txt

