====
Data
====
• *Terje Stanium's* plasmids which are called:

1. 1-F5-96_S1_L001_R1_001.fastq.gz  > Library **F5**, *Forward*
2. 1-F5-96_S1_L001_R2_001.fastq.gz  > Library **F5**, *Reverse*
3. 4-F20-96_S2_L001_R1_001.fastq.gz > Library **F20**, *Forward*
4. 4-F20-96_S2_L001_R2_001.fastq.gz > Library **F20**, *Reverse*

Used in this analysis.


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


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Aligning reads with best candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
• We choose two plasmids considering the breadth of coverage results. Accession number of these plasmids are NC_025138.1 and
NC_025175.1, respectively.

------------------
F5-Best candidates
------------------

• For NC_025138.1

.. code:: bash

    bwa index NC_025138.1.fasta
    bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' NC_025138.1.fasta 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-NC_025138.1.sam
    samtools fixmate -O bam F5-NC_025138.1.sam F5-NC_025138.1fixmate.bam
    samtools sort -O bam -o F5-NC_025138.1fixmatesorted.bam F5-NC_025138.1fixmate.bam

• We can generalize the code like this;

.. code:: bash

    bwa index {plasmid accession}.fasta
    bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' {plasmid accession}.fasta 1-F5-96_S1_L001_R1_001.fastq.gz 1-F5-96_S1_L001_R2_001.fastq.gz > F5-{plasmid accession}.sam
    samtools fixmate -O bam F5-{plasmid accession}.sam F5-{plasmid accession}fixmate.bam
    samtools sort -O bam -o F5-{plasmid accession}fixmatesorted.bam F5-{plasmid accession}fixmate.bam


-------------------
F20-Best candidates
-------------------

• For NC_025138.1

.. code:: bash

    bwa index NC_025138.1.fasta
    bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' NC_025138.1.fasta 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-NC_025138.1.sam
    samtools fixmate -O bam F20-NC_025138.1.sam F20-NC_025138.1fixmate.bam
    samtools sort -O bam -o F20-NC_025138.1fixmatesorted.bam F20-NC_025138.1fixmate.bam

• We can generalize the code like this;

.. code:: bash

    bwa index {plasmid accession}.fasta
    bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' {plasmid accession}.fasta 4-F20-96_S2_L001_R1_001.fastq.gz 4-F20-96_S2_L001_R2_001.fastq.gz > F20-{plasmid accession}.sam
    samtools fixmate -O bam F20-{plasmid accession}.sam F20-{plasmid accession}fixmate.bam
    samtools sort -O bam -o F20-{plasmid accession}fixmatesorted.bam F20-{plasmid accession}fixmate.bam


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Creating new fastq files contain only mapped reads with best candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

• Getting read id's for creating new fastq files for NC_025138.1

.. code:: bash
	
	samtools view F5-NC_025138.1fixmatesorted.bam |cut -f 1|uniq |head -300000|sort -u > F5-NC_025138.1fastqid.txt
	samtools view F20-NC_025138.1fixmatesorted.bam |cut -f 1|uniq |head -205000|sort -u > F20-NC_025138.1fastqid.txt


• Getting read id's for creating new fastq files for NC_025175.1

.. code:: bash

	samtools view F5-NC_025175.1fixmatesorted.bam |cut -f 1|uniq |head -200000|sort -u > F5-NC_025175.1fastqid.txt
	samtools view F20-NC_025175.1fixmatesorted.bam |cut -f 1|uniq |head -200000|sort -u > F20-NC_025175.1fastqid.txt

---------------------------------------------
Creating new fastq files with given read id's
---------------------------------------------

• This python script creates new fastq files with given forward read id's.

.. code-block:: python
   :linenos:

   from Bio import SeqIO

   input_file = "{Forward fastq file}"
   id_file = "{Library}-{plasmid accession}fastqid.txt"
   output_file = "{Library}-{plasmid accession}.fastq"
   wanted = set(line.rstrip("\n").split(None, 1)[0] for line in open(id_file))
   print("Found %i unique identifiers in %s" % (len(wanted), id_file))
   records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
   count = SeqIO.write(records, output_file, "fastq")
   print("Saved %i records from %s to %s" % (count, input_file, output_file))
   if count < len(wanted):
       print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))


• This python script creates new fastq files with given reverse read id's.

.. code-block:: python
   :linenos:
   
   from Bio import SeqIO

   input_file = "{Reverse fastq file}"
   id_file = "{Library}-{plasmid accession}fastqid.txt"
   output_file = "{Library}-{plasmid accession}.fastq"
   wanted = set(line.rstrip("\n").split(None, 1)[0] for line in open(id_file))
   print("Found %i unique identifiers in %s" % (len(wanted), id_file))
   records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
   count = SeqIO.write(records, output_file, "fastq")
   print("Saved %i records from %s to %s" % (count, input_file, output_file))
   if count < len(wanted):
       print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))


^^^^^^^^^^^^^^^^^^^^^^
Assembly with Geneious
^^^^^^^^^^^^^^^^^^^^^^

• After creating new fastq files, we should do assembly. You can do assembly easily with following geneious de novo assembly tutorial. You can access the tutorial here `De novo Assembly Tutorial`_.

.. _De novo Assembly Tutorial: https://www.geneious.com/tutorials/de-novo-assembly/


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Filtering reads with given contigs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

• After the assembly process we have ace files for each contig. Some regions in the contigs have lower depth for reference bases. We do not want those reads and we should delete them. The following python script will do the work.

.. code-block:: python
   :linenos:

   import sys
   from Bio.Sequencing import Ace
   from Bio import SeqIO
   import numpy as np

   if len(sys.argv)<4:
       print(sys.argv[0],"not enough arguments")
       print("Usage:", sys.argv[0],"ACE input.fastq output.fastq")
       exit(1)

   cmd_name = sys.argv.pop(0)
   input_file = sys.argv.pop(0)
   output_file = sys.argv.pop(0)

   def parse_af(ace_file):
       ans = {}
   for line in open(ace_file):
       if line.startswith("AF"):
           _, read_id, _, pos = line.strip().split()
           ans[read_id]=int(pos)
   return(ans)

   bad_reads = set()

   for ace_file in sys.argv:
       assembly = Ace.read(open(ace_file))
       contig = assembly.contigs[0]
  	   print("%s: %d reads" % (ace_file, contig.nreads))
   if len(contig.af)==0:
       af = parse_af(ace_file)
       all_reads = [(contig.reads[i].rd.name,
        			af[contig.reads[i].rd.name],
        			contig.reads[i].rd.padded_bases) for i in range(contig.nreads)]
   else:
       all_reads = [(contig.reads[i].rd.name,
        			contig.af[i].padded_start,
        			contig.reads[i].rd.padded_bases) for i in range(contig.nreads)]

   depth = np.zeros(contig.nbases+1, dtype=int)

   for name, start, length in all_reads:
       for j in range(start, start + length):
           depth[j] +=1

   bad_places = (depth < (depth.mean()-3*depth.std())) | (depth>= (depth.mean()+3*depth.std()))

   for name, start, length in all_reads:
       bad_bp_in_read = np.sum(bad_places[start:(start+length)])
       if bad_bp_in_read > length/5: # if over 20% of bp are "bad"...
           bad_reads.add(name[:name.index("_")]) # then remember the fragment name

   print("Now filtering %d bad fragments" % (len(bad_reads)))
   records = [r for r in SeqIO.parse(input_file, "fastq") if r.id not in bad_reads]
   count = SeqIO.write(records, output_file, "fastq")
   print("Saved %i records from %s to %s" % (count, input_file, output_file))

^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Assembly with filtered reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

• For NC_025138.1

.. code:: bash

	## Trimming
	
	sickle pe -f R1-filtered.fastq -r R2-filtered.fastq -t sanger -o R1-filtered-trimmed.fastq -p R2-filtered-trimmed.fastq -s singles.fastq -q 30 -l 45
	
	## Assembly
	
	### Spades
	
	spades.py --careful -o SPAdes_out -1 R1-filtered-trimmed.fastq -2 R2-filtered-trimmed.fastq -s singles.fastq
	
	### Abyss 
	
	abyss-pe  k=31 name=Kmer31 in='/home/fustunel/Plasmid/geneious-2019-02-21/Sickle-2019-03-19/F20-NC_025138.1/R1-filtered-trimmed.fastq /home/fustunel/Plasmid/geneious-2019-02-21/Sickle-2019-03-19/F20-NC_025138.1/R2-filtered-trimmed.fastq' se='/home/fustunel/Plasmid/geneious-2019-02-21/Sickle-2019-03-19/F20-NC_025138.1/singles.fastq'
	
	## Statistics
	
	quast.py SPAdes_out/scaffolds.fasta -o quast-SPAdes
	quast.py /home/fustunel/Plasmid/geneious-2019-02-21/Abyss-2019-03-19/F20-NC_025138.1/Kmer31-scaffolds.fa -o quast-Abyss


• We can generalize the code like this;

.. code:: bash

	## Trimming
	
	sickle pe -f R1-filtered.fastq -r R2-filtered.fastq -t sanger -o R1-filtered-trimmed.fastq -p R2-filtered-trimmed.fastq -s singles.fastq -q 30 -l 45
	
	## Assembly
	
	### Spades
	
	spades.py --careful -o SPAdes_out -1 $plasmid folder path/R1-filtered-trimmed.fastq -2 $plasmid folder path/R2-filtered-trimmed.fastq -s singles.fastq
	
	### Abyss 
	
	abyss-pe  k=31 name=Kmer31 in='$plasmid folder path/R1-filtered-trimmed.fastq $plasmid folder path/R2-filtered-trimmed.fastq' se='$plasmid folder path/singles.fastq'
	
	## Statistics
	
	quast.py SPAdes_out/scaffolds.fasta -o quast-SPAdes
	quast.py Kmer31-scaffolds.fa -o quast-Abyss


