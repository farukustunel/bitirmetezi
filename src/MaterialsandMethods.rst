====
Data
====
Samples were taken from plasmids isolated in *Escherichia coli* and *Citrobacter freundii*. These samples were sequenced into two libraries,  called F5 and F20, respectively. Each library contains forward and reverse reads, in two separated files. The file names of these libraries are: 

1. 1-F5-96_S1_L001_R1_001.fastq.gz  > Library **F5**, *Forward Reads*
2. 1-F5-96_S1_L001_R2_001.fastq.gz  > Library **F5**, *Reverse Reads*
3. 4-F20-96_S2_L001_R1_001.fastq.gz > Library **F20**, *Forward Reads*
4. 4-F20-96_S2_L001_R2_001.fastq.gz > Library **F20**, *Reverse Reads*

We aim to find the if there is a plasmid that matches with our reads in the NCBI *RefSeq database*. We dowloaded all the 15076 sequences that were avilable at https://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid on YY/YYYY, packed in 5 files.

=========
Protocols
=========

The workflow for bacterial genome assembly used in this study is shown in the figure below.

.. figure:: ../_static/protocol.svg
   :width: 300px
   :align: center
   :height: 800px
   :figclass: align-center 

   Workflow

--------
Trimming
--------

The first crucial step for the genome assembly process is trimming adaptors from both 5' and 3' ends of your sequences. Also, getting rid of low-quality reads improves the standard of your analyzes. We used ``trimmomatic``  by using the following bash code.

.. code-block:: bash
   :linenos:

   java -jar trimmomatic-0.36.jar PE -phred33 [Forward fastq file] [Reverse fastq file] \
   [Forward paired output] [Forward unpaired output] [Reverse paired output] [Reverse unpaired output] \
   ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold> \ 
   LEADING:<quality> TRAILING:<quality> SLIDINGWINDOW:<windowSize>:<requiredQuality> MINLEN:<length>


As you can see from the above code, ``trimmomatic`` is a java based application. In the first line, we used ``java -jar trimmomatic-0.36.jar`` command to run the tool. The ``PE`` parameter indicates that input is paired-end data. The ``-phred33`` parameter indicates the quality encoding of your sequences. There are different quality encodings, such as ``-phred64``. 

As a read input, we used our raw fastq files. In the second line, you can specify the output names. In the following line, ``ILLUMINACLIP`` parameter takes four inputs. With ``<fastaWithAdaptersEtc>``, you indicate a fasta file containing all adaptors and PCR sequences, etc. We used ``TruSeq3-PE.fa`` file for this option. The ``<seed mismatches>`` option, allows us to choose the maximum mismatch count which will still allow a full match to be performed. We used ``2`` for this option. The ``<palindrome clip threshold>`` specifies how accurate the match between the two 'adapter ligated' reads. We used ``30`` for this option. With ``<simple clip threshold>``, you can specify how accurate the match between any adapter etc. We used ``10`` for this option. In the last line, the ``LEADING`` parameter allow us to specify the minimum quality required to keep a base. We used ``3`` for this option. The ``TRAILING`` parameter acts as the previous one. We also used ``3`` for this option. The ``SLIDINGWINDOW`` parameter takes two inputs. With ``<windowSize>`` option, you can specify the number of bases to average across. With ``<requiredQuality>`` option, you can specify the average quality required. We used ``4 and 15`` for this option. The ``MINLEN`` parameter allows us to specify the minimum length of reads to be kept. We used ``36`` for this option.

Trimmomatic produces four fastq files: two for paired forward and reverse reads, and two for unpaired forward and reverse reads.

---------------
Quality Control
---------------

After the trimming process, we want to compare the files we get with our raw data. We used ``FastQC`` quality control tool for this purpose. This program can be run in two ways. Either as a graphical user program that you can dynamically upload your files with ``Fastq`` format and you can see the results easily. Or, you can use command-line as a non-graphical program. In both ways, ``FastQC`` produces ``HTML`` files as an output that contains basic statistics, per base and sequence quality, per base and sequence GC content, sequence length distribution, sequence duplication levels, over-represented sequences, K-mer content, *et cetera*.

We used ``FastQC`` as a graphical user program on windows operating system. You can run ``FastQC`` with the following steps, easily.

1. Click ``run_fastqc.bat`` file. This will open GUI (Graphical User Interface) for ``FastQC``.
2. Click ``File -> Open...`` or just simply push ``Ctrl + O`` from your keyboard. Choose your file and the program starts automatically.
3. Click ``File -> Save report...`` or push ``Ctrl + S`` from your keyboard. This will save your results as an ``HTML`` output.

-------
Mapping
-------

The following step in our protocol is to map or align the reads into reference plasmids. Hence, we compared each library with plasmid references using Burrows-Wheeler Aligner alignment tool. This tool has different algorithms. Two of them are ``Bwa-mem`` and ``Bwa-aln``. We used both of them and analyzed the results.

^^^^^^^
Bwa-mem
^^^^^^^

The following *bash* code gives us the result of the alignment, which contains only the reference plasmids having more than 1000 mapped reads.

In the first line, the program creates an index for the alignment. The ``-p`` parameter provides to title of output database. In the second line, *bwa* carries out the alignment process using the ``mem`` algorithm and creates an output in **Sequence Alignment Map** format (a.k.a *SAM*). The ``-o`` parameter indicates the output file.

In the third line, *samtools* sorts the *sam* file and converts it to a *bam* file. The ``-O`` parameter indicates the type of the output file. In the following line, *samtools* indexes the *bam* file. In the last line, *idxstats* gives us statistics about the results. It is a *TAB-delimited* file and each line consists of *reference sequence name*, *sequence length*, *mapped reads number*, *unmapped reads number*. Therefore, we use the *awk* command to select the third column, which represents *mapped reads number*, and we filtered with a given threshold. After that, we sorted the output with the *sort* command where ``N`` in ``-kN`` is the number of the key, and ``n`` means numeric sort. 

.. code-block:: bash
   :linenos:

   bwa index -p [Database name] [Reference_file_path]
   bwa mem [Database name] [Forward fastq file] [Reverse fastq file] -o [Output file].sam
   samtools sort -O BAM -o [Output file].bam [Output file].sam
   samtools index [Output file].bam
   samtools idxstats [Output file].bam |awk '$3>1000'|sort -k3n > [Output file].stats

^^^^^^^
Bwa-aln
^^^^^^^

Using ``bwa-aln`` is more or less the same as ``bwa-mem``, except minor differences. As it is seen in the previous explanation, the differences occur in the second, third, fourth and fifth lines. In the second line of code, we took up the database and we carried out the alignment process with our reads. The output is a ``sai`` file which is produced by ``bwa-aln``. In the third line, we follow the same idea for reverse reads. In the following line, ``sampe`` sub-command allow us to create a sam file with combining forward reads and reverse reads alignment. In the fifth line, we used ``samtools`` with ``view`` sub-command, ``-bS`` parameter allows us to create bam file from the *sam* file where ``-b`` parameter for the type of output file and the ``-S`` parameter for ignoring compatibility for previous versions of ``samtools``.

.. code-block:: bash
   :linenos:

   bwa index -p [Database name] [Reference_file_path]
   bwa aln [Database name] [Forward fastq file] > [Forward Output].sai
   bwa aln [Database name] [Reverse fastq file] > [Reverse Output].sai
   bwa sampe [Database name] [Forward Output].sai [Reverse Output].sai [Forward fastq file] [Reverse fastq file] > [Output file].sam
   samtools view -bS [Output file].sam > [Output file].bam
   samtools sort -O bam -o [Sorted output file].bam  [Output file].bam
   samtools index [Sorted output file].bam
   samtools idxstats [Sorted output file].bam |awk '$3>1000'|sort -k3n > [Output file].stats

---------------------------
Finding Breadth of Coverage
---------------------------

To determine the most appropriate reference plasmid, we should know the breadth of coverage of each reference plasmid. For this purpose, we used ``bedtools`` with the ``genomeCoverageBed`` sub-command. The ``genomeCoverageBed`` computes a histogram of coverage within a given genome. In the first line of code, we took the ``samtools idxstats`` output and we select columns one and two. As we mentioned previously, ``idxstats`` output contains *reference sequence name*, *sequence length*, *mapped reads number*, *unmapped reads number*, respectively. We need only first and second column for the make ``genomeCoverageBed`` work.

In the second line, the ``-ibam`` parameter allows us to use the *bam* file as input. With the ``-g`` parameter, the tool will report the depth of coverage at each base on each reference in the genome file. If we look at the output file, we can see that each line consists of *reference sequence name*, *depth of coverage*, *number of bases with given depth*, *length of the reference*, and *fraction of bases on reference with given depth*, respectively. In the third line, thanks to the ``awk`` command, we can calculate the breadth of coverage, easily.

.. warning::

   A **0** depth means that there are unmapped regions in the reference. We must pay attention, if we want the calculate the breadth of coverage, properly.

The breadth of coverage can be calculated with ``1-$5`` in the rows where the second column equals *0*, that it, the fraction of the reference with any depths, except 0.

.. note::

   The condition ``0.2>$5`` allow us to choose references with higher coverage.

In the following lines, we combine the results in one file with the ``cat`` command and we sort the file numerically.
   
.. code-block:: bash
   :linenos:

   samtools idxstats [Output file].bam |awk -v OFS='\t' '{print $1, $2}' > [Output file].txt
   genomeCoverageBed -ibam [Output file].bam -g [Output file].txt > [Coverage output file].txt
   awk -v OFS='\t' '$2==0 && 0.2>$5 {print $1,1-$5}' [Coverage output file].txt > [Coverage summary output file].txt
   cat [Coverage summary output file].txt > [All coverage summary file].txt
   sort -k2nr [All coverage summary file].txt > [Sorted all coverage summary file].txt


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Aligning Reads with Best Candidates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We determine the two most likely references. We choose one of them according to the breadth of coverage results. The other one was suggested to us by the author of the reads, Dr. Terje Steinum. You can access the references with accession number `NC_025175.1`_ and `NC_025138.1`_, respectively.

.. _NC_025175.1: https://www.ncbi.nlm.nih.gov/nuccore/NC_025175.1
.. _NC_025138.1: https://www.ncbi.nlm.nih.gov/nuccore/NC_025138.1

----------------------------------------------
Libraries - Choosing which reads to assesemble
----------------------------------------------

After choosing the reference plasmids to be used in the rest of this study, we want to see which reads mapped to them. This will allow us to filter the reads, and keep only the ones that will have a better chance to be assembled properly. The idea is to filter out the reads that may be part of the host DNA.

In the first line, ``bwa`` creates an index for the reference file. In the second line, the alignment occurs. The ``-R`` parameter allows us to configure the header line of the output ``sam`` file. In the following line, the ``fixmate`` sub-command allows us to fill in the mate coordinates from a name sorted alignment. In the last line, ``fixmate`` output sorted with ``sort`` sub-command.

.. code-block:: bash
   :linenos:

   bwa index [Reference plasmid].fasta
   bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' [Reference plasmid].fasta [Forward fastq file] [Reverse fastq file] > [Output file].sam
   samtools fixmate -O bam [Output file].sam [Fixmate output file].bam
   samtools sort -O bam -o [Sorted fixmate output file].bam [Fixmate output file].bam

---------------
Filtering Reads
---------------

Filtering reads is an important process if you want to assemble your reads successfully. Firstly, we merge the three ``fixmatesorted.bam`` files that we got from the alignment process, using the ``samtools merge`` command. You can see the code below.

..  code-block:: bash
    :linenos:

    samtools merge [merged].bam [fixmatesorted1].bam [fixmatesorted2].bam [fixmatesorted3].bam

..  warning::
    
    You should indicate the output file in the first place. Otherwise, you will get an error.

..  note::

    We can use up to three different bam files. But ``trimmomatic`` produced four different outputs. In the ``Mapping`` process, two of them are combined as a paired file. Hence, we got 3 bam files named ``pairedfixmatesorted``, ``forward-unpairedfixmatesorted``, ``reverse-unpairedfixmatesorted``.


..  note::

    We merged ``BAM`` files in order to keep all reads.


After merging the ``BAM`` files, we want to see the depth of coverage for each position in the plasmid genome. ``samtools mpileup`` gives detailed output for this. You can look at the ``mpileup`` output in the list given below. Each line consists of 5 *tab-separated* columns. Column 6 is optional.

1. Sequence name
2. Position (starting from 1)
3. Reference nucleotide at that position
4. Depth of coverage
5. Bases at that position from aligned reads
6. Phred Quality of those bases (OPTIONAL).

We need only columns ``1,2,4``. These columns contain sequence name, position and depth of coverage, respectively. We combine ``samtools mpileup`` with ``awk`` command to carry out this purpose. 


.. code-block:: bash
   :linenos:

   samtools mpileup [merged].bam | awk '{print $1"\t"$2"\t"$4}' > [depth].txt


The ``[depth].txt`` file allows us to filter high coverage regions on the plasmid genome. We parsed the file using following python script and we create a new file that contains reads to assemble.


.. code-block:: python
   :linenos:

   from Bio import SeqIO
   import pandas as pd
   import matplotlib as plt
   import numpy as np

   records = [x for x in SeqIO.parse("[reference].gb", "genbank")]

   df=pd.read_csv('[depth].txt',sep='\t', header=None, names=["Ref","Position","Depth"])

   def before_or_after(x,low, high):
       if x < low:
           return(-1)
       elif x > high:
           return(1)
       else:
           return(0)

   def inside(regions, x0, x1):
       for low, high in regions:
           a = before_or_after(x0, low, high)
           b = before_or_after(x1, low, high)
           if a != b or (a==0 and b==0):
               return((low,high))
       return None

   zero_region = [(y,x) for x,y in zip(df.Position[1:],df.Position[:-1]) if (x-y)>1]
   
   l = df.Depth.quantile(0.25)
   u = df.Depth.quantile(0.75)
   iqr = u-l
   lower_limit = l-1.5*iqr
   upper_limit = u+1.5*iqr

   a = df.Depth < lower_limit
   b = [df.Position[i] for i in range(1,len(a)) if a[i-1]!=a[i]]
   b.insert(0,1)

   low_cover_region=[(b[i],b[i+1]) for i in range(len(b)-1) if i%2==0 and b[i+1]-b[i]> 50]

   c = df.Depth > upper_limit # u+1.5*iqr
   d = [df.Position[i] for i in range(1,len(c)) if c[i-1]!=c[i]]

   high_cover_region=[(d[i],d[i+1]) for i in range(len(d)-1) if i%2==0 and d[i+1]-d[i]> 50]

   with open("reads_to_keep.txt", "w") as outfile:
       with open("[reference].sam","rt") as sam:
           for line in sam:
               if line[0]=="@":
                   continue
               cols = line.strip().split()
               read_name = cols[0]
               start = int(cols[3])
               if start==0:
                   continue
               end = start + len(cols[9])
               loc = inside(high_cover_region, start, end)
               if loc is None:
                   print(read_name, file=outfile)

With the filtered reads, we create new fastq file by using the code below.

.. code-block:: python
   :linenos:
   :caption: sam2fastq.py
   :name: sam2fastq.py

   input_file = sys.stdin
   id_file = sys.argv[1]
   output_file = sys.stdout
   wanted = set(line.rstrip("\n").split(None, 1)[0] for line in open(id_file))
   print("Found %i unique identifiers in %s" % (len(wanted), id_file),file=sys.stderr)
   records = (r for r in SeqIO.parse(input_file, "fastq") if r.id in wanted)
   count = SeqIO.write(records, output_file, "fastq")
   print("Saved %i records from %s to %s" % (count, "input_file", "output_file"), file=sys.stderr)
   if count < len(wanted):
       print("Warning %i IDs not found in %s" % (len(wanted) - count, "input_file"), file=sys.stderr)


--------------
Phrap Assembly 
--------------

^^^^^^^^^^^^^^^^^^^^
Assembly Preparation
^^^^^^^^^^^^^^^^^^^^

Phrap assembly tool takes a fasta file and a quality file as input. So, we need to create these files from the fastq file which we filtered in the previous step. The following bash code carries out this process. The ``Fastq`` files consist of 4 lines per sequence.

1. First line starts with a ``@`` character and this is followed by sequence identifier.
2. Second line contains the sequence of the read.
3. Third line contains only a ``+`` sign to separate raw sequences and quality scores of each read.
4. Fourth line consists of quality scores of reads corresponding to the sequences in the second line.

In the first line of code, we select first and second column to create a ``fasta`` file. Also, we replaced ``@`` symbol with ``>``. In the following line, we select only first and fourth column to create a ``quality file``.

.. code:: bash
	 :linenos:
   
   cat [filtered].fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > [filtered].fasta
   cat [filtered].fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-4 | tr '\t' '\n' > [filtered].qual

^^^^^^^^^^^^^^^^
Assembly Process
^^^^^^^^^^^^^^^^

To start the assembly with  ``Phrap``, we need to run the following code. With the ``-ace`` parameter, we will get ``ace`` files for the output.

.. code:: bash
   :linenos:

   phrap -ace [Library name]


.. warning::
   
   Fasta file and quality file names should be the same. Thereby, ``Phrap`` detects the file automatically.


.. _Quast:

^^^^^^^^^^^^^^^^^^^
Assembly Statistics
^^^^^^^^^^^^^^^^^^^

For the general statistics we want to see, we used ``Quast``. It is stand for Quality Assessment Tool for Genome Assemblies. The code we used for this is given below.

.. code:: bash
   :linenos:

   quast.py [contigs].fasta -o [Output folder]


---------------
SPAdes Assembly
---------------

^^^^^^^^^^^^^^^^^^^^
Assembly Preparation
^^^^^^^^^^^^^^^^^^^^

After the first assembly with ``Phrap``, we want to extend our contigs with the reads that we did not use in the first. For accomplish this we used ``samtools view`` with ``f`` parameter. You can look at the general code in the below.

.. code:: bash
   :linenos:

   samtools view -f 4 [Reference fixmatesorted].bam | cut -f 1 > [no align read ids].txt
   zcat [trimmed fastq].fq.gz | python3 sam2fastq.py [no align read ids].txt | gzip -c > [ready2assemble].fastq.gz

See :ref:`sam2fastq.py`.


^^^^^^^^^^^^^^^^
Assembly Process
^^^^^^^^^^^^^^^^

We used ``SPAdes`` as a second assembler. With the ``--trusted-contig`` parameter, we used contigs that we got in the first assembly as a base and we try to extend this contigs with we did not use before.

.. code:: bash
   :linenos:

   spades.py -o [Output folder] --only-assembler -1 [Forward ready2assemble fastq file] -2 [Reverse ready2assemble fastq file] --s1 [Singles ready2assemble 1] --s2 [Singles ready2assemble 2] --trusted-contigs [contigs].fasta


^^^^^^^^^^^^^^^^^^^
Assembly Statistics
^^^^^^^^^^^^^^^^^^^

We used ``Quast`` again for the statistics about assembly. See :ref:`Quast`.


------
Blastn
------

For scaffolding contigs that we get, we need to see *nucleotide blast* or *blastn* results. For this, we use following code in below. If we look at the parameters, we can see that ``query`` parameter for our contigs. The ``subject`` parameter for the reference plasmid. The ``outfmt`` parameter for the alignment view option and we used *6* for this which is *tabular*. The ``out`` parameter for the output file and the ``evalue`` parameter for the expectation value (E) threshold for saving hits.

.. code-block:: bash
   :linenos:

   blastn -query [spades].contigs -subject [Reference].fasta -outfmt 6 -out [assembly-reference.blastn].txt -evalue 1e-15


Also, you can get the mapped contigs with reference using the following bash code.

.. code-block:: bash
   :linenos:

   cut -f 1 [assembly-reference.blastn].txt |uniq > [contigsid].txt

   awk 'BEGIN{while((getline<"contigsid.txt")>0)l[">"$1]=1}/^>/{f=l[$1]}f' contigs.fasta > [mapped contigs].fasta


The following code gives us the order of contigs by positions, ascendingly.

.. code-block:: bash
   :linenos: 

   awk 'BEGIN {OFS="\t"} $9>$10 {x=$9; $9=$10; $10=x; print $0,"-"; next} {print $0,"+"}' [assembly-reference.blastn].txt |sort -k 9n |less -S > [Ordered Contigs].txt