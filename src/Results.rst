=======
Results
=======

As you can see the below table, we have **343434** reads for *Escherichia coli* and *251128* reads for *Citrobacter freundii*. If you can compare each library in proportion to the total number of bases, you will see that *Escherichia coli* library has more nucleotides and more reads.

.. list-table:: Raw Data Statistics
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - File Name
     - Number of Reads
     - Number of Bases
   * - 1-F5-96_S1_L001_R1_001.fastq.gz
     - 171717
     - 38124416
   * - 1-F5-96_S1_L001_R2_001.fastq.gz
     - 171717
     - 38253111
   * - 4-F20-96_S2_L001_R1_001.fastq.gz
     - 125564
     - 27633180
   * - 4-F20-96_S2_L001_R2_001.fastq.gz
     - 125564
     - 27727284


After trimming processes, we get 4 different output for each library. For each them, we have forward and reverse reads. For those reads, we have paired and unpaired reads. For *Escherichia coli* library, ``trimmomatic`` drops **443** reads for forward strand and drops **8394** reads for reverse strand. If we have a look *Citrobacter freundii*, **321** reads are dropped by ``trimmomatic`` for forward strand and **5926** reads are dropped for reverse strand. You can see the detail numbers in the table below.


.. list-table:: Statistics After Trimming
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - File Name
     - Number of Reads
     - Number of Bases
   * - F5_forward_paired.fq.gz
     - 163188
     - 35749448
   * - F5_forward_unpaired.fq.gz
     - 8086
     - 1514545
   * - F5_reverse_paired.fq.gz
     - 163188
     - 34123007
   * - F5_reverse_unpaired.fq.gz
     - 135
     - 21546
   * - F20_forward_paired.fq.gz
     - 119552
     - 25955781
   * - F20_forward_unpaired.fq.gz
     - 5691
     - 1059660
   * - F20_reverse_paired.fq.gz
     - 119552
     - 24575680
   * - F20_reverse_unpaired.fq.gz
     - 86
     - 11870


Also, you can see the comparasion of the raw reads and the remainder of reads after trimming process in the histogram below.

.. image:: ../_static/numberofreads.png
   :align: left


We aligned our read with ``NCBI RefSeq Database`` plasmids. You can see the table below, the total number of records of ``NCBI RefSeq Database`` are **15076**.

.. list-table:: Number of NCBI RefSeq Records
   :widths: 25 25 
   :header-rows: 1
   :align: left

   * - File Name
     - Number of Recors
   * - plasmid.1.1.genomic.fna
     - 4597
   * - plasmid.2.1.genomic.fna
     - 3234
   * - plasmid.3.1.genomic.fna
     - 2524
   * - plasmid.4.1.genomic.fna
     - 3023
   * - plasmid.5.1.genomic.fna
     - 1698


The summary of the results is listed in the below table. We put a threshold for the mapping reads number which is ``1000```. We use this number because we want to clear results. The threshold can be changed from the person who wants to follow this protocol. Also, we decided that at least half of the reference genome should be covered by our reads. Hence, the threshold for the breadth of coverage is ``50%``.

.. list-table:: Coverage Statistics
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - Library
     - Mapping Reads (>1000)
     - Breadth of Coverage (>50%)
   * - *Escherichia coli* (F5)
     - 173
     - 45
   * - *Citrobacter freundii* (F20)
     - 135    
     - 36 

But, these results can be interpreted wrongly. Because we do not know the answers of the following questions. 

1. Are the mapped reads separated the whole genome? 
2. Are the reads only mapped with particular regions of reference plasmid?

The breadth of coverage results can give us the answers of questions listed above.

.. list-table:: Top 5 Best Candidates
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - Library
     - Accession
     - Breadth of Coverage (%)
   * - *Escherichia coli* (F5)
     - NC_025175.1  
     - 90.2253%
   * - *Escherichia coli* (F5)
     - NC_024956.1  
     - 88.3060%
   * - *Escherichia coli* (F5)
     - NC_025139.1  
     - 79.7261%
   * - *Citrobacter freundii* (F20)
     - NC_019049.1  
     - 77.6670%
   * - *Escherichia coli* (F5)
     - NC_004998.1  
     - 60.1557%
      

The author of the original paper suggests that a reference with ``NC_025138.1`` NCBI accession number.


We find mapping reads number with using ``samtools view -F 4`` command. As you can see in the table below, ``NC_025175.1`` plasmid has more mapping reads than ``NC_025138.1`` plasmid comparing to both libraries. The difference between the number of mapping reads is **243443** for *Escherichia coli* (F5) and **165817** for *Citrobacter freundii* (F20).

.. list-table:: Number of Mapping Reads with Best Candidates
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - Library
     - Accession
     - Number of Mapping Reads
   * - *Escherichia coli* (F5)
     - NC_025175.1  
     - 277167
   * - *Citrobacter freundii* (F20)
     - NC_025175.1  
     - 187823
   * - *Escherichia coli* (F5)
     - NC_025138.1  
     - 33724
   * - *Citrobacter freundii* (F20)
     - NC_025138.1  
     - 22006
  


  