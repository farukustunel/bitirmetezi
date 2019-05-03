=======
Results
=======

As you can see the below table, we have **343434** reads for *Escherichia coli* and *251128* reads for *Citrobacter freundii*. If you can compare each library in proportion to the total number of bases, you will see that *Escherichia coli* library has more nucleotides, thus has more reads.

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


After trimming processes, we get 4 different output for each library. For each library, we have forward and reverse reads. For forward and reverse reads, we have paired and unpaired reads. For *Escherichia coli* library, ``trimmomatic`` drops **443** reads for forward strand and drops **8394** reads for reverse strand. If we have a look *Citrobacter freundii*, **321** reads are dropped by ``trimmomatic`` for forward strand and **5926** reads are dropped for reverse strand. You can see the detail numbers in the table below.


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
