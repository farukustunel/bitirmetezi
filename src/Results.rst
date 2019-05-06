=======
Results
=======

As you can see in the table below, we have **343434** reads in the sample taken from *Escherichia coli* and *251128* reads for *Citrobacter freundii*. If you  compare each library in proportion to the total number of bases, you will see that *Escherichia coli* library has more nucleotides and more reads.

.. list-table:: Raw Data Statistics
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - File Name
     - Number of Reads
     - Number of Bases
   * - 1-F5-96_S1_L001_R1_001
     - 171717
     - 38124416
   * - 1-F5-96_S1_L001_R2_001
     - 171717
     - 38253111
   * - 4-F20-96_S2_L001_R1_001
     - 125564
     - 27633180
   * - 4-F20-96_S2_L001_R2_001
     - 125564
     - 27727284


After trimming, we get 4 different outputs for each library. For each one of them, we have forward and reverse reads. For those reads, we have paired and unpaired reads. 

For *Escherichia coli* library, ``trimmomatic`` drops **443** reads for forward strand and drops **8394** reads for reverse strand. If we look at *Citrobacter freundii*, **321** reads are dropped by ``trimmomatic`` for forward strand and **5926** reads are dropped for reverse strand. You can see the detail numbers in the table below.


.. list-table:: Statistics After Trimming
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - File Name
     - Number of Reads
     - Number of Bases
   * - F5_forward_paired
     - 163188
     - 35749448
   * - F5_forward_unpaired
     - 8086
     - 1514545
   * - F5 forward dropped
     - 443
     - NA
   * - F5_reverse_paired
     - 163188
     - 34123007
   * - F5_reverse_unpaired
     - 135
     - 21546
   * - F20_forward_paired
     - 119552
     - 25955781
   * - F20_forward_unpaired
     - 5691
     - 1059660
   * - F20_reverse_paired
     - 119552
     - 24575680
   * - F20_reverse_unpaired
     - 86
     - 11870


Also, you can see the comparasion of the raw reads and the trimmed reads in the histogram below.

.. image:: ../_static/numberofreads.png
   :align: left


We aligned our reads with ``NCBI RefSeq Database`` plasmids. You can see in the table below the total number of records of ``NCBI RefSeq Database`` are **15076**.

.. list-table:: Number of NCBI RefSeq Records
   :widths: 25 25 
   :header-rows: 1
   :align: right

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
   * - Total
     - 15076


The summary of the results is listed in the below table. We put a threshold for the mapping reads number which is ``1000```. We use this number because we want to clear results. The threshold can be changed from the person who wants to follow this protocol. Also, we decided that at least half of the reference genome should be covered by our reads. Hence, the threshold for the breadth of coverage is ``50%``.

.. list-table:: Coverage Statistics
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - Library
     - Plasmids with >1K Mapping Reads
     - Plasmids where Breadth of Coverage >50%
   * - *Escherichia coli* (F5)
     - 173
     - 45
   * - *Citrobacter freundii* (F20)
     - 135    
     - 36 

These results can be interpreted wrongly, because we do not know the answers to the following questions:

1. Are the mapped reads separated in the  plasmid? 
2. Are the reads only mapped to particular regions of the reference plasmid?

The breadth of coverage results give us the answers of questions listed above. The following table shows the plasmids in the database with higher breadth of coverage.

.. list-table:: Top 5 Best Candidates
   :widths: 25 25 25
   :header-rows: 1
   :align: left

   * - Library
     - Accession
     - Breadth of Coverage (%)
   * - *Escherichia coli* (F5)
     - NC_025175.1  
     - 90.2%
   * - *Escherichia coli* (F5)
     - NC_024956.1  
     - 88.3%
   * - *Escherichia coli* (F5)
     - NC_025139.1  
     - 79.7%
   * - *Citrobacter freundii* (F20)
     - NC_019049.1  
     - 77.7%
   * - *Escherichia coli* (F5)
     - NC_004998.1  
     - 60.2%
      
We choose the plasmid with higher breadth of coverage (``NC_025175.1``) for further study. We also included in this study the plasmid ``NC_025138.1``, that was suggested by one of the authors of the original paper, who also isolated the samples we study.


We find the number of mapping reads using ``samtools view -F 4`` command. As you can see in the table below, more reads map to the plasmid ``NC_025175.1``  than to ``NC_025138.1``, in both libraries. The difference between the number of mapping reads is **243443** for *Escherichia coli* (F5) and **165817** for *Citrobacter freundii* (F20).

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
   * - *Escherichia coli* (F5)
     - NC_025138.1  
     - 33724
   * - *Citrobacter freundii* (F20)
     - NC_025175.1  
     - 187823
   * - *Citrobacter freundii* (F20)
     - NC_025138.1  
     - 22006
  


  
