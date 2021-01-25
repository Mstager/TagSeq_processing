# TagSeq_processing

These scripts trim, filter, and map TagSeq reads to a reference, then generate a read counts table for downstream analysis. The read counts table is then used to test for differential expression (after filtering and normalization) using a GLM-based method with the package edgeR in R. They were used in:

<i/> Stager & Cheviron 2020. Is there a role for sarcolipin in avian facultative thermogenesis in extreme cold? Biology Letters 16: 20200078. https://royalsocietypublishing.org/doi/full/10.1098/rsbl.2020.0078

Raw reads are available from NCBI: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA612334

Junco reads were mapped to the annotated Zonotrichia albicollis genome v. 1.0.1 available here: https://www.ncbi.nlm.nih.gov/genome/?term=txid44394[Organism:noexp]
