#!/bin/bash

SAMPLEID=`(ls reads/* | sed 's/reads\///g' | sed 's/\.fastq//g')` #name sample ID as filename of raw reads without .fastq extension

PERL=./rnaseq_clipper.pl #the location of the rnaseq_clipper.pl script from UT Austin
REF=/data/raw/maria/Zonotrichia_ref/GCF_000385455.1_Zonotrichia_albicollis-1.0.1_genomic.fna #the location of the reference genome
INPUT_DIR=/scratch/1/maria/TagSeq #the input dir should have a /reads /trim and /sam folder

for i in $SAMPLEID
do echo "trimming reads for  $i" 
perl ${PERL} ${INPUT_DIR}/reads/${i}.fastq [ATGC]?[ATGC][AC][AT]GGG+ | fastx_clipper -a AAAAAAAA -l 20 -Q33 | fastx_clipper -a AGATCGGAAG -l 20 -Q33 | fastq_quality_filter -Q33 -q 20 -p 90 > ${INPUT_DIR}/trim2/${i}.trim
done

for i in $SAMPLEID
do echo "bwa for $i"
bwa mem -t 10 ${REF} ${INPUT_DIR}/trim2/${i}.trim > ${INPUT_DIR}/sam2/${i}.sam
done

featureCounts -T 10 -F SAF -Q 20 -g 'gene_id' -a /scratch/1/maria/Zonotrichia_ref/GCF_000385455.1_Zonotrichia_albicollis-1.0.1_genomic.saf -o Junco_longacc_tagseq_9Feb2020.txt `(ls ${INPUT_DIR}/sam2/*.sam)`

exit 0;
