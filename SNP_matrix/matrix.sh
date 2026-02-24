#!/bin/bash

# Creating SNP matrix. 
# Splitting the combined VCF file into individual files.
mkdir -p VCFs
bcftools query -l snps.vcf.gz | parallel -j 40 'vcf-subset --exclude-ref -c {} snps.vcf.gz > VCFs/${}.vcf'

# Create matrix
git clone https://github.com/broadinstitute/broad-fungalgroup
cd broad-fungalgroup/scripts/SNPs
ls ../../../VCFs/*.vcf > vcf.txt

python3 vcfSnpsToFasta.py vcf.txt > tropicalis.fa
perl fasta2snpcounts.pl tropicalis.fa > tropicalis_matrix.tsv
