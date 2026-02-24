#!/usr/bin/env bash
source config.sh

${GATK} VariantFiltration \
  -R ${REF} \
  -V ${OUTDIR}/vcf/raw_variants.vcf.gz \
  -O ${OUTDIR}/vcf/filtered_variants.vcf.gz \
  --filter-name "QD_filter" --filter-expression "QD < 2.0" \
  --filter-name "FS_filter" --filter-expression "FS > 60.0" \
  --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
  --filter-name "DP_filter" --filter-expression "DP < 10"

# Extracting SNPs 
bcftools filter -i 'TYPE="snp" && (FMT/GQ)>20' \
   -Oz --threads ${THREADS} \
   -o ${OUTDIR}/vcf/snps.vcf.gz \
   ${OUTDIR}/vcf/raw_variants.vcf.gz 
