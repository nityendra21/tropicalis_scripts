#!/usr/bin/env bash
source config.sh

ls ${OUTDIR}/gvcf/*.g.vcf.gz > gvcf.list

${GATK} CombineGVCFs \
  -R ${REF} \
  -V gvcf.list \
  -O ${OUTDIR}/vcf/combined.g.vcf.gz

${GATK} GenotypeGVCFs \
  -R ${REF} \
  -V ${OUTDIR}/vcf/combined.g.vcf.gz \
  -O ${OUTDIR}/vcf/raw_variants.vcf.gz
