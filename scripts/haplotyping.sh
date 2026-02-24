#!/usr/bin/env bash
source config.sh
SAMPLE=$1

${GATK} HaplotypeCaller \
  -R ${REF} \
  -I ${OUTDIR}/bam/${SAMPLE}.dedup.bam \
  -O ${OUTDIR}/gvcf/${SAMPLE}.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads ${THREADS} \
  -ploidy ${PLOIDY} \
  -stand-call-conf 30 \
  --pcr-indel-model NONE
