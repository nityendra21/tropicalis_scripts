#!/usr/bin/env bash
source config.sh
SAMPLE=$1

${GATK} MarkDuplicates \
  -I ${OUTDIR}/bam/${SAMPLE}.sorted.bam \
  -O ${OUTDIR}/bam/${SAMPLE}.dedup.bam \
  -M ${OUTDIR}/logs/${SAMPLE}.metrics.txt \
  --CREATE_INDEX true
