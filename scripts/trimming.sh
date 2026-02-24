#!/usr/bin/env bash
source config.sh
SAMPLE=$1

${FASTP} \
  -i ${SAMPLE}_R1.fastq.gz \
  -I ${SAMPLE}_R2.fastq.gz \
  -o ${OUTDIR}/fastp/${SAMPLE}_R1.trimmed.fastq.gz \
  -O ${OUTDIR}/fastp/${SAMPLE}_R2.trimmed.fastq.gz \
  --thread ${THREADS} \
  --cut_front --cut_tail \ 
  --cut_window_size 4 --detect_adapter_for_pe \ 
  --qualified_quality_phred 30 --length_required 40 \
  --html ${OUTDIR}/fastp/${SAMPLE}.html \
  --json ${OUTDIR}/fastp/${SAMPLE}.json \
