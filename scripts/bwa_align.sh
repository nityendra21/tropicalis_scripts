#!/usr/bin/env bash

source config.sh

if [ ! -f ${REF}.bwt ]; then
    ${BWA} index ${REF}
fi

if [ ! -f ${REF}.fai ]; then
    ${SAMTOOLS} faidx ${REF}
fi

if [ ! -f ${REF%.fa}.dict ]; then
    ${GATK} CreateSequenceDictionary -R ${REF}
fi


SAMPLE=$1

${BWA} mem -t ${THREADS} -R '@RG\tID:'{SAMPLE}'\tSM:'${SAMPLE}'\tPL:illumina' \
    ${REF} \
    ${OUTDIR}/fastp/${SAMPLE}_R1.trimmed.fastq.gz \
    ${OUTDIR}/fastp/${SAMPLE}_R2.trimmed.fastq.gz | \
${SAMTOOLS} sort -@ ${THREADS} -o ${OUTDIR}/bam/${SAMPLE}.sorted.bam

${SAMTOOLS} index ${OUTDIR}/bam/${SAMPLE}.sorted.bam
