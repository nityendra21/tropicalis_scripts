#!/bin/bash

source config.sh
SAMPLE=$1

samtools depth \
    -aa -t ${THREADS} \
    -Q ${MIN_BASEQ} \
    -q ${MIN_MAPQ} \
    ../${OUTDIR}/bam/${SAMPLE}.dedup.bam \
    > ${CNV_OUTDIR}/depth/${SAMPLE}.depth
