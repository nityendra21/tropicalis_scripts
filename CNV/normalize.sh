#!/usr/bin/env bash
set -e

INPUT=$1
SAMPLE=$2

mkdir -p results/{normalized,desmoothed,byplot,islands,plot}

echo "Step 1: Normalizing coverage"
samtools view -@ ${THREADS} -h ${INPUT} | \
    python coverageFrame.py -f ${WINDOW} > results/normalized/${SAMPLE}.normalized.tsv

echo "Step 2: Removing smiley pattern"
Rscript splint.R \
    input=results/normalized/${SAMPLE}.normalized.tsv \
    tab=results/desmoothed/${SAMPLE}.desmoothed.tsv \
    byplot=results/byplot/${SAMPLE}.byplot.pdf \
    islands=results/islands/${SAMPLE}.islands.tsv \
    plot=results/plot/${SAMPLE}.pdf

echo "Step 3: Final CNV analysis"
Rscript plot_CNV.R \
    -i results/desmoothed/${SAMPLE}.desmoothed.tsv \
    -o results/final/${SAMPLE}.png \
    -t ${SAMPLE}
