#!/usr/bin/env bash
set -e
source config.sh

mkdir -p ${OUTDIR}/{fastp,bam,gvcf,vcf,logs}

bash scripts/02_bwa_index.sh

while read SAMPLE; do
    bash scripts/trimming.sh ${SAMPLE}
    bash scripts/bwa_align.sh ${SAMPLE}
    bash scripts/dedup.sh ${SAMPLE}
    bash scripts/haplotyping.sh ${SAMPLE}
done < ${SAMPLES}

bash scripts/genotyping.sh
bash scripts/filtration.sh

echo "Pipeline completed successfully."
